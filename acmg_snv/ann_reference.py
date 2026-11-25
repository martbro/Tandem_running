"""Direct base access to an indexed, uncompressed local GRCh38 FASTA."""

import os
from pathlib import Path
import tempfile
from threading import Lock

from .model import CHROMOSOMES, REFSEQ


def primary_chromosome(name):
    """Map a primary chromosome or GRCh38 RefSeq accession to a label."""
    name = name.upper()
    chrom = name.removeprefix("CHR")
    if chrom in CHROMOSOMES:
        return chrom
    return REFSEQ.get(name)


def _fai_entries(path):
    with open(path, encoding="utf-8") as index:
        for row in index:
            fields = row.rstrip("\r\n").split("\t")
            if len(fields) < 5:
                raise ValueError("Ligne FASTA .fai mal formée.")
            chrom = primary_chromosome(fields[0])
            if chrom is None:
                continue
            length, offset, bases, width = map(int, fields[1:5])
            if length < 1 or offset < 0 or bases < 1 or width < bases:
                raise ValueError(f"Index FASTA .fai invalide ({fields[0]}).")
            yield chrom, (length, offset, bases, width)


class Reference:
    """Read reference alleles without loading the FASTA into memory."""

    def __init__(self, path):
        self.path = Path(path)
        entries = list(_fai_entries(f"{self.path}.fai"))
        self.entries = dict(entries)
        if not self.entries:
            raise ValueError("Aucun chromosome nucléaire primaire du GRCh38.")
        if len(entries) != len(self.entries):
            raise ValueError("Un chromosome primaire apparaît deux fois.")
        self.lock = Lock()
        self.handle = self.path.open("rb")
        if self.handle.read(2) == b"\x1f\x8b":
            self.handle.close()
            raise ValueError("Décompressez le FASTA avant de l'utiliser.")

    def base(self, chrom, pos):
        """Return the forward-strand base at a 1-based position."""
        if chrom not in self.entries:
            raise ValueError(f"Chromosome {chrom} absent de l'index FASTA.")
        length, start, bases, width = self.entries[chrom]
        if not 1 <= pos <= length:
            raise ValueError(f"Position {chrom}:{pos} hors du FASTA.")
        line, column = divmod(pos - 1, bases)
        with self.lock:
            self.handle.seek(start + line * width + column)
            base = self.handle.read(1).upper()
        if base not in {b"A", b"C", b"G", b"T", b"N"}:
            raise ValueError(f"Base FASTA illisible à {chrom}:{pos}.")
        return base.decode("ascii")

    def check(self, variant):
        """Raise ValueError if the supplied reference allele mismatches."""
        observed = self.base(variant.chrom, variant.pos)
        if observed != variant.ref:
            raise ValueError(
                f"REF {variant.ref} incompatible avec le FASTA ({observed}) "
                f"pour {variant.key}."
            )

    def close(self):
        with self.lock:
            self.handle.close()


class _FastaEntry:
    def __init__(self, name, offset):
        self.name = name
        self.offset = offset
        self.length = 0
        self.bases = self.width = self.last_bases = self.last_width = 0

    def append(self, line):
        bases = len(line.rstrip(b"\r\n"))
        width = len(line)
        if bases == 0:
            raise ValueError(f"Ligne vide dans le FASTA ({self.name}).")
        if self.length == 0:
            self.bases, self.width = bases, width
        regular = (self.last_bases, self.last_width) == (self.bases, self.width)
        too_large = bases > self.bases or width > self.width
        if (self.length and not regular) or too_large:
            raise ValueError(f"Lignes FASTA irrégulières ({self.name}).")
        self.length += bases
        self.last_bases, self.last_width = bases, width

    def to_fai(self):
        if self.length == 0:
            raise ValueError(f"Séquence FASTA vide ({self.name}).")
        fields = (self.name, self.length, self.offset, self.bases, self.width)
        return "\t".join(map(str, fields)) + "\n"


def _fasta_lines(path):
    entry = None
    with open(path, "rb") as fasta:
        if fasta.read(2) == b"\x1f\x8b":
            raise ValueError("Décompressez le FASTA avant de l'indexer.")
        fasta.seek(0)
        for line in iter(fasta.readline, b""):
            if line.startswith(b">"):
                if entry:
                    yield entry.to_fai()
                fields = line[1:].split(maxsplit=1)
                if not fields:
                    raise ValueError("En-tête FASTA sans identifiant.")
                entry = _FastaEntry(fields[0].decode("ascii"), fasta.tell())
            elif entry is None:
                raise ValueError("L'en-tête FASTA initial manque.")
            else:
                entry.append(line)
        if not entry:
            raise ValueError("Fichier FASTA vide.")
        yield entry.to_fai()


def build_fasta_index(path):
    """Create a samtools-compatible FASTA index if missing."""
    output = Path(f"{path}.fai")
    if output.exists():
        raise FileExistsError(f"L'index existe déjà ({output}).")
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="\n", delete=False,
        dir=output.parent,
    ) as index:
        scratch = Path(index.name)
        try:
            index.writelines(_fasta_lines(path))
            index.close()
            os.replace(scratch, output)
        finally:
            index.close()
            scratch.unlink(missing_ok=True)
    return output
