"""Shared SNV and report values for the local classifier."""

import re
from dataclasses import dataclass


INHERITANCES = frozenset({"AD", "AR", "XL", "XLD", "XLR", "YL", "OTHER"})
CLASSES = ("B", "LB", "VUS", "LP", "P")
STRENGTHS = {"supporting": 1, "moderate": 2, "strong": 4, "very_strong": 8}
CHROMOSOMES = {str(number) for number in range(1, 23)} | {"X", "Y"}
ACCESSIONS = (
    11, 12, 12, 12, 10, 12, 14, 11, 12, 11, 10, 12,
    11, 9, 10, 10, 11, 10, 10, 11, 9, 11, 11, 10,
)
REFSEQ = {
    f"NC_{index:06d}.{version}": chrom
    for index, (chrom, version) in enumerate(
        zip([str(n) for n in range(1, 23)] + ["X", "Y"], ACCESSIONS), 1
    )
}
BARE = re.compile(r"(?:chr)?([\dXY]+)\s*[:\t, ]\s*(\d+)", re.I)
HGVS = re.compile(r"(NC_\d+\.\d+):g\.(\d+)([ACGT])>([ACGT])", re.I)


@dataclass(frozen=True, slots=True)
class Variant:
    """A one-based, forward-strand substitution on a GRCh38 chromosome."""

    chrom: str
    pos: int
    ref: str
    alt: str

    def __post_init__(self):
        chrom = str(self.chrom).upper().removeprefix("CHR")
        ref, alt = self.ref.upper(), self.alt.upper()
        if chrom not in CHROMOSOMES:
            raise ValueError("Chromosome nucléaire primaire GRCh38 requis.")
        if type(self.pos) is not int or self.pos < 1:
            raise ValueError("Position GRCh38 entière, à partir de 1.")
        if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
            raise ValueError("REF et ALT doivent être distincts et dans ACGT.")
        if len(ref) != 1 or len(alt) != 1:
            raise ValueError("REF et ALT doivent mesurer un nucléotide.")
        object.__setattr__(self, "chrom", chrom)
        object.__setattr__(self, "ref", ref)
        object.__setattr__(self, "alt", alt)

    @classmethod
    def parse(cls, text):
        """Accept a structured chromosome coordinate or genomic HGVS."""
        value = text.strip()
        match = HGVS.fullmatch(value)
        if match:
            accession, pos, ref, alt = match.groups()
            chrom = REFSEQ.get(accession.upper())
            if not chrom:
                raise ValueError("Accession RefSeq GRCh38 non reconnue.")
            return cls(chrom, int(pos), ref, alt)
        match = BARE.match(value)
        if match:
            chrom, pos = match.groups()
            alleles = re.fullmatch(
                r"[\s:,]+([ACGT])[\s:,>]+([ACGT])", value[match.end():], re.I
            )
            if alleles:
                return cls(chrom, int(pos), *alleles.groups())
        raise ValueError("Format attendu : chr1:123:A:G ou HGVS génomique.")

    @property
    def key(self):
        return f"chr{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    def as_dict(self):
        return {
            "chrom": self.chrom, "pos": self.pos,
            "ref": self.ref, "alt": self.alt,
            "key": self.key, "assembly": "GRCh38",
        }
