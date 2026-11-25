"""Derive SNV consequences from the RefSeq MANE genomic GFF3 and FASTA."""

import itertools
import json
from pathlib import Path
from urllib.parse import unquote

from .ann_reference import primary_chromosome
from .annotations import new_db, preferred_annotations, readonly_db, text_input


_BIN_SIZE = 100_000
_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODE = {"".join(codon): aa
         for codon, aa in zip(itertools.product(_BASES, repeat=3), _AA)}
_LONG_AA = dict(zip("ACDEFGHIKLMNPQRSTVWY*", (
    "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
    "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr", "Ter",
)))
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
_MANE_TAGS = ("MANE Select", "MANE Plus Clinical")


def _attributes(text):
    attrs = {}
    for entry in text.split(";"):
        key, sep, value = entry.partition("=")
        if sep:
            attrs[key] = unquote(value)
    return attrs


def _feature_rows(path):
    with text_input(path) as source:
        for line in source:
            if line.startswith("#!genome-build ") and "GRCh38" not in line:
                raise ValueError("Le GFF3 MANE doit utiliser GRCh38.")
            if line.startswith("##FASTA"):
                return
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError("Ligne GFF3 mal formée.")
            fields[0] = primary_chromosome(fields[0])
            if fields[0] and fields[2] in {"mRNA", "exon", "CDS"}:
                yield fields, _attributes(fields[-1])


def _transcript_record(fields, attrs):
    chrom, _, _, first, last, _, strand, _, _ = fields
    tags = attrs.get("tag", "").split(",")
    mane = next((tag for tag in _MANE_TAGS if tag in tags), None)
    transcript = attrs.get("transcript_id", attrs.get("Name", ""))
    if not mane or not transcript.startswith(("NM_", "XM_")):
        return None
    if strand not in {"+", "-"} or not attrs.get("gene"):
        raise ValueError(f"GFF3 MANE incomplet pour {transcript}.")
    return {
        "chrom": chrom, "start": int(first), "end": int(last),
        "strand": strand, "gene": attrs["gene"].upper(),
        "transcript": transcript, "mane": mane,
        "exon": [], "CDS": [],
        "partial": attrs.get("partial") == "true",
    }


def _assign_segments(transcripts, pending):
    for fields, attrs in pending:
        chrom, _, feature, start, end, _, strand, phase, _ = fields
        want = int(start), int(end)
        segment = (*want, int(phase)) if feature == "CDS" else want
        for parent in attrs.get("Parent", "").split(","):
            record = transcripts.get(parent)
            if record and record["chrom"] == chrom:
                if record["strand"] != strand:
                    raise ValueError("Brins contradictoires dans le GFF3.")
                record[feature].append(segment)
                if _atypical_segment(attrs):
                    record["partial"] = True


def _atypical_segment(attrs):
    return (
        attrs.get("partial") == "true"
        or bool(attrs.get("exception") or attrs.get("transl_except"))
    )


def _read_transcripts(path):
    transcripts = {}
    pending = []
    for fields, attrs in _feature_rows(path):
        if fields[2] == "mRNA":
            record = _transcript_record(fields, attrs)
            if record and attrs.get("ID"):
                key = attrs["ID"]
                if key in transcripts:
                    raise ValueError(f"Identifiant MANE dupliqué ({key}).")
                transcripts[key] = record
        else:
            pending.append((fields[:8] + [""], attrs))
    if not transcripts:
        raise ValueError("Aucun transcrit du GFF3 RefSeq MANE n'a été reconnu.")
    _assign_segments(transcripts, pending)
    return transcripts.values()


def _prepare_transcript(record, release):
    record["release"] = release
    record["exon"] = sorted(set(record["exon"]))
    cds = sorted(set(record["CDS"]), reverse=record["strand"] == "-")
    record["CDS"] = cds
    if not record["exon"] or not cds:
        raise ValueError(f"Exons ou CDS absents pour {record['transcript']}.")
    valid = not record["partial"]
    length = 0
    for start, end, phase in cds:
        valid &= phase == (-length % 3)
        valid &= any(a <= start <= end <= b for a, b in record["exon"])
        length += end - start + 1
    record["complete_cds"] = bool(valid and length % 3 == 0 and length >= 6)
    return record


def _bins(record):
    start = (record["start"] - 1) // _BIN_SIZE
    end = (record["end"] - 1) // _BIN_SIZE
    return range(start, end + 1)


def _insert_transcript(database, transcript, release):
    record = _prepare_transcript(transcript, release)
    payload = json.dumps(record, ensure_ascii=False)
    cursor = database.execute(
        "INSERT INTO transcript (start, end, payload) VALUES (?, ?, ?)",
        (record["start"], record["end"], payload),
    )
    database.executemany(
        "INSERT INTO bin VALUES (?, ?, ?)",
        ((record["chrom"], slot, cursor.lastrowid) for slot in _bins(record)),
    )


def build_mane_db(input_path, output_path):
    """Create a local transcript and splice-junction index from MANE GFF3."""
    transcripts = _read_transcripts(input_path)
    with new_db(output_path) as database:
        database.execute(
            "CREATE TABLE transcript (id INTEGER PRIMARY KEY, start INTEGER,"
            " end INTEGER, payload TEXT)"
        )
        database.execute(
            "CREATE TABLE bin (chrom TEXT, slot INTEGER, transcript INTEGER,"
            " PRIMARY KEY (chrom, slot, transcript)) WITHOUT ROWID"
        )
        for transcript in transcripts:
            _insert_transcript(database, transcript, Path(input_path).name)


def _coding_index(record, pos):
    index = 0
    for start, end, _ in record["CDS"]:
        if start <= pos <= end:
            offset = pos - start if record["strand"] == "+" else end - pos
            return index + offset
        index += end - start + 1
    return None


def _genomic_position(record, coding_index):
    for start, end, _ in record["CDS"]:
        length = end - start + 1
        if coding_index < length:
            if record["strand"] == "+":
                return start + coding_index
            return end - coding_index
        coding_index -= length
    raise ValueError("La position codante n'est pas présente dans le CDS.")


def _orient(base, record):
    return base if record["strand"] == "+" else base.translate(_COMPLEMENT)


def _read_codon(reference, record, index):
    first = (index // 3) * 3
    return "".join(
        _orient(reference.base(record["chrom"], _genomic_position(record, i)),
                record)
        for i in range(first, first + 3)
    )


def _protein_effect(index, codon, ref_aa, alt_aa):
    number = index // 3 + 1
    if number == 1 and codon == "ATG":
        return "start_lost", "p.(Met1?)"
    ref = _LONG_AA[ref_aa]
    alt = _LONG_AA[alt_aa]
    if ref_aa == alt_aa == "*":
        return "stop_retained_variant", f"p.(Ter{number}=)"
    if ref_aa == alt_aa:
        return "synonymous_variant", f"p.({ref}{number}=)"
    if alt_aa == "*":
        return "stop_gained", f"p.({ref}{number}Ter)"
    if ref_aa == "*":
        return "stop_lost", ""
    return "missense_variant", f"p.({ref}{number}{alt})"


def _coding_effect(reference, record, variant):
    index = _coding_index(record, variant.pos)
    if index is None:
        return {}
    ref, alt = (_orient(base, record) for base in (variant.ref, variant.alt))
    if not record["complete_cds"]:
        return {"consequence": "coding_sequence_variant"}
    hgvsc = f"{record['transcript']}:c.{index + 1}{ref}>{alt}"
    codon = _read_codon(reference, record, index)
    offset = index % 3
    changed = codon[:offset] + alt + codon[offset + 1:]
    consequence, hgvsp = "coding_sequence_variant", ""
    if codon in _CODE and changed in _CODE and codon[offset] == ref:
        consequence, hgvsp = _protein_effect(
            index, codon, _CODE[codon], _CODE[changed],
        )
    return {"consequence": consequence, "hgvsc": hgvsc, "hgvsp": hgvsp}


def _intronic_sites(record, pos):
    for (_, left_end), (right_start, _) in itertools.pairwise(record["exon"]):
        if left_end < pos < right_start:
            left, right = (left_end, right_start)
            if record["strand"] == "-":
                left, right = right, left
            yield "splice_donor_variant", abs(pos - left), left
            yield "splice_acceptor_variant", -abs(pos - right), right


def _intronic_effect(record, variant, sites):
    close = sorted(sites, key=lambda item: abs(item[1]))
    _, offset, anchor = close[0]
    categories = [t for t, n, _ in close if abs(n) <= 2]
    if abs(offset) > 8:
        return {"consequence": "intron_variant"}
    if not categories:
        categories.append("splice_region_variant")
    effect = {"consequence": "&".join(categories)}
    if len(close) == 1 or abs(close[1][1]) != abs(offset):
        effect["splice_offset"] = offset
        coding = _coding_index(record, anchor)
        if coding is not None and record["complete_cds"]:
            ref = _orient(variant.ref, record)
            alt = _orient(variant.alt, record)
            effect["hgvsc"] = (
                f"{record['transcript']}:c.{coding + 1}{offset:+}{ref}>{alt}"
            )
    return effect


def _exonic_splice_region(record, pos):
    for left, right in itertools.pairwise(record["exon"]):
        left_start, left_end = left
        right_start, right_end = right
        if left_end + 1 >= right_start:
            continue
        near_left = left_start <= pos <= left_end and left_end - pos <= 2
        near_right = right_start <= pos <= right_end and pos - right_start <= 2
        if near_left or near_right:
            return True
    return False


def _transcript_effect(reference, record, variant):
    sites = list(_intronic_sites(record, variant.pos))
    if sites:
        effect = _intronic_effect(record, variant, sites)
    else:
        effect = _coding_effect(reference, record, variant)
        if _exonic_splice_region(record, variant.pos):
            terms = [effect.get("consequence"), "splice_region_variant"]
            effect["consequence"] = "&".join(filter(None, terms))
    if not effect:
        effect = {"consequence": _utr_consequence(record, variant.pos)}
    return effect | {
        "gene": record["gene"], "transcript": record["transcript"],
        "mane": record["mane"], "assembly": "GRCh38",
        "source": "NCBI MANE RefSeq GFF3 + FASTA local",
        "version": record["release"],
    }


def _utr_consequence(record, pos):
    starts, ends, _ = zip(*record["CDS"])
    before = pos < min(starts)
    after = pos > max(ends)
    if before or after:
        return "5_prime_UTR_variant" if (
            before == (record["strand"] == "+")
        ) else "3_prime_UTR_variant"
    return "non_coding_transcript_exon_variant"


class ManeStore:
    """Query transcript effects using local FASTA bases and MANE annotations."""

    def __init__(self, path=None, reference=None):
        if path and reference is None:
            raise ValueError("Un FASTA indexé est requis avec la base MANE.")
        self.connection = readonly_db(path) if path else None
        self.reference = reference

    def _records(self, variant):
        pos = variant.pos
        slot = (pos - 1) // _BIN_SIZE
        rows = self.connection.execute(
            "SELECT payload FROM bin JOIN transcript t ON t.id=bin.transcript"
            " WHERE chrom=? AND slot=? AND t.start<=? AND t.end>=?",
            (variant.chrom, slot, pos, pos),
        )
        return (json.loads(row[0]) for row in rows)

    def lookup(self, variant):
        if self.connection is None:
            return []
        self.reference.check(variant)
        effects = (_transcript_effect(self.reference, record, variant)
                   for record in self._records(variant))
        return preferred_annotations(effects)

    def close(self):
        if self.connection is not None:
            self.connection.close()
