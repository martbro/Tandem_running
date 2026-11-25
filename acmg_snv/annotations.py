"""Normalized, local variant annotations."""

from contextlib import contextmanager
import csv
import gzip
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import tempfile

from .ann_reference import Reference, build_fasta_index, primary_chromosome
from .ann_vcf import csq_columns, parse_vcf_annotation
from .model import Variant


_COORDS = ("chrom", "pos", "ref", "alt")
_TEXT = ("gene", "transcript", "mane", "consequence", "source", "version")
_FLOATS = {
    "revel", "bayesdel", "alphamissense", "spliceai", "phylop",
    "ba1_threshold", "bs1_threshold", "population_af", "popmax_faf",
}
_INTS = {
    "population_ac", "population_an", "popmax_an",
    "protein_position", "splice_offset",
}
_BOOLEAN = {
    "population_pass", "population_eligible", "population_absence_verified",
    "bp7_rna_used", "protein_splice_impact_excluded",
}
_JSON = {"lof", "comparators", "functional_ps3"}


def text_input(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, encoding="utf-8-sig", newline="")


def readonly_db(path):
    uri = Path(path).resolve().as_uri()
    return sqlite3.connect(f"{uri}?mode=ro", uri=True, check_same_thread=False)


def _reject_constant(value):
    raise ValueError(f"Constante JSON non prise en charge ({value}).")


def _finite_float(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Nombre non fini ({value}).")
    return number


def load_json(text):
    return json.loads(
        text, parse_constant=_reject_constant, parse_float=_finite_float,
    )


def _coordinate_match(row, variant):
    expected = (variant.chrom, str(variant.pos), variant.ref, variant.alt)
    provided = [row.get(key) for key in _COORDS]
    for index, (actual, target) in enumerate(zip(provided, expected)):
        if actual in (None, ""):
            continue
        normalized = str(actual).strip().upper()
        if index == 0:
            normalized = primary_chromosome(normalized)
        if normalized != target:
            raise ValueError(
                f"Coordonnée {_COORDS[index]} discordante pour {variant.key}."
            )


def _mane_status(value):
    status = value.lower().replace("_", " ").replace("-", " ").strip()
    if not status:
        return ""
    if status in {"select", "mane select"}:
        return "MANE Select"
    if status in {"plus clinical", "mane plus clinical"}:
        return "MANE Plus Clinical"
    raise ValueError(f"Statut MANE inconnu ({value}).")


def normalize_annotation(row, variant):
    """Validate identity and normalize an annotation's core fields."""
    if not isinstance(row, dict):
        raise ValueError("Chaque annotation doit être un objet JSON.")
    _coordinate_match(row, variant)
    build = str(row.get("assembly", "GRCh38"))
    if not re.fullmatch(r"GRCh38(?:\.p\d+)?|hg38", build, flags=re.I):
        raise ValueError("L'annotation doit être repérée en GRCh38.")
    annotation = dict(row)
    for field in _TEXT:
        value = annotation.pop(field, None)
        if value is None:
            continue
        if not isinstance(value, str):
            raise ValueError(f"Le champ {field} doit être du texte.")
        if value.strip():
            annotation[field] = value.strip()
    if "gene" in annotation:
        annotation["gene"] = annotation["gene"].upper()
    if "mane" in annotation:
        annotation["mane"] = _mane_status(annotation["mane"])
    annotation["assembly"] = "GRCh38"
    annotation.update(zip(_COORDS, (variant.chrom, variant.pos,
                                    variant.ref, variant.alt)))
    return annotation


def preferred_annotations(rows):
    """Present MANE options first while retaining all transcript records."""
    return sorted(rows, key=_annotation_order)


def _annotation_order(row):
    ranks = {"MANE Select": 0, "MANE Plus Clinical": 1}
    return (ranks.get(row.get("mane"), 2), row.get("gene", ""),
            row.get("transcript", ""))


def _paste_vcf(text, variant):
    columns = ()
    results = []
    for line in text.splitlines():
        if line.startswith("##INFO=<ID=CSQ,"):
            columns = csq_columns(line)
        elif line and not line.startswith("#"):
            if len(line.split("\t")) < 8:
                raise ValueError("Une ligne VCF requiert au moins 8 champs.")
            if not columns:
                raise ValueError("Collez également l'en-tête INFO/CSQ de VEP.")
            results.extend(parse_vcf_annotation(variant, line, columns))
    if not results:
        raise ValueError(f"Aucune annotation VCF pour {variant.key}.")
    return results


def parse_annotation(text, variant):
    """Parse pasted JSON, normalized TSV or VEP CSQ and VCF records."""
    text = text.strip().lstrip("\ufeff")
    if not text:
        return []
    if text.startswith(("{", "[")):
        payload = load_json(text)
        rows = payload if isinstance(payload, list) else [payload]
    elif set(_COORDS) <= set(text.splitlines()[0].split("\t")):
        tsv = csv.DictReader(text.splitlines(), delimiter="\t")
        rows = map(_tsv_payload, tsv)
    else:
        rows = _paste_vcf(text, variant)
    return preferred_annotations(
        normalize_annotation(row, variant) for row in rows
    )


@contextmanager
def new_db(output_path):
    """Build a SQLite database atomically; retain existing outputs."""
    output = Path(output_path)
    if output.exists():
        raise FileExistsError(f"Le fichier existe déjà ({output}).")
    with tempfile.NamedTemporaryFile(dir=output.parent, delete=False) as handle:
        scratch = Path(handle.name)
    try:
        connection = sqlite3.connect(scratch)
        try:
            with connection:
                yield connection
        finally:
            connection.close()
        if output.exists():
            raise FileExistsError(f"Le fichier existe déjà ({output}).")
        os.replace(scratch, output)
    finally:
        scratch.unlink(missing_ok=True)


class AnnotationStore:
    """Lookup local annotations by canonical genome position and alleles."""

    def __init__(self, path=None):
        self.connection = readonly_db(path) if path else None

    def lookup(self, variant):
        if self.connection is None:
            return []
        rows = self.connection.execute(
            "SELECT payload FROM annotation WHERE variant_key = ?",
            (variant.key,),
        )
        return preferred_annotations(load_json(row[0]) for row in rows)

    def close(self):
        if self.connection is not None:
            self.connection.close()


def _tsv_value(field, value):
    if field in _INTS:
        return int(value)
    if field in _FLOATS:
        return _finite_float(value)
    if field in _JSON:
        return load_json(value)
    if field in _BOOLEAN or field.endswith("_applicable"):
        parsed = load_json(value.lower())
        if not isinstance(parsed, bool):
            raise ValueError(f"{field} doit valoir true ou false.")
        return parsed
    return value


def _tsv_payload(row):
    if None in row:
        raise ValueError("Nombre de colonnes TSV supérieur à l'en-tête.")
    nested = row.pop("annotation_json", "")
    extras = load_json(nested) if nested else {}
    if not isinstance(extras, dict):
        raise ValueError("annotation_json doit contenir un objet JSON.")
    return extras | {
        key: _tsv_value(key, val)
        for key, val in row.items() if val not in ("", None, ".")
    }


def _source_rows(path):
    with text_input(path) as source:
        extension = str(path).removesuffix(".gz").lower()
        if extension.endswith(".tsv"):
            yield from map(_tsv_payload, csv.DictReader(source, delimiter="\t"))
        elif extension.endswith(".json"):
            data = load_json(source.read())
            yield from data if isinstance(data, list) else [data]
        elif extension.endswith((".jsonl", ".ndjson")):
            yield from (load_json(line) for line in source if line.strip())
        else:
            raise ValueError("Utilisez un fichier TSV, JSON ou JSONL.")


def _index_rows(path):
    for row_number, row in enumerate(_source_rows(path), 1):
        try:
            chrom, pos, ref, alt = (row[field] for field in _COORDS)
            label = primary_chromosome(str(chrom))
            variant = Variant.parse(f"chr{label}:{pos}:{ref}:{alt}")
            annotation = normalize_annotation(row, variant)
            annotation.setdefault("source", Path(path).name)
            payload = json.dumps(
                annotation, ensure_ascii=False, allow_nan=False,
            )
        except (KeyError, ValueError, TypeError) as exc:
            error = f"Annotation {row_number} invalide : {exc}"
            raise ValueError(error) from exc
        yield variant.key, payload


def build_annotation_db(input_path, output_path):
    """Create a bulk-import annotation index from normalized data."""
    with new_db(output_path) as database:
        database.execute(
            "CREATE TABLE annotation (variant_key TEXT, payload TEXT,"
            " PRIMARY KEY(variant_key, payload)) WITHOUT ROWID"
        )
        database.executemany(
            "INSERT OR IGNORE INTO annotation VALUES (?, ?)",
            _index_rows(input_path),
        )
