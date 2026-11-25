"""Read an allele-specific, laboratory-curated PS3 table."""

import csv
import io
from pathlib import Path

from .model import STRENGTHS, Variant


def load_ps3(path=None, *, text=None):
    """Read a UTF-8 TSV/CSV table of laboratory evidence."""
    if path is not None and text is not None:
        raise ValueError("Indiquez une seule source pour le fichier PS3.")
    if path is not None:
        text = Path(path).read_text(encoding="utf-8-sig")
    if not text or not text.strip():
        return {}
    content = text.lstrip("\ufeff")
    delimiter = "\t" if "\t" in content.splitlines()[0] else ","
    reader = csv.DictReader(io.StringIO(content), delimiter=delimiter)
    if not {"variant", "strength", "reason"} <= set(reader.fieldnames or []):
        raise ValueError("PS3 exige les colonnes variant, strength, reason.")
    source = Path(path).name if path is not None else "fichier PS3 importé"
    result = {}
    for raw in reader:
        key, entry = ps3_entry(raw, reader.line_num, source)
        if key in result:
            raise ValueError(f"Plusieurs lignes PS3 pour {key}.")
        result[key] = entry
    return result


def ps3_entry(raw, line, source):
    if None in raw or any(value is None for value in raw.values()):
        raise ValueError(f"Ligne PS3 mal formée : {line}.")
    variant = Variant.parse(raw["variant"])
    strength = raw["strength"].strip().lower()
    reason = raw["reason"].strip()
    if strength not in STRENGTHS or not reason:
        raise ValueError(f"Force ou justification PS3 invalide, ligne {line}.")
    return variant.key, {
        "status": "applied", "strength": strength, "reason": reason,
        "source": raw.get("source", "").strip() or source,
        "evidence_ids": [raw.get("evidence_id", "").strip()
                         or f"laboratoire:PS3:{variant.key}"],
        "origin": "laboratory_ps3",
    }


def functional_evidence(annotation, laboratory):
    if laboratory:
        return dict(laboratory)
    record = annotation.get("functional_ps3")
    if not isinstance(record, dict) or record.get("validated") is not True:
        return None
    strength = record.get("strength")
    required = ("reason", "source", "evidence_id")
    complete = all(record.get(key) for key in required)
    if strength not in STRENGTHS or not complete:
        return None
    return {
        "status": "applied", "strength": strength,
        "reason": record["reason"], "source": record["source"],
        "evidence_ids": [str(record["evidence_id"])], "origin": "curated",
    }
