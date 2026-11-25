"""Extract allele-specific VEP annotations using the VCF-declared format."""

import math
import re
from urllib.parse import unquote

from .ann_reference import primary_chromosome


def csq_columns(header_line):
    """Read ordered VEP fields from the INFO/CSQ header declaration."""
    match = re.search(r'\bFormat:\s*([^"]+)"', header_line)
    if not header_line.startswith("##INFO=<ID=CSQ,") or not match:
        raise ValueError("En-tête VCF INFO/CSQ invalide.")
    columns = tuple(name.strip() for name in match[1].strip().split("|"))
    if not {"Allele", "Consequence"} <= set(columns):
        raise ValueError("Champs Allele et Consequence absents du CSQ.")
    if any(not name for name in columns) or len(set(columns)) != len(columns):
        raise ValueError("Colonnes CSQ invalides ou répétées.")
    return columns


def _allele_index(variant, fields):
    chrom = primary_chromosome(fields[0])
    if chrom != variant.chrom or int(fields[1]) != variant.pos:
        return None
    alts = fields[4].upper().split(",")
    if fields[3].upper() != variant.ref or variant.alt not in alts:
        return None
    if alts.count(variant.alt) != 1:
        raise ValueError("Allèle alternatif VCF dupliqué.")
    return alts.index(variant.alt) + 1


def _csq_value(info):
    for field in info.split(";"):
        if field.startswith("CSQ="):
            return field.removeprefix("CSQ=")
    return ""


def _csq_match(row, variant, allele_index):
    number = row.get("ALLELE_NUM")
    if number:
        return int(number) == allele_index
    return row["Allele"].upper() == variant.alt


def _probability(value, label):
    if not value or value == ".":
        return None
    number = float(value)
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise ValueError(f"Score VEP {label} hors de [0, 1].")
    return number


def _score(row, names):
    for name in names:
        if row.get(name):
            value = _probability(row[name], name)
            if value is not None:
                return value
    return None


def _add_predictions(annotation, row):
    scores = {
        "revel": _score(row, ("REVEL", "REVEL_score")),
        "alphamissense": _score(row, ("am_pathogenicity", "AlphaMissense")),
    }
    delta_scores = []
    symbol = row.get("SpliceAI_pred_SYMBOL")
    if not symbol or symbol.upper() == annotation.get("gene"):
        for kind in ("AG", "AL", "DG", "DL"):
            names = (f"SpliceAI_pred_DS_{kind}", f"SpliceAI_DS_{kind}")
            score = _score(row, names)
            if score is not None:
                delta_scores.append(score)
    if delta_scores:
        scores["spliceai"] = max(delta_scores)
    for predictor, score in scores.items():
        if score is not None:
            annotation[predictor] = score
            annotation[f"{predictor}_source"] = "VEP INFO/CSQ"


def _mane(row):
    if row.get("MANE_SELECT"):
        return "MANE Select"
    if row.get("MANE_PLUS_CLINICAL"):
        return "MANE Plus Clinical"
    return ""


def _to_annotation(row):
    if row.get("MANE_SELECT") or row.get("MANE_PLUS_CLINICAL"):
        transcript = row.get("MANE_SELECT") or row.get("MANE_PLUS_CLINICAL")
    else:
        transcript = row.get("Feature", "")
    annotation = {
        "gene": row.get("SYMBOL", "").upper(),
        "transcript": transcript,
        "mane": _mane(row),
        "consequence": row.get("Consequence", ""),
        "hgvsc": row.get("HGVSc", ""),
        "hgvsp": row.get("HGVSp", ""),
        "source": "VEP INFO/CSQ",
        "vep": {name: value for name, value in row.items() if value},
    }
    annotation = {key: value for key, value in annotation.items() if value}
    _add_predictions(annotation, row)
    return annotation


def parse_vcf_annotation(variant, record_line, csq_fields):
    """Read CSQ annotations for one REF/ALT pair in a VCF record."""
    if not csq_fields:
        return []
    fields = record_line.rstrip("\r\n").split("\t")
    if len(fields) < 8:
        raise ValueError("Une ligne VCF requiert au moins 8 champs.")
    allele_index = _allele_index(variant, fields)
    if allele_index is None:
        return []
    result = []
    for value in filter(None, _csq_value(fields[7]).split(",")):
        values = [unquote(item) for item in value.split("|")]
        if len(values) > len(csq_fields):
            raise ValueError("Nombre de valeurs CSQ supérieur à l'en-tête.")
        row = {name: "" for name in csq_fields}
        row.update(zip(csq_fields, values))
        if _csq_match(row, variant, allele_index):
            result.append(_to_annotation(row))
    return result
