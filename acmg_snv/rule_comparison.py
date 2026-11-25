"""Reason from independently curated comparator substitutions."""

import re

from .model import Variant
from .rule_prediction import MISSENSE, splice_score
from .rule_results import number, result, token


AMINO_ACIDS = dict(zip(
    "ala arg asn asp cys gln glu gly his ile leu lys met phe pro "
    "ser thr trp tyr val".split(), "ARNDCQEGHILKMFPSTWYV",
))
PROTEIN_HGVS = re.compile(r"p\.\(?([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})\)?$")
CLASSES = {"p", "lp"}


def amino_acid(name):
    """Normalize one- and three-letter amino acid names."""
    key = token(name)
    if key in AMINO_ACIDS:
        return AMINO_ACIDS[key]
    value = key.upper()
    return value if value in AMINO_ACIDS.values() else ""


def protein(row):
    """Distinguish residue identity, position and alternative amino acid."""
    reference = row.get("protein_reference")
    position = row.get("protein_position")
    alternate = row.get("protein_alternate")
    if not reference or not position or not alternate:
        match = PROTEIN_HGVS.search(str(row.get("hgvsp") or ""))
        if match is None:
            return None
        reference, raw_position, alternate = match.groups()
        position = int(raw_position)
    ref, alt = amino_acid(reference), amino_acid(alternate)
    if not ref or not alt or ref == alt or type(position) is not int:
        return None
    if position < 1:
        return None
    return ref, position, alt


def variant_id(row):
    """Normalize genomic identity so a variant cannot use itself as evidence."""
    try:
        if token(row.get("variant_id")):
            return Variant.parse(row["variant_id"]).key
        return Variant(row["chrom"], row["pos"], row["ref"], row["alt"]).key
    except (AttributeError, KeyError, TypeError, ValueError):
        return None


def comparable(annotation, row, context):
    """Require a valid context and independent primary evidence."""
    if not isinstance(row, dict):
        return False
    if row.get("independently_curated") is not True:
        return False
    if row.get("protein_splice_impact_excluded") is not True:
        return False
    if token(row.get("classification")) not in CLASSES:
        return False
    if not token(row.get("source")):
        return False
    target = variant_id(annotation)
    comparator = variant_id(row)
    if target is None or comparator is None or target == comparator:
        return False
    for field in ("gene", "disease", "transcript"):
        expected = token(context.get(field) or annotation.get(field))
        if not expected or token(row.get(field)) != expected:
            return False
    return True


def effect_comparison(row, annotation):
    """Require a protein impact at least as large as the comparator's."""
    method = token(row.get("comparison_method"))
    target = number(row.get("target_score"))
    other = number(row.get("comparator_score"))
    if target is None or other is None:
        return False
    if method == "revel":
        score = number(annotation.get("revel"), 0, 1)
        return 0 <= other <= target <= 1 and (score is None or target == score)
    if method == "grantham":
        return 0 < other <= target <= 215
    if method == "blosum62":
        return -4 <= target <= other <= 11
    return False


def from_comparator(row, code, annotation):
    """Describe the comparator and its supporting observation."""
    classification = token(row["classification"])
    identity = variant_id(row)
    evidence_ids = [f"comparator:{identity}"]
    why = f"SNV comparateur {identity}, classe {classification.upper()}."
    strength = "strong" if classification == "p" else "moderate"
    if code == "PM5":
        strength = "moderate" if classification == "p" else "supporting"
        why += f" {row['comparison_method']} {row['target_score']} contre "
        why += f"{row['comparator_score']} pour le comparateur."
        if token(row["comparison_method"]) == "revel":
            transcript = annotation.get("transcript", "")
            evidence_ids.append(f"in_silico:revel:{transcript}")
    return result(why, strength, row["source"], evidence_ids)


def comparator_code(annotation, row, context, observed):
    """Use a comparator only at the same residue and with appropriate effect."""
    if not comparable(annotation, row, context):
        return None
    changed = protein(row)
    if changed is None or changed[:2] != observed[:2]:
        return None
    if changed == observed:
        return "PS1"
    mechanism = token(context.get("mechanism"))
    if mechanism in {"gof", "gain_of_function", "activating"}:
        if row.get("mechanism_concordant") is not True:
            return None
    return "PM5" if effect_comparison(row, annotation) else None


def select_rows(annotation, context, comparators, observed):
    """Keep the best supported comparison for each of PS1 and PM5."""
    rows = {}
    strengths = ("supporting", "moderate", "strong")
    for comparator in comparators:
        code = comparator_code(annotation, comparator, context, observed)
        if code is None:
            continue
        candidate = from_comparator(comparator, code, annotation)
        previous = rows.get(code)
        if previous is None:
            rows[code] = candidate
            continue
        rank = strengths.index(candidate["strength"])
        if rank > strengths.index(previous["strength"]):
            rows[code] = candidate
    return rows


def assess_comparison(annotation, context):
    """Automate the missense branches of PS1 and PM5 from curated evidence."""
    why = "SNV comparateur compatible et preuves primaires à documenter."
    missing = {code: result(why) for code in ("PS1", "PM5")}
    if annotation.get("consequence") != MISSENSE:
        return missing
    score = splice_score(annotation)
    if score is not None and score >= 0.2:
        why = "Effet d'épissage prédit. Revoir la comparaison protéique."
        return {code: result(why) for code in missing}
    if annotation.get("protein_splice_impact_excluded") is not True:
        return missing
    observed = protein(annotation)
    comparators = annotation.get("comparators")
    if observed is None or not isinstance(comparators, list):
        return missing
    return missing | select_rows(annotation, context, comparators, observed)
