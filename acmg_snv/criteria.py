"""The 2026 points, eligibility decisions, and logged curator changes."""

from .model import STRENGTHS


SOURCE = "ACGS 2026, table S4, doi:10.1136/jmg-2026-111835"
CODES = (
    "PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5",
    "PM6", "PP1", "PP2", "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3",
    "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
)
INACTIVE = {
    "PP5": "Une assertion de classification seule ne constitue pas une preuve.",
    "BP6": "Une assertion de classification seule ne constitue pas une preuve.",
    "BP3": "Le critère vise les insertions et délétions, hors périmètre SNV.",
}
MANUAL = {
    "PS3": "Essai fonctionnel protéique spécifique du variant à documenter.",
    "PS4": "Enrichissement en cas ou données cas-témoins à documenter.",
    "PM1": "Domaine critique ou point chaud propre au mécanisme à examiner.",
    "PM4": "Conséquence sur la longueur et alternatives à examiner.",
    "PP2": "Mécanisme faux-sens et contrainte régionale à documenter.",
    "BS2": "Génotype, âge et phénotypage des témoins à documenter.",
    "BS3": "Essai fonctionnel protéique validé à documenter.",
    "BS4": "Non-ségrégation informative à documenter.",
    "BP1": "Mécanisme de la maladie et localisation du variant à examiner.",
    "BP2": "Phase, pénétrance et phénotype de l'autre variant à documenter.",
    "BP5": "Autre diagnostic moléculaire explicatif à documenter.",
}
MAX_GENERIC = {
    "PVS1": "very_strong", "PS1": "strong", "PS2": "very_strong",
    "PS3": "strong", "PS4": "strong", "PM1": "moderate",
    "PM2": "moderate", "PM3": "very_strong", "PM4": "moderate",
    "PM5": "strong", "PM6": "very_strong", "PP1": "strong",
    "PP2": "supporting", "PP3": "supporting", "PP4": "strong",
    "BS1": "strong", "BS2": "strong", "BS3": "strong", "BS4": "strong",
    "BP1": "supporting", "BP2": "supporting", "BP4": "supporting",
    "BP5": "supporting", "BP7": "strong",
}


def base_record(code):
    return {
        "code": code,
        "status": "not_applicable" if code in INACTIVE else "not_evaluable",
        "strength": None,
        "points": None,
        "reason": INACTIVE.get(code) or MANUAL.get(code) or (
            "Éléments nécessaires à l’évaluation non renseignés."
        ),
        "source": SOURCE,
        "evidence_ids": [],
    }


def fill_records(assessment):
    all_rows = {code: base_record(code) for code in CODES}
    for code, record in assessment.items():
        if code not in INACTIVE:
            all_rows[code].update(record)
    return all_rows


def note(row, message, status="not_applied"):
    row["previous_strength"] = row["strength"]
    row["status"] = status
    row["reason"] = f"{row['reason']} {message}"
    row["points"] = None


def nominal_points(code, strength):
    if code == "BA1" and strength == "standalone":
        return None
    weight = STRENGTHS[strength]
    return -weight if code.startswith("B") else weight


def is_generic_strength(code, strength):
    if code == "BA1":
        return strength == "standalone"
    maximum = MAX_GENERIC.get(code)
    return strength in STRENGTHS and bool(
        maximum and STRENGTHS[strength] <= STRENGTHS[maximum]
    )


def apply_overrides(rows, overrides, allowed):
    seen = set()
    for update in overrides or []:
        if not isinstance(update, dict):
            raise ValueError("Modification de critère mal formée.")
        code = update.get("code")
        if code in seen or code not in rows or code in INACTIVE:
            raise ValueError(f"Modification invalide pour {code}.")
        seen.add(code)
        apply_override(rows[code], update, allowed)


def apply_override(row, update, allowed):
    reason = str(update.get("reason") or "").strip()
    if not reason:
        raise ValueError(f"Justification requise pour {row['code']}.")
    previous = dict(row)
    if update.get("action", "apply") == "remove":
        note(row, f"Retrait par le biologiste : {reason}")
    elif update.get("action", "apply") == "apply":
        row.update(override_record(row["code"], update, reason, allowed))
    else:
        raise ValueError("Action manuelle inconnue.")
    row.update(origin="manual", previous=previous)


def override_record(code, update, reason, allowed):
    strength = update.get("strength")
    valid = strength in (*STRENGTHS, "standalone")
    if not valid or not allowed(code, strength):
        raise ValueError(f"Force {strength} indisponible pour {code}.")
    source = str(update.get("source") or "").strip()
    evidence = update.get("evidence_ids") or []
    if not source or not isinstance(evidence, list):
        raise ValueError(f"Source et liste des preuves requises pour {code}.")
    if any(not isinstance(item, str) or not item.strip() for item in evidence):
        raise ValueError(f"Identifiant de preuve invalide pour {code}.")
    return {
        "status": "applied", "strength": strength, "reason": reason,
        "source": source,
        "evidence_ids": evidence or [f"source:{source.casefold()}"],
        "standalone": strength == "standalone",
    }


def score_class(score):
    if score >= 10:
        return "P"
    if score >= 6:
        return "LP"
    if score >= 0:
        return "VUS"
    if score >= -6:
        return "LB"
    return "B"


def signed_rows(rows):
    return [row for row in rows.values() if row["status"] == "applied"]


def assign_points(rows):
    for row in signed_rows(rows):
        row["points"] = nominal_points(row["code"], row["strength"])
