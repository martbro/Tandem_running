"""Évaluation prudente des preuves cliniques renseignées par le biologiste."""

from fractions import Fraction
from math import prod

from .model import Variant


S4 = "ACGS 2026, supplément S4, doi:10.1136/jmg-2026-111835"
S2 = f"{S4} et tableau S2"
S3 = f"{S4} et tableau S3, ClinGen SVI de novo v1.1"
PM3 = f"{S4}, ClinGen SVI PM3 v1.0"
PP1 = f"{S4}, Jarvik et Browning 2016, tableau 1"
LEVELS = ("supporting", "moderate", "strong", "very_strong")
DENOVO_POINTS = {"specific": 8, "consistent": 4, "heterogeneous": 2}
RECESSIVE = {"AR", "XLR"}


def row(reason, source, *, strength=None, ids=(), inapplicable=False):
    status = "not_applicable" if inapplicable else "not_evaluable"
    if strength:
        status = "applied"
    return {
        "status": status,
        "strength": strength,
        "reason": reason,
        "source": source,
        "evidence_ids": sorted(set(ids)),
    }


def quarter_strength(points):
    for threshold, level in zip((16, 8, 4, 2), reversed(LEVELS)):
        if points >= threshold:
            return level
    return None


def text_field(data, name, default=""):
    value = data.get(name, default)
    return value.strip() if isinstance(value, str) else ""


def has_checks(data, *fields):
    return all(data.get(field) is True for field in fields)


def countable_cases(value):
    if not isinstance(value, list):
        return []
    return [part for part in value if isinstance(part, dict)]


def case_details(cases):
    details = []
    for case in cases:
        reason = text_field(case, "reason")
        source = text_field(case, "source")
        details.append(f"{reason} ({source})" if source else reason)
    return " | ".join(details)


def denovo_case(case):
    status = text_field(case, "parental_relationships")
    specification = text_field(case, "phenotype_fit")
    case_id = text_field(case, "case_id", "proband")
    if status not in {"confirmed", "assumed"} or not case_id:
        return None
    if not has_checks(case, "parents_tested", "denovo_established"):
        return None
    if not text_field(case, "reason"):
        return None
    if case.get("germline_mosaicism") is True and status != "confirmed":
        return None
    points = DENOVO_POINTS.get(specification, 0)
    if status == "assumed":
        points //= 2
    if not points:
        return None
    return (case_id, points, status, specification, case)


def denovo_strength(points, cases, inheritance):
    strength = quarter_strength(points)
    if inheritance != "AR" or strength is None:
        return strength, ""
    other_alleles = {case.get("additional_plp") is True for case in cases}
    if other_alleles == {True}:
        return strength, ""
    if len(other_alleles) > 1:
        return None, " Cas AR hétérogènes, à évaluer individuellement."
    rank = LEVELS.index(strength)
    if rank == 0:
        return None, " Seuil insuffisant après la réduction AR d’un niveau."
    return LEVELS[rank - 1], " Force réduite d’un niveau en AR sans autre P/LP."


def denovo_rows(context):
    cases = countable_cases(context.get("denovo"))
    valid = [item for case in cases if (item := denovo_case(case))]
    missing = "Résultats parentaux et concordance des observations requis."
    if not valid:
        return {code: row(missing, S3) for code in ("PS2", "PM6")}
    ids = [item[0] for item in valid]
    if len(set(ids)) != len(ids) or not unrelated(context, valid):
        reason = "Vérifier les doublons et l’indépendance des cas de novo."
        return {code: row(reason, S3) for code in ("PS2", "PM6")}
    return score_denovo(valid, context)


def unrelated(context, cases):
    return len(cases) == 1 or context.get("denovo_unrelated") is True


def score_denovo(valid, context):
    mixed = sum(item[1] for item in valid if item[3] == "heterogeneous")
    points = min(mixed, 4)
    points += sum(item[1] for item in valid if item[3] != "heterogeneous")
    strength, note = denovo_strength(
        points, [item[4] for item in valid], context.get("inheritance")
    )
    confirmed = any(item[2] == "confirmed" for item in valid)
    code, other = ("PS2", "PM6") if confirmed else ("PM6", "PS2")
    summary = f"{points / 4:g} point(s) ClinGen de novo, {len(valid)} cas."
    details = case_details(item[4] for item in valid)
    if not strength:
        details = f"Force non atteinte. {details}"
    ids = [f"phenotype:{item[0]}" for item in valid]
    ids += [f"denovo:{item[0]}" for item in valid]
    return {
        code: row(f"{summary}{note} {details}", S3, strength=strength, ids=ids),
        other: row(f"Observations évaluées sous {code}.", S3,
                   inapplicable=True),
    }


def partner_key(text):
    try:
        return Variant.parse(text).key.upper()
    except ValueError:
        return text.upper()


def candidate_phase(case, variant_key):
    relationship = text_field(case, "relationship").lower()
    if relationship not in {"homozygous", "trans", "unknown"}:
        return None
    if not has_checks(case, "affected", "rare"):
        return None
    if not text_field(case, "reason"):
        return None
    if relationship == "homozygous" and has_checks(case, "phase_confirmed"):
        return (2, "homozygous")
    if not has_checks(case, "partner_rare", "independent_classification",
                      "same_condition"):
        return None
    partner = partner_key(text_field(case, "variant"))
    if not partner or partner == variant_key:
        return None
    classification = text_field(case, "classification").upper()
    if classification not in {"P", "LP"}:
        return None
    if relationship == "trans" and has_checks(case, "phase_confirmed"):
        return (4, "trans")
    if relationship == "unknown" and case.get("predominantly_cis") is False:
        return (2 if classification == "P" else 1, "unknown")
    return None


def phase_cases(context, annotation):
    observations = {}
    variant_key = partner_key(text_field(annotation, "key"))
    for case in countable_cases(context.get("phase")):
        weight = candidate_phase(case, variant_key)
        identifier = text_field(case, "case_id", "proband")
        group = text_field(case, "family_id") or identifier
        if not weight or not group:
            continue
        if group not in observations or observations[group][0] < weight[0]:
            observations[group] = (*weight, case)
    return observations


def phase_points(observations):
    unknown, homozygous, definite = {}, 0, 0
    for points, relation, case in observations.values():
        if relation == "unknown":
            partner = partner_key(text_field(case, "variant"))
            unknown[partner] = max(unknown.get(partner, 0), points)
        elif relation == "homozygous":
            homozygous += points
        else:
            definite += points
    return definite + min(homozygous, 4) + sum(unknown.values())


def pm3_row(context, annotation):
    inheritance = context.get("inheritance")
    if inheritance not in RECESSIVE:
        return row("Critère de transmission récessive.", PM3, inapplicable=True)
    if inheritance == "XLR":
        return row("Vérifier manuellement les cas récessifs liés à l’X.", PM3)
    observations = phase_cases(context, annotation)
    if not observations:
        reason = "Vérifier phase, rareté et validité du partenaire pour PM3."
        return row(reason, PM3)
    if len(observations) > 1 and context.get("phase_unrelated") is not True:
        return row("Indépendance des cas PM3 à confirmer.", PM3)
    points = phase_points(observations)
    strength = quarter_strength(points)
    summary = f"{points / 4:g} point(s) de phase ClinGen."
    if not strength:
        summary += " Seuil minimal de 0,5 non atteint."
    details = case_details(case for _, _, case in observations.values())
    ids = [f"biallelic:{identifier}" for identifier in observations]
    return row(f"{summary} {details}", PM3, strength=strength, ids=ids,
               inapplicable=strength is None)


def pp4_strength(evidence):
    finding = text_field(evidence, "specificity")
    if finding == "specific":
        return "supporting"
    if not has_checks(evidence, "mdt_confirmed"):
        return None
    return {"additional_features": "moderate", "pathognomonic": "strong"}.get(
        finding
    )


def pp4_row(context, denovo):
    evidence = context.get("phenotype_evidence")
    if not isinstance(evidence, dict):
        return row("Spécificité et étendue du test à valider.", S2)
    reason = text_field(evidence, "reason")
    strength = pp4_strength(evidence)
    case_id = f"phenotype:{text_field(evidence, 'case_id', 'proband')}"
    if not reason or not has_checks(evidence, "majority_causes_tested"):
        return row("Justifier la spécificité et la couverture du test.", S2)
    if not strength:
        return row("Confirmer la spécificité et la discussion en RCP.", S2)
    if any(case_id in item["evidence_ids"] for item in denovo.values()
           if item["status"] == "applied"):
        return row("Spécificité déjà comptée dans PS2/PM6.", S2,
                   inapplicable=True)
    return row(reason, S2, strength=strength, ids=[case_id])


def family_probability(family, context):
    value = family.get("probability")
    if value is not None and not isinstance(value, bool):
        try:
            probability = Fraction(str(value))
        except (ValueError, ZeroDivisionError):
            return None
        return probability if 0 < probability < 1 else None
    meioses = family.get("informative_meioses")
    simple = has_checks(family, "jarvik_simple_assumptions")
    if simple and context.get("inheritance") == "AD":
        if type(meioses) is int and 0 < meioses <= 100:
            return Fraction(1, 2 ** meioses)
    return None


def family_strength(probability, count):
    thresholds = (32, 16, 8) if count == 1 else (16, 8, 4)
    for denominator, strength in zip(thresholds, reversed(LEVELS[:3])):
        if probability <= Fraction(1, denominator):
            return strength
    return None


def pp1_row(context):
    evidence = context.get("segregation_evidence")
    if not isinstance(evidence, dict):
        return row("Méioses informatives ou probabilités documentées requises.",
                   PP1)
    if not has_checks(evidence, "definitive_gene_disease", "cosegregation"):
        return row("Valider le lien gène-maladie et la coségrégation.", PP1)
    families = countable_cases(evidence.get("families"))
    ids = [text_field(family, "family_id") for family in families]
    if not ids or not all(ids) or len(set(ids)) != len(ids):
        return row("Identifier chaque famille sans doublon.", PP1)
    if len(ids) > 1 and evidence.get("independent_families") is not True:
        return row("Vérifier que les familles sont indépendantes.", PP1)
    return score_pp1(families, ids, evidence, context)


def score_pp1(families, ids, evidence, context):
    values = [family_probability(family, context) for family in families]
    reason = text_field(evidence, "reason")
    if not reason or any(value is None for value in values):
        return row("Documenter les méioses et le calcul par famille.", PP1)
    probability = prod(values)
    strength = family_strength(probability, len(ids))
    note = f"Probabilité sous l’hypothèse nulle {probability}."
    note += f" Nombre de familles {len(ids)}."
    if not strength:
        note += " Seuil de coségrégation non atteint."
    ids = [f"segregation:{identifier}" for identifier in ids]
    return row(f"{note} {reason}", PP1, strength=strength, ids=ids,
               inapplicable=strength is None)


def assess_family(annotation, context):
    """Propose les critères à partir de faits cliniques structurés."""
    denovo = denovo_rows(context)
    return {
        **denovo,
        "PM3": pm3_row(context, annotation),
        "PP1": pp1_row(context),
        "PP4": pp4_row(context, denovo),
    }
