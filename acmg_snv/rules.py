"""Conservative ACGS 2026 SNV assessments from locally supplied evidence."""

from .rule_comparison import assess_comparison
from .rule_lof import assess_lof
from .rule_population import assess_population
from .rule_prediction import assess_prediction
from .rule_results import result


def assess(annotation: dict, context: dict) -> dict[str, dict]:
    """Return generic assessments for the supported criteria."""
    rows = assess_population(annotation, context)
    rows.update(assess_prediction(annotation, context))
    rows.update(assess_lof(annotation, context))
    rows.update(assess_comparison(annotation, context))
    if rows["PVS1"]["status"] == "applied":
        why = "La prédiction PP3 ne se cumule pas avec PVS1."
        rows["PP3"] = result(why, applicable=False)
    stop_revel_double_count(rows)
    return rows


def stop_revel_double_count(rows):
    """Prefer direct variant comparison when PM5 and PP3 reuse REVEL."""
    pp3 = rows["PP3"]
    pm5 = rows["PM5"]
    if pp3["status"] != "applied" or pm5["status"] != "applied":
        return
    if set(pp3["evidence_ids"]).intersection(pm5["evidence_ids"]):
        why = "Le score REVEL est déjà utilisé pour la comparaison PM5."
        rows["PP3"] = result(why, applicable=False)
