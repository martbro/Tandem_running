"""Shared classification workflow for the browser and the batch command."""

from . import __version__
from .ann_mane import ManeStore
from .annotations import AnnotationStore, Reference
from .combination import (
    apply_caps, opposing_evidence, remove_duplicates, remove_incompatible,
    suppress_pm2, total_and_class,
)
from .criteria import (
    INACTIVE, apply_overrides, assign_points, fill_records, is_generic_strength,
    note, signed_rows,
)
from .family import assess_family
from .functional import functional_evidence, load_ps3
from .rules import assess
from .selection import (
    brief_annotation, choose_annotation, combine_annotations, prepare_context,
    select_consequence,
)
from .vcep import Catalog, apply_profile, vcep_level


class Classifier:
    """Resources are opened locally and reused across variants."""

    def __init__(
        self, annotation_db=None, reference=None, vcep_dir=None, ps3=None,
        mane_db=None,
    ):
        self.reference = Reference(reference) if reference else None
        self.store = AnnotationStore(annotation_db)
        self.mane = ManeStore(mane_db, self.reference)
        self.catalog = Catalog(vcep_dir)
        self.ps3 = load_ps3(ps3)

    def close(self):
        self.mane.close()
        self.store.close()
        if self.reference:
            self.reference.close()

    def classify(
        self, variant, context, annotation=None, overrides=None, ps3=None,
    ):
        context = prepare_context(context)
        report = base_report(variant, context)
        if not context["inheritance"]:
            return stop(report, "needs_context", "Renseignez la transmission.")
        if self.reference:
            self.reference.check(variant)
            report["reference_checked"] = True
        rows = self.annotations(variant, annotation)
        if not rows:
            return stop(report, "not_evaluable", "Annotation GRCh38 absente.")
        report["candidates"] = [brief_annotation(row) for row in rows]
        selected, message = choose_annotation(rows, context)
        if message:
            return stop(report, "needs_context", message)
        selected.setdefault("key", variant.key)
        selected.setdefault("variant_id", variant.key)
        if not apply_annotation(report, selected):
            return report
        internal = (ps3 or {}).get(variant.key) or self.ps3.get(variant.key)
        return self.evaluate(report, selected, overrides, internal)

    def annotations(self, variant, annotation):
        database = self.store.lookup(variant)
        mane = self.mane.lookup(variant)
        pasted = [annotation] if isinstance(annotation, dict) else annotation
        pasted = pasted or []
        return combine_annotations(combine_annotations(mane, database), pasted)

    def evaluate(self, report, annotation, overrides, internal):
        context = report["context"]
        resolution = self.catalog.resolve(
            context.get("gene"), context.get("disease"), context["inheritance"],
        )
        report["vcep"] = resolution
        report["alerts"].extend(resolution["alerts"])
        if self.vcep_exclusion_needed(report, resolution):
            return report
        generic = not resolution["matches"]
        rows = self.criteria(annotation, context, resolution)
        allowed = strength_check(resolution, context)
        functional = functional_evidence(annotation, internal)
        if functional:
            add_functional(rows, functional, allowed)
        apply_overrides(rows, overrides, allowed)
        hard_alerts = combination_alerts(rows, annotation, generic)
        report["alerts"].extend(hard_alerts)
        report["criteria"] = list(rows.values())
        return finish(report, rows, hard_alerts, generic)

    def vcep_exclusion_needed(self, report, resolution):
        context = report["context"]
        if resolution["matches"] or not self.catalog.for_gene(context["gene"]):
            return False
        gene = context["gene"]
        resolution["other_scopes"] = self.catalog.for_gene(gene)
        disease = context["disease"]
        if self.catalog.for_disease(gene, disease):
            stop(report, "needs_context",
                 "Transmission incompatible avec le VCEP de cette maladie.")
            return True
        if str(context.get("vcep_exclusion_reason") or "").strip():
            return False
        stop(report, "review_required",
             "Un VCEP existe pour un autre contexte. Vérifiez son périmètre "
             "ou documentez la raison de son exclusion.")
        return True

    def criteria(self, annotation, context, resolution):
        if not resolution["matches"]:
            return fill_records(assess(annotation, context) | assess_family(
                annotation, context,
            ))
        rows = fill_records({})
        for row in rows.values():
            if row["code"] not in INACTIVE:
                row["reason"] = "Évaluer selon la spécification VCEP affichée."
                row["source"] = "ClinGen CSpec"
        profile = resolution.get("profile")
        if profile:
            for code, assessment in apply_profile(
                profile, annotation, context, rows,
            ).items():
                rows[code].update(assessment)
        return rows


def base_report(variant, context):
    return {
        "schema_version": 1, "software_version": __version__,
        "framework": "ACGS 2026", "variant": variant.as_dict(),
        "context": context,
        "annotation": None, "status": "not_evaluable", "classification": None,
        "score": None, "alerts": [], "criteria": [], "reference_checked": False,
    }


def stop(report, status, message):
    report["status"] = status
    report["alerts"].append(message)
    return report


def apply_annotation(report, selected):
    report["annotation"] = selected
    context = report["context"]
    context["disease"] = context.get("disease") or selected.get("disease", "")
    context["gene"] = context.get("gene") or selected.get("gene", "")
    transcript = selected.get("transcript") or context.get("transcript")
    context["transcript"] = transcript
    consequence = select_consequence(selected)
    if not selected.get("consequence"):
        stop(report, "not_evaluable", "Conséquence moléculaire non renseignée.")
    elif not consequence:
        stop(report, "out_of_scope", "Hors périmètre codant ou épissage.")
    elif not context["gene"] or not context["transcript"]:
        stop(report, "needs_context", "Gène et transcrit clinique requis.")
    else:
        selected["all_consequences"] = selected["consequence"]
        selected["consequence"] = consequence
        return True
    return False


def strength_check(resolution, context):
    reviewed = context.get("vcep_specification")

    def allowed(code, strength):
        if strength == "standalone" and code != "BA1":
            return False
        level = vcep_level(resolution, code, strength, reviewed)
        if resolution["matches"]:
            return level is True
        return is_generic_strength(code, strength)

    return allowed


def add_functional(rows, functional, allowed):
    row = rows["PS3"]
    row.update(functional)
    if not allowed("PS3", functional["strength"]):
        note(row, "Force du fichier hors référentiel applicable.")


def combination_alerts(rows, annotation, generic):
    alerts = remove_incompatible(rows, [], generic)
    alerts.extend(remove_duplicates(rows))
    assign_points(rows)
    alerts.extend(apply_caps(rows, annotation, generic))
    if generic:
        suppress_pm2(rows)
    return alerts


def vcep_confirmed(report):
    resolution = report["vcep"]
    if not resolution["requires_review"]:
        return True
    context = report["context"]
    if context.get("vcep_reviewed") is not True or not context.get("disease"):
        return False
    selected = context.get("vcep_specification")
    matches = resolution["matches"]
    exact = {f"{row['id']}@{row['version']}" for row in matches if row["exact"]}
    if selected:
        return selected in exact
    return len(exact) == 1


def finish(report, rows, hard_alerts, generic):
    context = report["context"]
    score, classification, reason = total_and_class(rows, context, generic)
    report.update(score=score, score_reason=reason)
    has_conflict = opposing_evidence(rows)
    if has_conflict:
        message = "Preuves pathogènes et bénignes : revue experte requise."
        report["alerts"].append(message)
    if not vcep_confirmed(report):
        return stop(report, "review_required", "Spécification VCEP à revoir.")
    if not signed_rows(rows) and not hard_alerts:
        report["score"] = None
        return stop(report, "not_evaluable", "Aucun critère coté à ce stade.")
    return finish_review(
        report, rows, hard_alerts, classification, has_conflict,
    )


def finish_review(report, rows, hard_alerts, classification, has_conflict):
    context = report["context"]
    resolved = context.get("reviewed") is True and bool(
        str(context.get("conflict_resolution") or "").strip(),
    )
    if hard_alerts or has_conflict and (not resolved or classification == "B"
                                        and rows["BA1"]["status"] == "applied"):
        report["status"] = "review_required"
    else:
        report["classification"] = classification
        approved = context.get("reviewed") is True and context.get("disease")
        report["status"] = "approved" if approved else "provisional"
        if context.get("reviewed") is True and not approved:
            report["alerts"].append("Renseignez la maladie pour valider.")
    return report
