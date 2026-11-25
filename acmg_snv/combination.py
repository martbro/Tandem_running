"""Eligibility and combination constraints for the 2026 sequence framework."""

from collections import defaultdict

from .criteria import note, score_class, signed_rows


EXCLUSIVE = (
    ("PVS1", "PM1"), ("PVS1", "PM4"), ("PVS1", "PP2"), ("PVS1", "PP3"),
    ("PS1", "PM4"),
)


def active(rows, *codes):
    return all(rows[code]["status"] == "applied" for code in codes)


def remove_incompatible(rows, alerts, generic):
    for priority, other in EXCLUSIVE if generic else ():
        if active(rows, priority, other):
            note(rows[other], f"Combinaison exclue avec {priority} (S4).")
            alerts.append(f"{priority} et {other} sont incompatibles (S4).")
    if generic and active(rows, "PS2", "PM6"):
        note(rows["PM6"], "Intégrer les observations au critère de novo.")
        alerts.append("PS2 et PM6 doivent être évalués ensemble (S4).")
    if active(rows, "BP4", "BP7") and uses_rna(rows["BP7"]):
        note(rows["BP4"], "La preuve d’épissage est évaluée sous BP7(RNA).")
        alerts.append("BP4 écarté en présence de BP7(RNA).")
    population = [code for code in ("BA1", "BS1", "PM2") if active(rows, code)]
    if len(population) > 1:
        for code in population:
            note(rows[code], "Conflit entre critères de fréquence.")
        alerts.append("Fréquences incompatibles : " + ", ".join(population))
    return alerts


def uses_rna(row):
    text = f"{row.get('source', '')} {row.get('reason', '')}".lower()
    return row.get("rna") is True or "bp7(rna)" in text


def remove_duplicates(rows):
    ids = defaultdict(set)
    for row in signed_rows(rows):
        for evidence in row["evidence_ids"]:
            ids[evidence].add(row["code"])
    groups = {tuple(sorted(group)) for group in ids.values() if len(group) > 1}
    groups = sorted(groups)
    for codes in groups:
        message = "Preuve partagée avec " + ", ".join(codes) + "."
        for code in codes:
            note(rows[code], message)
    return ["Preuve partagée : " + ", ".join(codes) for codes in groups]


def apply_caps(rows, annotation, generic):
    if not generic:
        return []
    alerts = []
    if active(rows, "PVS1", "PS3"):
        lof = annotation.get("lof") or {}
        permitted = isinstance(lof, dict) and (
            lof.get("predicted_nmd") is False
            or lof.get("leaky_splicing") is True
        )
        if permitted:
            cap_pair(rows, "PVS1", "PS3", 8)
        else:
            note(rows["PS3"], "Association PVS1/PS3 à documenter selon S4.")
            alerts.append("Association PVS1/PS3 à revoir selon la table S4.")
    if active(rows, "PP1", "PP4"):
        cap_pair(rows, "PP1", "PP4", 5)
    return alerts


def cap_pair(rows, first, second, cap):
    total = rows[first]["points"] + rows[second]["points"]
    if total <= cap:
        return
    for code, counterpart in ((second, first), (first, second)):
        row = rows[code]
        reduction = min(row["points"], total - cap)
        if reduction:
            row["nominal_points"] = row["points"]
            row["points"] -= reduction
            row["reason"] += f" Plafond de {cap} points avec {counterpart}."
            total -= reduction


def suppress_pm2(rows):
    others = [r for r in signed_rows(rows) if r["code"] != "PM2"]
    if not active(rows, "PM2") or len(others) < 2:
        return
    scores = [r["points"] for r in others if r["points"] is not None]
    if scores and all(points < 0 for points in scores):
        total = sum(scores)
        if score_class(total) == "LB" and total + rows["PM2"]["points"] >= 0:
            note(rows["PM2"], "La rareté seule déplacerait LB vers VUS (S4).")


def total_and_class(rows, context, generic):
    applied = signed_rows(rows)
    points = [r["points"] for r in applied if r["points"] is not None]
    total = sum(points)
    if active(rows, "BA1"):
        return total if points else None, "B", "BA1 autonome"
    result = score_class(total)
    strong_bs1 = active(rows, "BS1") and rows["BS1"]["strength"] == "strong"
    if result == "LB" and strong_bs1:
        standalone = rows["BS1"].get("standalone") is True or (
            generic and context["inheritance"] == "AD"
            and context.get("high_penetrance") is True
        )
        if standalone:
            return total, result, "BS1 autonome"
    agreeing = sum((p > 0 if total > 0 else p < 0) for p in points)
    if result != "VUS" and agreeing < 2:
        return total, "VUS", "Deux critères corroborants requis (ACGS R6)."
    return total, result, "Seuils SNV ACGS 2026 (figure 2)."


def opposing_evidence(rows):
    positive = any(r["code"].startswith("P") for r in signed_rows(rows))
    negative = any(r["code"].startswith("B") for r in signed_rows(rows))
    return positive and negative
