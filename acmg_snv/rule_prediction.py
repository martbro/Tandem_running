"""Assess computational protein and splice evidence."""

from .rule_results import number, result, token


MISSENSE = "missense_variant"
CANONICAL = {"splice_acceptor_variant", "splice_donor_variant"}
UTR = {"5_prime_UTR_variant", "3_prime_UTR_variant"}
INTRONIC = {
    "intron_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_donor_5th_base_variant",
    "splice_polypyrimidine_tract_variant",
}
SPLICE_TESTABLE = INTRONIC | UTR | {
    MISSENSE,
    "synonymous_variant",
    "non_coding_transcript_exon_variant",
}
MISSENSE_LIMITS = {"revel": 0.644, "bayesdel": 0.130, "alphamissense": 0.792}
DISPLAY_NAMES = {
    "revel": "REVEL",
    "bayesdel": "BayesDel",
    "alphamissense": "AlphaMissense",
    "spliceai": "SpliceAI",
    "phylop": "PhyloP",
}


def read_score(annotation, key):
    """Validate the expected value range for a predictor."""
    if key == "phylop":
        return number(annotation.get(key))
    lower = -1 if key == "bayesdel" else 0
    return number(annotation.get(key), lower, 1)


def measured(annotation, key, score, condition, strength="supporting"):
    """Record the chosen prediction and comparison."""
    display = DISPLAY_NAMES[key]
    reason = f"{display} {score:g} ({condition})."
    source = annotation.get(f"{key}_source") or display
    evidence = f"in_silico:{key}:{annotation.get('transcript', '')}"
    return result(reason, strength, source, (evidence,))


def missense_tool(annotation, context):
    """Select one missense predictor and screen known limitations of REVEL."""
    chosen = context.get("missense_predictor")
    if chosen is None:
        chosen = annotation.get("missense_predictor", "revel")
    name = token(chosen)
    if name not in MISSENSE_LIMITS:
        return None
    if annotation.get(f"{name}_applicable") is False:
        return None
    if name != "revel" or annotation.get("revel_applicable") is True:
        return name
    gene = token(context.get("gene") or annotation.get("gene"))
    mechanism = token(context.get("mechanism"))
    if gene in {"pkd1", "dsg2"} or mechanism in {"gof", "gain_of_function"}:
        return None
    return name


def protein_result(annotation, context):
    """Choose a supporting PP3 or BP4 row for a missense consequence."""
    name = missense_tool(annotation, context)
    score = read_score(annotation, name) if name else None
    if score is None:
        why = "Prédicteur faux-sens applicable ou score non documenté."
        return {code: result(why) for code in ("PP3", "BP4")}
    pp3 = result("Seuil délétère non atteint.", applicable=False)
    bp4 = result("Pas de prédiction bénigne recevable.", applicable=False)
    limit = MISSENSE_LIMITS[name]
    if score >= limit:
        pp3 = measured(annotation, name, score, f"≥ {limit:.3f}")
    if name == "revel" and score <= 0.290:
        bp4 = measured(annotation, name, score, "≤ 0,290")
    return {"PP3": pp3, "BP4": bp4}


def splice_score(annotation):
    """Read the delta score only when the tool is eligible for the site."""
    if annotation.get("spliceai_applicable") is False:
        return None
    return read_score(annotation, "spliceai")


def splice_result(annotation, consequence, context):
    """Interpret one splice prediction or retain the missense interpretation."""
    if annotation.get("bp7_rna_used") is True:
        if consequence == MISSENSE:
            return protein_result(annotation, context)
        why = "Données d'ARN déjà retenues en BP7."
        return {code: result(why, applicable=False) for code in ("PP3", "BP4")}
    score = splice_score(annotation)
    if score is not None and score >= 0.2:
        pp3 = measured(annotation, "spliceai", score, "≥ 0,20")
        bp4 = result("Effet d'épissage prédit.", applicable=False)
        return {"PP3": pp3, "BP4": bp4}
    if consequence == MISSENSE:
        return protein_result(annotation, context)
    if score is None:
        why = "Score SpliceAI applicable non disponible."
        return {code: result(why) for code in ("PP3", "BP4")}
    pp3 = result("Seuil délétère SpliceAI non atteint.", applicable=False)
    bp4 = result("Seuil bénin SpliceAI non atteint.", applicable=False)
    if score <= 0.1:
        bp4 = measured(annotation, "spliceai", score, "≤ 0,10")
    return {"PP3": pp3, "BP4": bp4}


def bp7_eligible(annotation, consequence):
    """Check that the consequence is in the BP7 computation scope."""
    if consequence in UTR or consequence == "synonymous_variant":
        return True
    offset = annotation.get("splice_offset")
    if consequence not in INTRONIC or type(offset) is not int:
        return False
    return offset >= 7 or offset <= -21


def bp7(annotation, consequence):
    """Pair reassuring splice prediction with lack of conservation."""
    if annotation.get("bp7_rna_used") is True:
        return result("Rapporter les données d'ARN utilisées en BP7(RNA).")
    if consequence in INTRONIC and annotation.get("splice_offset") is None:
        return result("Distance intronique au site d'épissage non documentée.")
    if not bp7_eligible(annotation, consequence):
        why = "Conséquence hors du champ du BP7 computationnel."
        return result(why, applicable=False)
    score = splice_score(annotation)
    conservation = read_score(annotation, "phylop")
    if score is None or conservation is None:
        return result("Prédiction d'épissage ou conservation non documentée.")
    if score > 0.1 or conservation >= 0.1:
        return result("Conditions SpliceAI et PhyloP non réunies.",
                      applicable=False)
    splice_source = annotation.get("spliceai_source") or "SpliceAI"
    row = measured(annotation, "phylop", conservation, "< 0,10")
    row["reason"] += f" SpliceAI {score:g} (≤ 0,10)."
    row["source"] += f" | {splice_source}"
    return row


def assess_prediction(annotation, context):
    """Apply generic computational rules only to covered consequences."""
    consequence = annotation.get("consequence")
    if not isinstance(consequence, str):
        consequence = ""
    if consequence in SPLICE_TESTABLE:
        rows = splice_result(annotation, consequence, context)
    else:
        why = "Conséquence hors du champ des prédictions automatisées."
        if consequence in CANONICAL:
            why = "Site canonique d'épissage à évaluer avec PVS1."
        rows = {code: result(why) for code in ("PP3", "BP4")}
    rows["BP7"] = bp7(annotation, consequence)
    return rows
