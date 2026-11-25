"""Implement documented branches of the ACGS PVS1 and PM4 decision tree."""

import re

from .rule_prediction import CANONICAL, SPLICE_TESTABLE
from .rule_results import number, result, token


FIGURE = "ACGS 2026, figure 1"
TRUNCATING = {"stop_gained", "start_lost", "stop_lost"}
VALIDITIES = {"moderate", "strong", "definitive"}
PLUS_TWO = re.compile(r"\+2[ACGT]>[ACGT]$", re.I)


def recorded(data, reason, strength):
    """Attribute the decision tree result to the local assessment."""
    source = f"{FIGURE} | {data['source']}"
    return result(reason, strength, source, (data["evidence_id"],))


def affected_region(data):
    """Score a region deleted or skipped without triggering NMD."""
    critical = data.get("critical_domain")
    fraction = number(data.get("protein_fraction_lost"), 0, 1)
    if critical is True:
        return recorded(data, "Perte d'un domaine critique.", "very_strong")
    if critical is not False or fraction is None:
        return result("Domaine critique ou proportion protéique non évalués.")
    strength = "strong" if fraction > 0.10 else "moderate"
    return recorded(data, f"Perte de {fraction:.1%} de la protéine.", strength)


def five_prime(data):
    """Assess initiation loss and premature stops in the 5′ escape region."""
    if data.get("alternative_functional_transcript") is True:
        why = "Transcrit fonctionnel avec un codon initiateur préservé."
        return result(why, applicable=False)
    if data.get("alternative_functional_transcript") is not False:
        return result("Transcrits fonctionnels alternatifs non évalués.")
    alternative = data.get("alternative_start_before_critical")
    if alternative is False:
        why = "Aucun initiateur alternatif en amont des domaines critiques."
        return recorded(data, why, "very_strong")
    if alternative is not True:
        return result("Codons d'initiation alternatifs non évalués.")
    upstream = data.get("pathogenic_upstream")
    if type(upstream) is not bool:
        return result("Variants pathogènes en amont non évalués.")
    if upstream:
        why = "Initiateur alternatif et variant pathogène en amont."
        return recorded(data, why, "moderate")
    why = "Initiateur alternatif sans variant pathogène connu en amont."
    return recorded(data, why, "supporting")


def premature_stop(data):
    """Assess a premature termination event with known transcript position."""
    region = data.get("nmd_region")
    if region == "5prime_escape":
        return five_prime(data)
    if region == "3prime_escape" or data.get("predicted_nmd") is False:
        return affected_region(data)
    if region == "standard" and data.get("predicted_nmd") is True:
        why = "Codon d'arrêt prématuré avec dégradation de l'ARN prévue."
        return recorded(data, why, "very_strong")
    return result("Effet sur la dégradation de l'ARN à évaluer.")


def lost_stop(data):
    """Use the nonstop decay path when the following stop codon is absent."""
    next_stop = data.get("stop_before_polya")
    if next_stop is True:
        why = "Arrêt alternatif en phase avec allongement de la protéine."
        return result(why, applicable=False)
    if next_stop is not False:
        return result("Arrêt en phase avant polyadénylation non évalué.")
    why = "Absence d'arrêt en phase avant polyadénylation, NSD probable."
    strength = "strong"
    if data.get("nsd_confirmed") is True:
        why = "Dégradation de l'ARN sans codon d'arrêt (NSD) démontrée."
        strength = "very_strong"
    return recorded(data, why, strength)


def splice_lof(data, annotation):
    """Apply the loss-of-function tree after assessing splice exceptions."""
    change = splice_change(data, annotation)
    if change == "+2C>T":
        why = "Substitution +2C>T sans motif de perte du site canonique."
        return result(why, applicable=False)
    if change == "+2T>C" and data.get("rna_confirmed") is not True:
        loss = number(data.get("spliceai_donor_loss"), 0, 1)
        if loss is None or loss < 0.8:
            return result("Site +2T>C à examiner pour un donneur GC résiduel.")
    if data.get("splice_cautions_reviewed") is not True:
        return result("Site d'épissage et sites cryptiques à évaluer.")
    if data.get("splice_out_of_frame") is True:
        return premature_stop(data)
    if data.get("splice_out_of_frame") is False:
        return affected_region(data)
    return result("Cadre de lecture de l'événement d'épissage non évalué.")


def splice_change(data, annotation):
    """Read a +2 donor substitution from the selected transcript."""
    match = PLUS_TWO.search(str(annotation.get("hgvsc") or ""))
    if match:
        return match[0].upper()
    return str(data.get("splice_dinucleotide_change") or "").upper()


def transcript_assessed(data):
    """Require traceable assessment of transcript impact and rescue."""
    if not isinstance(data, dict):
        return False
    traceable = token(data.get("source")) and token(data.get("evidence_id"))
    return bool(traceable) and data.get("transcript_relevant") is True


def pvs1(annotation, context, consequence):
    """Gate the decision tree on disease mechanism and transcript relevance."""
    data = annotation.get("lof") or {}
    rna = isinstance(data, dict) and data.get("rna_confirmed") is True
    rna = rna and consequence in CANONICAL | SPLICE_TESTABLE
    if consequence not in TRUNCATING | CANONICAL and not rna:
        why = "Conséquence sans événement nul démontré ou prédit."
        return result(why, applicable=False)
    mechanism = token(context.get("mechanism"))
    if mechanism in {"gof", "gain_of_function", "dominant_negative"}:
        why = "Mécanisme causal hors perte de fonction."
        return result(why, applicable=False)
    validity = token(context.get("gene_disease_validity"))
    valid_lof = mechanism in {"lof", "loss_of_function"}
    if not valid_lof or validity not in VALIDITIES:
        why = "Perte de fonction ou validité gène-maladie à documenter."
        return result(why)
    if not transcript_assessed(data):
        return result("Impact sur le transcrit pertinent à vérifier.")
    if type(data.get("rescue_possible")) is not bool:
        return result("Mécanisme de sauvetage à vérifier.")
    if data["rescue_possible"]:
        return result("Mécanisme de sauvetage plausible.", applicable=False)
    return pvs1_event(data, annotation, consequence, rna)


def pvs1_event(data, annotation, consequence, rna):
    """Select the relevant Figure 1 branch."""
    if consequence == "start_lost":
        return five_prime(data)
    if consequence == "stop_lost":
        return lost_stop(data)
    if consequence == "stop_gained":
        return premature_stop(data)
    row = splice_lof(data, annotation)
    if rna and row["status"] == "applied":
        row["reason"] = "PVS1(RNA). " + row["reason"]
    return row


def pm4(annotation, consequence):
    """Use a characterized in-frame downstream stop for a stop-loss SNV."""
    if consequence != "stop_lost":
        return result("Allongement protéique à évaluer si pertinent.")
    data = annotation.get("lof")
    if not transcript_assessed(data):
        return result("Conséquence sur le transcrit pertinent à documenter.")
    if data.get("stop_before_polya") is True:
        why = "Arrêt alternatif en phase avec allongement de la protéine."
        return recorded(data, why, "moderate")
    if data.get("stop_before_polya") is False:
        return result("Arrêt alternatif en phase absent.", applicable=False)
    return result("Arrêt en phase avant polyadénylation non évalué.")


def assess_lof(annotation, context):
    """Assess PVS1 and the supported loss-of-stop branch for PM4."""
    consequence = annotation.get("consequence")
    if not isinstance(consequence, str):
        why = "Conséquence transcriptomique non renseignée."
        return {code: result(why) for code in ("PVS1", "PM4")}
    rows = {"PVS1": pvs1(annotation, context, consequence)}
    rows["PM4"] = pm4(annotation, consequence)
    return rows
