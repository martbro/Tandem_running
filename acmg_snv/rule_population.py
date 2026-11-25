"""Assess the generic population evidence using qualified local observations."""

from fractions import Fraction

from .rule_results import number, result, token


POPULATION_CODES = ("BA1", "BS1", "PM2")


def counts(annotation):
    """Require an eligible group and a quality assessment at this site."""
    source = annotation.get("population_source")
    group = annotation.get("population_group")
    count = annotation.get("population_ac")
    total = annotation.get("population_an")
    if not token(source) or not token(group):
        return None
    if annotation.get("population_pass") is not True:
        return None
    if annotation.get("population_eligible") is not True:
        return None
    if type(count) is not int or type(total) is not int:
        return None
    if total <= 0 or not 0 <= count <= total:
        return None
    return count, total, source, group


def threshold(context, field, default=None):
    """Accept a disease threshold when its derivation has been recorded."""
    if field not in context:
        return default
    score = number(context[field], 0, 1)
    if not score or not token(context.get("population_threshold_source")):
        return None
    return Fraction(str(score))


def benign_thresholds(context):
    """Load BA1 and BS1 cutoffs and reject inconsistent pairs."""
    ba1 = threshold(context, "ba1_threshold", Fraction(1, 20))
    bs1 = threshold(context, "bs1_threshold")
    if context.get("ba1_exception") is True:
        ba1 = threshold(context, "ba1_threshold")
    if ba1 is not None and bs1 is not None and ba1 <= bs1:
        return None, None
    return ba1, bs1


def compare_frequency(freq, limit, strength, observation):
    """Determine whether an eligible frequency exceeds its cutoff."""
    source, evidence_id, summary = observation
    if limit is None:
        return result("Seuil de fréquence applicable non documenté.")
    if freq > limit:
        why = f"{summary}. Seuil {float(limit):.8g} dépassé."
        return result(why, strength, source, (evidence_id,))
    relation = f"{summary}. Seuil {float(limit):.8g} non dépassé."
    return result(relation, source=source, applicable=False)


def rarity(annotation, context, values, bs1, observation):
    """Apply PM2 only for resolved absence or qualified extreme rarity."""
    count, total, _, _ = values
    source, evidence_id, summary = observation
    if total < 50_000:
        why = "Couverture insuffisante pour PM2 (moins de 50 000 allèles)."
        return result(why, source=source)
    if count == 0 and annotation.get("population_absence_verified") is True:
        why = f"{summary}. Absence vérifiée dans les données pertinentes."
        return result(why, "moderate", source, (evidence_id,))
    if count == 0:
        return result("Absence dans les jeux pertinents non vérifiée.")
    if token(context.get("inheritance")) not in {"ad", "ar", "xlr"}:
        why = "Rareté à interpréter pour ce mode de transmission."
        return result(why)
    if bs1 is not None and Fraction(count, total) <= bs1 / 10:
        if context.get("pm2_rare_compatible") is True:
            why = f"{summary}. Rareté compatible, fréquence ≤ BS1/10."
            return result(why, "moderate", source, (evidence_id,))
        return result("Compatibilité de la rareté avec la maladie à valider.")
    if bs1 is None:
        return result("Seuil BS1 absent pour interpréter la rareté.")
    return result("Fréquence supérieure à BS1/10.", applicable=False)


def assess_population(annotation, context):
    """Prefer BA1 over BS1 and count each frequency observation once."""
    values = counts(annotation)
    if values is None:
        why = "Données de fréquence, groupe ou qualité non documentés."
        return {code: result(why) for code in POPULATION_CODES}
    freq, observation = describe(values, context)
    ba1, bs1 = benign_thresholds(context)
    rows = {"BA1": compare_frequency(freq, ba1, "standalone", observation)}
    rows["BS1"] = compare_frequency(freq, bs1, "strong", observation)
    if rows["BA1"]["status"] == "applied":
        why = "Fréquence déjà retenue en BA1."
        rows["BS1"] = result(why, applicable=False)
    if rows["BS1"]["status"] == "applied":
        high_penetrance = context.get("high_penetrance") is True
        rows["BS1"]["standalone"] = high_penetrance and (
            token(context.get("inheritance")) == "ad"
        )
    if any(row["status"] == "applied" for row in rows.values()):
        rows["PM2"] = result("Fréquence non compatible.", applicable=False)
    else:
        rows["PM2"] = rarity(annotation, context, values, bs1, observation)
    return rows


def describe(values, context):
    """Record numerator, denominator, group, and their provenance."""
    count, total, source, group = values
    provenance = context.get("population_threshold_source")
    evidence_id = f"population:{source}:{group}"
    if token(provenance):
        source = f"{source} | seuils {provenance}"
    freq = Fraction(count, total)
    summary = f"{count}/{total} ({float(freq):.8g}), groupe {group}"
    return freq, (source, evidence_id, summary)
