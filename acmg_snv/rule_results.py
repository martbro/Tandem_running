"""Represent assessed criteria consistently."""

from math import isfinite


SUPPLEMENT = "ACGS 2026, tableau supplémentaire S4"


def result(reason, strength=None, source=None, ids=(), applicable=True):
    """Return a criterion result before scoring."""
    status = "not_evaluable" if applicable else "not_applicable"
    row = {
        "status": "applied" if strength else status,
        "reason": reason,
        "source": f"{SUPPLEMENT} | {source}" if source else SUPPLEMENT,
        "evidence_ids": list(ids),
    }
    if strength:
        row["strength"] = strength
    return row


def number(value, minimum=None, maximum=None):
    """Accept finite numeric annotation values in the stated range."""
    if type(value) not in (int, float):
        return None
    try:
        finite = isfinite(value)
    except OverflowError:
        return None
    if not finite:
        return None
    if minimum is not None and value < minimum:
        return None
    if maximum is not None and value > maximum:
        return None
    return value


def token(value):
    """Normalize a name or return an empty string."""
    return value.strip().casefold() if isinstance(value, str) else ""
