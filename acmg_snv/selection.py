"""Pick a clinically relevant annotation without guessing between genes."""

from .model import INHERITANCES


CODING = {
    "missense_variant", "synonymous_variant", "stop_gained", "stop_lost",
    "start_lost", "stop_retained_variant", "start_retained_variant",
    "coding_sequence_variant", "protein_altering_variant", "initiator_codon",
    "incomplete_terminal_codon_variant",
}
SPLICE = {
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "splice_polypyrimidine_tract_variant", "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
}
PRIORITY = (
    "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
    "start_lost", "stop_lost", "missense_variant", "synonymous_variant",
)


def prepare_context(context):
    value = dict(context or {})
    for key in ("inheritance", "gene", "mechanism", "disease", "transcript"):
        if key in value:
            value[key] = str(value[key] or "").strip()
    for key in ("inheritance", "gene", "mechanism"):
        value[key] = value.get(key, "").upper()
    inheritance = value["inheritance"]
    if inheritance and inheritance not in INHERITANCES:
        raise ValueError("Transmission attendue non reconnue.")
    if value.get("phase") is not None and not isinstance(value["phase"], list):
        raise ValueError("La phase doit être une liste d'observations.")
    return value


def transcript_matches(first, second):
    if not first or not second:
        return False
    if "." in first and "." in second:
        return first == second
    return first.split(".")[0] == second.split(".")[0]


def same_context(first, second):
    return first.get("gene") == second.get("gene") and transcript_matches(
        first.get("transcript", ""), second.get("transcript", "")
    )


def combine_annotations(earlier, later):
    rows = []
    site_rows = []
    for annotation in (*earlier, *later):
        row = dict(annotation)
        if not row.get("gene") and not row.get("transcript"):
            site_rows.append(row)
            continue
        match = next((item for item in rows if same_context(item, row)), None)
        if match is None:
            rows.append(row)
        else:
            match.update({k: v for k, v in row.items() if v not in (None, "")})
    return overlay_site_data(rows, site_rows)


def overlay_site_data(rows, site_rows):
    if not rows:
        return site_rows
    fields = {
        "phylop", "phylop_source", "popmax_faf", "popmax_an", "clinvar_id",
        "clinvar_significance", "clinvar_review_status",
    }
    for site in site_rows:
        site_data = {key: val for key, val in site.items() if (
            key.startswith("population_") or key in fields
        )}
        rows = [site_data | row for row in rows]
    return rows


def mane_rank(annotation):
    mane = str(annotation.get("mane", "")).lower().replace("_", " ")
    if "select" in mane or mane in {"true", "1"}:
        return 0
    if "plus" in mane:
        return 1
    return 2


def consequence_set(annotation):
    consequence = annotation.get("consequence") or ""
    terms = str(consequence).replace(",", "&").split("&")
    return {item.strip() for item in terms}


def select_consequence(annotation):
    terms = consequence_set(annotation)
    ordered = list(PRIORITY) + sorted((CODING | SPLICE) - set(PRIORITY))
    return next((term for term in ordered if term in terms), None)


def annotation_context(rows, context):
    target = context.get("target_genes") or []
    if isinstance(target, str):
        target = target.replace(",", " ").split()
    targets = {str(gene).strip().upper() for gene in target}
    gene, transcript = context.get("gene"), context.get("transcript")
    return [
        row for row in rows
        if (not gene or row.get("gene", "").upper() == gene)
        and (not targets or row.get("gene", "").upper() in targets)
        and (not transcript or transcript_matches(row.get("transcript"),
                                                  transcript))
    ]


def choose_annotation(rows, context):
    candidates = annotation_context(rows, context)
    if not candidates:
        return None, "Aucune annotation pour le gène ou transcrit demandé."
    scoped = [row for row in candidates if select_consequence(row)]
    candidates = scoped or candidates
    if len({row.get("gene") for row in candidates}) > 1:
        return None, "Plusieurs gènes candidats : choisissez le gène clinique."
    rank = min(map(mane_rank, candidates))
    best = [row for row in candidates if mane_rank(row) == rank]
    if len(best) != 1:
        return None, "Plusieurs annotations : choisissez le transcrit."
    return dict(best[0]), None


def brief_annotation(annotation):
    fields = ("gene", "transcript", "mane", "consequence", "hgvsc", "hgvsp")
    return {field: annotation[field] for field in fields if field in annotation}
