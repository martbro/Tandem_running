"""Resolve versioned ClinGen specifications from an offline CSpec copy."""

import gzip
import json
from pathlib import Path
import re
from urllib.parse import urlencode

from .rule_results import number, token


_DATA = Path(__file__).with_name("data")
_BUNDLE = "clingen_cspec.json.gz"
_OVERLAYS = "vcep_profiles.json"
_STRENGTHS = {"standalone", "verystrong", "strong", "moderate", "supporting"}
_FREQUENCY_STRENGTH = {"BA1": "standalone", "BS1": "strong"}
_MODES = {
    "autosomal dominant inheritance": "AD",
    "autosomal recessive inheritance": "AR",
    "x-linked inheritance": "XL",
    "x-linked dominant inheritance": "XLD",
    "x-linked recessive inheritance": "XLR",
    "semidominant inheritance": "SEMIDOMINANT",
    "undetermined mode of inheritance": "",
    "other": "",
}
_MYH7 = {
    "id": "GN002", "version": "2.0", "gene": "MYH7",
    "diseases": ["MONDO:0005021", "MONDO:0005045"], "inheritance": "AD",
    "coverage": "partial", "minimum_an": 10_000,
    "thresholds": {"BA1": 0.001, "BS1": 0.0001},
}


def _identifier(document):
    return document["@id"].rsplit("/", 1)[-1]


def _key(identifier, version):
    return f"{identifier}@{version}"


def _mode(value):
    name = token(value)
    return _MODES.get(name, name.upper())


def _strength(value):
    return token(value).replace(" ", "").replace("_", "").replace("-", "")


def _compatible_mode(wanted, available):
    if not wanted or wanted == "OTHER" or not available or wanted == available:
        return True
    x_linked = {"XL", "XLD", "XLR"}
    if "XL" in (wanted, available):
        return wanted in x_linked and available in x_linked
    return available == "SEMIDOMINANT" and wanted in {"AD", "AR"}


def _scope(document, rule_set, gene, disease, mode, names):
    identifier = _identifier(document)
    version = document["version"]
    return {
        "id": identifier, "version": version,
        "key": _key(identifier, version),
        "source": f'{document["url"]}?{urlencode({"version": version})}',
        "panel": document["affiliation"]["label"], "gene": gene,
        "disease": disease, "disease_name": names.get(disease, disease),
        "inheritance": mode, "status": document["status"],
        "rule_set": _identifier(rule_set),
    }


def _scopes(document, names):
    for rule_set in document.get("ruleSets", ()):
        for gene in rule_set.get("genes", ()):
            for disease in gene.get("diseases") or ({"label": ""},):
                modes = disease.get("modeOfInheritance") or ({},)
                for entry in modes:
                    mode = _mode(entry.get("@label", ""))
                    yield _scope(document, rule_set, gene["label"],
                                 disease["label"], mode, names)


def _match(scope, diseases, requested_mode):
    disease = scope["disease"]
    mode = scope["inheritance"]
    if diseases and disease and disease not in diseases:
        return None
    if not _compatible_mode(requested_mode, mode):
        return None
    exact = bool(diseases and disease in diseases and requested_mode)
    return dict(scope, exact=exact)


def _read_bundle(directory):
    local = Path(directory) / _BUNDLE if directory else None
    path = local if local and local.is_file() else _DATA / _BUNDLE
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def _read_overlays(directory):
    if not directory:
        return []
    path = Path(directory) / _OVERLAYS
    if not path.is_file():
        return []
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, list):
        raise ValueError("vcep_profiles.json doit contenir une liste.")
    return [_validate_overlay(profile) for profile in payload]


def _validate_overlay(profile):
    required = ("id", "version", "gene", "inheritance", "validation")
    if not isinstance(profile, dict) or not all(
        token(profile.get(field)) for field in required
    ):
        raise ValueError("Identité, transmission et validation VCEP requises.")
    diseases = profile.get("diseases")
    if not isinstance(diseases, list) or not diseases or not all(
        re.fullmatch(r"mondo:\d{7}", token(disease)) for disease in diseases
    ):
        raise ValueError("Un profil VCEP requiert des identifiants MONDO.")
    allowed = profile.get("allowed", {})
    manual = profile.get("manual", [])
    if not _valid_policy(allowed, manual):
        raise ValueError("Codes ou forces invalides dans un profil VCEP.")
    _validate_thresholds(profile)
    return dict(profile, coverage="partial", id=profile["id"].strip().upper(),
                version=profile["version"].strip(),
                gene=profile["gene"].strip().upper(),
                diseases=[disease.strip().upper() for disease in diseases],
                inheritance=_mode(profile["inheritance"]))


def _valid_code(code):
    return bool(re.fullmatch(r"(?:PVS|PS|PM|PP|BA|BS|BP)\d", code))


def _valid_policy(allowed, manual):
    return (isinstance(allowed, dict) and isinstance(manual, list)
            and all(isinstance(code, str) and _valid_code(code)
                    for code in manual)
            and all(_valid_code(code) and isinstance(levels, list)
                    and all(_strength(level) in _STRENGTHS for level in levels)
                    for code, levels in allowed.items()))


def _validate_thresholds(profile):
    thresholds = profile.get("thresholds", {})
    if not isinstance(thresholds, dict) or any(
        key not in _FREQUENCY_STRENGTH or not number(value, 0, 1)
        for key, value in thresholds.items()
    ):
        raise ValueError("Les seuils VCEP BA1 et BS1 doivent être dans ]0, 1].")
    if thresholds and (profile.get("metric") != "gnomad_popmax_faf"
            or type(profile.get("minimum_an")) is not int
            or profile["minimum_an"] < 1):
        raise ValueError("FAF popmax gnomAD et AN minimal requis.")
    if len(thresholds) == 2 and thresholds["BS1"] >= thresholds["BA1"]:
        raise ValueError("Le seuil VCEP BS1 doit être inférieur au seuil BA1.")


def _official_documents(bundle):
    return {
        _key(_identifier(doc), doc["version"]): doc
        for doc in bundle["specifications"] if doc["cspecStatus"] == "Released"
    }


def _active_documents(bundle, docs):
    for doc in bundle["registry"]["data"]:
        if doc["status"] != "Released":
            continue
        key = _key(_identifier(doc), doc["version"])
        if not docs.get(key, {}).get("legacyFullySuperseded", False):
            yield doc


def _profile_matches(profile, match):
    return (match["exact"] and match["id"] == profile["id"]
            and match["version"] == profile["version"]
            and match["gene"] == profile["gene"]
            and match["disease"] in profile["diseases"]
            and match["inheritance"] == profile["inheritance"])


def _auto_profile(matches, overlays):
    if len(matches) != 1:
        return None
    match = matches[0]
    local = [p for p in overlays if _profile_matches(p, match)]
    if len(local) > 1:
        raise ValueError("Plusieurs profils locaux s'appliquent au même VCEP.")
    if local:
        return dict(local[0], source=match["source"])
    if _profile_matches(_MYH7, match):
        return dict(_MYH7, source=match["source"])
    return None


def _guides(documents, checksums):
    return {
        (key, _identifier(rule_set)): {
            "title": doc["label"],
            "detail_sha256": checksums.get(_detail_url(doc), ""),
            "updated": doc.get("lastUpdated"),
            "criteria": _criteria(rule_set.get("criteriaCodes", ())),
        }
        for key, doc in documents.items()
        for rule_set in doc.get("ruleSets", ())
    }


def _detail_url(document):
    params = urlencode({"version": document["version"]})
    return f'{document["@id"]}?{params}'


def _guidance(match, guides):
    details = guides.get((match["key"], match["rule_set"]))
    return dict(match, **(details or {"title": match["panel"], "criteria": {}}))


def _criteria(codes):
    return {code["label"].upper(): _criterion(code) for code in codes}


def _criterion(code):
    return {
        "description": code.get("description", ""),
        "strengths": {
            _strength(row["label"]): {
                "applicability": row.get("applicability", ""),
                "description": row.get("description", ""),
            } for row in code.get("evidenceStrengths", ())
        },
    }


def _scope_alerts(matches, profile):
    if not matches:
        return []
    if len({match["key"] for match in matches}) > 1:
        return ["Plusieurs VCEP possibles. Choisir une spécification."]
    if not any(match["exact"] for match in matches):
        return ["Préciser la maladie MONDO et la transmission VCEP."]
    if profile:
        return ["Profil VCEP partiel. Appliquer les autres règles en revue."]
    return ["Appliquer les règles du VCEP en revue avant de classer."]


class Catalog:
    """Match released specifications without network access."""

    def __init__(self, directory=None):
        bundle = _read_bundle(directory)
        self.overlays = _read_overlays(directory)
        self.catalog = {key: bundle[key] for key in
                        ("source", "retrieved", "index_sha256")}
        self.documents = _official_documents(bundle)
        checksums = bundle.get("specification_sha256", {})
        self.guides = _guides(self.documents, checksums)
        self.names = bundle.get("disease_labels", {})
        self.scopes = {}
        for document in _active_documents(bundle, self.documents):
            for scope in _scopes(document, self.names):
                self.scopes.setdefault(scope["gene"], []).append(scope)

    def _diseases(self, disease):
        name = token(disease)
        identifier = re.fullmatch(r"(?:.*/)?mondo[:_](\d+)", name)
        if identifier:
            return {f"MONDO:{identifier[1]}"}
        return {key for key, label in self.names.items()
                if token(label) == name}

    def for_gene(self, gene):
        return self.scopes.get(token(gene).upper(), ())

    def for_disease(self, gene, disease):
        identifiers = self._diseases(disease)
        return [row for row in self.for_gene(gene)
                if row["disease"] in identifiers]

    def resolve(self, gene, disease, inheritance):
        requested = self._diseases(disease)
        mode = _mode(inheritance)
        matches = [
            match for scope in self.for_gene(gene)
            if (match := _match(scope, requested, mode))
        ]
        exact = [match for match in matches if match["exact"]]
        profile = _auto_profile(exact, self.overlays)
        guidance = [_guidance(match, self.guides) for match in matches]
        _check_automatic_levels(profile, guidance)
        return {"matches": matches, "profile": profile, "guidance": guidance,
                "alerts": _scope_alerts(matches, profile),
                "requires_review": bool(matches), "catalog": self.catalog}


def _reviewed_guidance(resolution, reviewed):
    guidance = resolution.get("guidance", ())
    if isinstance(reviewed, str):
        return [row for row in guidance if row["key"] == reviewed]
    exact = [row for row in guidance if row["exact"]]
    return exact or guidance


def _applicability(row, code, strength):
    criterion = row.get("criteria", {}).get(code, {})
    selected = criterion.get("strengths", {}).get(_strength(strength), {})
    status = token(selected.get("applicability"))
    if status.startswith("not applicable"):
        return False
    if status.startswith("applicable"):
        return True
    return None


def _check_automatic_levels(profile, guidance):
    if not profile:
        return
    key = _key(profile["id"], profile["version"])
    rows = [row for row in guidance if row["key"] == key and row["exact"]]
    for code in profile.get("thresholds", {}):
        strength = _FREQUENCY_STRENGTH[code]
        if not rows or not all(_applicability(r, code, strength) for r in rows):
            raise ValueError(f"Force du profil VCEP incompatible ({code}).")


def vcep_level(resolution, code, strength, reviewed=None):
    """Use explicit CSpec labels, preserving unresolved or mixed entries."""
    rows = _reviewed_guidance(resolution, reviewed)
    code = code.upper()
    profile = resolution.get("profile") or {}
    if not _level_allowed(profile, code, strength):
        return False
    values = {_applicability(row, code, strength) for row in rows}
    return values.pop() if len(values) == 1 else None


def _level_allowed(profile, code, strength):
    allowed = profile.get("allowed", {}).get(code)
    return allowed is None or _strength(strength) in map(_strength, allowed)


def _population_profile_values(profile, annotation, context):
    faf = number(annotation.get("popmax_faf"), 0, 1)
    total = annotation.get("popmax_an")
    source = annotation.get("population_source")
    if (faf is None or type(total) is not int
            or total < profile["minimum_an"]):
        return None
    if "gnomad" not in token(source):
        return None
    if annotation.get("population_pass") is not True:
        return None
    if context.get("vcep_population_confirmed") is not True:
        return None
    return faf, total, source


def _population_row(code, profile, values):
    faf, total, source = values
    limit = profile["thresholds"][code]
    reason = (f"FAF popmax {faf:.8g}, AN {total}, seuil {limit:.8g}. "
              "Comparabilité de la population confirmée.")
    provenance = f'{profile["source"]} | {source}'
    if profile.get("validation"):
        provenance += f' | {profile["validation"]}'
    return {
        "status": "applied", "strength": _FREQUENCY_STRENGTH[code],
        "reason": reason,
        "source": provenance,
        "evidence_ids": [f"vcep:population:{source}:FAF_popmax"],
    }


def apply_profile(profile, annotation, context, criteria):
    """Assess the supported frequency subset under qualified inputs."""
    if not profile or not profile.get("thresholds"):
        return {}
    values = _population_profile_values(profile, annotation, context)
    if not values:
        return {}
    for code in ("BA1", "BS1"):
        limit = profile["thresholds"].get(code)
        if limit is None or values[0] < limit:
            continue
        if code in profile.get("manual", []):
            return {}
        if criteria.get(code, {}).get("status") == "applied":
            return {}
        row = _population_row(code, profile, values)
        if _level_allowed(profile, code, row["strength"]):
            return {code: row}
        return {}
    return {}
