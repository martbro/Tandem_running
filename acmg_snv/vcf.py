"""Streaming annotation of GRCh38 VCF files, with one result per ALT."""

import json
import re
import struct
import zlib
from dataclasses import dataclass

from .ann_reference import primary_chromosome
from .annotations import csq_columns, parse_vcf_annotation
from .model import Variant


_INFO = {
    "ACMG_CLASS": ("String", "ACMG class code"),
    "ACMG_SCORE": ("Integer", "Bayesian point total"),
    "ACMG_STATUS": ("String", "ACMG review status"),
}
_REPORT_FIELDS = ("classification", "score", "status")
_HEADER_ID = re.compile(r"^##INFO=<ID=([^,>]+)")
_META_PAIR = re.compile(r'(?:^|,)(\w+)=("(?:[^"\\]|\\.)*"|[^,]*)')
_FILE_FORMAT = re.compile(r"##fileformat=VCFv4\.[0-9]+")
_OLD_BUILD = re.compile(r"\b(?:GRCh37|hg19|hs37d5|b37)\b", re.IGNORECASE)
_COLUMNS = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
_BLOCK_HEADER = bytes.fromhex("1f8b08040000000000ff060042430200")
_BLOCK_PAYLOAD = 64 * 1024 - 256


class BgzfWriter:
    """Write UTF-8 VCF text as concatenated BGZF blocks."""

    def __init__(self, stream):
        self.stream = stream
        self.pending = bytearray()

    def write(self, text):
        self.pending.extend(text.encode("utf-8"))
        while len(self.pending) >= _BLOCK_PAYLOAD:
            self._block(self.pending[:_BLOCK_PAYLOAD])
            del self.pending[:_BLOCK_PAYLOAD]

    def _block(self, data):
        compressor = zlib.compressobj(wbits=-15)
        compressed = compressor.compress(data) + compressor.flush()
        bsize = len(_BLOCK_HEADER) + 2 + len(compressed) + 8 - 1
        header = _BLOCK_HEADER + struct.pack("<H", bsize)
        tail = struct.pack("<II", zlib.crc32(data), len(data))
        self.stream.write(header + compressed + tail)

    def finish(self):
        if self.pending:
            self._block(self.pending)
            self.pending.clear()
        self._block(b"")


@dataclass
class Batch:
    """Shared classification settings for one VCF stream."""

    classifier: object
    context: dict
    contexts: dict

    def classify(self, variant, row, columns):
        """Evaluate one ALT with its applicable context and annotations."""
        context = self.context | self.contexts.get(variant.key, {})
        annotation = None
        if columns:
            annotation = parse_vcf_annotation(variant, row, columns) or None
        return self.classifier.classify(variant, context, annotation)


def _definition(key):
    kind, description = _INFO[key]
    return f'##INFO=<ID={key},Number=A,Type={kind},Description="{description}">'


def _line_error(number, message):
    return ValueError(f"VCF, ligne {number} — {message}")


def _check_metadata(line, number):
    if line.startswith(("##reference=", "##contig=")) and _OLD_BUILD.search(
        line
    ):
        raise _line_error(number, "GRCh37 déclaré. GRCh38 est requis.")
    match = _HEADER_ID.match(line)
    if match and match[1] in _INFO:
        key = match[1]
        fields = dict(_META_PAIR.findall(line.removeprefix("##INFO=<")[:-1]))
        if fields.get("Number") != "A" or fields.get("Type") != _INFO[key][0]:
            raise _line_error(number, f"Définition incompatible de {key}.")
        return False
    return True


def _check_columns(line, number):
    keys = line.split("\t")
    valid = tuple(keys[:8]) == _COLUMNS
    if len(keys) > 8:
        valid = valid and len(keys) > 9 and keys[8] == "FORMAT"
    samples = keys[9:]
    if not valid or len(samples) != len(set(samples)) or "" in samples:
        raise _line_error(number, "Colonnes de l’en-tête invalides.")
    return len(keys)


def _header(numbered, output):
    csq = ()
    for number, row in numbered:
        line = row.rstrip("\r\n")
        if number == 1 and not _FILE_FORMAT.fullmatch(line):
            raise _line_error(number, "En-tête ##fileformat VCFv4.x requis.")
        if line.startswith("##"):
            if line.startswith("##INFO=<ID=CSQ,"):
                csq = csq_columns(line)
            if _check_metadata(line, number):
                output.write(line + "\n")
            continue
        columns = _check_columns(line, number)
        for key in _INFO:
            output.write(_definition(key) + "\n")
        output.write(line + "\n")
        return columns, csq
    raise ValueError("VCF sans en-tête #CHROM.")


def _fields(line, columns, number):
    record = line.rstrip("\r\n").split("\t")
    if len(record) != columns or line.startswith("#"):
        raise _line_error(number, f"{columns} colonnes attendues.")
    position = record[1]
    if not position.isascii() or not position.isdecimal() or int(position) < 1:
        raise _line_error(number, "Position entière positive attendue.")
    if not all(record[:8]):
        raise _line_error(number, "Une colonne VCF obligatoire est vide.")
    return record


def _out_of_scope(chrom, pos, ref, alt, reason):
    return {
        "variant": {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt},
        "classification": None,
        "score": None,
        "status": "out_of_scope",
        "criteria": [],
        "alerts": [reason],
    }


def _evaluate(record, alt, line, batch, csq):
    chrom, pos, _, ref = record[:4]
    try:
        canonical = primary_chromosome(chrom) or chrom
        variant = Variant(canonical, int(pos), ref, alt)
    except ValueError as exc:
        return _out_of_scope(chrom, int(pos), ref, alt, str(exc))
    return batch.classify(variant, line, csq)


def _info_fields(existing):
    if existing == ".":
        return []
    return [
        field
        for field in existing.split(";")
        if field.partition("=")[0] not in _INFO
    ]


def _entry(result, key):
    value = result.get(key)
    return "." if value is None else str(value)


def _update_info(record, results):
    info = _info_fields(record[7])
    for key, field in zip(_INFO, _REPORT_FIELDS):
        values = ",".join(_entry(result, field) for result in results)
        info.append(f"{key}={values}")
    record[7] = ";".join(info)


def _write_details(details, result, number, index, identifier):
    report = {**result, "vcf": {"line": number, "alt_index": index,
                               "id": identifier}}
    if "vcep" in report:
        report["vcep"] = {
            key: value for key, value in report["vcep"].items()
            if key != "guidance"
        }
    details.write(json.dumps(report, ensure_ascii=False, allow_nan=False))
    details.write("\n")


def _record(line, columns, csq, batch, number, details):
    record = _fields(line, columns, number)
    if record[4] == ".":
        record[7] = ";".join(_info_fields(record[7])) or "."
        return record
    results = []
    for index, alt in enumerate(record[4].split(","), 1):
        if not alt or alt == ".":
            raise _line_error(number, "ALT manquant au site multiallélique.")
        try:
            result = _evaluate(record, alt, line, batch, csq)
        except ValueError as exc:
            raise _line_error(number, str(exc)) from exc
        _write_details(details, result, number, index, record[2])
        results.append(result)
    _update_info(record, results)
    return record


def annotate_vcf(source, output, details, batch):
    """Write a VCF with ACMG INFO values and stream detailed JSONL results."""
    numbered = enumerate(source, 1)
    columns, csq = _header(numbered, output)
    for number, line in numbered:
        record = _record(line, columns, csq, batch, number, details)
        output.write("\t".join(record) + "\n")
