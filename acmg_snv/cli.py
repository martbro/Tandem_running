"""Commands for the local browser interface and Linux VCF pipelines."""

import argparse
import gzip
import json
import os
import sqlite3
import sys
import tempfile
from contextlib import contextmanager
from itertools import combinations
from pathlib import Path

from .annotations import (
    build_annotation_db,
    build_fasta_index,
    parse_annotation,
)
from .ann_mane import build_mane_db
from .engine import Classifier
from .model import INHERITANCES, Variant
from .vcep import Catalog
from .vcf import Batch, BgzfWriter, annotate_vcf


@contextmanager
def _classifier(args):
    classifier = Classifier(
        annotation_db=args.annotation_db,
        reference=args.reference,
        mane_db=args.mane_db,
        vcep_dir=args.vcep_dir,
        ps3=args.ps3,
    )
    try:
        yield classifier
    finally:
        classifier.close()


def _read_text(path):
    if path == "-":
        return sys.stdin.read()
    return Path(path).read_text(encoding="utf-8-sig")


def _object(path):
    if path is None:
        return {}
    data = json.loads(_read_text(path))
    if not isinstance(data, dict):
        raise ValueError(f"{path} doit contenir un objet JSON.")
    return data


def _context(args):
    context = _object(args.context)
    if args.inheritance is not None:
        context["inheritance"] = args.inheritance
    return context


def _contexts(path):
    contexts = {}
    for key, context in _object(path).items():
        canonical = Variant.parse(key).key
        if not isinstance(context, dict):
            raise ValueError(f"Le contexte de {canonical} doit être un objet.")
        if canonical in contexts:
            raise ValueError(f"Contexte défini en double pour {canonical}.")
        contexts[canonical] = context
    return contexts


def _classify(args):
    _check_stdin(args.annotation, args.context, args.overrides)
    variant = Variant.parse(args.variant)
    context = _context(args)
    annotation = None
    if args.annotation is not None:
        annotation = parse_annotation(_read_text(args.annotation), variant)
    with _classifier(args) as classifier:
        result = classifier.classify(
            variant, context, annotation, _object(args.overrides)
        )
    print(json.dumps(result, indent=2, ensure_ascii=True, allow_nan=False))


def _web(args):
    from .web import serve

    with _classifier(args) as classifier:
        serve(
            classifier, port=args.port, open_browser=not args.no_browser
        )


def _vcep(args):
    catalog = Catalog(args.vcep_dir)
    result = catalog.resolve(args.gene, args.disease, args.inheritance)
    print(json.dumps(result, indent=2, ensure_ascii=True, allow_nan=False))


def _open_text(path, mode, compressed):
    opener = gzip.open if compressed else open
    return opener(path, mode, encoding="utf-8", newline="")


def _compressed(path):
    return path.lower().endswith((".gz", ".bgz"))


@contextmanager
def _input(path):
    if path == "-":
        yield sys.stdin
        return
    with _open_text(path, "rt", _compressed(path)) as source:
        yield source


@contextmanager
def _staged(path):
    if path == "-":
        yield path
        return
    destination = Path(path)
    fd, name = tempfile.mkstemp(prefix=".acmg-", dir=destination.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        yield temporary
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def _output(temporary, destination, bgzf=False):
    if destination == "-":
        yield sys.stdout
        sys.stdout.flush()
        return
    compressed = _compressed(destination)
    if compressed and bgzf:
        with open(temporary, "wb") as stream:
            writer = BgzfWriter(stream)
            yield writer
            writer.finish()
        return
    with _open_text(temporary, "wt", compressed) as stream:
        yield stream


@contextmanager
def _outputs(output, report):
    with _staged(output) as staged_output, _staged(report) as staged_report:
        with (
            _output(staged_output, output, bgzf=True) as target,
            _output(staged_report, report) as details,
        ):
            yield target, details
        if output != "-":
            staged_output.replace(output)
        if report != "-":
            staged_report.replace(report)


def _same_path(first, second):
    if first.resolve() == second.resolve():
        return True
    return first.exists() and second.exists() and first.samefile(second)


def _check_paths(args, report):
    if args.output == report == "-":
        raise ValueError("Le VCF et le rapport exigent des sorties distinctes.")
    sources = (
        args.input, args.context, args.contexts, args.reference,
        args.annotation_db, args.mane_db, args.ps3,
    )
    outputs = [Path(path) for path in (args.output, report) if path != "-"]
    inputs = [Path(path) for path in sources if path not in (None, "-")]
    checks = list(combinations(outputs, 2))
    checks += [(output, source) for output in outputs for source in inputs]
    for first, second in checks:
        if _same_path(first, second):
            raise ValueError(f"Les fichiers doivent être distincts ({first}).")
    for output in outputs:
        if output.is_dir():
            raise ValueError(f"Choisissez un fichier de sortie ({output}).")


def _check_stdin(*paths):
    if paths.count("-") > 1:
        raise ValueError("Une seule entrée peut utiliser stdin.")


def _vcf(args):
    _check_stdin(args.input, args.context, args.contexts)
    report = args.report
    if report is None:
        if args.output == "-":
            raise ValueError("Indiquez --report si le VCF va sur stdout.")
        report = args.output + ".jsonl"
    _check_paths(args, report)
    context = _context(args)
    contexts = _contexts(args.contexts)
    with (
        _classifier(args) as classifier,
        _input(args.input) as source,
        _outputs(args.output, report) as (output, details),
    ):
        batch = Batch(classifier, context, contexts)
        annotate_vcf(source, output, details, batch)


def _resources(parser):
    parser.add_argument("--annotation-db", help="Base locale d’annotations")
    parser.add_argument("--reference", help="FASTA GRCh38 indexé")
    parser.add_argument("--mane-db", help="Base locale des transcrits MANE")
    parser.add_argument("--vcep-dir", help="Répertoire local des profils VCEP")
    parser.add_argument("--ps3", help="Fichier des preuves fonctionnelles PS3")


def _common_context(parser):
    parser.add_argument("--context", help="Contexte clinique commun, en JSON")
    parser.add_argument("--inheritance", choices=sorted(INHERITANCES))


def _classify_parser(commands):
    command = commands.add_parser("classify", help="Évaluer un SNV")
    _resources(command)
    _common_context(command)
    command.add_argument("--variant", required=True, help="chr1:123:A:G")
    command.add_argument("--annotation", help="Fichier d’annotation ou -")
    command.add_argument("--overrides", help="Critères modifiés au format JSON")
    command.set_defaults(run=_classify)


def _vcf_parser(commands):
    command = commands.add_parser("vcf", help="Évaluer les ALT d’un flux VCF")
    _resources(command)
    _common_context(command)
    command.add_argument("input", help="VCF, VCF.gz ou - pour stdin")
    command.add_argument("-o", "--output", default="-", help="VCF, VCF.gz ou -")
    command.add_argument("--report", help="Rapport détaillé JSONL ou JSONL.gz")
    command.add_argument("--contexts", help="Contextes par variant, objet JSON")
    command.set_defaults(run=_vcf)


def _port(value):
    port = int(value)
    if not 0 <= port <= 65535:
        raise argparse.ArgumentTypeError("Port hors de l’intervalle 0 à 65535.")
    return port


def _web_parser(commands):
    command = commands.add_parser("web", help="Ouvrir l’interface locale")
    _resources(command)
    command.add_argument("--port", type=_port, default=8765, help="0 pour auto")
    command.add_argument("--no-browser", action="store_true")
    command.set_defaults(run=_web)


def _vcep_parser(commands):
    command = commands.add_parser("vcep", help="Consulter les VCEP locaux")
    command.add_argument("--gene", required=True, help="Symbole HGNC")
    command.add_argument("--disease", default="", help="Identifiant MONDO")
    command.add_argument("--inheritance", choices=sorted(INHERITANCES))
    command.add_argument("--vcep-dir", help="Répertoire local des profils VCEP")
    command.set_defaults(run=_vcep)


def _index_annotations(args):
    build_annotation_db(args.input, args.output)


def _index_mane(args):
    build_mane_db(args.input, args.output)


def _index_fasta(args):
    build_fasta_index(args.input)


def _index_parsers(commands):
    command = commands.add_parser("index", help="Indexer les annotations")
    command.add_argument("input", help="Annotations normalisées TSV ou JSONL")
    command.add_argument("output", help="Destination SQLite")
    command.set_defaults(run=_index_annotations)
    command = commands.add_parser("index-mane", help="Indexer les données MANE")
    command.add_argument("input", help="Transcrits MANE normalisés")
    command.add_argument("output", help="Destination SQLite")
    command.set_defaults(run=_index_mane)
    command = commands.add_parser("index-fasta", help="Créer un index FASTA")
    command.add_argument("input", help="FASTA GRCh38 sans compression")
    command.set_defaults(run=_index_fasta)


def _parser():
    parser = argparse.ArgumentParser(
        prog="acmg_snv",
        description="Aide à la classification des SNV GRCh38, ACGS 2026.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    _classify_parser(commands)
    _vcf_parser(commands)
    _web_parser(commands)
    _vcep_parser(commands)
    _index_parsers(commands)
    return parser


def main(argv=None):
    """Run one offline command and return its exit code."""
    args = _parser().parse_args(argv)
    try:
        args.run(args)
    except BrokenPipeError:
        return 1
    except KeyboardInterrupt:
        return 130
    except (ValueError, OSError, sqlite3.Error) as exc:
        print(f"Erreur ACMG — {exc}", file=sys.stderr)
        return 2
    return 0
