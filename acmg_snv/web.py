"""Interface web locale, servie avec la bibliothèque standard de Python."""

import json
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Lock

from .annotations import parse_annotation
from .functional import load_ps3
from .model import Variant

ASSETS = Path(__file__).with_name("web_assets")
ROUTES = {
    "/": ("index.html", "text/html"),
    "/style.css": ("style.css", "text/css"),
    "/app.js": ("app.js", "text/javascript"),
}
MAX_REQUEST_BYTES = 32 * 1024 * 1024
CONTENT_POLICY = (
    "default-src 'none'; script-src 'self'; style-src 'self'; "
    "connect-src 'self'; base-uri 'none'; form-action 'none'; "
    "frame-ancestors 'none'"
)


def _read_input(payload, key, expected_type, default):
    value = payload.get(key, default)
    if not isinstance(value, expected_type):
        raise ValueError(f"Type de donnée invalide pour {key}.")
    return value


def _parse_variant(text):
    value = text.strip().lstrip("\ufeff")
    fields = value.split("\t")
    if len(fields) >= 8:
        if len(value.splitlines()) != 1:
            raise ValueError("Collez une seule ligne VCF sans en-tête.")
        chrom, pos, _, ref, alt = fields[:5]
        return Variant(chrom, int(pos), ref, alt)
    return Variant.parse(value)


def _read_ps3(payload):
    ps3_text = _read_input(payload, "ps3", str, "")
    if not ps3_text.strip():
        return None
    rows = load_ps3(text=ps3_text)
    filename = _read_input(payload, "ps3_name", str, "")
    if filename:
        for entry in rows.values():
            if entry["source"] == "fichier PS3 importé":
                entry["source"] = filename
    return rows


def _classify(classifier, payload):
    if not isinstance(payload, dict):
        raise ValueError("La saisie doit être un objet JSON.")
    variant_text = _read_input(payload, "variant", str, "")
    context = _read_input(payload, "context", dict, {})
    annotation_text = _read_input(payload, "annotation", str, "")
    overrides = _read_input(payload, "overrides", list, [])
    variant = _parse_variant(variant_text)
    annotation = None
    if annotation_text.strip():
        annotation = parse_annotation(annotation_text, variant)
    ps3 = _read_ps3(payload)
    return classifier.classify(variant, context, annotation, overrides, ps3)


class LocalServer(ThreadingHTTPServer):
    """Lie le service au poste local et fournit le classificateur."""

    def __init__(self, port, classifier):
        super().__init__(("127.0.0.1", port), LocalHandler)
        self.classifier = classifier
        self.classifier_lock = Lock()
        self.authority = f"127.0.0.1:{self.server_port}"
        self.url = f"http://{self.authority}"


class LocalHandler(BaseHTTPRequestHandler):
    """Sert les ressources statiques et l'API de classification."""

    def _send(self, status, content_type, body):
        self.send_response(status)
        self.send_header("Content-Type", f"{content_type}; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Security-Policy", CONTENT_POLICY)
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, status, data):
        body = json.dumps(data, ensure_ascii=False).encode("utf-8")
        self._send(status, "application/json", body)

    def _local_request(self):
        host = self.headers.get("Host")
        origin = self.headers.get("Origin")
        if host == self.server.authority and origin in {None, self.server.url}:
            return True
        self._json(HTTPStatus.FORBIDDEN, {"error": "Accès local requis."})
        return False

    def do_GET(self):
        if not self._local_request():
            return
        resource = ROUTES.get(self.path)
        if resource is None:
            self._json(HTTPStatus.NOT_FOUND, {"error": "Page introuvable."})
            return
        filename, content_type = resource
        body = (ASSETS / filename).read_bytes()
        self._send(HTTPStatus.OK, content_type, body)

    def _payload(self):
        if self.headers.get_content_type() != "application/json":
            raise ValueError("Le contenu de la requête doit être en JSON.")
        length = int(self.headers.get("Content-Length", "0"))
        if not 0 < length <= MAX_REQUEST_BYTES:
            raise ValueError("La saisie est vide ou dépasse 32 Mio.")
        return json.loads(self.rfile.read(length))

    def do_POST(self):
        if not self._local_request():
            return
        if self.path != "/classify":
            self._json(HTTPStatus.NOT_FOUND, {"error": "Page introuvable."})
            return
        try:
            payload = self._payload()
            with self.server.classifier_lock:
                result = _classify(self.server.classifier, payload)
        except ValueError as error:
            self._json(HTTPStatus.BAD_REQUEST, {"error": str(error)})
            return
        self._json(HTTPStatus.OK, result)

    def log_message(self, *_args):
        """Évite d'écrire les requêtes HTTP dans la console."""


def serve(classifier, port=8765, open_browser=True):
    """Démarre l'interface locale et attend Ctrl+C."""
    with LocalServer(port, classifier) as server:
        print(f"Interface locale : {server.url}", flush=True)
        print("Appuyez sur Ctrl+C pour arrêter le service.", flush=True)
        if open_browser:
            webbrowser.open(server.url)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
