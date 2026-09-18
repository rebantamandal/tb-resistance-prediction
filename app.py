"""Local browser UI for amr.py, using only Python's standard-library web server.

Binds to 127.0.0.1. No cloud service, account, or external API is used.
This development UI is NOT a production or multi-user medical application.
"""
from __future__ import annotations

import argparse
import io
import json
import re
import secrets
import shutil
import threading
import traceback
import webbrowser
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from amr import Config, NOTICE, load_bundle, predict_rows, read_csv, train_models

ROOT = Path(__file__).resolve().parent
MAX_BODY = 32 * 1024 * 1024
TRAIN_LOCK = threading.Lock()


def _run_path(root: Path, run_id: str) -> Path:
    if not re.fullmatch(r"run_[0-9]{8}_[0-9]{6}_[a-f0-9]{8}", run_id):
        raise ValueError("Invalid run identifier.")
    result = root / run_id
    if result.is_symlink() or not result.is_dir():
        raise ValueError("Model run not found.")
    return result


def _runs(root: Path) -> list[dict[str, Any]]:
    result = []
    for path in sorted(root.glob("run_*"), reverse=True):
        try:
            valid = _run_path(root, path.name)
            report = json.loads((valid / "report.json").read_text(encoding="utf-8"))
            result.append({"run_id": path.name, "organism": report["config"]["organism"],
                           "antibiotic": report["config"]["antibiotic"],
                           "created_utc": report["created_utc"], "variants": list(report["models"])})
        except (ValueError, OSError, KeyError):
            continue
    return result


def create_server(port: int = 8765, model_root: Path | None = None) -> ThreadingHTTPServer:
    models = model_root or ROOT / "models"
    models.mkdir(parents=True, exist_ok=True)
    token = secrets.token_urlsafe(32)

    class Handler(BaseHTTPRequestHandler):
        server_version = "AMRLocal/1.0"

        def log_message(self, fmt: str, *args: Any) -> None:
            # Do not log request bodies, uploaded records, or model predictions.
            pass

        def _reply(self, status: int, body: bytes, kind: str = "application/json") -> None:
            self.send_response(status)
            self.send_header("Content-Type", kind + "; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("X-Frame-Options", "DENY")
            self.send_header("Referrer-Policy", "no-referrer")
            self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; connect-src 'self'; base-uri 'none'; frame-ancestors 'none'")
            self.end_headers()
            try:
                self.wfile.write(body)
            except (BrokenPipeError, ConnectionResetError):
                pass

        def _json(self, status: int, payload: Any) -> None:
            self._reply(status, json.dumps(payload, allow_nan=False).encode("utf-8"))

        def _allowed(self, api: bool) -> bool:
            actual_port = self.server.server_address[1]
            hosts = {f"127.0.0.1:{actual_port}", f"localhost:{actual_port}"}
            if self.headers.get("Host") not in hosts:
                self._json(403, {"error": "Only local browser access is allowed."})
                return False
            origin = self.headers.get("Origin")
            if origin and origin not in {"http://" + host for host in hosts}:
                self._json(403, {"error": "Cross-origin requests are not allowed."})
                return False
            if api and not secrets.compare_digest(self.headers.get("X-AMR-Token", ""), token):
                self._json(403, {"error": "Refresh the local app page to establish a valid session."})
                return False
            return True

        def do_GET(self) -> None:
            if not self._allowed(api=self.path.startswith("/api/")):
                return
            if self.path == "/":
                html = (ROOT / "web" / "index.html").read_text(encoding="utf-8")
                html = html.replace("__SESSION_TOKEN__", token)
                self._reply(200, html.encode("utf-8"), "text/html")
            elif self.path == "/api/models":
                self._json(200, {"runs": _runs(models)})
            else:
                self._json(404, {"error": "Not found."})

        def do_POST(self) -> None:
            if not self._allowed(api=True):
                return
            try:
                length = int(self.headers.get("Content-Length", "0"))
                if not 0 < length <= MAX_BODY:
                    self._json(413, {"error": "Request is empty or exceeds the 32 MiB local prototype limit. Use the CLI for larger matrices."})
                    return
                if self.headers.get("Content-Type", "").split(";")[0] != "application/json":
                    raise ValueError("Expected JSON request.")
                payload = json.loads(self.rfile.read(length))
                if not isinstance(payload, dict):
                    raise ValueError("Expected a JSON object.")
                if self.path == "/api/inspect":
                    frame = read_csv(io.StringIO(payload["csv"]))
                    if frame.empty:
                        raise ValueError("The CSV has no observations. Blank templates show the required structure; add real data first.")
                    columns = []
                    for col in frame:
                        distinct = frame[col].dropna().unique()
                        numeric = pd.to_numeric(frame[col].dropna(), errors="coerce")
                        columns.append({"name": col, "example": str(frame[col].dropna().iloc[0]) if frame[col].notna().any() else "",
                            "numeric": bool(len(numeric) and numeric.notna().all()), "unique_n": len(distinct),
                            "values": sorted(str(v) for v in distinct)[:500]})
                    self._json(200, {"row_count": len(frame), "columns": columns})
                elif self.path == "/api/train":
                    if not TRAIN_LOCK.acquire(blocking=False):
                        self._json(409, {"error": "Another training run is in progress."})
                        return
                    run_id = "run_" + datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S") + "_" + secrets.token_hex(4)
                    temporary = models / (".working_" + run_id)
                    try:
                        cfg = Config(**payload["config"])
                        report = train_models(read_csv(io.StringIO(payload["csv"])), cfg, temporary)
                        temporary.rename(models / run_id)
                        self._json(200, {"run_id": run_id, "report": report})
                    finally:
                        if temporary.exists():
                            shutil.rmtree(temporary)
                        TRAIN_LOCK.release()
                elif self.path == "/api/model":
                    run = _run_path(models, payload["run_id"])
                    report = json.loads((run / "report.json").read_text(encoding="utf-8"))
                    self._json(200, {"run_id": run.name, "report": report})
                elif self.path == "/api/predict":
                    run = _run_path(models, payload["run_id"])
                    variant = payload["variant"]
                    if variant not in {"genomic_only", "full"}:
                        raise ValueError("Unknown model variant.")
                    path = run / (variant + ".joblib")
                    if path.is_symlink():
                        raise ValueError("Linked model files are not supported.")
                    # Only local app-created runs are selectable; no model upload is provided.
                    if not payload.get("trust_local_model"):
                        raise ValueError("Confirm that this is a model created locally by you.")
                    bundle = load_bundle(path, trusted=True)
                    if "csv" in payload:
                        frame = read_csv(io.StringIO(payload["csv"]))
                    elif "records" in payload:
                        frame = pd.DataFrame(payload["records"]).replace("", np.nan)
                    else:
                        raise ValueError("Provide a prediction CSV or a single record.")
                    results = predict_rows(bundle, frame)
                    self._json(200, {"records": json.loads(results.head(1000).to_json(orient="records")),
                        "total_rows": len(results), "csv": results.to_csv(index=False), "notice": NOTICE})
                else:
                    self._json(404, {"error": "Not found."})
            except (ValueError, TypeError, OSError, KeyError) as exc:
                self._json(400, {"error": str(exc)})
            except Exception:
                traceback.print_exc()
                self._json(500, {"error": "Unexpected local error. Check the terminal for details; no records were sent to an external service."})

    server = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    server.daemon_threads = True
    return server


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()
    try:
        server = create_server(args.port)
    except OSError as exc:
        raise SystemExit(f"Could not start local app: {exc}. Try --port 8766.") from exc
    url = f"http://127.0.0.1:{server.server_address[1]}"
    print(f"AMR Prediction: {url}\n{NOTICE}\nPress Ctrl+C to stop.")
    if not args.no_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
