"""A local browser interface for the resistance predictor.

Run it, a page opens, you drop a CSV in, you get a table back. No command line
and no options to get wrong. It is the same code path as
`predict_resistance.py` — this only wraps it in a page.

    python predict_app.py

The server binds to localhost only and is never reachable from the network. It
uses the Python standard library, so there is nothing extra to install.

RESEARCH USE ONLY. The models are not clinically validated and the scores are
not probabilities. Nothing here should inform a treatment decision.
"""
from __future__ import annotations

import argparse
import html
import io
import json
import secrets
import sys
import tempfile
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pandas as pd

import predict_resistance as pr

MAX_UPLOAD_BYTES = 64 * 1024 * 1024      # a 100-patient file is ~15 MB
UNCERTAIN_LOW, UNCERTAIN_HIGH = 0.40, 0.60

PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TB Resistance Predictor</title>
<style>
  :root {
    color-scheme: light dark;
    --bg: #fbfbf8; --surface: #ffffff; --line: #e1e5ea;
    --ink: #14213d; --ink-soft: #4a5568; --ink-faint: #6a7179;
    --accent: #2a78d6; --danger: #c1442a; --danger-bg: #fdf0ec;
    --ok: #17724f; --ok-bg: #edf7f2; --warn: #8a5a00; --warn-bg: #fdf6e6;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      --bg: #14161a; --surface: #1c1f25; --line: #2e333c;
      --ink: #f2f4f8; --ink-soft: #b8c0cc; --ink-faint: #8a93a1;
      --accent: #5a9ae8; --danger: #ef8b72; --danger-bg: #2c1d19;
      --ok: #62c79b; --ok-bg: #16281f; --warn: #e0b355; --warn-bg: #2a2317;
    }
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; background: var(--bg); color: var(--ink);
    font: 16px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif;
  }
  .wrap { max-width: 1100px; margin: 0 auto; padding: 40px 16px 80px; }
  h1 { font-size: 30px; margin: 0 0 6px; letter-spacing: -0.2px; }
  .sub { color: var(--ink-soft); margin: 0 0 28px; }
  .notice {
    background: var(--warn-bg); border: 1px solid var(--line);
    border-left: 4px solid var(--warn); border-radius: 10px;
    padding: 14px 18px; margin-bottom: 28px; color: var(--ink);
  }
  .notice b { color: var(--warn); }
  .card {
    background: var(--surface); border: 1px solid var(--line);
    border-radius: 14px; padding: 28px; margin-bottom: 24px;
  }
  #drop {
    border: 2px dashed var(--line); border-radius: 12px; padding: 44px 20px;
    text-align: center; cursor: pointer; transition: border-color .15s, background .15s;
  }
  #drop.over, #drop:hover { border-color: var(--accent); background: var(--bg); }
  #drop p { margin: 6px 0; }
  .hint { color: var(--ink-faint); font-size: 14px; }
  button {
    background: var(--accent); color: #fff; border: 0; border-radius: 9px;
    padding: 11px 20px; font-size: 15px; font-weight: 600; cursor: pointer;
    font-family: inherit;
  }
  button:disabled { opacity: .5; cursor: default; }
  button.secondary { background: transparent; color: var(--accent); border: 1px solid var(--line); }
  table { border-collapse: collapse; width: 100%; font-size: 14px; }
  th, td { padding: 9px 11px; text-align: left; border-bottom: 1px solid var(--line); }
  th { font-weight: 600; color: var(--ink-soft); font-size: 13px;
       text-transform: uppercase; letter-spacing: .4px; }
  td.id { font-weight: 600; white-space: nowrap; }
  .pill {
    display: inline-block; padding: 2px 9px; border-radius: 999px;
    font-size: 13px; font-weight: 600; white-space: nowrap;
  }
  .r { background: var(--danger-bg); color: var(--danger); }
  .s { background: var(--ok-bg); color: var(--ok); }
  .u { background: var(--warn-bg); color: var(--warn); }
  .score { color: var(--ink-faint); font-variant-numeric: tabular-nums; font-size: 13px; }
  .flag { background: var(--warn-bg); border-left: 4px solid var(--warn);
          padding: 12px 16px; border-radius: 8px; margin-top: 18px; font-size: 14px; }
  .err { background: var(--danger-bg); border-left: 4px solid var(--danger);
         padding: 14px 18px; border-radius: 8px; color: var(--danger); }
  .legend { display: flex; gap: 22px; flex-wrap: wrap; margin: 18px 0 0;
            font-size: 14px; color: var(--ink-soft); }
  code { background: var(--bg); border: 1px solid var(--line); border-radius: 5px;
         padding: 1px 6px; font-size: 13px; }
  .row { display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
  .scroll { overflow-x: auto; }
  a { color: var(--accent); }
</style>
</head>
<body>
<div class="wrap">
  <h1>TB Resistance Predictor</h1>
  <p class="sub">Upload a variant CSV. Get a resistance prediction for every sample, for every drug.</p>

  <div class="notice">
    <b>Research use only.</b> These models are not medically approved and the
    scores are not probabilities. Never use this output to decide anyone's
    treatment.
  </div>

  <div class="card">
    <div id="drop">
      <p><b>Drop a CSV here</b>, or click to choose a file</p>
      <p class="hint">Columns: SAMPLE, CHROM, POS, REF, ALT &middot; positions must use the NC_000962.3 reference</p>
      <p class="hint" id="chosen"></p>
    </div>
    <div class="row" style="margin-top:18px">
      <button id="go" disabled>Predict</button>
      <button class="secondary" id="demo">Use the example file</button>
      <span class="hint" id="status"></span>
    </div>
  </div>

  <div id="out"></div>
</div>

<script>
const drop = document.getElementById('drop');
const go = document.getElementById('go');
const status = document.getElementById('status');
const chosen = document.getElementById('chosen');
const out = document.getElementById('out');
const picker = document.createElement('input');
picker.type = 'file'; picker.accept = '.csv,text/csv';
let file = null;

function setFile(f) {
  file = f;
  chosen.textContent = f ? f.name + '  (' + (f.size/1048576).toFixed(1) + ' MB)' : '';
  go.disabled = !f;
}
drop.onclick = () => picker.click();
picker.onchange = () => setFile(picker.files[0]);
drop.ondragover = e => { e.preventDefault(); drop.classList.add('over'); };
drop.ondragleave = () => drop.classList.remove('over');
drop.ondrop = e => {
  e.preventDefault(); drop.classList.remove('over');
  if (e.dataTransfer.files.length) setFile(e.dataTransfer.files[0]);
};

function busy(on, msg) {
  go.disabled = on || !file;
  document.getElementById('demo').disabled = on;
  status.textContent = msg || '';
}

async function send(body, headers) {
  busy(true, 'Reading the file and scoring every sample...');
  out.innerHTML = '';
  try {
    const res = await fetch('/predict', {method: 'POST', headers, body});
    const data = await res.json();
    if (data.error) {
      out.innerHTML = '<div class="card err">' + data.error + '</div>';
    } else {
      render(data);
    }
  } catch (err) {
    out.innerHTML = '<div class="card err">Could not reach the local server: ' + err + '</div>';
  }
  busy(false, '');
}

go.onclick = () => send(file, {'Content-Type': 'text/csv'});
document.getElementById('demo').onclick = () => send('', {'X-Use-Example': '1'});

function render(d) {
  let h = '<div class="card">';
  h += '<div class="row" style="justify-content:space-between">';
  h += '<h2 style="margin:0;font-size:20px">' + d.rows.length + ' sample' +
       (d.rows.length === 1 ? '' : 's') + ' scored</h2>';
  h += '<a href="' + d.download + '" download="predictions.csv"><button class="secondary">Download CSV</button></a>';
  h += '</div>';
  h += '<div class="legend"><span><span class="pill r">Resistant</span> drug will probably fail</span>' +
       '<span><span class="pill s">Susceptible</span> drug will probably work</span>' +
       '<span><span class="pill u">Uncertain</span> send for lab testing</span></div>';
  h += '<div class="scroll" style="margin-top:18px"><table><thead><tr><th>Sample</th><th>Recognised changes</th>';
  d.drugs.forEach(x => h += '<th>' + x + '</th>');
  h += '</tr></thead><tbody>';
  d.rows.forEach(r => {
    h += '<tr><td class="id">' + r.id + '</td><td class="score">' + r.known +
         ' of ' + r.supplied + '</td>';
    d.drugs.forEach(x => {
      const c = r.calls[x];
      h += '<td><span class="pill ' + c.cls + '">' + c.label + '</span> ' +
           '<span class="score">' + c.score.toFixed(2) + '</span></td>';
    });
    h += '</tr>';
  });
  h += '</tbody></table></div>';
  if (d.flagged.length) {
    h += '<div class="flag"><b>' + d.flagged.length + ' sample' +
         (d.flagged.length === 1 ? '' : 's') +
         ' had no recognisable DNA changes:</b> ' + d.flagged.join(', ') +
         '. Do not read these as Susceptible &mdash; it usually means the positions ' +
         'were measured against a different reference than <code>NC_000962.3</code>.</div>';
  }
  h += '</div>';
  h += '<div class="card"><h2 style="margin-top:0;font-size:18px">How to read this</h2>' +
       '<p>The number beside each result is confidence, not a probability. ' +
       '<b>0.80 does not mean an 80% chance.</b> Use it as a ranking: higher means ' +
       'more confident. Anything between 0.40 and 0.60 is the model saying it is ' +
       'genuinely unsure.</p>' +
       '<p class="hint">Expected accuracy on samples from a collection the model ' +
       'never trained on: ' + d.accuracy + '</p></div>';
  out.innerHTML = h;
}
</script>
</body>
</html>
"""


class Predictor:
    """Loads every registered model once, then scores uploads against them."""

    def __init__(self, dictionary_path: Path) -> None:
        self.models: dict[str, dict] = {}
        for drug, entry in pr.DEFAULT_REGISTRY.items():
            model_path = Path(entry["run"]) / "genomic_only.joblib"
            if not model_path.is_file():
                continue
            bundle = pr.load_bundle(model_path)
            self.models[drug] = {
                "bundle": bundle,
                "columns": pr.bundle_feature_columns(bundle, model_path),
                "study_auc": entry.get("auc_study_grouped"),
            }
        if not self.models:
            raise SystemExit("No trained models found. Run from the project folder.")

        reference = next(iter(self.models))
        self.columns = self.models[reference]["columns"]
        for drug, info in self.models.items():
            if info["columns"] != self.columns:
                raise SystemExit(f"{drug} was trained on a different feature set. "
                                 "Retrain every model from one prepare_variants.py run.")

        if not dictionary_path.is_file():
            raise SystemExit(f"Feature dictionary not found: {dictionary_path}")
        dictionary = pd.read_csv(dictionary_path, dtype=str)
        self.key_of_column = pr.column_to_key(self.columns, dictionary)

    def accuracy_sentence(self) -> str:
        parts = [f"{d} {i['study_auc']:.2f}" for d, i in sorted(self.models.items())
                 if i["study_auc"] is not None]
        return " · ".join(parts) if parts else "see RESULTS.md"

    def score(self, csv_bytes: bytes) -> dict:
        # build_matrix streams from a path, so the upload goes to a temp file
        # that is removed as soon as the matrix is built.
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as handle:
            handle.write(csv_bytes)
            temp = Path(handle.name)
        try:
            matrix, supplied = pr.build_matrix(temp, self.key_of_column, progress=False)
        except SystemExit as error:
            # build_matrix names the file it was given, which here is a
            # meaningless temporary path. Rewrite it for the person reading.
            message = str(error).replace(str(temp), "Your file")
            raise SystemExit(message) from None
        finally:
            # On Windows the streaming reader can still hold the file when the
            # parse fails. Deleting is best-effort so the real error survives
            # instead of being replaced by a PermissionError.
            try:
                temp.unlink(missing_ok=True)
            except OSError:
                pass

        features = matrix[self.columns]
        known = features.sum(axis=1)

        result = pd.DataFrame({
            "isolate_id": matrix["isolate_id"],
            "variant_calls_supplied": supplied.to_numpy(),
            "known_features_present": known.to_numpy(),
        })
        rows = [{"id": str(i), "supplied": int(s), "known": int(k), "calls": {}}
                for i, s, k in zip(result["isolate_id"],
                                   result["variant_calls_supplied"],
                                   result["known_features_present"])]

        for drug in sorted(self.models):
            scores = self.models[drug]["bundle"]["pipeline"].predict_proba(features)[:, 1]
            result[f"{drug}_score"] = scores.round(4)
            result[f"{drug}_call"] = ["Resistant" if s >= 0.5 else "Susceptible"
                                      for s in scores]
            for row, value in zip(rows, scores):
                if UNCERTAIN_LOW <= value <= UNCERTAIN_HIGH:
                    label, cls = "Uncertain", "u"
                elif value >= 0.5:
                    label, cls = "Resistant", "r"
                else:
                    label, cls = "Susceptible", "s"
                row["calls"][drug] = {"label": label, "cls": cls, "score": float(value)}

        flagged = [r["id"] for r in rows if r["known"] == 0]
        result["warning"] = ["No known resistance feature detected; not evidence of "
                             "susceptibility." if k == 0 else ""
                             for k in result["known_features_present"]]

        buffer = io.StringIO()
        result.to_csv(buffer, index=False)
        return {
            "drugs": sorted(self.models),
            "rows": rows,
            "flagged": flagged,
            "accuracy": self.accuracy_sentence(),
            "csv": buffer.getvalue(),
        }


def make_handler(predictor: Predictor, token: str, example: Path):
    class Handler(BaseHTTPRequestHandler):
        last_csv = ""

        def _send(self, code: int, body: bytes, content_type: str,
                  extra: dict[str, str] | None = None) -> None:
            self.send_response(code)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("X-Content-Type-Options", "nosniff")
            for key, value in (extra or {}).items():
                self.send_header(key, value)
            self.end_headers()
            self.wfile.write(body)

        def _json(self, payload: dict, code: int = 200) -> None:
            self._send(code, json.dumps(payload).encode("utf-8"), "application/json")

        def _local(self) -> bool:
            # Loopback only. The server never binds elsewhere, but check anyway.
            return self.client_address[0] in {"127.0.0.1", "::1"}

        def do_GET(self) -> None:  # noqa: N802
            if not self._local():
                self._send(403, b"local only", "text/plain")
                return
            if self.path in ("/", "/index.html"):
                self._send(200, PAGE.encode("utf-8"), "text/html; charset=utf-8")
            elif self.path == "/predictions.csv":
                self._send(200, Handler.last_csv.encode("utf-8"),
                           "text/csv; charset=utf-8",
                           {"Content-Disposition": 'attachment; filename="predictions.csv"'})
            else:
                self._send(404, b"not found", "text/plain")

        def do_POST(self) -> None:  # noqa: N802
            if not self._local() or self.path != "/predict":
                self._send(403, b"local only", "text/plain")
                return
            try:
                length = int(self.headers.get("Content-Length") or 0)
            except ValueError:
                self._json({"error": "Malformed request."}, 400)
                return
            if length > MAX_UPLOAD_BYTES:
                self._json({"error": f"That file is larger than "
                                     f"{MAX_UPLOAD_BYTES // (1024*1024)} MB. Use the "
                                     f"command line for very large batches."}, 413)
                return

            if self.headers.get("X-Use-Example"):
                if not example.is_file():
                    self._json({"error": "The example file is not in this folder."}, 404)
                    return
                payload = example.read_bytes()
            else:
                payload = self.rfile.read(length) if length else b""
            if not payload.strip():
                self._json({"error": "No file received."}, 400)
                return

            try:
                outcome = predictor.score(payload)
            except SystemExit as error:
                self._json({"error": html.escape(str(error))}, 400)
                return
            except Exception as error:                     # noqa: BLE001
                self._json({"error": "Could not read that file. It needs the columns "
                                     "SAMPLE, CHROM, POS, REF, ALT. "
                                     f"({html.escape(type(error).__name__)})"}, 400)
                return

            Handler.last_csv = outcome.pop("csv")
            outcome["download"] = "/predictions.csv"
            self._json(outcome)

        def log_message(self, *args) -> None:              # keep the console quiet
            return

    return Handler


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Browser interface for the predictor.")
    parser.add_argument("--port", type=int, default=8777)
    parser.add_argument("--feature-dictionary", type=Path,
                        default=Path("data/variant_panel/feature_dictionary.csv"))
    parser.add_argument("--example", type=Path,
                        default=Path("examples/sample_isolates.csv"))
    parser.add_argument("--no-browser", action="store_true",
                        help="Don't open a browser window automatically.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    print(pr.NOTICE)
    print("Loading models...")
    predictor = Predictor(args.feature_dictionary)
    print(f"Ready: {', '.join(sorted(predictor.models))}")

    token = secrets.token_urlsafe(16)
    server = ThreadingHTTPServer(("127.0.0.1", args.port),
                                 make_handler(predictor, token, args.example))
    url = f"http://127.0.0.1:{args.port}/"
    print(f"\n  Open {url}\n  Press Ctrl+C to stop.\n")
    if not args.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("Stopped.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
