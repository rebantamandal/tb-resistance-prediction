"""Local HTTP interface and request-isolation tests; no external services used."""
import json
import re
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from pathlib import Path

from app import create_server


class LocalAppTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.server = create_server(0, Path(cls.temp.name))
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.base = f"http://127.0.0.1:{cls.server.server_address[1]}"
        cls.html = urllib.request.urlopen(cls.base + "/").read().decode()
        cls.token = re.search(r"const TOKEN='([^']+)';", cls.html).group(1)

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join()
        cls.temp.cleanup()

    def request(self, path, payload=None, token=True, origin=None):
        headers = {"Content-Type": "application/json"}
        if token:
            headers["X-AMR-Token"] = self.token
        if origin:
            headers["Origin"] = origin
        data = json.dumps(payload).encode() if payload is not None else None
        req = urllib.request.Request(self.base + path, data=data, headers=headers)
        return urllib.request.urlopen(req)

    def test_page_renders_with_research_notice(self):
        self.assertIn("No pretrained model", self.html)
        self.assertNotIn("__SESSION_TOKEN__", self.html)

    def test_empty_model_list(self):
        result = json.load(self.request("/api/models"))
        self.assertEqual(result["runs"], [])

    def test_inspect_preserves_column_values(self):
        result = json.load(self.request("/api/inspect", {"csv": "isolate_id,value\n00001,1\n"}))
        self.assertEqual(result["columns"][0]["example"], "00001")

    def test_missing_token_rejected(self):
        with self.assertRaises(urllib.error.HTTPError) as err:
            self.request("/api/models", token=False)
        self.assertEqual(err.exception.code, 403)

    def test_cross_origin_rejected(self):
        with self.assertRaises(urllib.error.HTTPError) as err:
            self.request("/api/models", origin="https://untrusted.example")
        self.assertEqual(err.exception.code, 403)

    def test_header_only_csv_is_not_training_data(self):
        with self.assertRaises(urllib.error.HTTPError) as err:
            self.request("/api/inspect", {"csv": "a,b\n"})
        self.assertEqual(err.exception.code, 400)

    def test_model_path_traversal_rejected(self):
        with self.assertRaises(urllib.error.HTTPError) as err:
            self.request("/api/model", {"run_id": "../../etc"})
        self.assertEqual(err.exception.code, 400)


if __name__ == "__main__":
    unittest.main()
