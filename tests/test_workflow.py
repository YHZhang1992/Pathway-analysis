import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class PortableWorkflowTest(unittest.TestCase):
    def test_example_run(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temporary:
            result = subprocess.run([sys.executable, str(root / "workflow.py"), "--input", str(root / "examples" / "input.csv"), "--config", str(root / "config" / "workflow.json"), "--output", temporary], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            output = Path(temporary) / "results.csv"
            self.assertTrue(output.is_file())
            with output.open(encoding="utf-8", newline="") as handle:
                self.assertGreater(len(list(csv.DictReader(handle))), 0)
            manifest = json.loads((Path(temporary) / "run_manifest.json").read_text(encoding="utf-8"))
            with output.open(encoding="utf-8", newline="") as handle:
                output_rows = list(csv.DictReader(handle))
            self.assertEqual(manifest["row_count"], len(output_rows))


if __name__ == "__main__":
    unittest.main()
