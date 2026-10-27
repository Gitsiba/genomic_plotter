import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


class PlotterSmokeTests(unittest.TestCase):
    def test_all_demo_commands_create_pngs(self):
        with tempfile.TemporaryDirectory() as tmp:
            for script, name in [
                ("plot_metagene.py", "metagene.png"),
                ("plot_metaorf.py", "metaorf.png"),
                ("plot_circos_data.py", "circos.png"),
                ("plot_mutation_impact.py", "impact.png"),
            ]:
                out = Path(tmp) / name
                result = subprocess.run(
                    [sys.executable, str(SCRIPTS / script), "--demo", "--out", str(out)],
                    capture_output=True, text=True, timeout=60,
                )
                self.assertEqual(result.returncode, 0, msg=result.stderr)
                self.assertTrue(out.is_file(), msg=f"{script} did not create output")
                self.assertGreater(out.stat().st_size, 1000)

    def test_metagene_accepts_numeric_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            matrix = Path(tmp) / "matrix.tsv"
            matrix.write_text("1\t2\t3\n2\t3\t4\n", encoding="utf-8")
            out = Path(tmp) / "nested" / "plot.png"
            result = subprocess.run(
                [sys.executable, str(SCRIPTS / "plot_metagene.py"),
                 "--input", str(matrix), "--out", str(out)],
                capture_output=True, text=True, timeout=60,
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertTrue(out.is_file())

    def test_real_input_required_unless_demo_explicit(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPTS / "plot_mutation_impact.py")],
            capture_output=True, text=True, timeout=30,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--input", result.stderr)


if __name__ == "__main__":
    unittest.main()
