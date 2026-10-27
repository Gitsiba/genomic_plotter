import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pyBigWig

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "plot_bigwig_metagene.py"
SPEC = importlib.util.spec_from_file_location("plot_bigwig_metagene", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class BigWigMetageneInputTests(unittest.TestCase):
    def test_reads_valid_bed6(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "genes.bed"
            bed.write_text("chr1\t10\t110\tGENE1\t0\t+\nchr2\t20\t90\tGENE2\t0\t-\n")
            genes = MODULE.read_bed6(bed)
            self.assertEqual(len(genes), 2)
            self.assertEqual(genes[0], ("chr1", 10, 110, "GENE1", "+"))
            self.assertEqual(genes[1][4], "-")

    def test_rejects_invalid_bed_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "bad.bed"
            bed.write_text("chr1\t100\t10\tGENE1\t0\t+\n")
            with self.assertRaises(ValueError):
                MODULE.read_bed6(bed)

    def test_rejects_missing_bed_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "short.bed"
            bed.write_text("chr1\t10\t100\n")
            with self.assertRaises(ValueError):
                MODULE.read_bed6(bed)

    def test_bigwig_values_and_minus_strand_orientation(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            bw_path = tmp / "fixture.bw"
            bw = pyBigWig.open(str(bw_path), "w")
            bw.addHeader([("chr1", 20)])
            starts = list(range(20))
            bw.addEntries(["chr1"] * 20, starts, ends=[x + 1 for x in starts],
                          values=[float(x) for x in starts])
            bw.close()

            bed = tmp / "genes.bed"
            bed.write_text(
                "chr1\t0\t20\tPLUS\t0\t+\n"
                "chr1\t0\t20\tMINUS\t0\t-\n",
                encoding="utf-8",
            )
            matrix, labels = MODULE.summarize_bigwig(bw_path, bed, bins=2, min_coverage=1.0)
            rows = {name: row for name, row in zip(labels, matrix)}
            self.assertIn("PLUS", rows)
            self.assertIn("MINUS", rows)
            np.testing.assert_allclose(rows["PLUS"], [4.5, 14.5])
            np.testing.assert_allclose(rows["MINUS"], [14.5, 4.5])

    def test_end_to_end_cli_creates_plot_and_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            bw_path = tmp / "fixture.bw"
            bw = pyBigWig.open(str(bw_path), "w")
            bw.addHeader([("chr1", 20)])
            starts = list(range(20))
            bw.addEntries(["chr1"] * 20, starts, ends=[x + 1 for x in starts],
                          values=[float(x + 1) for x in starts])
            bw.close()
            bed = tmp / "genes.bed"
            bed.write_text("chr1\t0\t20\tGENE1\t0\t+\n", encoding="utf-8")
            out = tmp / "plots" / "metagene.png"
            matrix_out = tmp / "plots" / "matrix.tsv"
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "--bigwig", str(bw_path),
                 "--genes-bed", str(bed), "--bins", "4", "--min-coverage", "1.0",
                 "--out", str(out), "--matrix-out", str(matrix_out)],
                capture_output=True, text=True, timeout=30,
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertTrue(out.is_file())
            self.assertGreater(out.stat().st_size, 1000)
            self.assertTrue(matrix_out.is_file())
            self.assertIn("GENE1", matrix_out.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
