import importlib.util
import tempfile
import unittest
from pathlib import Path

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


if __name__ == "__main__":
    unittest.main()
