"""Regression tests for matching absolute weights and physical grids across the two formats."""
import csv
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1] / "chatgpt-python-unfold"
sys.path.insert(0, str(PYTHON_ROOT))

from dftb_unfold import COLUMNS
from spectrum import spectral_grid
from tools.compare_backends import compare_maps, read_map


class MapChecks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.prefix = Path(self.tmp.name) / "unfold"
        self.csv_file = self.prefix.with_suffix(".csv")
        self.q = np.array([[0., 0., 0.], [.2, 0., 0.]])
        with self.csv_file.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLUMNS)
            w.writeheader()
            for i in range(2):
                w.writerow(dict(point=i+1, distance_inv_ang=i*.25, spin=1, band=1, energy_ev=.8,
                                energy_rel_ev=.3, weight_full=.2, weight_refold=.75))
        for suffix, col in (("unfoldedBands", "weight_full"), ("refoldedBands", "weight_refold")):
            _, _, v = spectral_grid(self.csv_file, col, 1, -1., 1., .25, .25)
            self.write(suffix, v.ravel()*.25)

    def write(self, suffix, vals):
        rows = ["2 9 -1.000 1.000 0.500000"]
        for i, q in enumerate(self.q):
            rows += [" ".join(map(str, q)) + " 1", *map(str, vals[9*i:9*i+9])]
        self.prefix.with_suffix("." + suffix).write_text("\n".join(rows) + "\n")

    def compare(self):
        return compare_maps(self.prefix, self.csv_file, self.q, .5, 1, -1., 1., 8, .25, .003, .002)

    def test_bin_width_and_absolute_normalization(self):
        self.assertTrue(all(v["passed"] for v in self.compare().values()))
        q, grid, ef, bins = read_map(self.prefix.with_suffix(".refoldedBands"))
        self.write("refoldedBands", bins.ravel()*2)
        report = self.compare()
        self.assertFalse(report["weight_refold_spin1"]["passed"])
        self.assertAlmostEqual(report["weight_refold_spin1"]["relative_L1"], 1.)

    def test_zero_width_linear_bins(self):
        _, grid, v = spectral_grid(self.csv_file, "weight_full", 1, -1., 1., .25, 0.)
        # E=0.3 splits the weight .2 in the proportion 0.8:0.2 between .25 and .5 eV.
        self.assertAlmostEqual(v[0, np.where(grid == .25)[0][0]], .2*.8/.25)
        self.assertAlmostEqual(v[0, np.where(grid == .5)[0][0]], .2*.2/.25)
        self.assertAlmostEqual(v[0].sum()*.25, .2)

    def test_incompatible_q_and_incomplete_map(self):
        self.q[0, 0] = .9
        with self.assertRaisesRegex(ValueError, "path"):
            self.compare()
        path = self.prefix.with_suffix(".unfoldedBands")
        path.write_text("2 9 -1 1 .5\n0 0 0 1\n0.0\n")
        with self.assertRaisesRegex(ValueError, "incomplete"):
            read_map(path)
