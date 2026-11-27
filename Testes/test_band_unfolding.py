"""Regression tests for signs, non-orthogonal orbitals and reciprocal-space bookkeeping."""
import contextlib
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1] / "chatgpt-python-unfold"
sys.path.insert(0, str(PYTHON_ROOT))

from dftb_unfold import load_matrices, run, sampling, weights
from spectrum import spectral_grid
from unfold_fourier import (FourierBasis, bessel, eigenstates, group_momenta,
                            physical_k, reciprocal, refolding_vectors, tesseral)
from unfold_io import BOHR_ANG, HARTREE_EV, Geometry, RealMatrix, Shell, read_basis


CHAIN = PYTHON_ROOT / "examples" / "chain"


class UnfoldTests(unittest.TestCase):
    def test_native_xml_central_images_and_lattice(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"detailed.xml"
            xml = """<detailedout><geometry><typenames>"C" "B"</typenames>
<typesandcoordinates>2 1 2 3
1 5 2 3</typesandcoordinates><periodic>Yes</periodic><helical>No</helical>
<latticevectors>4 0 0 1 3 0 0 0 8</latticevectors>
<coordinateorigin>1 2 3</coordinateorigin></geometry></detailedout>"""
            path.write_text(xml)
            geom = Geometry.read(path)
            self.assertEqual(geom.species, ["B", "C"])
            np.testing.assert_array_equal(geom.positions, [[0,0,0], [4,0,0]])
            np.testing.assert_array_equal(geom.cell, [[4,0,0], [1,3,0], [0,0,8]])
            path.write_text(xml.replace("1 5 2 3", "1 6 2 3"))
            with self.assertRaisesRegex(ValueError, "ShowFoldedCoords"):
                Geometry.read(path)

    def test_pristine_chain_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            csv = Path(tmp)/"chain.csv"
            with contextlib.redirect_stderr(io.StringIO()):
                run(CHAIN/"unfold.json", csv)
            data = np.genfromtxt(csv, delimiter=",", names=True)
            self.assertEqual(len(data), 162)
            cosine = np.cos(2*np.pi*data["q1"])
            overlap = 0.07097404418511578  # independent two-centre integral in the example
            expected = (-.30-.16*cosine)/(1+2*overlap*cosine)*HARTREE_EV
            allowed = abs(data["energy_ev"]-expected) < 1e-9
            self.assertLess(max(data["weight_refold"][~allowed]), 1e-20)
            for point in set(data["point"]):
                visible = (data["point"] == point) & allowed
                self.assertGreater(data["weight_refold"][visible].sum(), .998)
                self.assertLess(data["weight_refold"][visible].sum(), 1.001)
            # Absolute weight and PDF/DOS integration, independent of path spacing.
            _, egrid, density = spectral_grid(csv, "weight_refold", 1, -13, -3, .01, .08)
            area = density.sum(axis=1)*(egrid[1]-egrid[0])
            np.testing.assert_allclose(area, [sum(data["weight_refold"][data["point"] == i])
                                              for i in range(1, 82)], atol=1e-9)

    def test_defect_and_state_normalization(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)/"case"
            shutil.copytree(CHAIN, target)
            hpath = target/"hamreal1.dat"
            text = hpath.read_text().splitlines()
            text[-1] = "-.18"
            hpath.write_text("\n".join(text)+"\n")
            hs, s = load_matrices(target)
            geom = Geometry.read(target/"chain.gen")
            shells = read_basis([target/"wfc.X.hsd"], geom.species, s.norb)
            fk = np.array([.18, 0, 0])
            k = physical_k(fk @ reciprocal(geom.cell/np.array([[2],[1],[1]])), geom.cell)
            h, sm = hs[0].bloch(k), s.bloch(k)
            e, c = eigenstates(h, sm)
            np.testing.assert_allclose(h@c, (sm@c)*e, atol=5e-15)
            np.testing.assert_allclose(c.conj().T@sm@c, np.eye(2), atol=5e-15)
            prim = np.linalg.solve(np.diag([2,1,1]), geom.cell)
            q = fk @ reciprocal(prim)
            g = refolding_vectors(prim, 36)
            basis = FourierBasis(geom, shells, s.norb, 8)
            _, w = weights(basis, c, q, g, np.arange(len(g)), 32)
            self.assertTrue(np.all(w > .01))  # translational symmetry is broken
            np.testing.assert_allclose(w.sum(), 1., atol=8e-4)
            partner = q + reciprocal(geom.cell)[0]
            _, other = weights(basis, c, partner, g, np.arange(len(g)), 32)
            np.testing.assert_allclose(w + other, np.ones(2), atol=8e-4)
            rotated = c * np.exp(1j*np.array([.76, -1.20]))
            np.testing.assert_allclose(abs(basis.factors(q[None]) @ c),
                                       abs(basis.factors(q[None]) @ rotated), atol=1e-13)

    def test_bessel_including_origin_against_angular_integral(self):
        t, dw = np.polynomial.legendre.leggauss(100)
        x = np.array([0., 1e-8, .8, 1., 2., 10., 31.])
        for l in range(4):
            pl = np.polynomial.legendre.Legendre.basis(l)(t)
            expected = ((np.exp(1j*x[:,None]*t) @ (dw*pl))/(2*(1j)**l)).real
            np.testing.assert_allclose(bessel(l, x), expected, atol=2e-13)
            direction = np.array([[.6, 0, .8], [-.2, .4, np.sqrt(.8)]])
            np.testing.assert_allclose((tesseral(l, direction)**2).sum(axis=1),
                                       (2*l+1)/(4*np.pi), atol=2e-14)

    def test_sp_interference_against_real_space_fourier_integral(self):
        cell = np.diag([12., 13., 14.])
        pos = np.array([[1.3, .8, 1.7]])
        geometry = Geometry(["X"], pos, cell)
        shells = [Shell(l, 5., np.array([1.]), np.array([[1.]])) for l in (0, 1)]
        basis = FourierBasis(geometry, {"X": shells}, [4], 2.)
        q = np.array([0., 0., 1.4])
        r, wr = np.polynomial.legendre.leggauss(70)
        r, wr = 2.5*(r+1), 2.5*wr
        t, wt = np.polynomial.legendre.leggauss(70)
        phase = np.exp(-1j*q[2]*r[:,None]*t)
        integ = wr*r*r*np.exp(-r)
        # Independent direct d^3 r integral for a 1*s + i*p_z combination.
        integrand = 1/np.sqrt(4*np.pi) + 1j*r[:,None]*np.sqrt(3/(4*np.pi))*t
        direct = 2*np.pi * (integ @ (phase*integrand) @ wt) / np.sqrt(np.linalg.det(cell))
        direct *= np.exp(-1j*float(q @ pos[0]))
        actual = basis.factors(q[None])[0] @ np.array([1,0,1j,0])
        np.testing.assert_allclose(actual, direct, atol=1e-12)
        self.assertEqual(basis.factors(np.array([[1.,.2,0.]]))[0,2], 0.)
        self.assertGreater(abs(basis.factors(q[None])[0,2]), 0.)

    def test_skew_and_incommensurate_refolding(self):
        ref = np.array([[4.,0,0],[2.,3.5,0],[.1,.2,6.]])
        mat = np.array([[2,1,0],[0,2,0],[0,0,1]])
        cell = mat @ ref
        shifts = refolding_vectors(ref, 10.)
        point = np.array([.13,.22,.01]) @ reciprocal(ref)
        groups, comm = group_momenta(point[None], shifts, cell, ref)
        self.assertTrue(comm)
        self.assertEqual(len(groups), 1)
        np.testing.assert_allclose(physical_k(point+shifts, cell),
                                   np.tile(next(iter(groups)), (len(shifts), 1)), atol=1e-9)
        skewed = cell.copy()
        skewed[0,0] += .23
        groups, comm = group_momenta(point[None], shifts, skewed, ref)
        self.assertFalse(comm)
        self.assertGreater(len(groups), 1)
        for key, jobs in groups.items():
            for _, ids in jobs:
                for g in shifts[ids]:
                    np.testing.assert_allclose(physical_k(point+g, skewed), key, atol=1e-9)
        self.assertEqual(sum(len(ids) for jobs in groups.values() for _,ids in jobs), len(shifts))
        self.assertEqual(tuple(physical_k(np.array([.5,-.5,1.5])@reciprocal(cell),cell)), (-.5,-.5,-.5))

    def test_reject_truncation_complex_and_indefinite_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)/"bad.dat"
            dest.write_text("1\n1 1 1\n1 0 1 0 0 0\n")
            with self.assertRaisesRegex(ValueError, "incomplete"):
                RealMatrix.read(dest)
            (Path(tmp)/"hamimag1.dat").touch()
            with self.assertRaisesRegex(ValueError, "collinear"):
                load_matrices(tmp)
            with self.assertRaisesRegex(ValueError, "positive definiteness"):
                eigenstates(np.eye(2), np.diag([1., -.2]))

    def test_onsite_image_pair_collinear_spin_and_cartesian_gen(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            # DFTB+ writes both translations for same-atom images.
            raw = "1\n1 3 1\n1 0 1 0 0 0\n-0.3\n1 1 1 1 0 0\n-0.1\n1 2 1 -1 0 0\n-0.1\n"
            (folder/"hamreal1.dat").write_text(raw)
            (folder/"hamreal2.dat").write_text(raw.replace("-0.3", "-0.28"))
            (folder/"overreal.dat").write_text(raw.replace("-0.3","1").replace("-0.1","0.1"))
            hs, s = load_matrices(folder)
            k = [.19,0,0]
            cos = np.cos(2*np.pi*k[0])
            np.testing.assert_allclose(hs[0].bloch(k), [[-.3-.2*cos]])
            np.testing.assert_allclose(hs[1].bloch(k), [[-.28-.2*cos]])
            np.testing.assert_allclose(s.bloch(k), [[1+.2*cos]])
            (folder/"s.gen").write_text("1 S\nX\n1 1 12 2 2\n1 1 1\n5 0 0\n0 5 0\n0 0 5\n")
            geom = Geometry.read(folder/"s.gen")
            np.testing.assert_allclose(geom.positions*BOHR_ANG, [[1,1,1]], atol=1e-12)
            cfg = json.loads((CHAIN/"unfold.json").read_text(encoding="utf-8"))
            ref, frac, points, _ = sampling(cfg, geom.cell)
            np.testing.assert_allclose(points, frac@reciprocal(ref))


if __name__ == "__main__":
    unittest.main()
