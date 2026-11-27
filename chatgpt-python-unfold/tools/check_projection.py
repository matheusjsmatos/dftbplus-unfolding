#!/usr/bin/env python3
"""Independent s,p carbon check of the tabulation/projection in the Fortran source.

This checker shares the DFTB+ H/S eigensolver with the Python implementation. It
evaluates Simpson integrals, orbital factors and absolute weights independently.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from dftb_unfold import load_matrices, sampling
from unfold_fourier import eigenstates, physical_k, refolding_vectors
from unfold_io import Geometry, HARTREE_EV, read_basis


def radial_table(shell, q):
    r = np.linspace(0., shell.cutoff, 6001)
    radial = np.zeros_like(r)
    for alpha, coeff in zip(shell.alpha, shell.coefficients):
        poly = sum(a * r**(shell.l + i) for i, a in enumerate(coeff))
        radial += np.exp(-alpha*r)*poly
    dr = np.ones(len(r))
    dr[1:-1:2], dr[2:-1:2] = 4., 2.
    dr *= shell.cutoff / (3*(len(r)-1))
    norm = float(dr @ ((r*radial)**2))
    integrand = dr * r**2 * radial
    vals = []
    for part in np.array_split(q, int(np.ceil(len(q)/24))):
        x = part[:, None] * r
        j = np.sinc(x/np.pi)
        if shell.l == 1:
            valid = abs(x) > .02
            j = np.divide(j-np.cos(x), x, out=np.zeros_like(x), where=valid)
            z = x[~valid]
            j[~valid] = z*(1/3-z*z*(1/30-z*z/840))
        vals.extend(np.sqrt(2/np.pi)*(j @ integrand))
    return np.array(vals), norm


def carbon_weights(qg, positions, table, volume, coeff, scale):
    length = np.linalg.norm(qg, axis=1)
    unit = np.divide(qg, length[:, None], out=np.zeros_like(qg), where=length[:, None] > 0)
    qgrid, s, p = table
    fs = np.interp(length, qgrid, s)*scale[0]
    fp = np.interp(length, qgrid, p)*scale[1]
    sp = np.empty((len(qg), 4), dtype=complex)
    sp[:, 0] = fs / np.sqrt(4*np.pi)
    sp[:, 1:] = -1j * np.sqrt(3/(4*np.pi)) * unit[:, [1, 2, 0]] * fp[:, None]
    phase = np.exp(-1j*(qg @ positions.T))
    pref = (2*np.pi)**3 / volume
    amplitude = (phase[:, :, None] * sp[:, None, :]).reshape(len(qg), -1) @ coeff
    weight = pref * np.abs(amplitude)**2
    return np.column_stack([weight[0], weight.sum(axis=0)])


def by_energy(energy, weight, tol=1e-7):
    # Rotations of a degenerate eigenbasis change individual weights but preserve their sum.
    starts = np.r_[0, np.flatnonzero(np.diff(energy) > tol) + 1]
    return np.add.reduceat(weight, starts, axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", type=Path, help="CSV from tools/check_graphene.py")
    ap.add_argument("--output", type=Path, default=Path("results/graphene/projection.json"))
    args = ap.parse_args()
    data = ROOT / "examples/graphene/reference"
    config = json.loads((data / "unfold.json").read_text())
    geo = Geometry.read(data / config["geometry"])
    ref, _, points, _ = sampling(config, geo.cell)
    shifts = refolding_vectors(ref, config["gcut_ry"])
    h, overlap = load_matrices(data)
    nc = int(overlap.norb.sum())
    files = [data / x for x in config["wfc_files"]]
    shells = read_basis(files, geo.species, overlap.norb)["C"]
    qgrid = .01 * np.arange(3+np.ceil((np.linalg.norm(points, axis=1).max()+np.sqrt(config["gcut_ry"]))*100))
    (rs, ns), (rp, np_) = [radial_table(shell, qgrid) for shell in shells]
    series = np.full((len(points), nc, 3), np.nan)
    with args.csv.open(newline="") as stream:
        for row in csv.DictReader(stream):
            series[int(row["point"])-1, int(row["band"])-1] = [float(row[k]) for k in ("energy_ev","weight_full","weight_refold")]
    if not np.isfinite(series).all():
        raise ValueError("need a full graphene CSV")
    errors, original, expected = [], [], []
    for idx, q in enumerate(points):
        k = physical_k(q, geo.cell)
        e, c = eigenstates(h[0].bloch(k), overlap.bloch(k))
        if abs(e*HARTREE_EV-series[idx, :, 0]).max() > 1e-7:
            raise ValueError(f"incompatible energies at point {idx+1}")
        args_w = (q + shifts, geo.positions, (qgrid, rs, rp), abs(np.linalg.det(geo.cell)), c)
        ref_w = by_energy(e*HARTREE_EV, series[idx, :, 1:])
        corrected = by_energy(e*HARTREE_EV, carbon_weights(*args_w, [1., 1.]))
        legacy = by_energy(e*HARTREE_EV, carbon_weights(*args_w, 1/np.sqrt([ns, np_])))
        errors.extend(corrected - ref_w)
        original.extend(legacy - ref_w)
        expected.extend(ref_w)
    expected = np.asarray(expected)
    def errors_for(values):
        values = np.abs(values)
        return {"max_abs": values.max(axis=0).tolist(),
                "relative_max": (values.max(axis=0)/expected.max(axis=0)).tolist(),
                "relative_L1": (values.sum(axis=0)/expected.sum(axis=0)).tolist()}
    report = {"method":"independent NumPy reproduction of Fortran equations (shared H/S eigenvectors)",
              "points":len(points), "g_vectors":len(shifts), "columns":["full","refold"],
              "carbon_radial_norms_s_p":[ns, np_], "without_rescaling":errors_for(errors),
              "with_original_rescaling":errors_for(original)}
    report["passed"] = max(report["without_rescaling"]["relative_max"]) < .001
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = json.dumps(report, indent=2)
    args.output.write_text(result + "\n")
    print(result)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
