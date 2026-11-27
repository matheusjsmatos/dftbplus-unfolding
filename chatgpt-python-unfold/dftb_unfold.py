#!/usr/bin/env python3
"""Electronic band unfolding/refolding from DFTB+ WriteRealHS output."""

import argparse
import csv
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

from unfold_fourier import FourierBasis, eigenstates, group_momenta, reciprocal, refolding_vectors
from unfold_io import BOHR_ANG, HARTREE_EV, Geometry, RealMatrix, read_basis


COLUMNS = ("point", "distance_inv_ang", "q1", "q2", "q3", "K1_sc", "K2_sc", "K3_sc",
           "spin", "band", "energy_ev", "energy_rel_ev", "weight_full", "weight_refold")


def sampling(config, cell):
    if ("supercell_matrix" in config) == ("reference_cell_angstrom" in config):
        raise ValueError("specify exactly one of supercell_matrix and reference_cell_angstrom")
    if "supercell_matrix" in config:
        m = np.asarray(config["supercell_matrix"], dtype=float)
        if m.shape != (3, 3) or not np.all(np.isfinite(m)) or not np.array_equal(m, np.rint(m)):
            raise ValueError("supercell_matrix must be an integer 3 x 3 matrix")
        if abs(np.linalg.det(m)) < 0.5:
            raise ValueError("singular supercell_matrix")
        reference = np.linalg.solve(m, cell)
    else:
        reference = np.asarray(config["reference_cell_angstrom"], dtype=float) / BOHR_ANG
    if (reference.shape != (3, 3) or not np.all(np.isfinite(reference))
            or abs(np.linalg.det(reference)) < 1e-12):
        raise ValueError("invalid reference cell")
    corners = np.asarray(config["path_frac"], dtype=float)
    if corners.ndim != 2 or corners.shape[1] != 3 or len(corners) == 0 or not np.isfinite(corners).all():
        raise ValueError("path_frac must contain at least one fractional reciprocal vector")
    intervals = config["intervals"]
    if len(intervals) != len(corners)-1 or any(type(n) != int or n < 1 for n in intervals):
        raise ValueError("one positive integer in intervals is required for each path segment")
    if "labels" in config and len(config["labels"]) != len(corners):
        raise ValueError("one label is required for each vertex")
    parts = [corners[:1]]
    for a, b, count in zip(corners[:-1], corners[1:], intervals):
        parts.append(np.linspace(a, b, count+1)[1:])
    frac = np.concatenate(parts)
    points = frac @ reciprocal(reference)
    distance = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))] / BOHR_ANG
    return reference, frac, points, distance


def load_matrices(directory):
    directory = Path(directory)
    if (any(directory.glob("hamimag*.dat"))
            or any(p.name not in {"hamreal1.dat", "hamreal2.dat"} for p in directory.glob("hamreal*.dat"))):
        raise ValueError("only a scalar Hamiltonian or two real collinear spin channels are supported")
    s = RealMatrix.read(directory / "overreal.dat")
    hs = [RealMatrix.read(directory / "hamreal1.dat")]
    if (directory / "hamreal2.dat").is_file():
        hs.append(RealMatrix.read(directory / "hamreal2.dat"))
    if any(not np.array_equal(s.norb, h.norb) for h in hs):
        raise ValueError("incompatible atom/orbital tables in H and S")
    return hs, s


def weights(basis, c, point, shifts, indices, batch):
    central = np.zeros(c.shape[1])
    total = np.zeros_like(central)
    for start in range(0, len(indices), batch):
        selected = indices[start:start+batch]
        amplitudes = basis.factors(point + shifts[selected]) @ c
        contribution = np.abs(amplitudes)**2
        total += np.sum(contribution, axis=0)
        if selected[0] == 0:  # indices in each group are in increasing order
            central = contribution[0].copy()
    return central, total


def run(config_path, output, plan_only=False):
    config_path = Path(config_path)
    base = config_path.parent
    cfg = json.loads(config_path.read_text(encoding="utf-8-sig"))
    geo = Geometry.read(base / cfg["geometry"])
    ref, frac, points, distance = sampling(cfg, geo.cell)
    gcut = float(cfg.get("gcut_ry", 0))
    shifts = refolding_vectors(ref, gcut)
    groups, commensurate = group_momenta(points, shifts, geo.cell, ref)
    print(f"{len(geo.species)} atoms | {len(points)} path points | {len(shifts)} refolding vectors"
          f" | {len(groups)} distinct K per spin | commensurate: {commensurate}", file=sys.stderr)
    if plan_only:
        return
    hs, s = load_matrices(base / cfg.get("matrix_dir", "."))
    if len(geo.species) != len(s.norb):
        raise ValueError("geometry and matrix have different numbers of atoms")
    wavefiles = cfg["wfc_files"]
    if not isinstance(wavefiles, list) or not wavefiles:
        raise ValueError("wfc_files must list the Waveplot radial file(s) of your SK set")
    shells = read_basis([base / f for f in wavefiles], geo.species, s.norb)
    qmax = np.linalg.norm(points, axis=1).max() + np.sqrt(gcut)
    basis = FourierBasis(geo, shells, s.norb, qmax)
    eref = float(cfg.get("reference_ev", 0))
    window = cfg.get("window_ev")
    if not np.isfinite(eref):
        raise ValueError("reference_ev must be finite")
    if window is not None and (len(window) != 2 or not np.isfinite(window).all() or window[0] > window[1]):
        raise ValueError("window_ev must contain ordered finite bounds relative to reference_ev")
    batch = cfg.get("q_batch", 64)
    if type(batch) != int or batch < 1:
        raise ValueError("q_batch must be a positive integer")
    dest = Path(output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    pending = tempfile.NamedTemporaryFile(mode="w", encoding="ascii", newline="", delete=False,
                                          prefix=dest.name+".", dir=dest.parent)
    nrows = 0
    try:
        with pending as stream:
            writer = csv.writer(stream)
            writer.writerow(COLUMNS)
            for k, jobs in groups.items():
                overlap = s.bloch(k)
                for spin, h in enumerate(hs, 1):
                    e, c = eigenstates(h.bloch(k), overlap)
                    e *= HARTREE_EV
                    relative = e - eref
                    bands = np.arange(len(e))
                    if window is not None:
                        bands = bands[(relative >= window[0]) & (relative <= window[1])]
                    if not len(bands):
                        continue
                    coeff = c[:, bands]
                    for i, selected in jobs:
                        full, refold = weights(basis, coeff, points[i], shifts, selected, batch)
                        prefix = (i+1, distance[i], *frac[i], *k, spin)
                        writer.writerows((*prefix, b+1, e[b], relative[b], full[j], refold[j])
                                         for j, b in enumerate(bands))
                        nrows += len(bands)
        if nrows == 0:
            raise ValueError("no eigenvalues in the requested energy window")
        Path(pending.name).replace(dest)
    finally:
        Path(pending.name).unlink(missing_ok=True)
    print(f"{nrows} spectral lines written to {dest}", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="JSON path, geometry, and Waveplot basis")
    parser.add_argument("-o", "--output", default="unfolded.csv")
    parser.add_argument("--plan", action="store_true", help="count q, g and unique K before diagonalization")
    args = parser.parse_args(argv)
    try:
        run(args.config, args.output, args.plan)
    except (ValueError, OSError, KeyError, np.linalg.LinAlgError) as exc:
        parser.exit(2, f"unfold: {exc}\n")


if __name__ == "__main__":
    main()
