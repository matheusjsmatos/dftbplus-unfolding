#!/usr/bin/env python3
"""Compare absolute full/refolded spectral densities on a common path and energy grid."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from dftb_unfold import load_matrices, sampling
from spectrum import spectral_grid
from unfold_io import Geometry, read_basis


def numbers(values):
    return " ".join(f"{value:.16g}" for value in values)


def hsd_path(path):
    value = Path(path).resolve().as_posix()
    if '"' in value or "\n" in value:
        raise ValueError(f"Cannot quote path in HSD: {path}")
    return '"' + value + '"'


def make_input(config_file, native, lo, hi, n, sigma, normalize=False):
    base = config_file.parent
    cfg = json.loads(config_file.read_text(encoding="utf-8-sig"))
    geo = Geometry.read(base / cfg["geometry"])
    reference, _, points, _ = sampling(cfg, geo.cell)
    # Serialize just the orbitals used by the Hamiltonian, even if a distributed WFC lists more shells.
    hs, overlap = load_matrices(base / cfg.get("matrix_dir", "."))
    atomic_numbers = {}
    basis = read_basis([base / name for name in cfg["wfc_files"]], geo.species, overlap.norb,
                       atomic_numbers=atomic_numbers)
    ef = float(cfg.get("reference_ev", 0.0))
    content = ["Action = Unfold", "Output = unfold", "KPointSign = -1",
               "NormalizeOrbitals = " + ("Yes" if normalize else "No"), "SumRuleCheck = No"]
    for key, name in (("DetailedXML", "detailed.xml"), ("BandOut", "band.out"), ("EigenvecBin", "eigenvec.bin")):
        content.append(f"{key} = {hsd_path(native / name)}")
    content += ["RefoldingLatticeVectors = {", *[numbers(v) for v in reference], "}",
                f"RefoldingGCutoff [Hartree] = {float(cfg.get('gcut_ry', 0))/2:.16g}", "QPoints = {"]
    for count, q in zip([1, *cfg["intervals"]], cfg["path_frac"]):
        content.append(f"{count} {numbers(q)}")
    content += ["}", "Energy = {", f"Window [eV] = {{ {lo:.16g} {hi:.16g} }}", f"NrOfIntervals = {n}",
                f"FermiLevel [eV] = {ef:.16g}", f"Broadening [eV] = {sigma:.16g}", "}", "Basis = {"]
    for species, shells in basis.items():
        content.extend([f"{species} = {{", f"AtomicNumber = {atomic_numbers[species]}"])
        for shell in shells:
            content += ["Orbital = {", f"AngularMomentum = {shell.l}", f"Cutoff = {shell.cutoff:.16g}",
                        "Exponents = { " + numbers(shell.alpha) + " }",
                        "Coefficients = { " + numbers(shell.coefficients.ravel()) + " }", "}"]
        content.append("}")
    content.append("}")
    return "\n".join(content) + "\n", points, ef, len(hs)


def read_map(path):
    with path.open(encoding="ascii") as stream:
        fields = stream.readline().split()
        if len(fields) != 5:
            raise ValueError(f"{path}: expected nq nE Emin Emax EF")
        nq, ne = map(int, fields[:2])
        lo, hi, ef = map(float, fields[2:])
        if nq < 1 or ne < 2 or hi <= lo:
            raise ValueError(f"{path}: invalid map dimensions")
        q, bins = np.empty((nq, 3)), np.empty((nq, ne))
        for i in range(nq):
            fields = stream.readline().split()
            if len(fields) != 4:
                raise ValueError(f"{path}: missing q record {i+1}")
            q[i] = np.array(fields[:3], dtype=float)
            values = []
            while len(values) < ne:
                row = stream.readline().split()
                if not row:
                    raise ValueError(f"{path}: incomplete energy block at q {i+1}")
                values.extend(float(x.replace("D", "E")) for x in row)
            if len(values) != ne:
                raise ValueError(f"{path}: extra values at q {i+1}")
            bins[i] = values
        if stream.read().strip() or not all(np.isfinite(a).all() for a in (q, bins, [lo, hi, ef])):
            raise ValueError(f"{path}: nonfinite values or trailing data")
    return q, np.linspace(lo, hi, ne), ef, bins


def compare_maps(prefix, csv_file, points, ef, spins, lo, hi, n, sigma, tol_peak, tol_l1):
    report = {}
    hard_window = np.zeros((spins, len(points), 2))
    with csv_file.open(encoding="ascii", newline="") as stream:
        rows = csv.DictReader(stream)
        seen_spins = set()
        for row in rows:
            ispin, iq = int(row["spin"])-1, int(row["point"])-1
            if not (0 <= ispin < spins and 0 <= iq < len(points)):
                raise ValueError("CSV spin or point outside the expected range")
            seen_spins.add(ispin+1)
            if lo <= float(row["energy_rel_ev"]) <= hi:
                hard_window[ispin, iq] += [float(row["weight_full"]), float(row["weight_refold"])]
            if abs(float(row["energy_ev"])-float(row["energy_rel_ev"])-ef) > 1e-7:
                raise ValueError("The spectrum CSV uses a different reference energy.")
        if seen_spins != set(range(1, spins+1)):
            raise ValueError("Inconsistent spin channels in CSV and matrices.")
    for spin in range(1, spins+1):
        suffix = "" if spins == 1 else (".spinUp" if spin == 1 else ".spinDn")
        for name, column in (("unfoldedBands", "weight_full"), ("refoldedBands", "weight_refold")):
            q, grid, fer, bins = read_map(Path(str(prefix) + "." + name + suffix))
            # Native format stores energy bounds with three decimals and q with eight.
            if len(grid) != n+1 or not np.allclose([grid[0], grid[-1]], [lo, hi], atol=1e-7, rtol=0):
                raise ValueError(f"{name}: energy grid differs; use bounds with at most three decimals")
            if abs(ef-fer) > 5.1e-7 or q.shape != points.shape or not np.allclose(q, points, atol=1.1e-8, rtol=0):
                raise ValueError(f"{name}: path or energy reference differs")
            step = (hi-lo)/n
            _, py_grid, density = spectral_grid(csv_file, column, spin, lo, hi, np.nextafter(step, np.inf), sigma)
            if density.shape != bins.shape or not np.allclose(grid, py_grid, atol=1e-10, rtol=0):
                raise ValueError("Python and Fortran grid dimensions differ")
            delta = np.abs(bins / step - density)
            peak = float(density.max())
            l1 = float(np.abs(density).sum())
            if min(peak, l1) <= 0:
                raise ValueError("Empty comparison window for " + column)
            py_mass = density.sum(axis=1) * step
            native_mass = bins.sum(axis=1)
            bright = py_mass > 1e-10 * py_mass.max()
            mass_ratio = (native_mass[bright]-py_mass[bright])/py_mass[bright]
            col = 0 if column == "weight_full" else 1
            report[f"{column}_spin{spin}"] = {
                "native_binned_weight": float(native_mass.sum()),
                "python_binned_weight": float(py_mass.sum()),
                "python_hard_window_weight": float(hard_window[spin-1, :, col].sum()),
                "relative_integrated_difference_qmin_qmax": [float(mass_ratio.min()), float(mass_ratio.max())],
                "max_abs_per_eV": float(delta.max()), "relative_max": float(delta.max()/peak),
                "relative_L1": float(delta.sum()/l1), "passed": bool(delta.max()/peak <= tol_peak and delta.sum()/l1 <= tol_l1)}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("csv", type=Path, help="unpruned Python spectrum")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--exe", help="native unfold executable, on PATH or absolute")
    group.add_argument("--maps", type=Path, help="prefix of existing maps, e.g. /run/unfold")
    parser.add_argument("--native", type=Path, help="same-run band.out, detailed.xml and eigenvec.bin")
    parser.add_argument("--normalize-native-orbitals", action="store_true",
                        help="use the original Fortran radial rescaling with --exe")
    parser.add_argument("--workdir", type=Path, default=Path("results"))
    parser.add_argument("--range", type=float, nargs=2, default=[-20., 10.], dest="bounds")
    parser.add_argument("--intervals", type=int, default=600)
    parser.add_argument("--sigma", type=float, default=.08)
    parser.add_argument("--tol-peak", type=float, default=.003)
    parser.add_argument("--tol-l1", type=float, default=.002)
    args = parser.parse_args()
    lo, hi = args.bounds
    if not np.isfinite([lo, hi, args.sigma, args.tol_peak, args.tol_l1]).all() or hi <= lo or args.intervals < 1 or args.sigma < 0 or min(args.tol_peak, args.tol_l1) < 0:
        parser.error("require ordered finite bounds, positive interval count, nonnegative sigma and tolerances")
    if max(abs(v-round(v, 3)) for v in (lo, hi)) > 1e-10:
        parser.error("the native map header requires bounds with at most three decimals")
    if args.exe and not args.native:
        parser.error("--exe requires --native")
    cfg = json.loads(args.config.read_text(encoding="utf-8-sig"))
    if cfg.get("window_ev") is not None:
        parser.error("remove window_ev and regenerate the CSV to include all Gaussian tails")
    args.workdir.mkdir(parents=True, exist_ok=True)
    out = Path(tempfile.mkdtemp(prefix="compare_", dir=args.workdir.resolve()))
    try:
        source, points, ef, nspin = make_input(args.config, args.native or Path("."), lo, hi, args.intervals, args.sigma,
                                                      args.normalize_native_orbitals)
        prefix = args.maps
        if args.exe:
            exe = shutil.which(args.exe)
            if exe is None:
                raise ValueError(f"executable unavailable: {args.exe}")
            (out / "unfold_in.hsd").write_text(source, encoding="utf-8")
            log = out / "native.log"
            with log.open("w") as stream:
                proc = subprocess.run([str(Path(exe).resolve())], cwd=out, stdout=stream, stderr=subprocess.STDOUT)
            if proc.returncode:
                raise ValueError(f"native unfold exited {proc.returncode}; see {log}")
            prefix = out / "unfold"
        result = compare_maps(prefix, args.csv, points, ef, nspin, lo, hi, args.intervals, args.sigma,
                              args.tol_peak, args.tol_l1)
        result["passed"] = all(v["passed"] for v in result.values())
        result["grid"] = {"lo_eV": lo, "hi_eV": hi, "intervals": args.intervals, "sigma_eV": args.sigma,
                          "reference_eV": ef, "points": len(points)}
        text = json.dumps(result, indent=2)
        (out / "comparison.json").write_text(text + "\n")
        print(text, "\n", out)
        return 0 if result["passed"] else 1
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"compare: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
