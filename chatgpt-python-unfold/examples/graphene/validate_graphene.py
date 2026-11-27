#!/usr/bin/env python3
"""Run a pristine graphene integration test with an installed DFTB+ and a matching C basis."""

import argparse
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

EXAMPLE = Path(__file__).resolve().parent
sys.path.insert(0, str(EXAMPLE.parents[1]))
from dftb_unfold import run, sampling
from unfold_fourier import physical_k
from unfold_io import Geometry


def kblock(points):
    rows = "\n".join("    " + " ".join(f"{x:.15g}" for x in point) + " 1" for point in points)
    return "{\n" + rows + "\n  }"


def calculate(exe, root, tag, gen, sk, kpoints, restart=None, dump=False):
    folder = root / tag
    folder.mkdir()
    shutil.copyfile(gen, folder / "geom.gen")
    shutil.copyfile(sk, folder / "C-C.skf")
    if restart:
        shutil.copyfile(restart / "charges.bin", folder / "charges.bin")
    text = f"""Geometry = GenFormat {{
  <<< "geom.gen"
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SCCTolerance = 1e-10
  MaxSCCIterations = {1 if restart else 150}
  ReadInitialCharges = {"Yes" if restart else "No"}
  ConvergentSCCOnly = {"No" if restart else "Yes"}
  RecomputeAfterDensity = {"No" if restart else "Yes"}
  SlaterKosterFiles = {{ C-C = "C-C.skf" }}
  MaxAngularMomentum = {{ C = "p" }}
  Filling = Fermi {{ Temperature [Kelvin] = 100 }}
  KPointsAndWeights = {kpoints}
}}
Options {{
  WriteCharges = Yes
  WriteDetailedXML = {"No" if restart else "Yes"}
  ShowFoldedCoords = Yes
  WriteRealHS = {"Yes" if dump else "No"}
}}
Analysis {{ WriteBandOut = Yes }}
"""
    (folder / "dftb_in.hsd").write_text(text, encoding="ascii")
    log = folder / "run.log"
    with log.open("w", encoding="utf-8") as stream:
        proc = subprocess.run([exe], cwd=folder, stdout=stream, stderr=subprocess.STDOUT)
    if dump:
        complete = "Hamilton/Overlap written, exiting program." in log.read_text(errors="replace")
        files = [folder / name for name in ("hamreal1.dat", "overreal.dat")]
        if complete and all(p.is_file() and p.stat().st_size for p in files):
            return folder
        raise RuntimeError(f"matrix export incomplete (exit {proc.returncode}); see {log}")
    if proc.returncode:
        raise RuntimeError(f"DFTB+ exited with status {proc.returncode}; see {log}")
    if restart is None:
        result = (folder / "detailed.out").read_text(errors="replace") + log.read_text(errors="replace")
        if not re.search(r"SCC\s+converged", result, re.I):
            raise RuntimeError(f"SCC convergence was not confirmed; see {log}")
    return folder


def native_bands(file):
    blocks = []
    for line in Path(file).read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "KPT":
            if int(fields[1]) != len(blocks)+1 or int(fields[3]) != 1:
                raise ValueError(f"{file}: unexpected K/spin order")
            blocks.append([])
        else:
            if not blocks or int(fields[0]) != len(blocks[-1])+1:
                raise ValueError(f"{file}: invalid band index")
            blocks[-1].append(float(fields[1]))
    return np.array(blocks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dftb", default="dftb+", help="DFTB+ executable on PATH or absolute path")
    parser.add_argument("--sk-dir", type=Path, required=True, help="directory containing C-C.skf")
    parser.add_argument("--wfc", type=Path, required=True, help="matching Waveplot wfc.*.hsd")
    parser.add_argument("--workdir", type=Path, default=Path("."))
    args = parser.parse_args()
    exe = shutil.which(args.dftb)
    if not exe:
        parser.error(f"DFTB+ executable unavailable: {args.dftb}")
    exe = str(Path(exe).resolve())
    sk = (args.sk_dir / "C-C.skf").resolve()
    if not sk.is_file() or not args.wfc.is_file():
        parser.error("both C-C.skf and a matching radial basis file are required")
    args.workdir.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="graphene_", dir=args.workdir.resolve()))
    print(f"Graphene run: {root}", flush=True)
    for name in ("unfold.json", "primitive.gen", "supercell.gen"):
        shutil.copyfile(EXAMPLE/name, root/name)
    shutil.copyfile(args.wfc, root/"wfc.hsd")
    cfg = json.loads((root/"unfold.json").read_text(encoding="utf-8"))
    sc = Geometry.read(root/"supercell.gen")
    pc = Geometry.read(root/"primitive.gen")
    _, _, q, _ = sampling(cfg, sc.cell)
    ks = physical_k(q, sc.cell)
    unique_sc, lookup = np.unique(ks, axis=0, return_inverse=True)
    lookup = lookup.reshape(-1)
    # For these 2x2 supercells the equivalent primitive SCF mesh is twice as dense.
    for name, geom, n, path in [
            ("pc", root/"primitive.gen", 12, -physical_k(q, pc.cell)),
            ("sc", root/"supercell.gen", 6, -unique_sc)]:
        mesh = f"SupercellFolding {{\n {n} 0 0\n 0 {n} 0\n 0 0 1\n 0 0 0\n }}"
        scf = calculate(exe, root, name+"_scf", geom, sk, mesh)
        calculate(exe, root, name+"_bands", geom, sk, kblock(path), restart=scf)
        if name == "sc":
            calculate(exe, root, "export", geom, sk, mesh, restart=scf, dump=True)
    out = root/"unfolded.csv"
    cfg["geometry"] = "sc_scf/detailed.xml"
    (root/"unfold.json").write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    run(root/"unfold.json", out)
    primitive = native_bands(root/"pc_bands"/"band.out")
    native = native_bands(root/"sc_bands"/"band.out")[lookup]
    energy = np.full_like(native, np.nan)
    weight = np.zeros_like(native)
    with out.open(newline="") as stream:
        for row in csv.DictReader(stream):
            i, j = int(row["point"])-1, int(row["band"])-1
            energy[i, j] = float(row["energy_ev"])
            weight[i, j] = float(row["weight_refold"])
    if not np.isfinite(energy).all() or not np.isfinite(weight).all():
        raise RuntimeError("missing band or nonfinite weight")
    delta = np.abs(energy[:, :, None] - primitive[:, None, :])
    tol_ev = 2e-4  # DFTB+ writes band.out with 4 decimal places in eV
    leak = np.where(delta.min(axis=2) > tol_ev, weight, 0).max()
    coverage = np.max(np.where(delta <= tol_ev, weight[:, :, None], 0), axis=1)
    report = {"points": len(q), "max_supercell_band_error_eV": float(abs(energy-native).max()),
              "max_weight_outside_primitive_bands": float(leak),
              "min_matching_weight_per_primitive_level": float(coverage.min()),
              "energy_tolerance_eV": tol_ev, "weight_tolerance": 1e-6}
    report["passed"] = bool(report["max_supercell_band_error_eV"] <= tol_ev
                            and leak < 1e-6 and coverage.min() > 1e-6)
    (root/"validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
