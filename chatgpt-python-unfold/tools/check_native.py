#!/usr/bin/env python3
"""Compare both programs on fresh DFTB+ 26.1 output for graphene or diamond."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
sys.path[:0] = [str(ROOT), str(ROOT / "examples/graphene")]
from validate_graphene import calculate, compare_native, kblock, native_bands
from check_native_states import compare_state_files
from dftb_unfold import run, sampling
from unfold_fourier import physical_k
from unfold_io import Geometry


def write_gen(path, positions, lattice):
    lines = [f"{len(positions)} F", "C"]
    for i, xyz in enumerate(positions, 1):
        lines.append(f"{i} 1 " + " ".join(map(str, xyz)))
    lines.append("0 0 0")
    lines.extend(" ".join(map(str, v)) for v in lattice)
    path.write_text("\n".join(lines)+"\n", encoding="ascii")


def prepare_diamond(work, dftb, sk, wfc):
    a = 3.567
    fc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    write_gen(work/"sc.gen", np.vstack([fc, fc+.25]), np.eye(3)*a)
    write_gen(work/"pc.gen", [[0, 0, 0], [.25, .25, .25]], fc[1:]*a)
    shutil.copyfile(wfc, work/"wfc.hsd")
    cfg = json.loads((PROJECT/"validation/diamond/unfold.json").read_text())
    cfg.update(geometry="sc_scf/detailed.xml", matrix_dir="export", wfc_files=["wfc.hsd"])
    config = work/"unfold.json"
    config.write_text(json.dumps(cfg, indent=2))
    sc, pc = Geometry.read(work/"sc.gen"), Geometry.read(work/"pc.gen")
    _, _, q, _ = sampling(cfg, sc.cell)
    k, lookup = np.unique(physical_k(q, sc.cell), axis=0, return_inverse=True)
    mesh = "SupercellFolding {\n 4 0 0\n 0 4 0\n 0 0 4\n .5 .5 .5\n}"
    # Use the same converged supercell charges for both band-path and H/S export.
    scf = calculate(dftb, work, "sc_scf", work/"sc.gen", sk, mesh)
    calculate(dftb, work, "sc_bands", work/"sc.gen", sk, kblock(-k), restart=scf, eigenvectors=True)
    calculate(dftb, work, "export", work/"sc.gen", sk, mesh, restart=scf, dump=True)
    # Unfold the conventional SCF mesh to an equivalent set for the primitive cell.
    sc_mesh = (np.indices((4, 4, 4)).reshape(3, -1).T + .5)/4.
    shifts = np.vstack([[0, 0, 0], np.eye(3)])
    primitive_k = ((sc_mesh[:, None, :] + shifts) @ np.linalg.inv(sc.cell).T @ pc.cell.T).reshape(-1, 3)
    primitive_k = primitive_k.round(12)
    primitive_k -= np.floor(primitive_k + .5)
    primitive_k = np.unique(primitive_k.round(10), axis=0)
    if len(primitive_k) != 256:
        raise ValueError("unexpected unfolded diamond SCF mesh")
    pscf = calculate(dftb, work, "pc_scf", work/"pc.gen", sk, kblock(primitive_k))
    calculate(dftb, work, "pc_bands", work/"pc.gen", sk, kblock(-physical_k(q, pc.cell)), restart=pscf)
    csv = work/"unfolded.csv"
    run(config, csv)
    report = compare_native(csv, work/"pc_bands/band.out", work/"sc_bands/band.out", lookup.ravel())
    (work/"validation.json").write_text(json.dumps(report, indent=2)+"\n")
    return work


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("case", choices=["graphene", "diamond"])
    ap.add_argument("--dftb", required=True, type=Path)
    ap.add_argument("--unfold", required=True, type=Path)
    ap.add_argument("--mio", required=True, type=Path, help="checkout with skfiles/C-C.skf and extras/wfc.hsd")
    ap.add_argument("--output", type=Path, default=Path("results/native"))
    args = ap.parse_args()
    dftb, unfold = args.dftb.resolve(), args.unfold.resolve()
    sk, wfc = args.mio/"skfiles/C-C.skf", args.mio/"extras/wfc.hsd"
    hashes = ("a0e2fa003c192f7ff3f518045bd6fca9f146eecb79d44c90fe2ed27a361c8996",
              "692f8d91392b8970c44ff1457224164078e344d0537b2a89e3c19fd2edc68d47")
    for path, expected in zip([sk, wfc], hashes):
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            ap.error(f"unexpected mio-1-1 parameter file: {path}")
    args.output.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=args.case+"_", dir=args.output.resolve()))
    print(f"Working directory: {work}", flush=True)
    if args.case == "diamond":
        calc = prepare_diamond(work, str(dftb), sk.resolve(), wfc.resolve())
        bounds, intervals = (-25, 15), 400
    else:
        subprocess.run([sys.executable, str(ROOT/"examples/graphene/validate_graphene.py"),
                        "--dftb", str(dftb), "--sk-dir", str(sk.resolve().parent),
                        "--wfc", str(wfc.resolve()), "--workdir", str(work), "--save-eigenvectors"], check=True)
        calc, = work.glob("graphene_*")
        bounds, intervals = (-20, 10), 600
    ok = bool(json.loads((calc/"validation.json").read_text())["passed"])
    state_report = compare_state_files(calc/"export", calc/"sc_bands", native_bands)
    (calc/"states.json").write_text(json.dumps(state_report, indent=2)+"\n")
    print(json.dumps(state_report, indent=2))
    ok = ok and state_report["passed"]
    # Match both common definitions of a binned spectrum, including the original diamond default.
    for mode, sigma in (("linear", 0), ("gaussian", .08)):
        status = subprocess.run([sys.executable, str(ROOT/"tools/compare_backends.py"), str(calc/"unfold.json"),
                                 str(calc/"unfolded.csv"), "--exe", str(unfold), "--native", str(calc/"sc_bands"),
                                 "--workdir", str(work/mode), "--range", *map(str, bounds), "--intervals",
                                 str(intervals), "--sigma", str(sigma)]).returncode
        ok = ok and status == 0
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
