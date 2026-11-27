#!/usr/bin/env python3
"""Recompute the 2x2 graphene result against two archived DFTB+ band structures."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "examples/graphene"))
from validate_graphene import compare_native
from dftb_unfold import run, sampling
from unfold_fourier import physical_k
from unfold_io import Geometry


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/graphene"))
    args = parser.parse_args()
    data = ROOT / "examples/graphene/reference"
    config = json.loads((data / "unfold.json").read_text())
    geo = Geometry.read(data / config["geometry"])
    _, _, q, _ = sampling(config, geo.cell)
    _, mapping = np.unique(physical_k(q, geo.cell), axis=0, return_inverse=True)
    args.output.mkdir(parents=True, exist_ok=True)
    csv_file = args.output / "unfolded.csv"
    run(data / "unfold.json", csv_file)
    report = compare_native(csv_file, data / "primitive.band.out", data / "supercell.band.out", mapping.ravel())
    text = json.dumps(report, indent=2)
    (args.output / "validation.json").write_text(text + "\n")
    print(text)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
