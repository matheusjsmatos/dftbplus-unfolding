#!/usr/bin/env python3
"""Plot an .unfoldedBands / .refoldedBands file written by `unfold` (SIESTA format).

Usage: plot_unfolded.py FILE [-o fig.png] [--vmax FRAC] [--ref primitive_band.out FERMI_EV]

The optional --ref overlays the bands of a primitive-cell DFTB+ calculation (band.out) in red
(energies shifted by FERMI_EV). The x axis is the point index along the q path.
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_unfold(fn):
    with open(fn) as f:
        head = f.readline().split()
        nq, ne = int(head[0]), int(head[1])
        emin, emax, efermi = float(head[2]), float(head[3]), float(head[4])
        w = np.zeros((nq, ne))
        for i in range(nq):
            f.readline()                      # q vector, path index
            for j in range(ne):
                w[i, j] = float(f.readline())
    return np.linspace(emin, emax, ne), w


def read_band_out(fn):
    """DFTB+ band.out: blocks 'KPT n SPIN s KWEIGHT w' followed by index/energy/occupation."""
    bands, cur = [], None
    for line in open(fn):
        p = line.split()
        if not p:
            continue
        if p[0] == "KPT":
            cur = []
            bands.append(cur)
        elif cur is not None and len(p) >= 2:
            cur.append(float(p[1]))
    return np.array(bands)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("-o", "--output", default="unfolded.png")
    ap.add_argument("--vmax", type=float, default=0.3, help="fraction of the maximum weight at saturation")
    ap.add_argument("--ref", nargs=2, metavar=("BAND_OUT", "FERMI_EV"))
    a = ap.parse_args()
    e, w = read_unfold(a.file)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.imshow(w.T, origin="lower", aspect="auto", cmap="Greys",
              extent=[0, w.shape[0] - 1, e[0], e[-1]], vmax=a.vmax * w.max())
    if a.ref:
        b = read_band_out(a.ref[0]) - float(a.ref[1])
        x = np.linspace(0, w.shape[0] - 1, len(b))
        ax.plot(x, b, "r-", lw=0.6, alpha=0.7)
    ax.set_xlabel("ponto do caminho q")
    ax.set_ylabel("E - E_F (eV)")
    fig.tight_layout()
    fig.savefig(a.output, dpi=150)
    print("escrito", a.output)
