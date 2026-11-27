#!/usr/bin/env python3
"""Broaden a CSV spectrum on an eV grid. Optionally draw a standalone PDF."""
import argparse
import csv
from pathlib import Path

import numpy as np


def spectral_grid(filename, column, spin, lo, hi, step, sigma):
    if not np.isfinite([lo, hi, step, sigma]).all() or hi <= lo or min(step, sigma) <= 0:
        raise ValueError("invalid energy limits, step or Gaussian width")
    grid = np.linspace(lo, hi, int(np.ceil((hi-lo)/step))+1)
    profiles, coordinates = {}, {}
    with open(filename, newline="", encoding="ascii") as f:
        for row in csv.DictReader(f):
            if int(row["spin"]) != spin:
                continue
            idx = int(row["point"])
            coordinates[idx] = float(row["distance_inv_ang"])
            if idx not in profiles:
                profiles[idx] = np.zeros_like(grid)
            energy = float(row["energy_rel_ev"])
            weight = float(row[column])
            left, right = np.searchsorted(grid, [energy-8*sigma, energy+8*sigma])
            if weight and right > left:
                t = (grid[left:right]-energy)/sigma
                profiles[idx][left:right] += weight*np.exp(-0.5*t*t)/(sigma*np.sqrt(2*np.pi))
    if not profiles:
        raise ValueError(f"no data for spin {spin}")
    ids = sorted(profiles)
    if ids != list(range(1, ids[-1]+1)):
        raise ValueError("some path points are absent; rerun unfolding with all bands for a map")
    return np.array([coordinates[i] for i in ids]), grid, np.array([profiles[i] for i in ids])


def draw_pdf(path, x, energy, values, column, sigma, spin):
    from PIL import Image
    from reportlab.lib.utils import ImageReader
    from reportlab.pdfgen import canvas

    if len(x) < 2 or not np.all(np.diff(x) > 0):
        raise ValueError("the PDF requires a path with positive length and distinct adjacent points")
    vmax = values.max()
    if vmax <= 0:
        raise ValueError("all plotted weights are zero within this energy window")
    # The data need not have equal path steps. Interpolate in physical distance.
    columns = np.linspace(x[0], x[-1], 860)
    pixels = np.array([np.interp(columns, x, row) for row in values.T])[::-1]
    image = Image.fromarray(np.uint8(255*(1-np.clip(pixels/vmax, 0, 1))))
    pdf = canvas.Canvas(str(path), pagesize=(620, 424))
    pdf.setTitle("DFTB+ unfolded electronic spectrum")
    left, bottom, width, height = 70, 70, 460, 285
    pdf.setFont("Helvetica-Bold", 17)
    pdf.drawString(left, 391, "DFTB+ electronic spectrum")
    pdf.setFont("Helvetica", 10)
    name = "Refolded" if column == "weight_refold" else "Fully unfolded"
    pdf.drawString(left, 372, f"{name}  |  spin {spin}  |  Gaussian width {sigma:g} eV")
    pdf.drawImage(ImageReader(image), left, bottom, width, height)
    pdf.setLineWidth(0.6)
    pdf.rect(left, bottom, width, height)
    pdf.setFont("Helvetica", 9)
    for v in np.linspace(energy[0], energy[-1], 6):
        y = bottom+height*(v-energy[0])/(energy[-1]-energy[0])
        pdf.line(left-4, y, left, y)
        pdf.drawRightString(left-7, y-3, f"{v:.2f}")
    for v in np.linspace(x[0], x[-1], 6):
        px = left+width*(v-x[0])/(x[-1]-x[0])
        pdf.line(px, bottom-4, px, bottom)
        pdf.drawCentredString(px, bottom-17, f"{v:.2f}")
    pdf.drawCentredString(left+width/2, bottom-37, "Distance along reciprocal-space path (1/Angstrom)")
    pdf.saveState()
    pdf.translate(19, bottom+height/2)
    pdf.rotate(90)
    pdf.drawCentredString(0, 0, "Energy relative to reference (eV)")
    pdf.restoreState()
    for i in range(100):
        pdf.setFillGray(1-i/99)
        pdf.rect(548, bottom+i*height/100, 10, height/100+0.1, fill=1, stroke=0)
    pdf.setFillGray(0)
    pdf.drawString(563, bottom-2, "0")
    pdf.drawString(563, bottom+height-4, f"{vmax:.2g}")
    pdf.drawString(544, bottom+height+9, "1/eV")
    pdf.showPage()
    pdf.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--range", nargs=2, type=float, required=True, metavar=("MIN", "MAX"))
    parser.add_argument("--sigma", type=float, default=.06)
    parser.add_argument("--step", type=float, default=.01)
    parser.add_argument("--spin", type=int, default=1)
    parser.add_argument("--weight", choices=["full", "refold"], default="refold")
    parser.add_argument("-o", "--output", default="spectrum.npz", type=Path)
    parser.add_argument("--pdf", type=Path, help="PDF figure (requires Pillow and ReportLab)")
    args = parser.parse_args()
    if args.output.suffix.lower() != ".npz":
        parser.error("output must have the .npz extension")
    try:
        column = "weight_"+args.weight
        x, energy, values = spectral_grid(args.csv, column, args.spin, *args.range, args.step, args.sigma)
        np.savez_compressed(args.output, distance_inv_ang=x, energy_rel_ev=energy, spectral_weight_per_ev=values,
                            spin=args.spin, sigma_ev=args.sigma, weight=column)
        if args.pdf:
            draw_pdf(args.pdf, x, energy, values, column, args.sigma, args.spin)
    except (ValueError, KeyError, ImportError, OSError) as exc:
        parser.exit(2, f"spectrum: {exc}\n")


if __name__ == "__main__":
    main()
