"""Check scalar DFTB+ eigenvectors against the exported H/S in the native integration tests."""
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

from dftb_unfold import load_matrices
from unfold_fourier import eigenstates
from unfold_io import HARTREE_EV


def compare_state_files(matrix_dir, band_dir, native_bands):
    band_dir = Path(band_dir)
    doc = ET.parse(band_dir/"detailed.xml").getroot()
    n = int(doc.findtext("nroforbitals"))
    nk = int(doc.findtext("nrofkpoints"))
    if int(doc.findtext("nrofspins")) != 1 or int(doc.findtext("nrofstates")) != n:
        raise ValueError("native state diagnostic requires a scalar full eigensystem")
    real = doc.findtext("real").strip().lower() in ("yes", "true", "1")
    dtype = np.dtype("=f8" if real else "=c16")
    kpoints = np.fromstring(doc.findtext("kpointsandweights").replace("D", "E"), sep=" ").reshape(nk, 4)
    energies = native_bands(band_dir/"band.out")
    if energies.shape != (nk, n):
        raise ValueError("the band.out dimensions differ from detailed.xml")
    h, overlap = load_matrices(matrix_dir)
    eigenfile = band_dir/"eigenvec.bin"
    if len(h) != 1 or eigenfile.stat().st_size != 4 + nk*n*n*dtype.itemsize:
        raise ValueError("incompatible spin, stream layout or dimensions")
    audit = dict(max_metric_error=0., max_offdiagonal_H_eV=0.,
                 max_energy_error_eV=0., max_weight_outside_degenerate_subspace=0.)
    with eigenfile.open("rb") as stream:
        if int(np.fromfile(stream, dtype="=i4", count=1)[0]) != int(doc.findtext("identity")):
            raise ValueError("eigenvec.bin and detailed.xml have different run IDs")
        for ik, row in enumerate(kpoints):
            c = np.fromfile(stream, dtype=dtype, count=n*n).reshape((n,n), order="F")
            k = -row[:3]
            s, hk = overlap.bloch(k), h[0].bloch(k)
            e, reconstructed = eigenstates(hk, s)
            e *= HARTREE_EV
            gram = c.conj().T @ s @ c
            hc = (c.conj().T @ hk @ c)*HARTREE_EV
            weights = abs(reconstructed.conj().T @ s @ c)**2
            far = abs(e[:,None]-energies[ik][None,:]) > 2e-4
            audit["max_metric_error"] = max(audit["max_metric_error"], float(abs(gram-np.eye(n)).max()))
            audit["max_offdiagonal_H_eV"] = max(audit["max_offdiagonal_H_eV"],
                                                 float(abs(hc-np.diag(np.diag(hc))).max()))
            audit["max_energy_error_eV"] = max(audit["max_energy_error_eV"], float(abs(e-energies[ik]).max()))
            audit["max_weight_outside_degenerate_subspace"] = max(audit["max_weight_outside_degenerate_subspace"],
                                                                  float(np.where(far, weights, 0.).sum(axis=0).max()))
    audit["passed"] = bool(audit["max_metric_error"] < 1e-6 and audit["max_offdiagonal_H_eV"] < 2e-4
                           and audit["max_energy_error_eV"] < 2e-4
                           and audit["max_weight_outside_degenerate_subspace"] < 1e-5)
    return audit
