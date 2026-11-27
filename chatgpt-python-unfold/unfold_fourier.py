"""Fourier projection of non-orthogonal DFTB eigenstates onto normalized plane waves."""

from itertools import product
from math import ceil, prod

import numpy as np


def reciprocal(cell):
    return 2 * np.pi * np.linalg.inv(cell).T


def refolding_vectors(cell, cutoff_ry):
    """Enumerate |g|^2 < cutoff_Ry in atomic units, including the origin first."""
    if cutoff_ry == 0:
        return np.zeros((1, 3))
    if not np.isfinite(cutoff_ry) or cutoff_ry < 0:
        raise ValueError("gcut_ry must be finite and nonnegative")
    # |g.a_i|/(2 pi) bounds each integer even for skew cells.
    limits = np.ceil(np.sqrt(cutoff_ry) * np.linalg.norm(cell, axis=1) / (2*np.pi)).astype(int)
    b = reciprocal(cell)
    vectors = [np.zeros(3)]
    plane = np.arange(-limits[2], limits[2] + 1)[:, None] * b[2]
    for i, j in product(range(-limits[0], limits[0]+1), range(-limits[1], limits[1]+1)):
        g = i*b[0] + j*b[1] + plane
        squared = np.einsum("ij,ij->i", g, g)
        vectors.extend(g[(squared < cutoff_ry) & (squared > 0)])
    return np.asarray(vectors)


def physical_k(qs, cell):
    """Represent q mod G_SC in the centered reciprocal parallelepiped."""
    f = np.asarray(qs) @ cell.T / (2*np.pi)
    # Round before wrapping, so floating error cannot split +/-0.5 into different keys.
    f = np.round(f, 10)
    return np.round(f - np.floor(f + 0.5), 10)


def bessel(l, x):
    """Spherical j_l for l <= 3. The series avoids cancellation at the origin."""
    out = np.empty_like(x, dtype=float)
    small = (np.abs(x) < 1)
    z = x[small]
    term = z**l / prod(range(1, 2*l+2, 2))
    series = term.copy()
    for m in range(1, 10):
        term = term * (-z*z / (2*m*(2*l + 2*m + 1)))
        series += term
    out[small] = series
    z = x[~small]
    j0 = np.sin(z)/z
    if l == 0:
        out[~small] = j0
    else:
        j1 = (j0-np.cos(z))/z
        for order in range(1, l):
            j0, j1 = j1, (2*order+1)*j1/z-j0
        out[~small] = j1
    return out


def tesseral(l, unit):
    """Normalized real harmonics in DFTB+/Waveplot order m=-l,...,+l."""
    x, y, z = unit.T
    if l == 0:
        return np.full((len(x), 1), 1/np.sqrt(4*np.pi))
    if l == 1:
        return np.column_stack([y, z, x]) * np.sqrt(3/(4*np.pi))
    if l == 2:
        return np.column_stack([np.sqrt(15)*x*y, np.sqrt(15)*y*z,
                                np.sqrt(5)/2*(2*z*z-x*x-y*y), np.sqrt(15)*x*z,
                                np.sqrt(15)/2*(x*x-y*y)]) / np.sqrt(4*np.pi)
    if l == 3:
        a, b = np.sqrt(35/(32*np.pi)), np.sqrt(21/(32*np.pi))
        return np.column_stack([
            a*y*(3*x*x-y*y), np.sqrt(105/(4*np.pi))*x*y*z,
            b*y*(4*z*z-x*x-y*y), np.sqrt(7/(16*np.pi))*z*(2*z*z-3*x*x-3*y*y),
            b*x*(4*z*z-x*x-y*y), np.sqrt(105/(16*np.pi))*z*(x*x-y*y),
            a*x*(x*x-3*y*y)])
    raise ValueError("only s, p, d and f shells are supported")


def radial_grid(shell, qmax):
    # 12-point composite Gauss-Legendre. Limit phase change and exponential
    # decay per panel. Coefficients are the Waveplot polynomial in ascending powers.
    nodes, weights = np.polynomial.legendre.leggauss(12)
    panels = max(4, ceil(shell.cutoff * max(qmax, max(shell.alpha)) / 4))
    half = shell.cutoff / (2*panels)
    r = (half * (2*np.arange(panels)[:, None] + 1 + nodes)).ravel()
    dr = np.tile(half*weights, panels)
    radial = np.zeros_like(r)
    for alpha, coeff in zip(shell.alpha, shell.coefficients):
        radial += np.polynomial.polynomial.polyval(r, coeff) * np.exp(-alpha*r)
    return r, dr*r**(shell.l+2)*radial


class FourierBasis:
    def __init__(self, geometry, basis, counts, qmax):
        self.volume = abs(np.linalg.det(geometry.cell))
        self.norb = int(sum(counts))
        self.groups = []
        offsets = np.zeros(len(counts), dtype=np.intp)
        offsets[1:] = np.cumsum(counts[:-1], dtype=np.intp)
        names = np.asarray(geometry.species)
        for name, shells in basis.items():
            atoms = np.flatnonzero(names == name)
            entries = []
            for shell in shells:
                l = shell.l
                cols = offsets[atoms, None] + l*l + np.arange(2*l+1)
                r, w = radial_grid(shell, qmax)
                entries.append((l, cols, r, w))
            self.groups.append((geometry.positions[atoms], entries))

    def factors(self, q):
        """<q|AO> in the cell-translation Bloch gauge, rows q and columns AO."""
        sizes = np.linalg.norm(q, axis=1)
        unit = np.divide(q, sizes[:, None], out=np.zeros_like(q), where=sizes[:, None] != 0)
        form = np.empty((len(q), self.norb), dtype=complex)
        common = 4*np.pi/np.sqrt(self.volume)
        for positions, shells in self.groups:
            phase = np.exp(-1j*(q @ positions.T))
            for l, cols, radius, weights in shells:
                radial = bessel(l, sizes[:, None]*radius) @ weights
                angular = tesseral(l, unit)
                shell = (common * (-1j)**l * radial[:, None]) * angular
                form[:, cols] = phase[:, :, None] * shell[:, None, :]
        return form


def eigenstates(h, s):
    """Solve H C = S C E, C^H S C = I; use SciPy's generalized LAPACK if installed."""
    try:
        from scipy.linalg import eigh
    except ImportError:
        eigh = None
    try:
        if eigh is not None:
            return eigh(h, s, check_finite=False)
        lower = np.linalg.cholesky(s)
        left = np.linalg.solve(lower, h)
        orthogonal = np.linalg.solve(lower, left.conj().T).conj().T
        e, u = np.linalg.eigh((orthogonal + orthogonal.conj().T) / 2)
        return e, np.linalg.solve(lower.conj().T, u)
    except np.linalg.LinAlgError as exc:
        raise ValueError("generalized eigensolver failed; check positive definiteness of S and the input pair") from exc


def group_momenta(points, shifts, cell, reference):
    """Group (path point, reciprocal-vector indices) by a common diagonalization."""
    ratio = cell @ np.linalg.inv(reference)
    commensurate = np.allclose(ratio, np.rint(ratio), atol=2e-9, rtol=0)
    groups = {}
    all_g = np.arange(len(shifts))
    for i, point in enumerate(points):
        if commensurate or len(shifts) == 1:
            k = tuple(physical_k(point, cell))
            groups.setdefault(k, []).append((i, all_g))
        else:
            folded = physical_k(point + shifts, cell)
            unique, inverse = np.unique(folded, axis=0, return_inverse=True)
            inverse = inverse.reshape(-1)
            order = np.argsort(inverse, kind="stable")
            bounds = np.r_[0, np.cumsum(np.bincount(inverse))]
            for j, k in enumerate(unique):
                groups.setdefault(tuple(k), []).append((i, order[bounds[j]:bounds[j+1]]))
    return groups, commensurate
