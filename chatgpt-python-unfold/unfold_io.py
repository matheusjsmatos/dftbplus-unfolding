"""Readers for DFTB+ real-space matrix dumps, periodic GEN, and Waveplot basis files."""

from dataclasses import dataclass
from pathlib import Path
import re
import shlex
from xml.etree import ElementTree

import numpy as np


BOHR_ANG = 0.529177210903
HARTREE_EV = 27.211386245988


def floats(tokens):
    values = np.array([float(t.replace("D", "E").replace("d", "e")) for t in tokens])
    if not np.all(np.isfinite(values)):
        raise ValueError("non-finite input value")
    return values


def data_lines(stream):
    for line in stream:
        text = line.partition("#")[0].strip()
        if text:
            yield text.split()


@dataclass
class Geometry:
    species: list
    positions: np.ndarray  # Cartesian, bohr, folded relative to the GEN origin
    cell: np.ndarray       # lattice vectors are rows, bohr

    @classmethod
    def read(cls, path):
        if Path(path).suffix.lower() == ".xml":
            return cls.from_detailed_xml(path)
        with Path(path).open(encoding="utf-8-sig") as stream:
            rows = list(data_lines(stream))
        if not rows or len(rows[0]) != 2 or rows[0][1].upper() not in {"S", "F"}:
            raise ValueError(f"{path}: a periodic GEN geometry (S or F) is required")
        n = int(rows[0][0])
        if n < 1 or len(rows) != n + 6:
            raise ValueError(f"{path}: invalid atom count or missing lattice vectors")
        names = rows[1]
        species, positions = [], []
        for number, row in enumerate(rows[2:2+n], 1):
            if len(row) != 5 or int(row[0]) != number or not 1 <= int(row[1]) <= len(names):
                raise ValueError(f"{path}: invalid atom index/species on atom {number}")
            species.append(names[int(row[1])-1])
            positions.append(floats(row[2:]))
        origin, *lattice = [floats(row) for row in rows[-4:]]
        cell = np.asarray(lattice) / BOHR_ANG
        if cell.shape != (3, 3) or origin.shape != (3,) or abs(np.linalg.det(cell)) < 1e-12:
            raise ValueError(f"{path}: singular or invalid lattice")
        if rows[0][1].upper() == "S":
            frac = np.linalg.solve(cell.T, ((positions - origin) / BOHR_ANG).T).T
        else:
            # The reference DFTB+ GEN reader first scales the origin to atomic
            # units, then applies the lattice to both fractional arrays.
            frac = np.asarray(positions) - origin / BOHR_ANG
        return cls(species, (frac % 1.0) @ cell, cell)

    @classmethod
    def from_detailed_xml(cls, path):
        """Keep the same central images selected by DFTB+ with ShowFoldedCoords."""
        try:
            geo = ElementTree.parse(path).getroot().find("geometry")
        except ElementTree.ParseError as exc:
            raise ValueError(f"{path}: invalid detailed.xml") from exc
        if geo is None:
            raise ValueError(f"{path}: no DFTB+ geometry block")
        yes = {"yes", "true", "t", ".true."}
        if (geo.findtext("periodic", "").strip().lower() not in yes
                or geo.findtext("helical", "no").strip().lower() in yes):
            raise ValueError(f"{path}: a periodic non-helical geometry is required")
        names = shlex.split(geo.findtext("typenames", ""))
        table = geo.findtext("typesandcoordinates", "").split()
        if not table or len(table) % 4:
            raise ValueError(f"{path}: invalid typesandcoordinates")
        indices = np.array([int(t) for t in table[::4]])
        if np.any((indices < 1) | (indices > len(names))):
            raise ValueError(f"{path}: invalid species index")
        coords = floats(table).reshape(-1, 4)[:, 1:]
        cell = floats(geo.findtext("latticevectors", "").split()).reshape(3, 3)
        origin = floats(geo.findtext("coordinateorigin", "").split()).reshape(3)
        if abs(np.linalg.det(cell)) < 1e-12:
            raise ValueError(f"{path}: singular lattice")
        coords = coords - origin
        frac = np.linalg.solve(cell.T, coords.T).T
        if np.any((frac < -1e-8) | (frac > 1+1e-8)):
            raise ValueError(f"{path}: regenerate with ShowFoldedCoords = Yes at the export geometry")
        return cls([names[i-1] for i in indices], coords, cell)


@dataclass
class RealMatrix:
    norb: np.ndarray
    rows: np.ndarray
    cols: np.ndarray
    shifts: np.ndarray
    shift_index: np.ndarray
    values: np.ndarray

    @classmethod
    def read(cls, path):
        with Path(path).open(encoding="ascii") as stream:
            rows = data_lines(stream)
            try:
                return cls._read(rows)
            except (ValueError, StopIteration, IndexError) as exc:
                raise ValueError(f"{path}: invalid or incomplete DFTB+ matrix ({exc})") from exc

    @classmethod
    def _read(cls, lines):
        def integers(count):
            fields = next(lines)
            if len(fields) != count:
                raise ValueError(f"expected {count} integers, found {fields}")
            return tuple(map(int, fields))

        (natom,) = integers(1)
        if natom < 1:
            raise ValueError("invalid NATOM")
        norb, nneigh = [], []
        for atom in range(1, natom + 1):
            index, nn, no = integers(3)
            if index != atom or min(nn, no) < 1:
                raise ValueError("invalid atom table")
            norb.append(no)
            nneigh.append(nn)
        offsets = np.cumsum([0] + norb)
        shift_ids = {}
        entries = []
        for i, nn in enumerate(nneigh):
            for neighbour in range(nn):
                ia, nb, ja, x, y, z = integers(6)
                j = ja - 1
                if (ia, nb) != (i + 1, neighbour) or not i <= j < natom:
                    raise ValueError("invalid neighbour block header")
                if nb == 0 and (i != j or (x, y, z) != (0, 0, 0)):
                    raise ValueError("first block must be onsite")
                shift_id = shift_ids.setdefault((x, y, z), len(shift_ids))
                block = np.array([floats(next(lines)) for _ in range(norb[i])])
                if block.shape != (norb[i], norb[j]):
                    raise ValueError("orbital block size disagrees with atom table")
                # writeSparse prints orbital j across each line. unpackHS
                # uses j as the row index. Only the global lower triangle
                # survives adjointLowerTriangle, including self-image blocks.
                col, row = np.indices(block.shape)
                col, row = col + offsets[i], row + offsets[j]
                keep = (row >= col) & (block != 0)
                entries.append((row[keep], col[keep], np.full(keep.sum(), shift_id), block[keep]))
        if next(lines, None) is not None:
            raise ValueError("extra entries after last block")
        row, col, ids, values = (np.concatenate(a) for a in zip(*entries))
        return cls(np.array(norb), row, col, np.array(list(shift_ids)), ids, values)

    def bloch(self, k):
        """Hermitian matrix for physical psi(r+R)=exp(+2 pi i k.R) psi(r)."""
        phase = np.exp(-2j * np.pi * (self.shifts @ k))
        n = int(self.norb.sum())
        lower = np.zeros((n, n), dtype=complex)
        np.add.at(lower, (self.rows, self.cols), self.values * phase[self.shift_index])
        if np.max(abs(lower.diagonal().imag)) > 2e-9 * max(1., np.max(abs(lower))):
            raise ValueError("complex diagonal in real-space dump; check self-image blocks")
        lower[np.diag_indices(n)] = lower.diagonal().real
        return lower + np.tril(lower, -1).conj().T


@dataclass
class Shell:
    l: int
    cutoff: float               # bohr
    alpha: np.ndarray
    coefficients: np.ndarray     # [exponent, polynomial power]


def _hsd_tree(text):
    """The numeric, literal HSD subset in distributed wfc.*.hsd files."""
    text = re.sub(r"#.*|!.*", "", text)
    tokens = re.findall(r'[^\s{}=;]+|[{}=]', text)
    pos = 0

    def group(nested):
        nonlocal pos
        nodes = {}
        while pos < len(tokens) and tokens[pos] != "}":
            key = tokens[pos].lower()
            pos += 1
            if pos < len(tokens) and tokens[pos] == "=":
                pos += 1
            if pos >= len(tokens):
                raise ValueError("missing value for " + key)
            if tokens[pos] == "{":
                pos += 1
                if key in {"exponents", "coefficients"}:
                    start = pos
                    while pos < len(tokens) and tokens[pos] != "}":
                        pos += 1
                    if pos == len(tokens):
                        raise ValueError("unclosed array " + key)
                    value = floats(tokens[start:pos])
                    pos += 1
                else:
                    value = group(True)
            else:
                value = tokens[pos]
                pos += 1
            nodes.setdefault(key, []).append(value)
        if nested:
            if pos == len(tokens):
                raise ValueError("unclosed HSD block")
            pos += 1
        elif pos != len(tokens):
            raise ValueError("unexpected closing brace")
        return nodes

    return group(False)


def read_basis(paths, species, orbital_counts):
    needed = dict(zip(species, orbital_counts))
    if any(needed[s] != n for s, n in zip(species, orbital_counts)):
        raise ValueError("inconsistent orbital counts for a single GEN species")
    basis = {}
    for path in paths:
        tree = _hsd_tree(Path(path).read_text(encoding="utf-8-sig"))
        if "basis" in tree:
            tree = tree["basis"][0]
        for name, target in needed.items():
            if name.lower() not in tree:
                continue
            if name in basis or len(tree[name.lower()]) != 1:
                raise ValueError(f"duplicate basis for {name}")
            shells = []
            count = 0
            for orb in tree[name.lower()][0].get("orbital", []):
                if count == target:
                    break
                l = int(orb["angularmomentum"][0])
                rc = float(orb["cutoff"][0])
                alpha = np.atleast_1d(orb["exponents"][0]).astype(float)
                coeff = np.atleast_1d(orb["coefficients"][0]).astype(float)
                if (l != len(shells) or l > 3 or not np.isfinite(rc) or rc <= 0
                        or not len(alpha) or np.any(alpha <= 0) or not len(coeff)
                        or len(coeff) % len(alpha)):
                    raise ValueError(f"unsupported or invalid Waveplot orbital for {name}")
                shells.append(Shell(l, rc, alpha, coeff.reshape(len(alpha), -1)))
                count += 2*l + 1
            if count != target:
                raise ValueError(f"{name}: basis supplies {count} orbitals, matrix requires {target}")
            basis[name] = shells
    missing = needed.keys() - basis.keys()
    if missing:
        raise ValueError("missing Waveplot radial basis for " + ", ".join(sorted(missing)))
    return basis
