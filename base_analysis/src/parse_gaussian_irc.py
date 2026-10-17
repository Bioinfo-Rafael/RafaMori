"""Minimal merge-aware Gaussian 16 IRC parser used by this project."""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ELEMENTS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe",
]
FLOAT = r"[-+]?(?:\d+\.?(?:\d*)?|\.\d+)(?:[DEde][-+]?\d+)?"
POINT_PATH = re.compile(r"Point\s+Number\s*:?\s*(\d+).*?Path\s+Number\s*:?\s*(\d+)", re.I)
POINT_DIRECTION = re.compile(r"Point\s+Number\s*:?\s*(\d+)\s+in\s+(FORWARD|REVERSE)", re.I)
NET_COORDINATE = re.compile(r"NET\s+REACTION\s+COORDINATE.*?=\s*(" + FLOAT + r")", re.I)
SCF_ENERGY = re.compile(r"SCF Done:.*?=\s*(" + FLOAT + r")", re.I)
CHECKPOINT_ENERGY = re.compile(r"Energy\s+From\s+Chk\s*=\s*(" + FLOAT + r")", re.I)


@dataclass
class Point:
    path_number: int
    point_number: int
    coordinate: float
    energy: float
    atom_numbers: np.ndarray
    atomic_numbers: np.ndarray
    xyz: np.ndarray
    forces: np.ndarray
    source_line: int
    coordinate_estimated: bool = False

    @property
    def completeness(self) -> int:
        return int(np.isfinite(self.coordinate)) + int(np.isfinite(self.energy)) + int(len(self.xyz) > 0) + int(np.isfinite(self.forces).any())


def _float(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def _orientation(lines: list[str], start: int) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    i, separators = start + 1, 0
    while i < min(len(lines), start + 16):
        if re.match(r"\s*-{5,}\s*$", lines[i]):
            separators += 1
            if separators == 2:
                i += 1
                break
        i += 1
    pattern = re.compile(r"^\s*(\d+)\s+(\d+)\s+-?\d+\s+(" + FLOAT + r")\s+(" + FLOAT + r")\s+(" + FLOAT + r")")
    atoms, atomic_numbers, xyz = [], [], []
    while i < len(lines):
        if re.match(r"\s*-{5,}\s*$", lines[i]):
            break
        match = pattern.match(lines[i])
        if not match:
            if atoms and lines[i].strip():
                break
        else:
            atoms.append(int(match.group(1)))
            atomic_numbers.append(int(match.group(2)))
            xyz.append([_float(match.group(3)), _float(match.group(4)), _float(match.group(5))])
        i += 1
    return i, np.asarray(atoms, int), np.asarray(atomic_numbers, int), np.asarray(xyz, float).reshape((-1, 3))


def _force_table(lines: list[str], start: int) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    i = start + 1
    while i < min(len(lines), start + 10):
        if re.match(r"\s*-{5,}\s*$", lines[i]):
            i += 1
            break
        i += 1
    pattern = re.compile(r"^\s*(\d+)\s+(\d+)\s+(" + FLOAT + r")\s+(" + FLOAT + r")\s+(" + FLOAT + r")")
    atoms, atomic_numbers, forces = [], [], []
    while i < len(lines):
        if re.match(r"\s*-{5,}\s*$", lines[i]):
            break
        match = pattern.match(lines[i])
        if not match:
            if atoms and lines[i].strip():
                break
        else:
            atoms.append(int(match.group(1)))
            atomic_numbers.append(int(match.group(2)))
            forces.append([_float(match.group(3)), _float(match.group(4)), _float(match.group(5))])
        i += 1
    return i, np.asarray(atoms, int), np.asarray(atomic_numbers, int), np.asarray(forces, float).reshape((-1, 3))


def _last_between(events: list[tuple], low: int, high: int):
    selected = [event for event in events if low <= event[0] < high]
    return selected[-1] if selected else None


def _scan(lines: list[str]):
    markers, geometries, forces, energies, coordinates = [], [], [], [], []
    i = 0
    while i < len(lines):
        line = lines[i]
        match = POINT_PATH.search(line)
        if match:
            markers.append({"line": i, "point": int(match.group(1)), "path": int(match.group(2))})
        else:
            match = POINT_DIRECTION.search(line)
            if match:
                markers.append({"line": i, "point": int(match.group(1)), "path": 1 if match.group(2).upper() == "FORWARD" else 2})
        match = SCF_ENERGY.search(line)
        if match:
            energies.append((i, _float(match.group(1))))
        match = NET_COORDINATE.search(line)
        if match:
            coordinates.append((i, abs(_float(match.group(1)))))
        if re.search(r"(?:Input|Standard|Z-Matrix) orientation:", line, re.I):
            end, atoms, atomic_numbers, xyz = _orientation(lines, i)
            if len(atoms):
                geometries.append((i, atoms, atomic_numbers, xyz))
            i = max(i, end - 1)
        elif re.search(r"Forces\s+\(Hartrees/Bohr\)", line, re.I):
            end, atoms, atomic_numbers, values = _force_table(lines, i)
            if len(atoms):
                forces.append((i, atoms, atomic_numbers, values))
            i = max(i, end - 1)
        i += 1
    return markers, geometries, forces, energies, coordinates


def parse_log(log_path: Path, series: str, direction_sign: int, expected_atoms: int, conversion: float) -> tuple[pd.DataFrame, dict]:
    lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    markers, geometries, force_blocks, energies, coordinates = _scan(lines)
    warnings: list[str] = []
    candidates: list[Point] = []
    for index, marker in enumerate(markers):
        low = marker["line"]
        high = markers[index + 1]["line"] if index + 1 < len(markers) else len(lines)
        previous = markers[index - 1]["line"] if index else 0
        geometry = _last_between(geometries, low, high)
        forces = _last_between(force_blocks, low, high)
        energy = _last_between(energies, low, high)
        coordinate = _last_between(coordinates, low, high)
        # Completion marker: coordinate follows immediately, optimized data precede it.
        if coordinate is not None and coordinate[0] - low <= 10 and index:
            geometry = _last_between(geometries, previous, low + 1)
            forces = _last_between(force_blocks, previous, low + 1)
            energy = _last_between(energies, previous, low + 1)
        elif geometry is None:
            geometry = _last_between(geometries, previous, low + 1)
            forces = forces or _last_between(force_blocks, previous, low + 1)
            energy = energy or _last_between(energies, previous, low + 1)
        if geometry is None:
            if marker["point"] != 0:
                warnings.append(f"line {low + 1}: path {marker['path']} point {marker['point']} has no geometry")
            continue
        _, atoms, atomic_numbers, xyz = geometry
        force_values = np.full_like(xyz, np.nan)
        if forces is not None and np.array_equal(forces[1], atoms):
            force_values = forces[3]
        elif forces is not None:
            warnings.append(f"line {low + 1}: force and geometry atom numbering differ")
        magnitude = coordinate[1] if coordinate else np.nan
        polarity = 1 if marker["path"] == 1 else -1
        signed = direction_sign * polarity * magnitude
        candidates.append(Point(marker["path"], marker["point"], signed, energy[1] if energy else np.nan, atoms, atomic_numbers, xyz, force_values, low + 1))

    unique: dict[tuple[int, int], Point] = {}
    duplicate_count = 0
    for point in candidates:
        key = (point.path_number, point.point_number)
        old = unique.get(key)
        if old is not None:
            duplicate_count += 1
        if old is None or (point.completeness, point.source_line) >= (old.completeness, old.source_line):
            unique[key] = point
    points = list(unique.values())

    # MoriPF contains an unconverged final reverse structure without NET coordinate.
    for path_number in (1, 2):
        path_points = [p for p in points if p.path_number == path_number]
        known = sorted((p.point_number, abs(p.coordinate)) for p in path_points if np.isfinite(p.coordinate))
        if len(known) >= 2:
            changes = np.diff([v for _, v in known]) / np.diff([n for n, _ in known])
            step = float(np.median(changes[-min(5, len(changes)):]))
            last_number, last_coordinate = known[-1]
            for point in path_points:
                if not np.isfinite(point.coordinate):
                    magnitude = last_coordinate + step * (point.point_number - last_number)
                    point.coordinate = direction_sign * (1 if path_number == 1 else -1) * magnitude
                    point.coordinate_estimated = True
                    warnings.append(f"path {path_number} point {point.point_number}: missing IRC coordinate extrapolated as {point.coordinate:.6g}")

    # These LOGs store exact TS energy in the checkpoint, but not TS Cartesian coordinates.
    if points and not any(abs(p.coordinate) < 1e-12 for p in points):
        first = {}
        for path_number in (1, 2):
            subset = [p for p in points if p.path_number == path_number]
            if subset:
                first[path_number] = min(subset, key=lambda p: p.point_number)
        ts_energy = np.nan
        first_marker = min((m["line"] for m in markers), default=len(lines))
        for line in lines[:first_marker + 1]:
            match = CHECKPOINT_ENERGY.search(line)
            if match:
                ts_energy = _float(match.group(1))
        if len(first) == 2 and np.array_equal(first[1].atomic_numbers, first[2].atomic_numbers):
            xyz = (first[1].xyz + first[2].xyz) / 2.0
            points.append(Point(0, 0, 0.0, ts_energy, first[1].atom_numbers, first[1].atomic_numbers, xyz, np.full_like(xyz, np.nan), min(first[1].source_line, first[2].source_line)))
            warnings.append("TS geometry is absent from LOG; Cartesian midpoint of first forward/reverse structures was used (TS energy is from Energy From Chk)")
        else:
            warnings.append("TS could not be reconstructed")
    points.sort(key=lambda point: point.coordinate)

    rows = []
    for point in points:
        if len(point.atom_numbers) != expected_atoms:
            warnings.append(f"path {point.path_number} point {point.point_number}: {len(point.atom_numbers)} atoms, expected {expected_atoms}")
        for position, atom_number in enumerate(point.atom_numbers):
            atomic_number = int(point.atomic_numbers[position])
            rows.append({
                "series": series, "path_number": point.path_number, "point_number": point.point_number,
                "signed_irc_coordinate": point.coordinate, "irc_coordinate_estimated": point.coordinate_estimated,
                "energy_hartree": point.energy, "atom_number": int(atom_number), "atomic_number": atomic_number,
                "element": ELEMENTS[atomic_number] if atomic_number < len(ELEMENTS) else f"Z{atomic_number}",
                "x": point.xyz[position, 0], "y": point.xyz[position, 1], "z": point.xyz[position, 2],
                "fx": point.forces[position, 0], "fy": point.forces[position, 1], "fz": point.forces[position, 2],
                "source_line": point.source_line,
            })
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError(f"No IRC points parsed from {log_path}")
    point_table = frame[["signed_irc_coordinate", "energy_hartree"]].drop_duplicates()
    ts_energy = float(point_table.loc[point_table.signed_irc_coordinate.abs().idxmin(), "energy_hartree"])
    frame["energy_relative_kcal"] = (frame.energy_hartree - ts_energy) * conversion
    counts = frame.groupby(["path_number", "point_number"]).size()
    if not (counts == expected_atoms).all():
        warnings.append("one or more parsed points do not contain the expected atom count")
    if (frame.groupby("atom_number").element.nunique() > 1).any():
        warnings.append("atom number/element identity changes along the path")
    metadata = {
        "series": series, "log_file": str(log_path), "point_count": int(len(counts)),
        "path_counts": {str(k): int(v) for k, v in frame[["path_number", "point_number"]].drop_duplicates().groupby("path_number").size().items()},
        "ts_energy_hartree": ts_energy, "duplicates_removed": duplicate_count,
        "force_missing_points": int(frame.groupby(["path_number", "point_number"])["fx"].apply(lambda values: values.isna().all()).sum()),
        "warnings": warnings,
    }
    return frame, metadata


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_or_parse(log_path: Path, series: str, direction_sign: int, expected_atoms: int, conversion: float, cache_dir: Path, force_reparse: bool = False) -> tuple[pd.DataFrame, dict, bool]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    csv_path = cache_dir / f"{series}_parsed.csv"
    metadata_path = cache_dir / f"{series}_parsed.meta.json"
    fingerprint = {"sha256": _sha256(log_path), "direction_sign": direction_sign, "parser_version": 2}
    if not force_reparse and csv_path.exists() and metadata_path.exists():
        saved = json.loads(metadata_path.read_text(encoding="utf-8"))
        if saved.get("fingerprint") == fingerprint:
            return pd.read_csv(csv_path), saved["metadata"], True
    frame, metadata = parse_log(log_path, series, direction_sign, expected_atoms, conversion)
    frame.to_csv(csv_path, index=False)
    metadata_path.write_text(json.dumps({"fingerprint": fingerprint, "metadata": metadata}, ensure_ascii=False, indent=2), encoding="utf-8")
    return frame, metadata, False
