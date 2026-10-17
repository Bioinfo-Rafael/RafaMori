"""Input parsing for optimized endpoints and optional IRC dipoles."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from parse_gaussian_irc import FLOAT, NET_COORDINATE, POINT_DIRECTION, POINT_PATH, _orientation


@dataclass
class OptimizedEndpoint:
    path: Path
    energy_hartree: float
    atom_numbers: np.ndarray
    atomic_numbers: np.ndarray
    xyz: np.ndarray
    optimization_completed: bool
    normal_termination: bool
    route: str
    method_basis: str
    solvent: str
    charge: int
    multiplicity: int
    stoichiometry: str
    energy_source_line: int
    geometry_source_line: int


SCF = re.compile(r"SCF Done:.*?=\s*(" + FLOAT + r")", re.I)
CHARGE = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", re.I)
STOICH = re.compile(r"Stoichiometry\s+(\S+)", re.I)
DIPOLE = re.compile(r"Dipole\s*=\s*(" + FLOAT + r")\s*(" + FLOAT + r")\s*(" + FLOAT + r")", re.I)


def _number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _route(lines: list[str]) -> str:
    for i, line in enumerate(lines):
        if re.match(r"\s*#", line):
            parts = [line.strip()]
            for following in lines[i + 1:i + 20]:
                if re.match(r"\s*-{5,}\s*$", following) or not following.strip():
                    break
                parts.append(following.strip())
            return " ".join(parts)
    return ""


def _method_basis(route: str) -> str:
    match = re.search(r"(?:R?M062X)\s*/\s*([^\s]+)", route, re.I)
    if match:
        return "M062X/" + match.group(1).upper()
    method = "M062X" if re.search(r"M062X", route, re.I) else "UNKNOWN"
    basis = "6-311++G(D,P)" if re.search(r"6-311\+\+G\(D,P\)", route, re.I) else "UNKNOWN"
    return f"{method}/{basis}"


def _solvent(route: str) -> str:
    compact = re.sub(r"\s+", "", route)
    if re.search(r"SCRF=.*SMD", compact, re.I) and re.search(r"SOLVENT=WATER", compact, re.I):
        return "SMD(water)"
    return "UNKNOWN"


def parse_optimized_endpoint(path: Path) -> OptimizedEndpoint:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    completed = [i for i, line in enumerate(lines) if "Optimization completed" in line]
    if not completed:
        raise ValueError(f"No 'Optimization completed' marker in {path}")
    completion_line = completed[-1]
    scf_events = []
    orientations = []
    for i, line in enumerate(lines[:completion_line + 1]):
        match = SCF.search(line)
        if match:
            scf_events.append((i, _number(match.group(1))))
        if re.search(r"(?:Input|Standard|Z-Matrix) orientation:", line, re.I):
            orientations.append(i)
    if not scf_events or not orientations:
        raise ValueError(f"Final optimized energy/geometry could not be located in {path}")
    energy_line, energy = scf_events[-1]
    geometry_line = orientations[-1]
    _, atom_numbers, atomic_numbers, xyz = _orientation(lines, geometry_line)
    if not len(atom_numbers):
        raise ValueError(f"Final orientation table is empty in {path}")
    charge, multiplicity = 0, 1
    for line in lines[:completion_line + 1]:
        match = CHARGE.search(line)
        if match:
            charge, multiplicity = int(match.group(1)), int(match.group(2))
            break
    stoichiometry = "UNKNOWN"
    for line in lines[:completion_line + 1]:
        match = STOICH.search(line)
        if match:
            stoichiometry = match.group(1)
            break
    route = _route(lines)
    return OptimizedEndpoint(
        path=path, energy_hartree=energy, atom_numbers=atom_numbers, atomic_numbers=atomic_numbers,
        xyz=xyz, optimization_completed=True,
        normal_termination=any("Normal termination of Gaussian" in line for line in lines[completion_line:]),
        route=route, method_basis=_method_basis(route), solvent=_solvent(route), charge=charge,
        multiplicity=multiplicity, stoichiometry=stoichiometry, energy_source_line=energy_line + 1,
        geometry_source_line=geometry_line + 1,
    )


def parse_irc_metadata(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    route = _route(lines)
    charge, multiplicity = 0, 1
    stoichiometry = "UNKNOWN"
    for line in lines[:1000]:
        match = CHARGE.search(line)
        if match:
            charge, multiplicity = int(match.group(1)), int(match.group(2))
        match = STOICH.search(line)
        if match:
            stoichiometry = match.group(1)
    return {"route": route, "method_basis": _method_basis(route), "solvent": _solvent(route), "charge": charge, "multiplicity": multiplicity, "stoichiometry": stoichiometry}


def parse_irc_dipoles(path: Path) -> dict[tuple[int, int], np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    markers, coordinates, dipoles = [], [], []
    for i, line in enumerate(lines):
        match = POINT_PATH.search(line)
        if match:
            markers.append({"line": i, "point": int(match.group(1)), "path": int(match.group(2))})
        else:
            match = POINT_DIRECTION.search(line)
            if match:
                markers.append({"line": i, "point": int(match.group(1)), "path": 1 if match.group(2).upper() == "FORWARD" else 2})
        match = NET_COORDINATE.search(line)
        if match:
            coordinates.append(i)
        match = DIPOLE.search(line)
        if match:
            dipoles.append((i, np.asarray([_number(match.group(1)), _number(match.group(2)), _number(match.group(3))])))
    result: dict[tuple[int, int], tuple[int, np.ndarray]] = {}
    for index, marker in enumerate(markers):
        low = marker["line"]
        high = markers[index + 1]["line"] if index + 1 < len(markers) else len(lines)
        previous = markers[index - 1]["line"] if index else 0
        immediate_coordinate = next((line for line in coordinates if low <= line < min(high, low + 11)), None)
        search_low, search_high = (previous, low + 1) if immediate_coordinate is not None and index else (low, high)
        candidates = [event for event in dipoles if search_low <= event[0] < search_high]
        if candidates:
            event = candidates[-1]
            key = (marker["path"], marker["point"])
            if key not in result or event[0] > result[key][0]:
                result[key] = event
    return {key: value for key, (_, value) in result.items()}
