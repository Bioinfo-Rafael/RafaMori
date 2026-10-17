"""Geometry-change and Cartesian atomic-force context distributions."""
from __future__ import annotations

import itertools
from collections import defaultdict

import numpy as np
import pandas as pd


COVALENT_RADII = {1: .31, 5: .84, 6: .76, 7: .71, 8: .66, 9: .57, 14: 1.11, 15: 1.07, 16: 1.05, 17: 1.02, 35: 1.20, 53: 1.39}


def _trajectory(frame: pd.DataFrame):
    keys = frame[["signed_irc_coordinate", "path_number", "point_number"]].drop_duplicates().sort_values("signed_irc_coordinate")
    atom_numbers = np.sort(frame.atom_number.unique()).astype(int)
    xyz, forces = [], []
    atomic_numbers = None
    for row in keys.itertuples(index=False):
        point = frame[(frame.path_number == row.path_number) & (frame.point_number == row.point_number)].set_index("atom_number").reindex(atom_numbers)
        xyz.append(point[["x", "y", "z"]].to_numpy(float))
        forces.append(point[["fx", "fy", "fz"]].to_numpy(float))
        if atomic_numbers is None:
            atomic_numbers = point.atomic_number.to_numpy(int)
    return keys.reset_index(drop=True), atom_numbers, np.asarray(atomic_numbers), np.asarray(xyz), np.asarray(forces)


def _distance_values(xyz: np.ndarray, pair: tuple[int, int]) -> np.ndarray:
    return np.linalg.norm(xyz[:, pair[0]] - xyz[:, pair[1]], axis=1)


def _angle_values(xyz: np.ndarray, triple: tuple[int, int, int]) -> np.ndarray:
    i, j, k = triple
    a, b = xyz[:, i] - xyz[:, j], xyz[:, k] - xyz[:, j]
    denominator = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
    cosine = np.divide(np.einsum("ij,ij->i", a, b), denominator, out=np.full(len(xyz), np.nan), where=denominator > 0)
    return np.degrees(np.arccos(np.clip(cosine, -1, 1)))


def _derivative(values: np.ndarray, coordinate: np.ndarray) -> np.ndarray:
    return np.gradient(values, coordinate, edge_order=2 if len(values) >= 3 else 1)


def _stats(matrix: np.ndarray, prefix: str) -> dict[str, np.ndarray]:
    return {
        f"{prefix}_mean": np.nanmean(matrix, axis=1), f"{prefix}_sd": np.nanstd(matrix, axis=1),
    }


def _percentile(matrix: np.ndarray, values: np.ndarray) -> np.ndarray:
    return np.asarray([100.0 * np.nanmean(row <= value) for row, value in zip(matrix, values)])


def _infer_bonds(xyz: np.ndarray, atomic_numbers: np.ndarray, factor: float, extra: float) -> list[tuple[int, int]]:
    bonds = []
    for i, j in itertools.combinations(range(len(atomic_numbers)), 2):
        cutoff = factor * (COVALENT_RADII.get(int(atomic_numbers[i]), .77) + COVALENT_RADII.get(int(atomic_numbers[j]), .77)) + extra
        if np.nanmin(_distance_values(xyz, (i, j))) <= cutoff:
            bonds.append((i, j))
    return bonds


def _bond_angles(bonds: list[tuple[int, int]]) -> list[tuple[int, int, int]]:
    neighbors: dict[int, set[int]] = defaultdict(set)
    for i, j in bonds:
        neighbors[i].add(j); neighbors[j].add(i)
    return [(i, j, k) for j in sorted(neighbors) for i, k in itertools.combinations(sorted(neighbors[j]), 2)]


def _target_indices(atom_numbers: np.ndarray, mapping: dict[str, int], aliases: list[str]) -> tuple[int, ...]:
    lookup = {int(number): index for index, number in enumerate(atom_numbers)}
    return tuple(lookup[mapping[alias]] for alias in aliases)


def analyze_geometry_context(system: str, frame: pd.DataFrame, atom_map: pd.DataFrame, targets: dict, config: dict) -> tuple[pd.DataFrame, dict]:
    keys, atom_numbers, atomic_numbers, xyz, _ = _trajectory(frame)
    coordinate = keys.signed_irc_coordinate.to_numpy(float)
    mapping = atom_map.dropna(subset=[system]).set_index("alias")[system].astype(int).to_dict()
    bonds = _infer_bonds(xyz, atomic_numbers, float(config["covalent_radius_factor"]), float(config["covalent_extra_tolerance_A"]))
    bond_rates_all = np.stack([np.abs(_derivative(_distance_values(xyz, pair), coordinate)) for pair in bonds], axis=1)
    angles = _bond_angles(bonds)
    angle_rates_all = np.stack([np.abs(_derivative(_angle_values(xyz, triple), coordinate)) for triple in angles], axis=1)

    fragments = targets["fragments"][system]
    lookup = {int(number): index for index, number in enumerate(atom_numbers)}
    fragment_pairs = set()
    pair_groups = [("amine_side", "coumarin_core"), ("amine_side", "substituent_side"), ("substituent_side", "coumarin_core")]
    cutoff = float(config["interfragment_cutoff_A"])
    for left, right in pair_groups:
        for a in fragments[left]:
            for b in fragments[right]:
                pair = tuple(sorted((lookup[a], lookup[b])))
                if pair[0] != pair[1] and np.nanmin(_distance_values(xyz, pair)) <= cutoff:
                    fragment_pairs.add(pair)
    fragment_pairs = sorted(fragment_pairs)
    contact_rates_all = np.stack([np.abs(_derivative(_distance_values(xyz, pair), coordinate)) for pair in fragment_pairs], axis=1)

    # Target traces are shown separately, so remove them from the reference
    # populations used for mean/SD and percentile calculations.
    reaction_pair = frozenset(_target_indices(atom_numbers, mapping, targets["distances"]["reaction_bond"]["atoms"]))
    target_contact_pairs = {
        frozenset(_target_indices(atom_numbers, mapping, targets["distances"][name]["atoms"]))
        for name in ("hb_substituent", "hb_coumarin", "interference")
    }
    target_angles = {
        (min(indices[0], indices[2]), indices[1], max(indices[0], indices[2]))
        for indices in (
            _target_indices(atom_numbers, mapping, spec["atoms"])
            for spec in targets["angles"].values()
        )
    }
    bond_keep = [frozenset(pair) != reaction_pair for pair in bonds]
    contact_keep = [frozenset(pair) not in target_contact_pairs for pair in fragment_pairs]
    angle_keep = [(min(triple[0], triple[2]), triple[1], max(triple[0], triple[2])) not in target_angles for triple in angles]
    bond_rates = bond_rates_all[:, bond_keep]
    contact_rates = contact_rates_all[:, contact_keep]
    angle_rates = angle_rates_all[:, angle_keep]
    if not bond_rates.shape[1] or not contact_rates.shape[1] or not angle_rates.shape[1]:
        raise ValueError(f"{system}: target exclusion left an empty geometry reference population")

    output = keys.copy()
    output.insert(0, "system", system)
    output["covalent_bond_reference_count"] = bond_rates.shape[1]
    output["interfragment_reference_count"] = contact_rates.shape[1]
    output["angle_reference_count"] = angle_rates.shape[1]
    for name, values in {**_stats(bond_rates, "covalent_bond_rate_A_per_R"), **_stats(contact_rates, "interfragment_rate_A_per_R"), **_stats(angle_rates, "angle_rate_degree_per_R")}.items():
        output[name] = values
    target_series = {"bond": {}, "contact": {}, "angle": {}}
    for target, spec in targets["distances"].items():
        pair = _target_indices(atom_numbers, mapping, spec["atoms"])
        values = np.abs(_derivative(_distance_values(xyz, pair), coordinate))
        population = bond_rates if target == "reaction_bond" else contact_rates
        category = "bond" if target == "reaction_bond" else "contact"
        output[f"target_{target}_rate_A_per_R"] = values
        output[f"target_{target}_percentile"] = _percentile(population, values)
        target_series[category][target] = values
    for target, spec in targets["angles"].items():
        triple = _target_indices(atom_numbers, mapping, spec["atoms"])
        values = np.abs(_derivative(_angle_values(xyz, triple), coordinate))
        output[f"target_{target}_rate_degree_per_R"] = values
        output[f"target_{target}_percentile"] = _percentile(angle_rates, values)
        target_series["angle"][target] = values
    details = {
        "bond_count": len(bonds), "bond_reference_count": bond_rates.shape[1],
        "interfragment_pair_count": len(fragment_pairs), "interfragment_reference_count": contact_rates.shape[1],
        "angle_count": len(angles), "angle_reference_count": angle_rates.shape[1],
        "targets": target_series,
    }
    return output, details


def analyze_atomic_forces(system: str, frame: pd.DataFrame, atom_map: pd.DataFrame, targets: dict) -> pd.DataFrame:
    keys, atom_numbers, _, _, forces = _trajectory(frame)
    magnitudes = np.linalg.norm(forces, axis=2)
    mapping = atom_map.dropna(subset=[system]).set_index("alias")[system].astype(int).to_dict()
    displayed_aliases = {
        alias
        for target_name in ("hb_substituent", "hb_coumarin", "interference")
        for alias in targets["distances"][target_name]["atoms"]
    }
    excluded_atoms = {mapping[alias] for alias in displayed_aliases}
    reference_mask = np.asarray([int(atom) not in excluded_atoms for atom in atom_numbers])
    if not reference_mask.any():
        raise ValueError(f"{system}: target exclusion left an empty atomic-force reference population")
    reverse_aliases: dict[int, list[str]] = defaultdict(list)
    for alias, atom in mapping.items(): reverse_aliases[atom].append(alias)
    rows = []
    for point_index, key in keys.iterrows():
        values = magnitudes[point_index]
        finite = values[reference_mask & np.isfinite(values)]
        mean = float(np.mean(finite)) if len(finite) else np.nan
        sd = float(np.std(finite)) if len(finite) else np.nan
        for atom_index, atom_number in enumerate(atom_numbers):
            rows.append({
                "system": system, "R": key.signed_irc_coordinate, "path_number": int(key.path_number), "point_number": int(key.point_number),
                "atom_number": int(atom_number), "force_magnitude_hartree_per_bohr": values[atom_index],
                "all_atom_mean": mean, "all_atom_sd": sd,
                "reference_atom_count": int(reference_mask.sum()),
                "excluded_from_reference": int(atom_number) in excluded_atoms,
                "target_aliases": ";".join(reverse_aliases.get(int(atom_number), [])), "is_target_atom": int(atom_number) in reverse_aliases,
            })
    return pd.DataFrame(rows)
