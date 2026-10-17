"""Cartesian force–motion diagnostics along the oriented IRC trajectory.

These quantities are useful chain-rule checks.  They are not advertised as a
rigorous atomwise decomposition of a mass-weighted Gaussian IRC coordinate.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd


def _trajectory(frame: pd.DataFrame):
    keys = frame[["signed_irc_coordinate", "path_number", "point_number"]].drop_duplicates().sort_values("signed_irc_coordinate").reset_index(drop=True)
    atoms = np.sort(frame.atom_number.unique()).astype(int)
    xyz, forces = [], []
    for key in keys.itertuples(index=False):
        point = frame[(frame.path_number == key.path_number) & (frame.point_number == key.point_number)].set_index("atom_number").reindex(atoms)
        xyz.append(point[["x", "y", "z"]].to_numpy(float))
        forces.append(point[["fx", "fy", "fz"]].to_numpy(float))
    return keys, atoms, np.asarray(xyz), np.asarray(forces)


def _align_to_ts(xyz: np.ndarray, forces: np.ndarray, ts_index: int) -> tuple[np.ndarray, np.ndarray]:
    """Remove rigid translation/rotation and rotate force vectors consistently."""
    reference = xyz[ts_index] - xyz[ts_index].mean(axis=0)
    aligned_xyz = np.empty_like(xyz)
    aligned_forces = np.empty_like(forces)
    for index, (coordinates, vectors) in enumerate(zip(xyz, forces)):
        centered = coordinates - coordinates.mean(axis=0)
        u, _, vt = np.linalg.svd(centered.T @ reference)
        correction = np.eye(3)
        correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
        rotation = u @ correction @ vt
        aligned_xyz[index] = centered @ rotation
        aligned_forces[index] = vectors @ rotation
    return aligned_xyz, aligned_forces


def _mapping(atom_map: pd.DataFrame, system: str) -> tuple[dict[str, int], dict[str, str]]:
    mapping = atom_map.dropna(subset=[system]).set_index("alias")[system].astype(int).to_dict()
    elements = atom_map.set_index("alias")["expected_element"].to_dict()
    return mapping, elements


def _aliases_by_atom(mapping: dict[str, int]) -> dict[int, str]:
    reverse: dict[int, list[str]] = defaultdict(list)
    for alias, atom in mapping.items():
        reverse[atom].append(alias)
    return {atom: ";".join(names) for atom, names in reverse.items()}


def _resolved_fragments(atoms: np.ndarray, mapping: dict[str, int], target_config: dict, config: dict, system: str) -> dict[str, list[int]]:
    resolved: dict[str, list[int]] = {}
    for name, spec in config["fragments"].items():
        if "shared_fragment" in spec:
            resolved[name] = sorted(int(atom) for atom in target_config["fragments"][system][spec["shared_fragment"]])
        elif "aliases" in spec:
            resolved[name] = sorted({mapping[alias] for alias in spec["aliases"]})
        elif "remaining_after" in spec:
            excluded = {atom for fragment in spec["remaining_after"] for atom in resolved[fragment]}
            resolved[name] = sorted(int(atom) for atom in atoms if int(atom) not in excluded)
        else:
            raise ValueError(f"Unsupported fragment definition for {name}: {spec}")
    return resolved


def analyze_force_motion(
    system: str,
    frame: pd.DataFrame,
    profile: pd.DataFrame,
    atom_map: pd.DataFrame,
    targets: dict,
    config: dict,
    conversion: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    keys, atoms, xyz, forces = _trajectory(frame)
    coordinate = keys.signed_irc_coordinate.to_numpy(float)
    ts_index = int(np.argmin(np.abs(coordinate)))
    xyz_aligned, forces_aligned = _align_to_ts(xyz, forces, ts_index)
    bohr_to_angstrom = float(config["bohr_to_angstrom"])
    velocity_bohr = np.gradient(xyz_aligned, coordinate, axis=0, edge_order=2) / bohr_to_angstrom
    contribution = np.einsum("pij,pij->pi", forces_aligned, velocity_bohr)
    speed = np.linalg.norm(velocity_bohr, axis=2)
    mapping, elements = _mapping(atom_map, system)
    reverse_aliases = _aliases_by_atom(mapping)

    force_by_r = profile.set_index("R")["F_smooth_hartree_per_R"].to_dict()
    atomic_rows, summary_rows = [], []
    for point_index, key in keys.iterrows():
        values = contribution[point_index]
        finite = np.isfinite(values)
        contribution_sum = float(np.sum(values[finite])) if finite.any() else np.nan
        reaction_force = float(force_by_r[key.signed_irc_coordinate])
        residual = contribution_sum - reaction_force if np.isfinite(contribution_sum) else np.nan
        relative_error = residual / max(abs(reaction_force), 1e-12) if np.isfinite(residual) else np.nan
        summary_rows.append({
            "system": system, "R": key.signed_irc_coordinate, "path_number": int(key.path_number), "point_number": int(key.point_number),
            "reaction_force_hartree_per_R": reaction_force,
            "sum_atomic_C_hartree_per_R": contribution_sum,
            "residual_hartree_per_R": residual, "relative_error": relative_error,
            "available_atom_count": int(finite.sum()), "total_atom_count": len(atoms),
        })
        for atom_index, atom in enumerate(atoms):
            atomic_rows.append({
                "system": system, "R": key.signed_irc_coordinate, "path_number": int(key.path_number), "point_number": int(key.point_number),
                "atom_number": int(atom), "target_aliases": reverse_aliases.get(int(atom), ""),
                "force_x_aligned_hartree_per_bohr": forces_aligned[point_index, atom_index, 0],
                "force_y_aligned_hartree_per_bohr": forces_aligned[point_index, atom_index, 1],
                "force_z_aligned_hartree_per_bohr": forces_aligned[point_index, atom_index, 2],
                "dx_dR_bohr_per_R": velocity_bohr[point_index, atom_index, 0],
                "dy_dR_bohr_per_R": velocity_bohr[point_index, atom_index, 1],
                "dz_dR_bohr_per_R": velocity_bohr[point_index, atom_index, 2],
                "cartesian_speed_bohr_per_R": speed[point_index, atom_index],
                "C_i_hartree_per_R": values[atom_index], "C_i_kcal_mol_per_R": values[atom_index] * conversion,
                "sum_atomic_C_hartree_per_R": contribution_sum, "reaction_force_hartree_per_R": reaction_force,
            })

    fragments = _resolved_fragments(atoms, mapping, targets, config, system)
    fragment_rows = []
    atom_lookup = {int(atom): index for index, atom in enumerate(atoms)}
    partition_names = set(config.get("disjoint_partition_fragments", []))
    for point_index, key in keys.iterrows():
        for fragment_name, fragment_atoms in fragments.items():
            indices = [atom_lookup[atom] for atom in fragment_atoms]
            values = contribution[point_index, indices]
            finite = np.isfinite(values)
            fragment_rows.append({
                "system": system, "R": key.signed_irc_coordinate, "path_number": int(key.path_number), "point_number": int(key.point_number),
                "fragment": fragment_name, "atom_numbers": " ".join(map(str, fragment_atoms)), "atom_count": len(fragment_atoms),
                "finite_atom_count": int(finite.sum()),
                "C_fragment_hartree_per_R": float(np.sum(values[finite])) if finite.any() else np.nan,
                "C_fragment_kcal_mol_per_R": float(np.sum(values[finite])) * conversion if finite.any() else np.nan,
                "is_disjoint_partition_fragment": fragment_name in partition_names,
            })

    alignment_rows = []
    for target_name, spec in targets["distances"].items():
        alias_a, alias_b = spec["atoms"]
        atom_a, atom_b = mapping[alias_a], mapping[alias_b]
        ia, ib = atom_lookup[atom_a], atom_lookup[atom_b]
        delta = xyz_aligned[:, ib] - xyz_aligned[:, ia]
        distance_angstrom = np.linalg.norm(delta, axis=1)
        unit = delta / distance_angstrom[:, None]
        signed_dq_dR_angstrom = np.gradient(distance_angstrom, coordinate, edge_order=2)
        signed_dq_dR_bohr = signed_dq_dR_angstrom / bohr_to_angstrom
        projected = .5 * np.einsum("pi,pi->p", forces_aligned[:, ib] - forces_aligned[:, ia], unit)
        alignment = projected * signed_dq_dR_bohr
        atom_text = f"{elements[alias_a]}{atom_a}–{elements[alias_b]}{atom_b}"
        for index, key in keys.iterrows():
            alignment_rows.append({
                "system": system, "R": key.signed_irc_coordinate, "path_number": int(key.path_number), "point_number": int(key.point_number),
                "target": target_name, "atom_numbers": f"{atom_a}-{atom_b}", "atom_label": atom_text,
                "q_angstrom": distance_angstrom[index], "signed_dq_dR_angstrom_per_R": signed_dq_dR_angstrom[index],
                "projected_generalized_force_hartree_per_bohr": projected[index],
                "A_q_hartree_per_R": alignment[index], "A_q_kcal_mol_per_R": alignment[index] * conversion,
            })

    return pd.DataFrame(atomic_rows), pd.DataFrame(summary_rows), pd.DataFrame(fragment_rows), pd.DataFrame(alignment_rows)
