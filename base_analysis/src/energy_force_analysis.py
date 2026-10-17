"""Energy derivatives and two-atom Cartesian-force projections."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter


def _point_energy(frame: pd.DataFrame) -> pd.DataFrame:
    columns = ["series", "signed_irc_coordinate", "path_number", "point_number", "energy_hartree", "energy_relative_kcal"]
    return frame[columns].drop_duplicates().sort_values("signed_irc_coordinate").reset_index(drop=True)


def _smooth(values: np.ndarray, config: dict) -> tuple[np.ndarray, int | None]:
    if not config.get("enabled", True):
        return values.copy(), None
    polyorder = int(config.get("polyorder", 3))
    requested = int(config.get("window_length", 7))
    window = min(requested, len(values) if len(values) % 2 else len(values) - 1)
    minimum = polyorder + 2 if (polyorder + 2) % 2 else polyorder + 3
    if window < minimum or not np.isfinite(values).all():
        return values.copy(), None
    return savgol_filter(values, window_length=window, polyorder=min(polyorder, window - 1), mode="interp"), window


def _derivative(values: np.ndarray, coordinate: np.ndarray) -> np.ndarray:
    return np.gradient(values, coordinate, edge_order=2 if len(values) >= 3 else 1)


def _overlay(frame: pd.DataFrame, value_column: str, path: Path, ylabel: str, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for series, group in frame.groupby("series", sort=False):
        group = group.sort_values("signed_irc_coordinate")
        ax.plot(group.signed_irc_coordinate, group[value_column], marker="o", markersize=2.5, linewidth=1.4, label=series)
    ax.axvline(0, color="black", linestyle="--", linewidth=.9, alpha=.7)
    ax.set(xlabel="Signed IRC coordinate", ylabel=ylabel, title=title); ax.grid(alpha=.22); ax.legend(); fig.tight_layout()
    fig.savefig(path, dpi=180); plt.close(fig)


def _raw_smoothed_overlay(frame: pd.DataFrame, raw: str, smoothed: str, path: Path, ylabel: str, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for index, (series, group) in enumerate(frame.groupby("series", sort=False)):
        group = group.sort_values("signed_irc_coordinate")
        color = colors[index % len(colors)]
        ax.plot(group.signed_irc_coordinate, group[raw], color=color, linestyle="--", linewidth=1.0, alpha=.55, label=f"{series} raw")
        ax.plot(group.signed_irc_coordinate, group[smoothed], color=color, marker="o", markersize=2.5, linewidth=1.5, label=f"{series} smoothed")
    ax.axvline(0, color="black", linestyle="--", linewidth=.9, alpha=.7)
    ax.set(xlabel="Signed IRC coordinate", ylabel=ylabel, title=title); ax.grid(alpha=.22); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(path, dpi=180); plt.close(fig)


def run_energy_analysis(frames: dict[str, pd.DataFrame], config: dict, tables_dir: Path, figures_dir: Path) -> dict:
    energy_pieces, force_pieces, constant_pieces = [], [], []
    smoothing = config["smoothing"]
    for series, frame in frames.items():
        energy = _point_energy(frame)
        coordinate = energy.signed_irc_coordinate.to_numpy(float)
        relative = energy.energy_relative_kcal.to_numpy(float)
        smoothed_energy, used_window = _smooth(relative, smoothing)
        force_raw = -_derivative(relative, coordinate)
        force_smoothed = -_derivative(smoothed_energy, coordinate)
        constant_raw = _derivative(force_raw, coordinate)
        constant_smoothed = _derivative(force_smoothed, coordinate)
        energy["energy_relative_kcal_mol"] = relative
        energy_pieces.append(energy.drop(columns="energy_relative_kcal"))
        common = energy[["series", "signed_irc_coordinate", "path_number", "point_number"]].copy()
        force = common.copy()
        force["reaction_force_raw_kcal_mol_per_irc"] = force_raw
        force["reaction_force_smoothed_kcal_mol_per_irc"] = force_smoothed
        force["savgol_window_length"] = used_window
        force["savgol_polyorder"] = int(smoothing.get("polyorder", 3)) if used_window else np.nan
        force_pieces.append(force)
        constant = common.copy()
        constant["reaction_force_constant_raw_kcal_mol_per_irc2"] = constant_raw
        constant["reaction_force_constant_smoothed_kcal_mol_per_irc2"] = constant_smoothed
        constant["savgol_window_length"] = used_window
        constant["savgol_polyorder"] = int(smoothing.get("polyorder", 3)) if used_window else np.nan
        constant_pieces.append(constant)
    energy_profile = pd.concat(energy_pieces, ignore_index=True)
    reaction_force = pd.concat(force_pieces, ignore_index=True)
    force_constant = pd.concat(constant_pieces, ignore_index=True)
    energy_profile.to_csv(tables_dir / "energy_profile.csv", index=False)
    reaction_force.to_csv(tables_dir / "reaction_force.csv", index=False)
    force_constant.to_csv(tables_dir / "reaction_force_constant.csv", index=False)
    _overlay(energy_profile, "energy_relative_kcal_mol", figures_dir / "energy" / "relative_energy.png", "Relative energy (kcal/mol)", "Relative IRC energy")
    _raw_smoothed_overlay(reaction_force, "reaction_force_raw_kcal_mol_per_irc", "reaction_force_smoothed_kcal_mol_per_irc", figures_dir / "energy" / "reaction_force.png", "Reaction force −dE/dξ (kcal mol⁻¹ ξ⁻¹)", "Reaction force")
    _raw_smoothed_overlay(force_constant, "reaction_force_constant_raw_kcal_mol_per_irc2", "reaction_force_constant_smoothed_kcal_mol_per_irc2", figures_dir / "energy" / "reaction_force_constant.png", "Reaction force constant dF/dξ", "Reaction force constant")
    return {"energy_profile": energy_profile, "reaction_force": reaction_force, "force_constant": force_constant}


def run_projected_forces(frames: dict[str, pd.DataFrame], atom_map: pd.DataFrame, targets: dict, tables_dir: Path, figures_dir: Path) -> tuple[pd.DataFrame, list[str]]:
    warnings, rows = [], []
    required = {"fx", "fy", "fz"}
    if any(not required.issubset(frame.columns) for frame in frames.values()):
        return pd.DataFrame(), ["Cartesian force columns are unavailable; projected-force analysis was skipped"]
    for series, frame in frames.items():
        mapping = atom_map.dropna(subset=[series]).set_index("alias")[series].astype(int).to_dict()
        for (coordinate, path_number, point_number), point in frame.groupby(["signed_irc_coordinate", "path_number", "point_number"], sort=True):
            point = point.set_index("atom_number")
            for target, spec in targets["distances"].items():
                atom_a, atom_b = [mapping[alias] for alias in spec["atoms"]]
                xyz_a = point.loc[atom_a, ["x", "y", "z"]].to_numpy(float)
                xyz_b = point.loc[atom_b, ["x", "y", "z"]].to_numpy(float)
                force_a = point.loc[atom_a, ["fx", "fy", "fz"]].to_numpy(float)
                force_b = point.loc[atom_b, ["fx", "fy", "fz"]].to_numpy(float)
                vector = xyz_b - xyz_a
                norm = np.linalg.norm(vector)
                value = np.nan
                if norm > 0 and np.isfinite(force_a).all() and np.isfinite(force_b).all():
                    unit = vector / norm
                    value = 0.5 * float(np.dot(force_b - force_a, unit))
                rows.append({
                    "series": series, "signed_irc_coordinate": coordinate, "path_number": int(path_number), "point_number": int(point_number),
                    "target": target, "label": spec["label"], "atom_aliases": "--".join(spec["atoms"]),
                    "atom_numbers": f"{atom_a}--{atom_b}", "projected_generalized_force_hartree_per_bohr": value,
                    "definition": "0.5 * (F_b - F_a) dot unit(r_b - r_a)",
                })
    result = pd.DataFrame(rows)
    if result.empty or not result.projected_generalized_force_hartree_per_bohr.notna().any():
        return result, ["No finite Cartesian forces were parsed; projected-force analysis was skipped"]
    result.to_csv(tables_dir / "projected_forces.csv", index=False)
    indexed = atom_map.set_index("alias")
    for target in targets["distances"]:
        subset = result[result.target == target]
        if not subset.projected_generalized_force_hartree_per_bohr.notna().any():
            warnings.append(f"{target}: no finite Cartesian force projections; figure skipped")
            continue
        spec = targets["distances"][target]
        number_text = " / ".join(
            f"{series} {'-'.join(str(int(indexed.loc[alias, series])) for alias in spec['atoms'])}"
            for series in ("MoriPF", "RAFA16")
        )
        title = f"Cartesian force projection: {spec['label']}\n{number_text}"
        _overlay(subset, "projected_generalized_force_hartree_per_bohr", figures_dir / "forces" / f"projected_force_{target}.png", "Projected generalized force (Hartree/Bohr)", title)
    missing = int(result.projected_generalized_force_hartree_per_bohr.isna().sum())
    if missing:
        warnings.append(f"Projected Cartesian force is unavailable for {missing} target-point rows (kept as NA)")
    return result, warnings
