"""Analysis restricted to configured distances, angles, and COM baselines."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from geometry import angle, com_distance, distance
from plotting import individual_plots, summary_plot


POINT_COLUMNS = ["series", "signed_irc_coordinate", "path_number", "point_number"]


def _series_mapping(atom_map: pd.DataFrame, series: str) -> dict[str, int]:
    mapping = atom_map.dropna(subset=[series]).set_index("alias")[series].astype(int).to_dict()
    return mapping


def validate_configuration(frames: dict[str, pd.DataFrame], atom_map: pd.DataFrame, targets: dict) -> tuple[pd.DataFrame, list[str]]:
    rows, warnings = [], []
    expected_by_alias = atom_map.set_index("alias").expected_element.to_dict()
    for series, frame in frames.items():
        mapping = _series_mapping(atom_map, series)
        available = set(frame.atom_number.unique())
        for alias, atom_number in mapping.items():
            elements = sorted(frame.loc[frame.atom_number == atom_number, "element"].unique())
            observed = "/".join(elements) if elements else "MISSING"
            expected = expected_by_alias.get(alias)
            status = "OK" if elements == [expected] else "MISMATCH"
            if status != "OK":
                warnings.append(f"{series}: {alias} atom {atom_number}, expected {expected}, observed {observed}")
            rows.append({"series": series, "type": "atom_mapping", "name": alias, "atom_numbers": str(atom_number), "observed_elements": observed, "status": status})
        fragments = targets["fragments"].get(series)
        if fragments is None:
            warnings.append(f"{series}: fragment definitions are missing")
            continue
        for name, atom_numbers in fragments.items():
            missing = sorted(set(atom_numbers) - available)
            duplicates = len(atom_numbers) - len(set(atom_numbers))
            status = "OK" if not missing and not duplicates else "INVALID"
            if missing: warnings.append(f"{series}: fragment {name} contains absent atoms {missing}")
            if duplicates: warnings.append(f"{series}: fragment {name} contains {duplicates} duplicate entries")
            elements = frame[frame.atom_number.isin(atom_numbers)][["atom_number", "element"]].drop_duplicates().element.value_counts().to_dict()
            rows.append({"series": series, "type": "fragment", "name": name, "atom_numbers": ",".join(map(str, atom_numbers)), "observed_elements": str(elements), "status": status})
        for baseline, spec in targets["com_baselines"].items():
            a, b = spec["fragments"]
            overlap = sorted(set(fragments[a]) & set(fragments[b]))
            if overlap: warnings.append(f"{series}: COM baseline {baseline} fragments overlap at {overlap}")
    return pd.DataFrame(rows), warnings


def _point_data(point: pd.DataFrame) -> tuple[dict[int, np.ndarray], dict[int, int]]:
    xyz = {int(row.atom_number): np.asarray([row.x, row.y, row.z], float) for row in point.itertuples(index=False)}
    atomic_numbers = {int(row.atom_number): int(row.atomic_number) for row in point.itertuples(index=False)}
    return xyz, atomic_numbers


def _add_delta(frame: pd.DataFrame, value_column: str, delta_column: str) -> pd.DataFrame:
    pieces = []
    for (_, target), group in frame.groupby(["series", "target"], sort=False):
        group = group.sort_values("signed_irc_coordinate").copy()
        ts_index = group.signed_irc_coordinate.abs().idxmin()
        group[delta_column] = group[value_column] - group.loc[ts_index, value_column]
        pieces.append(group)
    return pd.concat(pieces, ignore_index=True)


def _wide(frame: pd.DataFrame, value_column: str) -> pd.DataFrame:
    return frame.pivot_table(index=POINT_COLUMNS, columns="target", values=value_column, aggfunc="first").reset_index().rename_axis(columns=None)


def _atom_title(target_name: str, spec: dict, atom_map: pd.DataFrame) -> str:
    indexed = atom_map.set_index("alias")
    parts = []
    for series in ("RAFA16", "MoriPF"):
        atoms = [int(indexed.loc[alias, series]) for alias in spec["atoms"]]
        parts.append(f"{series} {'-'.join(map(str, atoms))}")
    return f"{spec['label']} | " + " / ".join(parts)


def analyze_targets(frames: dict[str, pd.DataFrame], atom_map: pd.DataFrame, targets: dict, tables_dir: Path, figures_dir: Path) -> dict:
    tables_dir.mkdir(parents=True, exist_ok=True)
    distance_rows, angle_rows, com_rows = [], [], []
    for series, frame in frames.items():
        mapping = _series_mapping(atom_map, series)
        fragments = targets["fragments"][series]
        point_groups = frame.groupby(["signed_irc_coordinate", "path_number", "point_number"], sort=True)
        for (coordinate, path_number, point_number), point in point_groups:
            xyz, atomic_numbers = _point_data(point)
            base = {"series": series, "signed_irc_coordinate": coordinate, "path_number": int(path_number), "point_number": int(point_number)}
            for target_name, spec in targets["distances"].items():
                atoms = [mapping[alias] for alias in spec["atoms"]]
                distance_rows.append({**base, "target": target_name, "label": spec["label"], "atom_aliases": "--".join(spec["atoms"]), "atom_numbers": "--".join(map(str, atoms)), "distance_A": distance(xyz, *atoms), "baseline_com_type": spec["com_baseline"]})
            for target_name, spec in targets["angles"].items():
                atoms = [mapping[alias] for alias in spec["atoms"]]
                angle_rows.append({**base, "target": target_name, "label": spec["label"], "atom_aliases": "-".join(spec["atoms"]), "atom_numbers": "-".join(map(str, atoms)), "angle_degree": angle(xyz, *atoms)})
            for target_name, spec in targets["com_baselines"].items():
                fragment_a, fragment_b = spec["fragments"]
                com_rows.append({**base, "target": target_name, "label": spec["label"], "fragment_a": fragment_a, "fragment_b": fragment_b, "com_distance_A": com_distance(xyz, atomic_numbers, fragments[fragment_a], fragments[fragment_b])})

    distances = _add_delta(pd.DataFrame(distance_rows), "distance_A", "delta_distance_from_ts_A")
    angles = _add_delta(pd.DataFrame(angle_rows), "angle_degree", "delta_angle_from_ts_degree")
    com = _add_delta(pd.DataFrame(com_rows), "com_distance_A", "delta_com_distance_from_ts_A")
    distances.to_csv(tables_dir / "distances_raw_long.csv", index=False)
    _wide(distances, "distance_A").to_csv(tables_dir / "distances_raw_wide.csv", index=False)
    angles.to_csv(tables_dir / "angles_raw_long.csv", index=False)
    _wide(angles, "angle_degree").to_csv(tables_dir / "angles_raw_wide.csv", index=False)
    com.to_csv(tables_dir / "com_distances_long.csv", index=False)
    _wide(com, "com_distance_A").to_csv(tables_dir / "com_distances_wide.csv", index=False)
    distances[POINT_COLUMNS + ["target", "label", "atom_numbers", "delta_distance_from_ts_A"]].to_csv(tables_dir / "distances_delta_from_ts_long.csv", index=False)
    angles[POINT_COLUMNS + ["target", "label", "atom_numbers", "delta_angle_from_ts_degree"]].to_csv(tables_dir / "angles_delta_from_ts_long.csv", index=False)

    corrected = distances.merge(
        com[[*POINT_COLUMNS, "target", "com_distance_A", "delta_com_distance_from_ts_A"]].rename(columns={"target": "baseline_com_type"}),
        on=[*POINT_COLUMNS, "baseline_com_type"], how="left", validate="many_to_one",
    )
    corrected = corrected.rename(columns={"com_distance_A": "baseline_com_distance_A", "delta_com_distance_from_ts_A": "delta_baseline_com_from_ts_A", "distance_A": "raw_distance_A"})
    corrected["residual_com_corrected_distance_A"] = corrected.delta_distance_from_ts_A - corrected.delta_baseline_com_from_ts_A
    corrected.to_csv(tables_dir / "distances_com_corrected_long.csv", index=False)
    _wide(corrected, "residual_com_corrected_distance_A").to_csv(tables_dir / "distances_com_corrected_wide.csv", index=False)

    distance_titles = {name: _atom_title(name, spec, atom_map) for name, spec in targets["distances"].items()}
    angle_titles = {name: _atom_title(name, spec, atom_map) for name, spec in targets["angles"].items()}
    com_titles = {name: spec["label"] for name, spec in targets["com_baselines"].items()}
    figure_paths = []
    figure_paths.append(summary_plot(distances, targets["distances"], "distance_A", "Distance (Å)", figures_dir / "summary" / "distances_all_raw.png", distance_titles))
    figure_paths.append(summary_plot(corrected, targets["distances"], "residual_com_corrected_distance_A", "Residual distance after COM-baseline correction (Å)", figures_dir / "summary" / "distances_all_com_corrected.png", distance_titles))
    figure_paths.extend(individual_plots(distances, targets["distances"], "distance_A", "Distance (Å)", figures_dir / "distances" / "raw", "raw", distance_titles))
    figure_paths.extend(individual_plots(corrected, targets["distances"], "residual_com_corrected_distance_A", "Residual distance after COM-baseline correction (Å)", figures_dir / "distances" / "com_corrected", "com_corrected", distance_titles))
    figure_paths.extend(individual_plots(distances, targets["distances"], "delta_distance_from_ts_A", "Δ distance from TS (Å)", figures_dir / "distances" / "delta_from_ts", "delta", distance_titles))
    figure_paths.append(summary_plot(angles, targets["angles"], "angle_degree", "Angle (degree)", figures_dir / "summary" / "angles_all_raw.png", angle_titles))
    figure_paths.extend(individual_plots(angles, targets["angles"], "angle_degree", "Angle (degree)", figures_dir / "angles" / "raw", "raw", angle_titles))
    figure_paths.extend(individual_plots(angles, targets["angles"], "delta_angle_from_ts_degree", "Δ angle from TS (degree)", figures_dir / "angles" / "delta_from_ts", "delta", angle_titles))
    figure_paths.append(summary_plot(com, targets["com_baselines"], "com_distance_A", "Center-of-mass distance (Å)", figures_dir / "summary" / "com_distances.png", com_titles, columns=3))
    figure_paths.extend(individual_plots(com, targets["com_baselines"], "com_distance_A", "Center-of-mass distance (Å)", figures_dir / "com", "com", com_titles))
    return {"distances": distances, "angles": angles, "com": com, "corrected": corrected, "figures": figure_paths}
