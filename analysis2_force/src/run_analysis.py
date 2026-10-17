"""Orchestrate the independent second-generation Reaction Force Analysis."""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

ANALYSIS_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = ANALYSIS_ROOT.parent
BASE_SRC = PROJECT_ROOT / "base_analysis" / "src"
sys.path.insert(0, str(BASE_SRC))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import yaml

from parse_gaussian_irc import load_or_parse
from geometry_context import analyze_atomic_forces, analyze_geometry_context
from force_motion import analyze_force_motion
from io_parsers import parse_irc_dipoles, parse_irc_metadata, parse_optimized_endpoint
from numerics import activation_decomposition, determine_orientation, orient_frame, reaction_force_profile
from plotting import activation_bars, atomic_force_context, atomic_force_motion_validation, dipole_comparison, energy_comparison, energy_force_kappa, energy_force_kappa_comparison, force_motion_alignment_plot, fragment_contribution_plot, geometry_context, target_geometry, target_geometry_comparison


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Reaction Force Analysis for configured Gaussian IRC systems")
    parser.add_argument("--config", type=Path, default=ANALYSIS_ROOT / "config" / "analysis2_config.yaml")
    parser.add_argument("--data-dir", type=Path, default=Path.home() / "Downloads", help="default directory for optimized endpoint LOGs")
    parser.add_argument("--reactant-log", action="append", default=[], metavar="SYSTEM=PATH")
    parser.add_argument("--product-log", action="append", default=[], metavar="SYSTEM=PATH")
    parser.add_argument("--force-reparse", action="store_true")
    return parser.parse_args()


def overrides(values: list[str]) -> dict[str, Path]:
    result = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Endpoint override must be SYSTEM=PATH, received {value!r}")
        system, path = value.split("=", 1); result[system] = Path(path).expanduser().resolve()
    return result


def resolve_config_path(config_path: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def endpoint_path(data_dir: Path, configured_name: str, override: Path | None) -> Path:
    if override is not None: return override
    direct = data_dir / configured_name
    if direct.exists(): return direct.resolve()
    matches = [path for path in data_dir.iterdir() if path.is_file() and path.name.lower() == configured_name.lower()]
    if len(matches) != 1:
        raise FileNotFoundError(f"Could not uniquely locate {configured_name} in {data_dir}")
    return matches[0].resolve()


def composition_consistent(endpoint, irc_meta: dict) -> bool:
    return endpoint.stoichiometry.upper() == irc_meta["stoichiometry"].upper() and len(endpoint.atom_numbers) == 57


def validate_calculation(system: str, reactant, product, irc_meta: dict) -> list[str]:
    warnings = []
    for role, endpoint in (("reactant", reactant), ("product", product)):
        if not endpoint.optimization_completed or not endpoint.normal_termination:
            warnings.append(f"{system} {role}: optimization/normal termination check failed")
        for key in ("method_basis", "solvent", "charge", "multiplicity"):
            if getattr(endpoint, key) != irc_meta[key]:
                warnings.append(f"{system} {role}: {key} mismatch ({getattr(endpoint,key)} vs IRC {irc_meta[key]})")
        if endpoint.stoichiometry.upper() != irc_meta["stoichiometry"].upper():
            warnings.append(f"{system} {role}: stoichiometry mismatch ({endpoint.stoichiometry}, {len(endpoint.atom_numbers)} atoms vs IRC {irc_meta['stoichiometry']}, 57 atoms)")
    return warnings


def main() -> int:
    args = arguments(); config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    atom_map = pd.read_csv(resolve_config_path(config_path, config["shared_base_config"]["atom_map"]))
    targets = yaml.safe_load(resolve_config_path(config_path, config["shared_base_config"]["targets"]).read_text(encoding="utf-8"))
    reactant_overrides, product_overrides = overrides(args.reactant_log), overrides(args.product_log)
    csv_dir, figure_dir, report_dir = [ANALYSIS_ROOT / "outputs" / name for name in ("csv", "figures", "reports")]
    cache_dir = ANALYSIS_ROOT / ".cache" / "parsed"
    for directory in (csv_dir, figure_dir, report_dir, cache_dir): directory.mkdir(parents=True, exist_ok=True)
    conversion = float(config["numerics"]["hartree_to_kcal_mol"])

    frames, profiles, landmarks, splines, endpoints, orientations = {}, {}, {}, {}, {}, {}
    orientation_records = {}
    geometry_results, atomic_results, dipole_pieces, endpoint_rows, population_rows, warning_messages = {}, {}, [], [], [], []
    atomic_contribution_pieces, force_validation_pieces, fragment_pieces, alignment_pieces = [], [], [], []
    decomposition_rows = []
    for system, spec in config["systems"].items():
        irc_path = resolve_config_path(config_path, spec["irc_log"])
        reactant_path = endpoint_path(args.data_dir.expanduser().resolve(), spec["reactant_log"], reactant_overrides.get(system))
        product_path = endpoint_path(args.data_dir.expanduser().resolve(), spec["product_log"], product_overrides.get(system))
        print(f"{system}: parsing optimized endpoints")
        reactant = parse_optimized_endpoint(reactant_path); product = parse_optimized_endpoint(product_path)
        irc_meta = parse_irc_metadata(irc_path)
        warning_messages.extend(validate_calculation(system, reactant, product, irc_meta))
        endpoint_rows.extend([
            {"system": system, "role": "reactant", "log_file": str(reactant.path), "energy_hartree": reactant.energy_hartree, "stoichiometry": reactant.stoichiometry, "atom_count": len(reactant.atom_numbers), "method_basis": reactant.method_basis, "solvent": reactant.solvent, "charge": reactant.charge, "multiplicity": reactant.multiplicity, "normal_termination": reactant.normal_termination},
            {"system": system, "role": "product", "log_file": str(product.path), "energy_hartree": product.energy_hartree, "stoichiometry": product.stoichiometry, "atom_count": len(product.atom_numbers), "method_basis": product.method_basis, "solvent": product.solvent, "charge": product.charge, "multiplicity": product.multiplicity, "normal_termination": product.normal_termination},
        ])
        frame, parse_meta, cached = load_or_parse(irc_path, system, 1, 57, conversion, cache_dir, args.force_reparse)
        warning_messages.extend(f"{system} IRC: {message}" for message in parse_meta["warnings"])
        orientation = determine_orientation(frame, reactant, product)
        if orientation.alternative_score_A - orientation.selected_score_A < 0.10:
            warning_messages.append(f"{system}: IRC orientation assignment is low-confidence; competing endpoint RMSD scores differ by only {orientation.alternative_score_A-orientation.selected_score_A:.4f} Å")
        coordinate_multiplier = float(spec.get("coordinate_multiplier", 1.0))
        if coordinate_multiplier not in (-1.0, 1.0):
            raise ValueError(f"{system}: coordinate_multiplier must be +1 or -1")
        # Explicit project convention takes precedence over the low-confidence
        # endpoint-RMSD diagnostic: MoriPF raw R, RAFA16 raw R multiplied by -1.
        oriented = orient_frame(frame, coordinate_multiplier == -1.0)
        frames[system], endpoints[system], orientations[system] = oriented, {"reactant": reactant, "product": product}, orientation
        effective_flipped = coordinate_multiplier == -1.0
        orientation_values = vars(orientation).copy()
        orientation_values.pop("flipped")
        orientation_records[system] = {
            "flipped": effective_flipped,
            "auto_flipped": bool(orientation.flipped),
            **orientation_values,
            "coordinate_multiplier": coordinate_multiplier,
            "effective_flipped": effective_flipped,
        }
        if bool(orientation.flipped) != effective_flipped:
            warning_messages.append(f"{system}: explicit coordinate convention overrides the endpoint-RMSD orientation diagnostic")
        profile, marks, spline = reaction_force_profile(oriented, conversion, int(config["numerics"]["dense_grid_points"]), config["numerics"]["smoothing_spline_lambda"])
        profile["E_relative_kcal"] = (profile.E_raw_hartree - reactant.energy_hartree) * conversion
        profile["E_smooth_relative_kcal"] = (profile.E_smooth_hartree - reactant.energy_hartree) * conversion
        profiles[system], landmarks[system], splines[system] = profile, marks, spline
        profile.to_csv(csv_dir / f"reaction_force_{system.upper()}.csv", index=False)
        consistent = composition_consistent(reactant, irc_meta) and composition_consistent(product, irc_meta)
        decomposition_rows.append(activation_decomposition(system, reactant, product, profile, marks, spline, conversion, consistent))

        spacing = float(np.median(np.diff(profile.R)))
        if not marks.alpha_continuous < 0 < marks.gamma_continuous: warning_messages.append(f"{system}: alpha < beta < gamma validation failed")
        if abs(marks.alpha_continuous - marks.alpha_raw) > float(config["numerics"]["extrema_raw_tolerance_steps"]) * spacing: warning_messages.append(f"{system}: smooth/raw alpha differ by more than configured tolerance")
        if abs(marks.gamma_continuous - marks.gamma_raw) > float(config["numerics"]["extrema_raw_tolerance_steps"]) * spacing: warning_messages.append(f"{system}: smooth/raw gamma differ by more than configured tolerance")
        scale = max(abs(profile.kappa_smooth_hartree_per_R2).max(), 1e-15)
        tolerance = float(config["numerics"]["kappa_zero_relative_tolerance"])
        if abs(marks.kappa_alpha) > tolerance * scale or abs(marks.kappa_gamma) > tolerance * scale: warning_messages.append(f"{system}: kappa at alpha/gamma is not sufficiently close to zero")

        geometry_data, geometry_detail = analyze_geometry_context(system, oriented, atom_map, targets, config["geometry_context"])
        population_rows.append({"system": system, **{key: value for key, value in geometry_detail.items() if key != "targets"}})
        warning_messages.append(f"{system}: geometry-change rates adjacent to TS use a reconstructed midpoint TS geometry and should be interpreted cautiously")
        atomic_data = analyze_atomic_forces(system, oriented, atom_map, targets)
        atomic_contribution, force_validation, fragment_contribution, force_alignment = analyze_force_motion(
            system, oriented, profile, atom_map, targets, config["force_motion"], conversion
        )
        atomic_contribution.to_csv(csv_dir / f"atomic_reaction_force_contribution_{system.upper()}.csv", index=False)
        atomic_contribution_pieces.append(atomic_contribution)
        force_validation_pieces.append(force_validation)
        fragment_pieces.append(fragment_contribution)
        alignment_pieces.append(force_alignment)
        geometry_data.to_csv(csv_dir / f"geometry_change_{system.upper()}.csv", index=False)
        atomic_data.to_csv(csv_dir / f"atomic_force_{system.upper()}.csv", index=False)
        geometry_results[system], atomic_results[system] = geometry_data, atomic_data

        dipole_map = parse_irc_dipoles(irc_path)
        point_keys = oriented[["signed_irc_coordinate", "path_number", "point_number"]].drop_duplicates()
        for row in point_keys.itertuples(index=False):
            vector = dipole_map.get((int(row.path_number), int(row.point_number)))
            dipole_pieces.append({"system": system, "R": row.signed_irc_coordinate, "path_number": int(row.path_number), "point_number": int(row.point_number), "dipole_x_debye": vector[0] if vector is not None else np.nan, "dipole_y_debye": vector[1] if vector is not None else np.nan, "dipole_z_debye": vector[2] if vector is not None else np.nan, "dipole_magnitude_debye": float(np.linalg.norm(vector)) if vector is not None else np.nan})
        effective = orientation_records[system]["effective_flipped"]
        print(f"{system}: {len(profile)} points, effective coordinate {'flipped' if effective else 'retained'} ({'cache' if cached else 'LOG'})")

    decomposition = pd.DataFrame(decomposition_rows)
    decomposition.to_csv(csv_dir / "activation_decomposition.csv", index=False)
    pd.DataFrame(endpoint_rows).to_csv(csv_dir / "optimized_endpoint_energies.csv", index=False)
    pd.DataFrame(population_rows).to_csv(csv_dir / "geometry_population_summary.csv", index=False)
    orientation_frame = pd.DataFrame([{"system": system, **orientation_records[system]} for system in orientations])
    orientation_frame.to_csv(csv_dir / "irc_orientation.csv", index=False)
    dipoles = pd.DataFrame(dipole_pieces).sort_values(["system", "R"])
    dipoles.to_csv(csv_dir / "dipole_profile.csv", index=False)
    force_validation = pd.concat(force_validation_pieces, ignore_index=True)
    fragment_contributions = pd.concat(fragment_pieces, ignore_index=True)
    force_alignments = pd.concat(alignment_pieces, ignore_index=True)
    force_validation.to_csv(csv_dir / "atomic_reaction_force_validation.csv", index=False)
    validation_summary_rows = []
    for system, group in force_validation.groupby("system", sort=False):
        finite = group.dropna(subset=["sum_atomic_C_hartree_per_R", "reaction_force_hartree_per_R"])
        residual = finite.sum_atomic_C_hartree_per_R - finite.reaction_force_hartree_per_R
        rmse = float(np.sqrt(np.mean(residual**2)))
        force_rms = float(np.sqrt(np.mean(finite.reaction_force_hartree_per_R**2)))
        record = {
            "system": system, "compared_point_count": len(finite),
            "pearson_correlation": float(finite[["sum_atomic_C_hartree_per_R", "reaction_force_hartree_per_R"]].corr().iloc[0, 1]),
            "mean_absolute_error_hartree_per_R": float(np.mean(np.abs(residual))),
            "root_mean_square_error_hartree_per_R": rmse,
            "normalized_rmse_by_reaction_force_rms": rmse / max(force_rms, 1e-15),
            "max_absolute_error_hartree_per_R": float(np.max(np.abs(residual))),
        }
        validation_summary_rows.append(record)
        if record["pearson_correlation"] < .95 or record["normalized_rmse_by_reaction_force_rms"] > .10:
            warning_messages.append(f"{system}: Cartesian force-motion sum does not closely reproduce the energy-derived Reaction Force")
    pd.DataFrame(validation_summary_rows).to_csv(csv_dir / "atomic_reaction_force_validation_summary.csv", index=False)
    fragment_contributions.to_csv(csv_dir / "fragment_force_motion_contribution.csv", index=False)
    force_alignments.to_csv(csv_dir / "force_motion_alignment.csv", index=False)

    for system in profiles:
        energy_force_kappa(system, profiles[system], landmarks[system], splines[system], endpoints[system]["reactant"].energy_hartree, conversion, figure_dir / f"01_{system.upper()}_energy_force_kappa.png")
        target_geometry(system, frames[system], atom_map, targets, landmarks[system], figure_dir / f"04_{system.upper()}_target_geometry_with_regions.png")
        geometry_context(system, geometry_results[system], landmarks[system], atom_map, targets, figure_dir / f"05_{system.upper()}_geometry_change_context.png")
        atomic_force_context(system, atomic_results[system], landmarks[system], atom_map, targets, figure_dir / f"06_{system.upper()}_atomic_force_context.png")
    target_geometry_comparison(frames, atom_map, targets, landmarks, figure_dir / "04_target_geometry_comparison.png")
    energy_force_kappa_comparison(profiles, landmarks, splines, endpoints, conversion, figure_dir / "01_energy_force_kappa_comparison.png")
    energy_comparison(profiles, landmarks, splines, endpoints, decomposition, conversion, figure_dir / "02_energy_profile_comparison_reactant_zero.png")
    activation_bars(decomposition, figure_dir / "03_activation_decomposition.png")
    atomic_force_motion_validation(force_validation, landmarks, figure_dir / "08_atomic_force_motion_validation.png")
    fragment_contribution_plot(fragment_contributions, landmarks, figure_dir / "09_fragment_force_motion_contribution.png")
    force_motion_alignment_plot(force_alignments, landmarks, atom_map, targets, figure_dir / "10_force_motion_alignment.png")
    if dipoles.dipole_magnitude_debye.notna().any(): dipole_comparison(dipoles.dropna(subset=["dipole_magnitude_debye"]), landmarks, figure_dir / "07_dipole_comparison.png")
    else: warning_messages.append("No dipole vectors could be associated with IRC points; Figure 07 skipped")

    summary = {
        "irc_logs": {system: str(resolve_config_path(config_path, spec["irc_log"])) for system, spec in config["systems"].items()},
        "endpoints": endpoint_rows, "orientation": orientation_frame.to_dict("records"),
        "activation": decomposition.to_dict("records"), "warnings": warning_messages,
        "figures": [str(path.relative_to(ANALYSIS_ROOT)) for path in sorted(figure_dir.glob("*.png"))],
        "csv": [str(path.relative_to(ANALYSIS_ROOT)) for path in sorted(csv_dir.glob("*.csv"))],
    }
    (report_dir / "analysis_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    lines = ["# Analysis 2 Execution Summary", "", "## Orientation", "", "```csv", orientation_frame.to_csv(index=False).strip(), "```", "", "## Activation decomposition", "", "```csv", decomposition.to_csv(index=False).strip(), "```", "", "## Warnings", ""]
    lines.extend([f"- {message}" for message in warning_messages] or ["- None"])
    (report_dir / "analysis_summary.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    implementation_files = [
        "analysis2_force/README.md", "analysis2_force/README_ANALYSIS.md", "analysis2_force/config/analysis2_config.yaml", "analysis2_force/run_analysis.py",
        "analysis2_force/src/force_motion.py", "analysis2_force/src/geometry_context.py", "analysis2_force/src/io_parsers.py", "analysis2_force/src/numerics.py",
        "analysis2_force/src/plotting.py", "analysis2_force/src/run_analysis.py", "analysis2_force/tests/test_analysis2.py",
        "base_analysis/README.md", "base_analysis/config/series_config.yaml", "base_analysis/src/energy_force_analysis.py",
    ]
    analysis2_outputs = [path for directory in (csv_dir, figure_dir, report_dir) for path in sorted(directory.glob("*")) if path.is_file() and path.name != "generated_files_manifest.txt"]
    base_outputs = sorted((PROJECT_ROOT / "base_analysis" / "outputs").rglob("*"))
    output_files = [str(path.relative_to(PROJECT_ROOT)) for path in [*analysis2_outputs, *base_outputs] if path.is_file() and not path.name.startswith(".")]
    cache_files = [str(path.relative_to(PROJECT_ROOT)) for path in sorted(cache_dir.glob("*")) if path.is_file()]
    manifest = [
        "Changed implementation/documentation files:", *[f"- {name}" for name in implementation_files],
        "", "Generated output files:", *[f"- {name}" for name in output_files],
        "- analysis2_force/outputs/reports/generated_files_manifest.txt",
        "", "Regenerated internal cache files (not analysis outputs):", *[f"- {name}" for name in cache_files],
    ]
    (report_dir / "generated_files_manifest.txt").write_text("\n".join(manifest)+"\n", encoding="utf-8")

    print("\nAnalysis 2 completed.")
    print(decomposition[["system","delta_E_act1_kcal_mol","delta_E_act2_kcal_mol","delta_E_act_kcal_mol","delta_E_rxn_kcal_mol","valid_for_physical_interpretation"]].to_string(index=False))
    print(f"Figures: {len(list(figure_dir.glob('*.png')))}; CSV: {len(list(csv_dir.glob('*.csv')))}")
    print("Warnings:")
    for message in warning_messages: print(f"- {message}")
    return 0


if __name__ == "__main__":
    try: raise SystemExit(main())
    except KeyboardInterrupt: raise
    except Exception as error:
        print(f"Analysis 2 failed: {error}", file=sys.stderr); traceback.print_exc(); raise SystemExit(1)
