"""Entry point for focused IRC target-coordinate analysis."""
from __future__ import annotations

import argparse
import shutil
import sys
import traceback
from pathlib import Path

import pandas as pd
import yaml

from energy_force_analysis import run_energy_analysis, run_projected_forces
from parse_gaussian_irc import load_or_parse
from targets_analysis import analyze_targets, validate_configuration


ROOT = Path(__file__).resolve().parents[1]


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze configured Gaussian IRC distances, angles, and COM baselines")
    parser.add_argument("--force-reparse", action="store_true", help="ignore parsed CSV caches")
    return parser.parse_args()


def _path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def _ensure_log(target: Path, warnings: list[str]) -> None:
    source = Path.home() / "Downloads" / target.name
    if target.exists():
        if source.exists() and source.resolve() != target.resolve():
            warnings.append(f"Both data and Downloads contain {target.name}; existing data file was retained without overwrite")
        return
    if source.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(target))
        print(f"Moved {source} -> {target}")


def _report(metadata: dict, atom_map: pd.DataFrame, targets: dict, warnings: list[str], table_paths: list[Path], figure_paths: list[Path], report_path: Path) -> None:
    lines = ["# Focused Gaussian IRC Analysis Summary", "", "## LOG files and parsed points", ""]
    for series, info in metadata.items():
        lines.extend([
            f"### {series}", "",
            f"- LOG: `{info['log_file']}`",
            f"- IRC points: {info['point_count']} (paths: {info['path_counts']})",
            f"- TS energy: {info['ts_energy_hartree']:.10f} Hartree",
            f"- Duplicate point records removed: {info['duplicates_removed']}",
            f"- Points without Cartesian force: {info.get('force_missing_points', 'not recorded')}", "",
        ])
    lines.extend(["## Target coordinates", ""])
    indexed = atom_map.set_index("alias")
    for category in ("distances", "angles"):
        lines.extend([f"### {category.title()}", ""])
        for name, spec in targets[category].items():
            mappings = []
            for series in metadata:
                numbers = [int(indexed.loc[alias, series]) for alias in spec["atoms"]]
                mappings.append(f"{series} {'-'.join(map(str, numbers))}")
            baseline = f"; COM baseline: `{spec['com_baseline']}`" if category == "distances" else ""
            lines.append(f"- `{name}`: {spec['label']} ({' / '.join(mappings)}){baseline}")
        lines.append("")
    lines.extend(["## Fragment definitions", ""])
    for series, fragments in targets["fragments"].items():
        lines.extend([f"### {series}", ""])
        for name, atoms in fragments.items():
            lines.append(f"- `{name}` ({len(atoms)} atoms): {', '.join(map(str, atoms))}")
        lines.append("")
    lines.extend([
        "## COM-baseline correction", "",
        "For each distance, the reported residual is:", "",
        "`d_residual(ξ) = [d(ξ) - d(TS)] - [D_COM(ξ) - D_COM(TS)]`", "",
        "COM is mass-weighted. No COM correction is applied to angles.", "",
        "| Distance target | COM baseline type | Fragment pair |", "|---|---|---|",
    ])
    for name, spec in targets["distances"].items():
        baseline = spec["com_baseline"]
        fragment_pair = " vs ".join(targets["com_baselines"][baseline]["fragments"])
        lines.append(f"| `{name}` | `{baseline}` | `{fragment_pair}` |")
    lines.extend([
        "", "## Energy and Cartesian-force definitions", "",
        "- Relative energy: `E(ξ) - E(TS)` in kcal/mol.",
        "- Reaction force: `F(ξ) = -dE/dξ`.",
        "- Reaction force constant: `κ(ξ) = dF/dξ`.",
        "- Smoothed derivatives use Savitzky–Golay filtering of energy with parameters from `series_config.yaml`; raw and smoothed values are both saved.",
        "- Two-atom projected generalized force: `0.5 × (F_b - F_a) · unit(r_b - r_a)` in Hartree/Bohr.",
        "", "## Output CSV files", "",
    ])
    lines.extend(f"- `{path.relative_to(ROOT)}`" for path in table_paths)
    lines.extend(["", "## Output figures", ""])
    lines.extend(f"- `{path.relative_to(ROOT)}`" for path in figure_paths)
    lines.extend(["", "## Warnings", ""])
    lines.extend([f"- {warning}" for warning in warnings] or ["- None"])
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = _arguments()
    series_config = yaml.safe_load((ROOT / "config" / "series_config.yaml").read_text(encoding="utf-8"))
    targets = yaml.safe_load((ROOT / "config" / "analysis_targets.yaml").read_text(encoding="utf-8"))
    atom_map = pd.read_csv(ROOT / "config" / "atom_map.csv")
    outputs = ROOT / "outputs"
    parsed_dir, tables_dir, figures_dir, reports_dir = [outputs / name for name in ("parsed", "tables", "figures", "reports")]
    for directory in (parsed_dir, tables_dir, figures_dir, reports_dir):
        directory.mkdir(parents=True, exist_ok=True)

    warnings, frames, metadata = [], {}, {}
    parser_config = series_config["parser"]
    missing = []
    for series, spec in series_config["series"].items():
        log_path = _path(spec["log_file"])
        _ensure_log(log_path, warnings)
        if not log_path.exists():
            missing.append(log_path)
            continue
        frame, info, cached = load_or_parse(
            log_path, series, int(spec.get("reaction_direction_sign", 1)),
            int(parser_config["expected_atom_count"]), float(parser_config["hartree_to_kcal_mol"]),
            parsed_dir, args.force_reparse,
        )
        frames[series], metadata[series] = frame, info
        warnings.extend(f"{series}: {warning}" for warning in info["warnings"])
        print(f"Parsed {series}: {info['point_count']} points ({'cache' if cached else 'LOG'})")
    if missing:
        print("Missing LOG files:", file=sys.stderr)
        for path in missing: print(f"- {path}", file=sys.stderr)
        return 2

    validation, validation_warnings = validate_configuration(frames, atom_map, targets)
    warnings.extend(validation_warnings)
    validation_path = tables_dir / "configuration_validation.csv"
    validation.to_csv(validation_path, index=False)
    results = analyze_targets(frames, atom_map, targets, tables_dir, figures_dir)
    energy_results = run_energy_analysis(frames, series_config["energy_analysis"], tables_dir, figures_dir)
    projected_forces, force_warnings = run_projected_forces(frames, atom_map, targets, tables_dir, figures_dir)
    warnings.extend(force_warnings)
    table_paths = sorted(tables_dir.glob("*.csv"))
    figure_paths = sorted(figures_dir.rglob("*.png"))
    report_path = reports_dir / "analysis_summary.md"
    _report(metadata, atom_map, targets, warnings, table_paths, figure_paths, report_path)

    print("\nAnalysis completed.\n")
    print("Parsed series:")
    for series, info in metadata.items(): print(f"- {series}: {info['point_count']} points; paths {info['path_counts']}")
    print("\nGenerated tables:")
    for path in table_paths: print(f"- {path.relative_to(ROOT)}")
    print("\nGenerated figures:")
    for path in figure_paths: print(f"- {path.relative_to(ROOT)}")
    print("\nGenerated reports:")
    print(f"- {report_path.relative_to(ROOT)}")
    print("\nWarnings:")
    if warnings:
        for warning in warnings: print(f"- {warning}")
    else:
        print("- None")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise
    except Exception as error:
        print(f"Analysis failed: {error}", file=sys.stderr)
        traceback.print_exc()
        raise SystemExit(1)
