"""Plots for raw, TS-delta, COM, and COM-corrected target coordinates."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def _draw(ax, frame: pd.DataFrame, value_column: str, ylabel: str, title: str) -> None:
    for series, group in frame.groupby("series", sort=False):
        group = group.sort_values("signed_irc_coordinate")
        ax.plot(group.signed_irc_coordinate, group[value_column], marker="o", markersize=2.5, linewidth=1.4, label=series)
    ax.axvline(0, color="black", linestyle="--", linewidth=.9, alpha=.7)
    ax.set(xlabel="Signed IRC coordinate", ylabel=ylabel, title=title)
    ax.grid(alpha=.22)
    ax.legend(fontsize=8)


def individual_plots(frame: pd.DataFrame, targets: dict, value_column: str, ylabel: str, output_dir: Path, suffix: str, titles: dict[str, str]) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for target in targets:
        path = output_dir / f"{target}_{suffix}.png"
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        _draw(ax, frame[frame.target == target], value_column, ylabel, titles[target])
        fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)
        paths.append(path)
    return paths


def summary_plot(frame: pd.DataFrame, targets: dict, value_column: str, ylabel: str, path: Path, titles: dict[str, str], columns: int = 2) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = len(targets)
    rows = (count + columns - 1) // columns
    fig, axes = plt.subplots(rows, columns, figsize=(7.2 * columns, 4.4 * rows), squeeze=False)
    for ax, target in zip(axes.flat, targets):
        _draw(ax, frame[frame.target == target], value_column, ylabel, titles[target])
    for ax in list(axes.flat)[count:]:
        ax.axis("off")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)
    return path
