"""Deterministic figures for second-generation Reaction Force Analysis."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline


SYSTEM_COLORS = {"MoriPF": "#4C78A8", "RAFA14": "#54A24B", "RAFA15": "#B279A2", "RAFA16": "#F58518"}
TARGET_DISPLAY_NAMES = {
    "reaction_bond": "Reaction N–C distance",
    "hb_substituent": "Substituent–amine contact",
    "hb_coumarin": "Coumarin–amine contact",
    "interference": "Coumarin–substituent interference",
    "hb_substituent_angle": "Substituent–amine angle",
    "hb_coumarin_angle": "Coumarin–amine angle",
}


def _boundaries(ax, landmarks, xmin: float, xmax: float) -> None:
    alpha, beta, gamma = landmarks.alpha_continuous, landmarks.beta, landmarks.gamma_continuous
    ax.axvline(alpha, color="#C44E52", linestyle="--", linewidth=.9, alpha=.75)
    ax.axvline(beta, color="black", linestyle="--", linewidth=.9, alpha=.75)
    ax.axvline(gamma, color="#C44E52", linestyle="--", linewidth=.9, alpha=.75)


def energy_force_kappa(system: str, profile: pd.DataFrame, landmarks, spline, reactant_energy: float, conversion: float, path: Path) -> None:
    x = profile.R.to_numpy(float); dense = np.linspace(x.min(), x.max(), 2001)
    color = SYSTEM_COLORS[system]
    fig, axes = plt.subplots(3, 1, figsize=(8.2, 10.5), sharex=True)
    y_raw = (profile.E_raw_hartree - reactant_energy) * conversion
    axes[0].plot(x, y_raw, "o-", ms=2.5, lw=.7, color=color, alpha=.45, label="Actual IRC points")
    axes[0].plot(dense, (spline(dense) - reactant_energy) * conversion, color=color, lw=1.7, label="Smoothing spline")
    axes[0].set_ylabel("V(R) − ETS (kcal/mol)")
    axes[1].plot(x, profile.F_raw_kcal_mol_per_R, "o-", ms=2.2, lw=.6, color=color, alpha=.3, label="Raw finite difference")
    axes[1].plot(dense, -spline.derivative(1)(dense) * conversion, color=color, lw=1.7, label="Smooth F(R)")
    axes[1].set_ylabel("F(R) (kcal mol⁻¹ R⁻¹)")
    axes[2].plot(dense, spline.derivative(2)(dense) * conversion, color=color, lw=1.7, label="Smooth κ(R)")
    axes[2].axhline(0, color="grey", lw=.7)
    axes[2].set_ylabel("κ = d²E/dR² (kcal mol⁻¹ R⁻²)"); axes[2].set_xlabel("Oriented IRC coordinate R")
    for ax in axes:
        _boundaries(ax, landmarks, x.min(), x.max()); ax.grid(alpha=.2); ax.legend(fontsize=8, ncol=2)
    fig.suptitle(f"{system}: energy, Reaction Force, and Reaction Force Constant")
    fig.tight_layout(); fig.savefig(path, dpi=190); plt.close(fig)


def energy_force_kappa_comparison(profiles: dict, landmarks: dict, splines: dict, endpoints: dict, conversion: float, path: Path) -> None:
    """Overlay both systems in the same energy/force/kappa three-panel layout."""
    colors = {"MoriPF": "#4C78A8", "RAFA14": "#54A24B", "RAFA15": "#B279A2", "RAFA16": "#F58518"}
    fig, axes = plt.subplots(3, 1, figsize=(8.5, 10.5), sharex=True)
    for system in landmarks:
        profile, spline, marks = profiles[system], splines[system], landmarks[system]
        x = profile.R.to_numpy(float)
        dense = np.linspace(x.min(), x.max(), 2001)
        reactant = endpoints[system]["reactant"].energy_hartree
        color = colors[system]
        axes[0].plot(x, (profile.E_raw_hartree-reactant)*conversion, "o", ms=2.3, alpha=.35, color=color)
        axes[0].plot(dense, (spline(dense)-reactant)*conversion, lw=1.6, color=color, label=system)
        axes[1].plot(x, profile.F_raw_kcal_mol_per_R, "o", ms=2.0, alpha=.22, color=color)
        axes[1].plot(dense, -spline.derivative(1)(dense)*conversion, lw=1.6, color=color, label=system)
        axes[2].plot(dense, spline.derivative(2)(dense)*conversion, lw=1.6, color=color, label=system)
        for ax in axes:
            ax.axvline(marks.alpha_continuous, color=color, ls="--", lw=.8, alpha=.55)
            ax.axvline(marks.gamma_continuous, color=color, ls="--", lw=.8, alpha=.55)
    for ax in axes:
        ax.axvline(0, color="black", ls="--", lw=.9, alpha=.75)
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    axes[0].set_ylabel("ΔE(R) − ΔE‡(MoriPF) (kcal/mol)")
    axes[1].set_ylabel("F(R) (kcal mol⁻¹ R⁻¹)")
    axes[2].axhline(0, color="grey", lw=.7)
    axes[2].set_ylabel("κ = d²E/dR² (kcal mol⁻¹ R⁻²)")
    axes[2].set_xlabel("IRC coordinate R")
    fig.suptitle("Four systems: electronic profile relative to MoriPF barrier; F and κ")
    fig.tight_layout(); fig.savefig(path, dpi=190); plt.close(fig)


def _connector(x0, y0, slope0, x1, y1, slope1):
    curve = CubicHermiteSpline([x0, x1], [y0, y1], [slope0, slope1])
    x = np.linspace(x0, x1, 80)
    return x, curve(x)


def energy_comparison(profiles: dict, landmarks: dict, splines: dict, endpoints: dict, decompositions: pd.DataFrame, conversion: float, path: Path) -> None:
    invalid = not decompositions.composition_consistent.all()
    if invalid:
        fig, (ax, zoom) = plt.subplots(2, 1, figsize=(10, 8.5), gridspec_kw={"height_ratios": [1.35, 1]})
    else:
        fig, ax = plt.subplots(figsize=(10, 6.2)); zoom = None
    colors = {"MoriPF": "#4C78A8", "RAFA14": "#54A24B", "RAFA15": "#B279A2", "RAFA16": "#F58518"}
    for system, profile in profiles.items():
        color = colors.get(system)
        x = profile.R.to_numpy(float); raw = profile.E_raw_hartree.to_numpy(float)
        er = endpoints[system]["reactant"].energy_hartree; ep = endpoints[system]["product"].energy_hartree
        y = (raw - er) * conversion; spacing = float(np.median(np.diff(x)))
        ax.plot(x, y, color=color, marker="o", ms=3, lw=1.5, label=f"{system} computed IRC")
        xr, xp = x.min() - spacing, x.max() + spacing
        ax.scatter([xr, xp], [0, (ep - er) * conversion], s=75, marker="D", color=color, edgecolor="black", zorder=5, label=f"{system} optimized minima")
        row = decompositions[decompositions.system == system].iloc[0]
        if bool(row.composition_consistent):
            left_slope = (y[1] - y[0]) / (x[1] - x[0]); right_slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
            cx, cy = _connector(xr, 0, 0, x.min(), y[0], left_slope)
            ax.plot(cx, cy, color=color, ls="--", alpha=.45, label=f"{system} visual connector (not computed IRC)")
            cx, cy = _connector(x.max(), y[-1], right_slope, xp, (ep - er) * conversion, 0)
            ax.plot(cx, cy, color=color, ls="--", alpha=.45)
        ax.axvline(landmarks[system].alpha_continuous, color=color, ls=":", lw=.8, alpha=.6)
        ax.axvline(landmarks[system].gamma_continuous, color=color, ls=":", lw=.8, alpha=.6)
        ax.text(.01, .88 if system == "MoriPF" else .76, f"{system}: ΔEact={row.delta_E_act_kcal_mol:.2f}; ΔErxn={row.delta_E_rxn_kcal_mol:.2f}", transform=ax.transAxes, fontsize=8, color=color)
        if zoom is not None:
            ts_energy = float(profile.loc[profile.nearest_to_beta, "E_raw_hartree"].iloc[0])
            zoom.plot(x, (raw-ts_energy)*conversion, color=color, marker="o", ms=3, lw=1.4, label=f"{system} computed IRC")
            zoom.axvline(landmarks[system].alpha_continuous, color=color, ls=":", lw=.9, alpha=.7)
            zoom.axvline(landmarks[system].gamma_continuous, color=color, ls=":", lw=.9, alpha=.7)
    ax.axvline(0, color="black", ls="--", lw=.9)
    if not decompositions.composition_consistent.all():
        ax.text(.01, .02, "WARNING: reactant/IRC stoichiometry mismatch; direct energies are not physically comparable and connectors are omitted.", transform=ax.transAxes, color="darkred", fontsize=8, bbox={"facecolor":"white","alpha":.8,"edgecolor":"darkred"})
    ax.set(xlabel="Oriented IRC coordinate R (optimized minima use pseudo-coordinate)", ylabel="Energy relative to each optimized reactant (kcal/mol)", title="Energy profile comparison: optimized reactant = 0")
    ax.grid(alpha=.2); ax.legend(fontsize=7, ncol=2)
    if zoom is not None:
        zoom.axvline(0, color="black", ls="--", lw=.9)
        zoom.set(xlabel="Oriented IRC coordinate R", ylabel="Computed IRC energy relative to TS (kcal/mol)", title="IRC-only zoom for visibility (diagnostic; TS=0)")
        zoom.grid(alpha=.2); zoom.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(path, dpi=190); plt.close(fig)


def activation_bars(decomposition: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8))
    systems = decomposition.system.tolist(); x = np.arange(len(systems)); width=.23
    fields = [("delta_E_act1_kcal_mol", "ΔEact,1"), ("delta_E_act2_kcal_mol", "ΔEact,2"), ("delta_E_act_kcal_mol", "ΔEact")]
    for offset, (field, label) in enumerate(fields):
        values = decomposition[field].to_numpy(float)
        bars = axes[0].bar(x + (offset-1)*width, values, width, label=label)
        axes[0].bar_label(bars, fmt="%.2f", fontsize=7, rotation=90, padding=2)
    axes[0].set_xticks(x, systems); axes[0].set_ylabel("Electronic energy (kcal/mol)"); axes[0].set_title("Activation-energy decomposition"); axes[0].legend(fontsize=8); axes[0].axhline(0,color="black",lw=.7)
    rxn = decomposition.delta_E_rxn_kcal_mol.to_numpy(float)
    bars = axes[1].bar(x, rxn, color=[SYSTEM_COLORS[s] for s in systems])
    axes[1].bar_label(bars, fmt="%.2f", fontsize=8); axes[1].set_xticks(x, systems); axes[1].set_ylabel("ΔErxn (kcal/mol)"); axes[1].set_title("Reaction energy"); axes[1].axhline(0,color="black",lw=.7)
    if not decomposition.composition_consistent.all():
        fig.text(.5, .01, "Stoichiometry mismatch: bars are arithmetic outputs, not physically interpretable barriers.", ha="center", color="darkred", fontsize=8)
    fig.tight_layout(rect=(0,.04,1,1)); fig.savefig(path, dpi=190); plt.close(fig)


def _point_xyz(frame: pd.DataFrame):
    keys = frame[["signed_irc_coordinate", "path_number", "point_number"]].drop_duplicates().sort_values("signed_irc_coordinate")
    atoms = np.sort(frame.atom_number.unique()).astype(int); xyz=[]
    for row in keys.itertuples(index=False):
        point=frame[(frame.path_number==row.path_number)&(frame.point_number==row.point_number)].set_index("atom_number").reindex(atoms)
        xyz.append(point[["x","y","z"]].to_numpy(float))
    return keys.signed_irc_coordinate.to_numpy(float), atoms, np.asarray(xyz)


def _system_mapping(atom_map: pd.DataFrame, system: str) -> dict[str, int]:
    return atom_map.dropna(subset=[system]).set_index("alias")[system].astype(int).to_dict()


def _target_atom_text(atom_map: pd.DataFrame, system: str, aliases: list[str]) -> str:
    mapping = _system_mapping(atom_map, system)
    elements = atom_map.set_index("alias")["expected_element"].to_dict()
    return "–".join(f"{elements[alias]}{mapping[alias]}" for alias in aliases)


def _target_number_text(atom_map: pd.DataFrame, system: str, aliases: list[str]) -> str:
    mapping = _system_mapping(atom_map, system)
    return "-".join(str(mapping[alias]) for alias in aliases)


def _target_values(frame: pd.DataFrame, atom_map: pd.DataFrame, system: str, kind: str, spec: dict) -> tuple[np.ndarray, np.ndarray]:
    r, atoms, xyz = _point_xyz(frame)
    lookup = {a: i for i, a in enumerate(atoms)}
    mapping = _system_mapping(atom_map, system)
    idx = [lookup[mapping[alias]] for alias in spec["atoms"]]
    if kind == "distance":
        values = np.linalg.norm(xyz[:, idx[0]] - xyz[:, idx[1]], axis=1)
    else:
        left = xyz[:, idx[0]] - xyz[:, idx[1]]
        right = xyz[:, idx[2]] - xyz[:, idx[1]]
        denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
        values = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", left, right) / denominator, -1, 1)))
    return r, values


def target_geometry(system: str, frame: pd.DataFrame, atom_map: pd.DataFrame, targets: dict, landmarks, path: Path) -> None:
    fig, axes=plt.subplots(3,2,figsize=(11,11),sharex=True)
    specs=[("distance",n,s) for n,s in targets["distances"].items()]+[("angle",n,s) for n,s in targets["angles"].items()]
    for ax,(kind,name,spec) in zip(axes.flat,specs):
        r, values = _target_values(frame, atom_map, system, kind, spec)
        unit = "Å" if kind == "distance" else "degree"
        atom_text = _target_atom_text(atom_map, system, spec["atoms"])
        ax.plot(r,values,marker="o",ms=2.5,lw=1.3,color=SYSTEM_COLORS[system])
        _boundaries(ax,landmarks,r.min(),r.max())
        ax.set(title=f"{TARGET_DISPLAY_NAMES[name]} | {system} {_target_number_text(atom_map, system, spec['atoms'])}",ylabel=unit); ax.grid(alpha=.2)
    axes[-1,0].set_xlabel("Oriented IRC coordinate R"); axes[-1,1].set_xlabel("Oriented IRC coordinate R")
    fig.suptitle(f"{system}: configured target geometry and Reaction Force regions"); fig.tight_layout(); fig.savefig(path,dpi=190); plt.close(fig)


def target_geometry_comparison(frames: dict, atom_map: pd.DataFrame, targets: dict, landmarks: dict, path: Path) -> None:
    """Overlay MoriPF and RAFA16 in the same 3x2 target-coordinate layout."""
    colors = {"MoriPF": "#4C78A8", "RAFA14": "#54A24B", "RAFA15": "#B279A2", "RAFA16": "#F58518"}
    fig, axes = plt.subplots(3, 2, figsize=(13, 12), sharex=True)
    specs = [("distance", name, spec) for name, spec in targets["distances"].items()] + [("angle", name, spec) for name, spec in targets["angles"].items()]
    for ax, (kind, _, spec) in zip(axes.flat, specs):
        for system in landmarks:
            r, values = _target_values(frames[system], atom_map, system, kind, spec)
            atom_text = _target_atom_text(atom_map, system, spec["atoms"])
            ax.plot(r, values, marker="o", ms=2.2, lw=1.3, color=colors[system], label=f"{system} {atom_text}")
            _boundaries(ax, landmarks[system], r.min(), r.max())
        numbers = "\n".join(f"{system} {_target_number_text(atom_map, system, spec['atoms'])}" for system in landmarks)
        ax.set(title=f"{TARGET_DISPLAY_NAMES[_]}\n{numbers}", ylabel="Å" if kind == "distance" else "degree")
        ax.grid(alpha=.2); ax.legend(fontsize=7)
    axes[-1,0].set_xlabel("IRC coordinate R"); axes[-1,1].set_xlabel("IRC coordinate R")
    fig.suptitle("Configured target geometry: four systems")
    fig.tight_layout(); fig.savefig(path, dpi=190); plt.close(fig)


def geometry_context(system: str, data: pd.DataFrame, landmarks, atom_map: pd.DataFrame, targets: dict, path: Path) -> None:
    r=data.signed_irc_coordinate.to_numpy(float); fig,axes=plt.subplots(3,1,figsize=(9,11),sharex=True)
    panels=[("covalent_bond_rate_A_per_R","target_reaction_bond_rate_A_per_R",["reaction_bond"],"|dr/dR| (Å/R)","Covalent bonds"),
            ("interfragment_rate_A_per_R",None,["hb_substituent","hb_coumarin","interference"],"|dr/dR| (Å/R)","Relevant interfragment contacts"),
            ("angle_rate_degree_per_R",None,["hb_substituent_angle","hb_coumarin_angle"],"|dθ/dR| (degree/R)","Bond-angle distribution vs target angles")]
    for ax,(prefix,_,target_names,ylabel,title) in zip(axes,panels):
        ax.fill_between(r,data[f"{prefix}_mean"]-data[f"{prefix}_sd"],data[f"{prefix}_mean"]+data[f"{prefix}_sd"],alpha=.14,label="Reference mean ± SD")
        ax.plot(r, data[f"{prefix}_mean"], color="#444444", lw=.8, alpha=.8, label="Reference mean")
        for target in target_names:
            unit="degree_per_R" if "angle" in prefix else "A_per_R"
            spec_group = "angles" if target.endswith("_angle") else "distances"
            atom_text = _target_atom_text(atom_map, system, targets[spec_group][target]["atoms"])
            ax.plot(r,data[f"target_{target}_rate_{unit}"],lw=1.3,label=f"{atom_text} ({system})")
        _boundaries(ax,landmarks,r.min(),r.max()); ax.set(ylabel=ylabel,title=title); ax.grid(alpha=.2); ax.legend(fontsize=7,ncol=3)
    axes[-1].set_xlabel("Oriented IRC coordinate R"); fig.suptitle(f"{system}: geometry-change context (separate coordinate classes)"); fig.tight_layout(); fig.savefig(path,dpi=190); plt.close(fig)


def atomic_force_context(system: str, data: pd.DataFrame, landmarks, atom_map: pd.DataFrame, targets: dict, path: Path) -> None:
    summary=data[["R","all_atom_mean","all_atom_sd"]].drop_duplicates().sort_values("R"); r=summary.R.to_numpy(float)
    fig,ax=plt.subplots(figsize=(9,5.5)); ax.fill_between(r,summary.all_atom_mean-summary.all_atom_sd,summary.all_atom_mean+summary.all_atom_sd,alpha=.14,label="Reference mean ± SD")
    ax.plot(r, summary.all_atom_mean, color="#444444", lw=.8, alpha=.8, label="Reference mean")
    colors = {"hb_substituent": "#4C78A8", "hb_coumarin": "#F58518", "interference": "#54A24B"}
    mapping = _system_mapping(atom_map, system)
    elements = atom_map.set_index("alias")["expected_element"].to_dict()
    for target_name in ("hb_substituent", "hb_coumarin", "interference"):
        aliases = targets["distances"][target_name]["atoms"]
        pair_text = _target_atom_text(atom_map, system, aliases)
        for index, alias in enumerate(aliases):
            atom = mapping[alias]
            group = data[data.atom_number == atom].sort_values("R")
            atom_text = f"{elements[alias]}{atom}"
            ax.plot(group.R, group.force_magnitude_hartree_per_bohr, color=colors[target_name], ls="-" if index == 0 else "--", lw=1.4, label=f"{atom_text} in {pair_text}")
    _boundaries(ax,landmarks,r.min(),r.max()); ax.set(xlabel="Oriented IRC coordinate R",ylabel="|Cartesian atomic force| (Hartree/Bohr)",title=f"{system}: Cartesian atomic force magnitude context (not Reaction Force)")
    ax.grid(alpha=.2); ax.legend(fontsize=6,ncol=3); fig.tight_layout(); fig.savefig(path,dpi=190); plt.close(fig)


def dipole_comparison(dipoles: pd.DataFrame, landmarks: dict, path: Path) -> None:
    fig,axes=plt.subplots(2,2,figsize=(12,9),sharey=False); axes=axes.ravel()
    for ax,(system,group) in zip(axes,dipoles.groupby("system",sort=False)):
        group=group.sort_values("R"); ax.plot(group.R,group.dipole_magnitude_debye,marker="o",ms=2.5,lw=1.3,color=SYSTEM_COLORS[system])
        _boundaries(ax,landmarks[system],group.R.min(),group.R.max()); ax.set(xlabel="Oriented IRC coordinate R",ylabel="Dipole magnitude (Debye)",title=system); ax.grid(alpha=.2)
    fig.suptitle("Dipole magnitude along IRC (supporting indicator)"); fig.tight_layout(); fig.savefig(path,dpi=190); plt.close(fig)


def atomic_force_motion_validation(data: pd.DataFrame, landmarks: dict, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12,9), sharey=False); axes=axes.ravel()
    for ax, system in zip(axes, landmarks):
        group = data[data.system == system].sort_values("R")
        ax.plot(group.R, group.reaction_force_hartree_per_R, color=SYSTEM_COLORS[system], lw=1.5, label=f"Reaction Force ({system})")
        ax.plot(group.R, group.sum_atomic_C_hartree_per_R, color="#555555", lw=1.2, ls="--", label="Σ atomic Cᵢ")
        _boundaries(ax, landmarks[system], group.R.min(), group.R.max())
        ax.axhline(0, color="grey", lw=.7)
        ax.set(title=system, xlabel="IRC coordinate R", ylabel="Force-motion quantity (Hartree/R)")
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle("Atomic force–motion chain-rule sanity check")
    fig.tight_layout(); fig.savefig(path, dpi=190); plt.close(fig)


def fragment_contribution_plot(data: pd.DataFrame, landmarks: dict, path: Path) -> None:
    labels = {"coumarin_side": "Coumarin side", "substituent_side": "Substituent side", "amine_side": "Amine side", "reaction_center": "Reaction-center atoms"}
    fig, axes = plt.subplots(2, 2, figsize=(12,9), sharey=True); axes=axes.ravel()
    for ax, system in zip(axes, landmarks):
        subset = data[data.system == system]
        for fragment, group in subset.groupby("fragment", sort=False):
            group = group.sort_values("R")
            ax.plot(group.R, group.C_fragment_hartree_per_R, lw=1.3, label=labels.get(fragment, fragment))
        _boundaries(ax, landmarks[system], subset.R.min(), subset.R.max())
        ax.axhline(0, color="grey", lw=.7)
        ax.set(title=system, xlabel="IRC coordinate R", ylabel="Σ Fᵢ·drᵢ/dR (Hartree/R)")
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle("Fragment force–motion contributions")
    fig.tight_layout(); fig.savefig(path, dpi=190); plt.close(fig)


def force_motion_alignment_plot(data: pd.DataFrame, landmarks: dict, atom_map: pd.DataFrame, targets: dict, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    target_order = ("reaction_bond", "hb_substituent", "hb_coumarin", "interference")
    for ax, target in zip(axes.flat, target_order):
        for system in landmarks:
            group = data[(data.system == system) & (data.target == target)].sort_values("R")
            ax.plot(group.R, group.A_q_hartree_per_R, color=SYSTEM_COLORS[system], marker="o", ms=2.0, lw=1.2, label=f"{group.atom_label.iloc[0]} ({system})")
            _boundaries(ax, landmarks[system], group.R.min(), group.R.max())
        numbers = "\n".join(f"{system} {_target_number_text(atom_map, system, targets['distances'][target]['atoms'])}" for system in landmarks)
        ax.set(title=f"{TARGET_DISPLAY_NAMES[target]}\n{numbers}", ylabel="Fproj,q × dq/dR (Hartree/R)")
        ax.axhline(0, color="grey", lw=.7); ax.grid(alpha=.2); ax.legend(fontsize=7)
    axes[-1, 0].set_xlabel("IRC coordinate R"); axes[-1, 1].set_xlabel("IRC coordinate R")
    fig.suptitle("Force–motion alignment of configured distances")
    fig.tight_layout(); fig.savefig(path, dpi=190); plt.close(fig)
