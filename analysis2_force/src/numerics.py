"""Orientation, smoothing-spline reaction force, and activation decomposition."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.interpolate import make_smoothing_spline
from scipy.optimize import minimize_scalar


@dataclass
class OrientationResult:
    flipped: bool
    negative_reactant_rmsd_A: float
    positive_reactant_rmsd_A: float
    negative_product_rmsd_A: float
    positive_product_rmsd_A: float
    selected_score_A: float
    alternative_score_A: float
    reactant_atoms_compared: int
    product_atoms_compared: int


@dataclass
class ForceLandmarks:
    alpha_continuous: float
    alpha_nearest: float
    beta: float
    gamma_continuous: float
    gamma_nearest: float
    alpha_raw: float
    gamma_raw: float
    kappa_alpha: float
    kappa_gamma: float


def kabsch_rmsd(mobile: np.ndarray, reference: np.ndarray) -> float:
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError(f"Kabsch inputs must have matching (N,3) shapes, got {mobile.shape}/{reference.shape}")
    mobile_centered = mobile - mobile.mean(axis=0)
    reference_centered = reference - reference.mean(axis=0)
    u, _, vt = np.linalg.svd(mobile_centered.T @ reference_centered)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    aligned = mobile_centered @ (u @ correction @ vt)
    return float(np.sqrt(np.mean(np.sum((aligned - reference_centered) ** 2, axis=1))))


def determine_orientation(frame: pd.DataFrame, reactant, product) -> OrientationResult:
    endpoints = frame[["signed_irc_coordinate", "path_number", "point_number"]].drop_duplicates().sort_values("signed_irc_coordinate")
    negative_key = endpoints.iloc[0]
    positive_key = endpoints.iloc[-1]

    def geometry(row) -> pd.DataFrame:
        return frame[(frame.path_number == row.path_number) & (frame.point_number == row.point_number)].sort_values("atom_number")

    negative, positive = geometry(negative_key), geometry(positive_key)

    def compare(endpoint: pd.DataFrame, optimized) -> tuple[float, int]:
        opt = {int(a): xyz for a, xyz in zip(optimized.atom_numbers, optimized.xyz)}
        common = [int(a) for a in endpoint.atom_number if int(a) in opt]
        endpoint_xyz = endpoint.set_index("atom_number").loc[common, ["x", "y", "z"]].to_numpy(float)
        optimized_xyz = np.asarray([opt[a] for a in common])
        return kabsch_rmsd(endpoint_xyz, optimized_xyz), len(common)

    nr, nr_count = compare(negative, reactant)
    pr, pr_count = compare(positive, reactant)
    np_, np_count = compare(negative, product)
    pp, pp_count = compare(positive, product)
    normal_score, flipped_score = nr + pp, pr + np_
    return OrientationResult(
        flipped=flipped_score < normal_score,
        negative_reactant_rmsd_A=nr, positive_reactant_rmsd_A=pr,
        negative_product_rmsd_A=np_, positive_product_rmsd_A=pp,
        selected_score_A=min(normal_score, flipped_score), alternative_score_A=max(normal_score, flipped_score),
        reactant_atoms_compared=min(nr_count, pr_count), product_atoms_compared=min(np_count, pp_count),
    )


def orient_frame(frame: pd.DataFrame, flip: bool) -> pd.DataFrame:
    result = frame.copy()
    if flip:
        result["signed_irc_coordinate"] *= -1
    return result.sort_values(["signed_irc_coordinate", "atom_number"]).reset_index(drop=True)


def reaction_force_profile(frame: pd.DataFrame, conversion: float, dense_grid_points: int, smoothing_lambda: float | None) -> tuple[pd.DataFrame, ForceLandmarks, object]:
    point = frame[["series", "signed_irc_coordinate", "path_number", "point_number", "energy_hartree"]].drop_duplicates().sort_values("signed_irc_coordinate").reset_index(drop=True)
    x = point.signed_irc_coordinate.to_numpy(float)
    energy = point.energy_hartree.to_numpy(float)
    if np.any(np.diff(x) <= 0):
        raise ValueError("IRC coordinates must be strictly increasing before spline fitting")
    spline = make_smoothing_spline(x, energy, lam=smoothing_lambda)
    first, second = spline.derivative(1), spline.derivative(2)
    e_smooth = spline(x)
    f_raw_h = -np.gradient(energy, x, edge_order=2)
    f_smooth_h = -first(x)
    # Project convention requested for Analysis 2: κ=dF/dR=-d²E/dR².
    kappa_smooth_h = -second(x)
    beta = float(x[np.argmin(np.abs(x))])
    epsilon = max(np.finfo(float).eps, np.median(np.diff(x)) * 0.05)
    # F=-V'; minimizing F is equivalent to maximizing V', hence minimize -V'.
    alpha_opt = minimize_scalar(lambda value: float(-first(value)), bounds=(float(x.min()), beta - epsilon), method="bounded")
    # gamma maximizes F=-V', so minimize V'.
    gamma_opt = minimize_scalar(lambda value: float(first(value)), bounds=(beta + epsilon, float(x.max())), method="bounded")
    alpha, gamma = float(alpha_opt.x), float(gamma_opt.x)
    negative = np.flatnonzero(x < beta)
    positive = np.flatnonzero(x > beta)
    alpha_raw_index = negative[np.argmin(f_raw_h[negative])]
    gamma_raw_index = positive[np.argmax(f_raw_h[positive])]
    alpha_nearest_index = int(np.argmin(np.abs(x - alpha)))
    gamma_nearest_index = int(np.argmin(np.abs(x - gamma)))
    landmarks = ForceLandmarks(
        alpha_continuous=alpha, alpha_nearest=float(x[alpha_nearest_index]), beta=beta,
        gamma_continuous=gamma, gamma_nearest=float(x[gamma_nearest_index]),
        alpha_raw=float(x[alpha_raw_index]), gamma_raw=float(x[gamma_raw_index]),
        kappa_alpha=float(-second(alpha)), kappa_gamma=float(-second(gamma)),
    )
    point["R"] = x
    point["E_raw_hartree"] = energy
    point["E_smooth_hartree"] = e_smooth
    point["F_raw_hartree_per_R"] = f_raw_h
    point["F_smooth_hartree_per_R"] = f_smooth_h
    point["kappa_smooth_hartree_per_R2"] = kappa_smooth_h
    point["F_raw_kcal_mol_per_R"] = f_raw_h * conversion
    point["F_smooth_kcal_mol_per_R"] = f_smooth_h * conversion
    point["kappa_smooth_kcal_mol_per_R2"] = kappa_smooth_h * conversion
    point["region"] = np.where(x < alpha, "structural_preparation", np.where(x <= gamma, "transition_region", "structural_relaxation"))
    point["nearest_to_alpha"] = np.arange(len(x)) == alpha_nearest_index
    point["nearest_to_beta"] = np.arange(len(x)) == int(np.argmin(np.abs(x - beta)))
    point["nearest_to_gamma"] = np.arange(len(x)) == gamma_nearest_index
    point["R_alpha_continuous"] = landmarks.alpha_continuous
    point["R_alpha_nearest"] = landmarks.alpha_nearest
    point["R_beta"] = landmarks.beta
    point["R_gamma_continuous"] = landmarks.gamma_continuous
    point["R_gamma_nearest"] = landmarks.gamma_nearest
    point["kappa_at_alpha_hartree_per_R2"] = landmarks.kappa_alpha
    point["kappa_at_gamma_hartree_per_R2"] = landmarks.kappa_gamma
    point["E_hartree"] = point["E_raw_hartree"]
    point["E_smooth"] = point["E_smooth_hartree"]
    point["F_raw"] = point["F_raw_hartree_per_R"]
    point["F_smooth"] = point["F_smooth_hartree_per_R"]
    point["kappa"] = point["kappa_smooth_hartree_per_R2"]
    return point, landmarks, spline


def activation_decomposition(system: str, reactant, product, profile: pd.DataFrame, landmarks: ForceLandmarks, spline, conversion: float, composition_consistent: bool) -> dict:
    e_reactant = float(reactant.energy_hartree)
    e_product = float(product.energy_hartree)
    e_alpha = float(spline(landmarks.alpha_continuous))
    nearest = profile.loc[profile.nearest_to_alpha].iloc[0]
    e_alpha_nearest = float(nearest.E_raw_hartree)
    beta_row = profile.loc[profile.nearest_to_beta].iloc[0]
    e_ts = float(beta_row.E_raw_hartree)
    act1, act2 = e_alpha - e_reactant, e_ts - e_alpha
    return {
        "system": system, "E_reactant_hartree": e_reactant, "E_alpha_continuous_hartree": e_alpha,
        "E_alpha_nearest_hartree": e_alpha_nearest, "E_TS_hartree": e_ts, "E_product_hartree": e_product,
        "delta_E_act1_hartree": act1, "delta_E_act2_hartree": act2,
        "delta_E_act_hartree": e_ts - e_reactant, "delta_E_rxn_hartree": e_product - e_reactant,
        "delta_E_act1_kcal_mol": act1 * conversion, "delta_E_act2_kcal_mol": act2 * conversion,
        "delta_E_act_kcal_mol": (e_ts - e_reactant) * conversion,
        "delta_E_rxn_kcal_mol": (e_product - e_reactant) * conversion,
        "R_alpha_continuous": landmarks.alpha_continuous, "R_alpha_nearest": landmarks.alpha_nearest,
        "R_beta": landmarks.beta, "R_gamma_continuous": landmarks.gamma_continuous,
        "R_gamma_nearest": landmarks.gamma_nearest,
        "composition_consistent": composition_consistent,
        "valid_for_physical_interpretation": composition_consistent,
    }
