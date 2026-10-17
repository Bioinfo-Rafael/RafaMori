# Reaction Force Analysis (Analysis Generation 2)

## Purpose

This analysis compares MoriPF and RAFA16 using a single, oriented reactant → transition state → product IRC coordinate. It asks whether differences in the electronic barrier can be separated into structural-preparation and transition contributions, and whether the configured reaction bond, hydrogen-bond/steric distances, and angles change unusually rapidly relative to appropriate whole-molecule coordinate populations.

The original focused-coordinate workflow is preserved under `../base_analysis/`. This directory is an independent second analysis generation and writes only under `analysis2_force/outputs/`.

## Critical input warning for the supplied files

The supplied optimized reactant LOGs and the IRC do **not** have the same stoichiometry:

| System | Reactant LOG | IRC / product LOG |
|---|---:|---:|
| MoriPF | C20H14F5NO4, 44 atoms | C23H23F5N2O4, 57 atoms |
| RAFA16 | C20H16F3NO4, 44 atoms | C23H25F3N2O4, 57 atoms |

The missing 13 atoms are the amine-side component. Consequently, direct subtraction of the specified reactant electronic energy from the 57-atom IRC energy yields approximately −109,000 kcal/mol and is **not a physically meaningful activation barrier or reaction energy**. The pipeline follows the requested arithmetic, writes `composition_consistent=false` and `valid_for_physical_interpretation=false`, omits visual endpoint connectors between unlike systems, and does not use these values to claim that either system has a lower structural-preparation cost. A same-stoichiometry optimized reactant complex, or the separately calculated missing-component electronic energy at matching theory/solvent conditions, is required for that comparison.

Endpoint orientation RMSDs are reported only as a diagnostic because the reactant comparison has just 44 common atoms and competing scores are close. The explicit project convention takes precedence: MoriPF retains the established reactant-left analysis coordinate, and RAFA16 uses `R_new=-R_original`. In the parser stream both are represented by `coordinate_multiplier: -1`. This is applied before every sort, derivative, statistic, CSV, and figure. See `README_ANALYSIS.md` for the complete reproducibility record.

## Running

From the project root:

```powershell
python analysis2_force\run_analysis.py
```

Available endpoint overrides do not hard-code absolute paths:

```powershell
python analysis2_force\run_analysis.py `
  --data-dir C:\path\to\logs `
  --reactant-log MoriPF=C:\path\MORIPF.LOG `
  --reactant-log RAFA16=C:\path\RAFA16.LOG `
  --product-log MoriPF=C:\path\MORIPF_LYS.LOG `
  --product-log RAFA16=C:\path\RAFA16_LYS.LOG `
  --force-reparse
```

Without overrides, endpoint LOG names come from `config/analysis2_config.yaml` and are searched case-insensitively in `--data-dir` (default `~/Downloads`). Shared atom mapping, target distances/angles, and fragments are read directly from `base_analysis/config`; no target atom indices are re-guessed.

## Theory

For an IRC path $\mathbf X(R)$ and the electronic potential energy evaluated along it, $V(R)=V(\mathbf X(R))$, the reaction force is

$$
F(R)=-\frac{dV(R)}{dR}.
$$

Equivalently,

$$
\frac{dV}{dR}=\nabla V\cdot\frac{d\mathbf X}{dR}.
$$

The Reaction Force is one scalar derivative of the system energy along the IRC. It is not a Cartesian force on an atom. The Reaction Force Constant is

$$
\kappa(R)=\frac{dF(R)}{dR}=-\frac{d^2V(R)}{dR^2}.
$$

The activation decomposition is

$$
\Delta E_{\mathrm{act},1}=V(\alpha)-V(R_{\mathrm{reactant}}),
$$

$$
\Delta E_{\mathrm{act},2}=V(\beta)-V(\alpha),
$$

$$
\Delta E_{\mathrm{act}}=\Delta E_{\mathrm{act},1}+\Delta E_{\mathrm{act},2},
$$

and the electronic reaction energy is

$$
\Delta E_{\mathrm{rxn}}=E_{\mathrm{product}}-E_{\mathrm{reactant}}.
$$

This is an electronic potential-energy analysis, not a Gibbs free-energy decomposition.

For a geometry coordinate $q_j(R)$, the context analysis uses its velocity magnitude

$$
v_j(R)=\left|\frac{dq_j}{dR}\right|.
$$

Covalent-bond distances, relevant interfragment distances, and angles have different units and are always analyzed as three separate populations.

## Alpha, beta, and gamma

- $\alpha$: continuous minimum of the smoothed $F(R)$ on the reactant side.
- $\beta$: energy maximum / TS, fixed at $R=0$.
- $\gamma$: continuous maximum of the smoothed $F(R)$ on the product side.

The intended order is $\alpha<\beta<\gamma$. At interior extrema, $\kappa(\alpha)\approx0$ and $\kappa(\gamma)\approx0$.

## Physical interpretation

Following Toro-Labbé and co-workers:

- Reactant → $\alpha$: structurally intensive preparation of the reactant.
- $\alpha$ → $\gamma$: chemically active / electronically intensive transition region.
- $\gamma$ → Product: structurally intensive relaxation.

The geometry distributions and atomic-force/dipole plots are additions made in this project. They are contextual diagnostics, not quantities defined as parts of the Reaction Force in the cited paper.

## Numerical procedure

1. The preserved base Gaussian parser reconstructs merged/restarted `(path_number, point_number)` records, energies, geometries, and Cartesian forces.
2. Each endpoint energy is the final `SCF Done` before the last `Optimization completed` marker; its geometry is the corresponding final orientation. A normal termination after that optimization is required.
3. Functional/basis, solvent, charge, multiplicity, stoichiometry, and atom count are compared with the IRC. Mismatches are warnings.
4. Each IRC terminal is Kabsch-aligned to both optimized endpoints as a diagnostic. The analysis coordinate follows the configured reactant-left convention; RAFA16 is explicitly multiplied by `-1` before all calculations.
5. `scipy.interpolate.make_smoothing_spline` fits only the computed IRC energy points. With `smoothing_spline_lambda: null`, SciPy selects smoothing strength by generalized cross-validation (GCV).
6. $F=-V'$ and $\kappa=F'=-V''$ are analytic derivatives of that spline. Raw nonuniform-grid finite differences are also calculated.
7. Continuous bounded optimization finds $\alpha$ and $\gamma$; nearest actual IRC points and raw-force extrema are recorded independently.
8. The configured tolerance checks ordering, raw/smooth extremum displacement, $\kappa$ near zero at extrema, and exact activation-decomposition closure.
9. The TS Cartesian geometry is absent from the supplied IRC LOGs and is reconstructed by the base parser as the midpoint of first forward/reverse points. Geometry velocities adjacent to TS are therefore less reliable; electronic energies are not reconstructed.

### Optimized minima and visual connectors

When same-stoichiometry endpoints are available, optimized minima are placed one median IRC step outside the terminal IRC coordinates and connected only across the uncomputed gap by a light dashed cubic Hermite curve. The stationary minimum slope is set to zero; the terminal slope follows the nearest IRC segment.

> Optimized stationary points and terminal IRC structures have no computed intermediate IRC points. The Hermite connection is for visualization only and is never used for activation energies, extrema, Reaction Force, or Reaction Force Constant.

For the currently supplied inconsistent reactants, connectors are omitted rather than implying that unlike compositions lie on one potential-energy surface.

## Geometry-change populations

### Covalent bonds

A bond is included if its distance falls below `factor × (r_cov,i + r_cov,j) + tolerance` at any IRC point. Before mean/SD and percentiles are calculated, the separately plotted target reaction bond is removed from the covalent reference population. The three plotted H-bond/interference pairs are likewise removed from the interfragment reference population; configured target angles are removed if present in the inferred angle population. At each $R$, reference mean±SD, target velocity, and target percentile are saved. Median/IQR output was removed from this analysis generation.

### Interfragment contacts

Candidate pairs are constructed only across configured fragment pairs (`amine_side/coumarin_core`, `amine_side/substituent_side`, and `substituent_side/coumarin_core`). A pair is retained if it is within the configurable 5 Å cutoff at any IRC point. H-bond and steric-contact velocities and percentiles are compared only with this population.

### Angles

All angles from two connected covalent bonds form a degree/R population. Configured H-bond angles are overlaid and percentiled against this angle-only population. Distances and angles are never mixed.

## Atomic Cartesian force and dipole

The auxiliary atomic-force CSV contains

$$
|\mathbf F_i|=\sqrt{F_{ix}^2+F_{iy}^2+F_{iz}^2}
$$

with reference-atom mean±SD. The six atoms belonging to the three plotted H-bond/interference pairs are excluded before this mean and SD are calculated. The plotted target traces are the two endpoint atoms of each pair; the two atoms share a color and use solid/dashed line styles. Median/IQR output was removed. This is explicitly a Cartesian atomic force magnitude, not the scalar Reaction Force. No claim of an atomic contribution to Reaction Force is made, and a mass-weighted IRC-tangent projection is intentionally not implemented because rigorous normalization requires additional care.

Dipole magnitude is

$$
|\boldsymbol\mu|=\sqrt{\mu_x^2+\mu_y^2+\mu_z^2}.
$$

It is shown only as a supporting indicator; dipole change is not identified with electronic rearrangement itself.

## Figures

- `01_*_energy_force_kappa.png`: energy, scalar Reaction Force, and $\kappa=dF/dR=-d^2E/dR^2$.
- `01_energy_force_kappa_comparison.png`: MoriPF and sign-reversed RAFA16 overlaid in the same three-panel layout.
- `02_energy_profile_comparison_reactant_zero.png`: requested reactant-zero arithmetic plus an IRC-only diagnostic zoom when composition is inconsistent.
- `03_activation_decomposition.png`: $\Delta E_{act,1}$, $\Delta E_{act,2}$, total, and reaction-energy comparison.
- `04_*_target_geometry_with_regions.png`: inherited 4 distances and 2 angles with atom-number labels and unlabeled $\alpha/\beta/\gamma$ boundary lines.
- `04_target_geometry_comparison.png`: MoriPF and sign-reversed RAFA16 overlaid in the same 3×2 target-coordinate layout.
- `05_*_geometry_change_context.png`: separate bond, interfragment-contact, and angle velocity populations.
- `06_*_atomic_force_context.png`: Cartesian atomic force magnitudes and target atoms.
- `07_dipole_comparison.png`: dipole magnitude supporting indicator.
- `08_atomic_force_motion_validation.png`: $\sum_i \mathbf F_i\cdot d\mathbf r_i/dR$ versus energy-derived Reaction Force.
- `09_fragment_force_motion_contribution.png`: configured fragment sums of the Cartesian force–motion product.
- `10_force_motion_alignment.png`: signed projected-force × signed distance-motion diagnostics.

## Important caveats

- $F(R)$ is not an atomic force; there is no meaningful “all-atom Reaction Force distribution.”
- $\kappa(R)$ is also one path scalar; an all-atom $\kappa$ distribution is not defined.
- The current reactant/IRC stoichiometry mismatch invalidates direct activation and reaction energies.
- Orientation assignment is low-confidence because the two endpoint RMSD scores are close and reactant comparison uses only 44 common atoms.
- Numerical differentiation and spline smoothing affect extrema. Raw/smoothed landmarks are both retained.
- The TS geometry is reconstructed, so geometry rates around $R=0$ may show artifacts.
- Endpoint connectors, when permitted, are visual only.
- This is electronic energy, not Gibbs energy.

## Tests

Run after analysis generation:

```powershell
python -m unittest discover -s analysis2_force\tests -v
python src\run_analysis.py
```

Tests cover optimized endpoint extraction, orientation, TS placement, derivative identities, extrema ordering/curvature, activation closure, endpoint relative energies, population separation, and base-analysis compatibility.

## Reference

Alejandro Toro-Labbé, Soledad Gutiérrez-Oliva, Jane S. Murray, and Peter Politzer, “The reaction force and the transition region of a reaction,” *Journal of Molecular Modeling* **15** (2009), 707–710. [https://doi.org/10.1007/s00894-008-0431-8](https://doi.org/10.1007/s00894-008-0431-8).
