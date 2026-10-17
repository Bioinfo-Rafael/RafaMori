# Focused Gaussian IRC Analysis Summary

## LOG files and parsed points

### MoriPF

- LOG: `/Users/cls-lab/Git/RafaMori/base_analysis/data/MORIPF_IRC1.LOG`
- IRC points: 65 (paths: {'0': 1, '1': 35, '2': 29})
- TS energy: -1799.6982910000 Hartree
- Duplicate point records removed: 63
- Points without Cartesian force: 1

### RAFA16

- LOG: `/Users/cls-lab/Git/RafaMori/base_analysis/data/RAFA16_IRC_COMPLETE_MERGED.LOG`
- IRC points: 71 (paths: {'0': 1, '1': 35, '2': 35})
- TS energy: -1601.2516332000 Hartree
- Duplicate point records removed: 70
- Points without Cartesian force: 1

## Target coordinates

### Distances

- `reaction_bond`: Reaction N--C distance (MoriPF 45-20 / RAFA16 43-20); COM baseline: `com_amine_coumarin_core`
- `hb_substituent`: Substituent--amine H-bond distance (MoriPF 33-46 / RAFA16 31-44); COM baseline: `com_amine_substituent`
- `hb_coumarin`: Coumarin--amine H-bond distance (MoriPF 19-47 / RAFA16 19-45); COM baseline: `com_amine_coumarin_core`
- `interference`: Coumarin--substituent interference distance (MoriPF 29-21 / RAFA16 29-21); COM baseline: `com_substituent_coumarin_core`

### Angles

- `hb_substituent_angle`: Substituent--amine H-bond angle (MoriPF 33-46-45 / RAFA16 31-44-43)
- `hb_coumarin_angle`: Coumarin--amine H-bond angle (MoriPF 19-47-45 / RAFA16 19-45-43)

## Fragment definitions

### RAFA16

- `amine_side` (13 atoms): 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55
- `coumarin_core` (22 atoms): 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
- `substituent_side` (9 atoms): 23, 24, 25, 26, 27, 28, 29, 30, 31

### MoriPF

- `amine_side` (13 atoms): 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57
- `coumarin_core` (22 atoms): 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
- `substituent_side` (11 atoms): 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33

## COM-baseline correction

For each distance, the reported residual is:

`d_residual(ξ) = [d(ξ) - d(TS)] - [D_COM(ξ) - D_COM(TS)]`

COM is mass-weighted. No COM correction is applied to angles.

| Distance target | COM baseline type | Fragment pair |
|---|---|---|
| `reaction_bond` | `com_amine_coumarin_core` | `amine_side vs coumarin_core` |
| `hb_substituent` | `com_amine_substituent` | `amine_side vs substituent_side` |
| `hb_coumarin` | `com_amine_coumarin_core` | `amine_side vs coumarin_core` |
| `interference` | `com_substituent_coumarin_core` | `substituent_side vs coumarin_core` |

## Energy and Cartesian-force definitions

- Relative energy: `E(ξ) - E(TS)` in kcal/mol.
- Reaction force: `F(ξ) = -dE/dξ`.
- Reaction force constant: `κ(ξ) = dF/dξ`.
- Smoothed derivatives use Savitzky–Golay filtering of energy with parameters from `series_config.yaml`; raw and smoothed values are both saved.
- Two-atom projected generalized force: `0.5 × (F_b - F_a) · unit(r_b - r_a)` in Hartree/Bohr.

## Output CSV files

- `outputs/tables/angles_delta_from_ts_long.csv`
- `outputs/tables/angles_raw_long.csv`
- `outputs/tables/angles_raw_wide.csv`
- `outputs/tables/com_distances_long.csv`
- `outputs/tables/com_distances_wide.csv`
- `outputs/tables/configuration_validation.csv`
- `outputs/tables/distances_com_corrected_long.csv`
- `outputs/tables/distances_com_corrected_wide.csv`
- `outputs/tables/distances_delta_from_ts_long.csv`
- `outputs/tables/distances_raw_long.csv`
- `outputs/tables/distances_raw_wide.csv`
- `outputs/tables/energy_profile.csv`
- `outputs/tables/projected_forces.csv`
- `outputs/tables/reaction_force.csv`
- `outputs/tables/reaction_force_constant.csv`

## Output figures

- `outputs/figures/angles/delta_from_ts/hb_coumarin_angle_delta.png`
- `outputs/figures/angles/delta_from_ts/hb_substituent_angle_delta.png`
- `outputs/figures/angles/raw/hb_coumarin_angle_raw.png`
- `outputs/figures/angles/raw/hb_substituent_angle_raw.png`
- `outputs/figures/com/com_amine_coumarin_core_com.png`
- `outputs/figures/com/com_amine_substituent_com.png`
- `outputs/figures/com/com_substituent_coumarin_core_com.png`
- `outputs/figures/distances/com_corrected/hb_coumarin_com_corrected.png`
- `outputs/figures/distances/com_corrected/hb_substituent_com_corrected.png`
- `outputs/figures/distances/com_corrected/interference_com_corrected.png`
- `outputs/figures/distances/com_corrected/reaction_bond_com_corrected.png`
- `outputs/figures/distances/delta_from_ts/hb_coumarin_delta.png`
- `outputs/figures/distances/delta_from_ts/hb_substituent_delta.png`
- `outputs/figures/distances/delta_from_ts/interference_delta.png`
- `outputs/figures/distances/delta_from_ts/reaction_bond_delta.png`
- `outputs/figures/distances/raw/hb_coumarin_raw.png`
- `outputs/figures/distances/raw/hb_substituent_raw.png`
- `outputs/figures/distances/raw/interference_raw.png`
- `outputs/figures/distances/raw/reaction_bond_raw.png`
- `outputs/figures/energy/reaction_force.png`
- `outputs/figures/energy/reaction_force_constant.png`
- `outputs/figures/energy/relative_energy.png`
- `outputs/figures/forces/projected_force_hb_coumarin.png`
- `outputs/figures/forces/projected_force_hb_substituent.png`
- `outputs/figures/forces/projected_force_interference.png`
- `outputs/figures/forces/projected_force_reaction_bond.png`
- `outputs/figures/summary/angles_all_raw.png`
- `outputs/figures/summary/com_distances.png`
- `outputs/figures/summary/distances_all_com_corrected.png`
- `outputs/figures/summary/distances_all_raw.png`

## Warnings

- MoriPF: path 2 point 29: missing IRC coordinate extrapolated as 3.55425
- MoriPF: TS geometry is absent from LOG; Cartesian midpoint of first forward/reverse structures was used (TS energy is from Energy From Chk)
- RAFA16: TS geometry is absent from LOG; Cartesian midpoint of first forward/reverse structures was used (TS energy is from Energy From Chk)
- Projected Cartesian force is unavailable for 8 target-point rows (kept as NA)
