from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "base_analysis" / "src"))
sys.path.insert(0, str(PROJECT / "analysis2_force" / "src"))

from io_parsers import parse_optimized_endpoint
from numerics import reaction_force_profile


CSV = PROJECT / "analysis2_force" / "outputs" / "csv"


class EndpointParsingTests(unittest.TestCase):
    def test_final_optimized_reactant_energies(self):
        expected = {"MORIPF.LOG": -1625.24459577, "RAFA16.LOG": -1426.79915417}
        for name, energy in expected.items():
            endpoint = parse_optimized_endpoint(Path.home() / "Downloads" / name)
            self.assertTrue(endpoint.optimization_completed)
            self.assertTrue(endpoint.normal_termination)
            self.assertAlmostEqual(endpoint.energy_hartree, energy, places=8)

    def test_final_optimized_product_energies(self):
        expected = {"MORIPF_LYS.LOG": -1799.70473558, "RAFA16_LYS.LOG": -1601.25653239}
        for name, energy in expected.items():
            endpoint = parse_optimized_endpoint(Path.home() / "Downloads" / name)
            self.assertEqual(len(endpoint.atom_numbers), 57)
            self.assertAlmostEqual(endpoint.energy_hartree, energy, places=8)


class ReactionForceOutputTests(unittest.TestCase):
    def test_orientation_and_ts(self):
        orientation = pd.read_csv(CSV / "irc_orientation.csv").set_index("system")
        self.assertTrue(bool(orientation.loc["MoriPF", "flipped"]))
        self.assertTrue(bool(orientation.loc["MoriPF", "auto_flipped"]))
        self.assertEqual(float(orientation.loc["MoriPF", "coordinate_multiplier"]), -1.0)
        self.assertTrue(bool(orientation.loc["RAFA16", "flipped"]))
        self.assertFalse(bool(orientation.loc["RAFA16", "auto_flipped"]))
        self.assertEqual(float(orientation.loc["RAFA16", "coordinate_multiplier"]), -1.0)
        self.assertTrue(bool(orientation.loc["RAFA16", "effective_flipped"]))
        for system in ("MORIPF", "RAFA16"):
            profile = pd.read_csv(CSV / f"reaction_force_{system}.csv")
            self.assertAlmostEqual(float(profile.loc[profile.nearest_to_beta, "R"].iloc[0]), 0.0, places=10)

    def test_derivative_signs_and_landmarks(self):
        for system in ("MORIPF", "RAFA16"):
            profile = pd.read_csv(CSV / f"reaction_force_{system}.csv")
            x = profile.R.to_numpy(float); energy = profile.E_raw_hartree.to_numpy(float)
            expected_force = -np.gradient(energy, x, edge_order=2)
            np.testing.assert_allclose(expected_force, profile.F_raw_hartree_per_R, rtol=1e-10, atol=1e-12)
            numerical_dF = np.gradient(profile.F_smooth, x, edge_order=2)
            # The supplied LOG omits the TS geometry and reconstructs its R=0
            # point, so exclude that immediate finite-difference neighborhood.
            mask = np.abs(x) > 1.1*np.median(np.diff(x))
            np.testing.assert_allclose(profile.kappa[mask], numerical_dF[mask], rtol=.08, atol=2e-4)
            alpha, beta, gamma = profile[["R_alpha_continuous", "R_beta", "R_gamma_continuous"]].iloc[0]
            self.assertLess(alpha, beta); self.assertLess(beta, gamma)
            scale = max(abs(profile.kappa).max(), 1e-15)
            self.assertLess(abs(profile.kappa_at_alpha_hartree_per_R2.iloc[0]), .12 * scale)
            self.assertLess(abs(profile.kappa_at_gamma_hartree_per_R2.iloc[0]), .12 * scale)

    def test_synthetic_force_and_kappa_identity(self):
        x = np.linspace(-3, 3, 61); energy = -x**2 + .05*x**4
        frame = pd.DataFrame({"series": "test", "signed_irc_coordinate": x, "path_number": np.where(x < 0, 2, 1), "point_number": np.arange(len(x)), "energy_hartree": energy})
        profile, marks, spline = reaction_force_profile(frame, 1.0, 1001, 1e-12)
        # Stay away from spline boundaries where natural-end conditions differ
        # slightly from the unconstrained analytic quartic derivatives.
        dense = np.linspace(-2.6, 2.6, 200)
        np.testing.assert_allclose(-spline.derivative(1)(dense), 2*dense-.2*dense**3, atol=3e-3)
        np.testing.assert_allclose(spline.derivative(2)(dense), -2+.6*dense**2, atol=2e-2)
        interior = np.abs(x) < 2.6
        np.testing.assert_allclose(profile.kappa[interior], (2-.6*x**2)[interior], atol=3e-2)
        self.assertLess(marks.alpha_continuous, 0); self.assertGreater(marks.gamma_continuous, 0)

    def test_activation_closure_and_endpoint_reference(self):
        data = pd.read_csv(CSV / "activation_decomposition.csv")
        np.testing.assert_allclose(data.delta_E_act1_hartree + data.delta_E_act2_hartree, data.delta_E_act_hartree, atol=1e-10)
        np.testing.assert_allclose(data.E_TS_hartree - data.E_reactant_hartree, data.delta_E_act_hartree, atol=1e-12)
        np.testing.assert_allclose(data.E_product_hartree - data.E_reactant_hartree, data.delta_E_rxn_hartree, atol=1e-12)
        self.assertFalse(data.valid_for_physical_interpretation.any())

    def test_coordinate_populations_are_separate(self):
        for system in ("MORIPF", "RAFA16"):
            data = pd.read_csv(CSV / f"geometry_change_{system}.csv")
            self.assertIn("covalent_bond_rate_A_per_R_mean", data)
            self.assertIn("interfragment_rate_A_per_R_mean", data)
            self.assertIn("angle_rate_degree_per_R_mean", data)
            self.assertNotEqual(data.covalent_bond_rate_A_per_R_mean.name, data.angle_rate_degree_per_R_mean.name)

    def test_plotted_targets_are_excluded_from_reference_statistics(self):
        populations = pd.read_csv(CSV / "geometry_population_summary.csv").set_index("system")
        for system in ("MoriPF", "RAFA16"):
            row = populations.loc[system]
            self.assertGreaterEqual(row.bond_count-row.bond_reference_count, 1)
            self.assertGreaterEqual(row.interfragment_pair_count-row.interfragment_reference_count, 3)
            atoms = pd.read_csv(CSV / f"atomic_force_{system.upper()}.csv")
            self.assertEqual(int(atoms.reference_atom_count.iloc[0]), 51)
            for _, point in atoms.groupby(["R", "path_number", "point_number"]):
                reference = point.loc[~point.excluded_from_reference.astype(bool), "force_magnitude_hartree_per_bohr"].dropna()
                if reference.empty:
                    self.assertTrue(pd.isna(point.all_atom_mean.iloc[0]))
                    self.assertTrue(pd.isna(point.all_atom_sd.iloc[0]))
                    continue
                self.assertAlmostEqual(float(point.all_atom_mean.iloc[0]), float(reference.mean()), places=12)
                self.assertAlmostEqual(float(point.all_atom_sd.iloc[0]), float(reference.std(ddof=0)), places=12)

    def test_visual_connector_not_used_in_numerics(self):
        data = pd.read_csv(CSV / "activation_decomposition.csv")
        self.assertFalse(any("connector" in column.lower() for column in data.columns))
        self.assertTrue(np.allclose(data.delta_E_act_hartree, data.E_TS_hartree-data.E_reactant_hartree))

    def test_coordinate_multiplier_is_applied_before_all_profiles(self):
        for system in ("MoriPF", "RAFA16"):
            raw = pd.read_csv(PROJECT / "analysis2_force" / ".cache" / "parsed" / f"{system}_parsed.csv")
            profile = pd.read_csv(CSV / f"reaction_force_{system.upper()}.csv")
            np.testing.assert_allclose(
                np.sort(profile.R.unique()),
                np.sort(-raw.signed_irc_coordinate.unique()),
                atol=1e-12,
            )

    def test_force_motion_chain_rule_and_fragment_partition(self):
        validation = pd.read_csv(CSV / "atomic_reaction_force_validation_summary.csv").set_index("system")
        self.assertTrue((validation.pearson_correlation > .99).all())
        self.assertTrue((validation.normalized_rmse_by_reaction_force_rms < .05).all())
        fragments = pd.read_csv(CSV / "fragment_force_motion_contribution.csv")
        partition = fragments[fragments.is_disjoint_partition_fragment.astype(bool)]
        for (system, coordinate), group in partition.groupby(["system", "R"]):
            if group.C_fragment_hartree_per_R.notna().all():
                summary = pd.read_csv(CSV / "atomic_reaction_force_validation.csv")
                expected = summary[(summary.system == system) & np.isclose(summary.R, coordinate)].sum_atomic_C_hartree_per_R.iloc[0]
                self.assertAlmostEqual(float(group.C_fragment_hartree_per_R.sum()), float(expected), places=12)

    def test_no_phase_shading_or_robust_distribution_outputs(self):
        plotting_source = (PROJECT / "analysis2_force" / "src" / "plotting.py").read_text()
        self.assertNotIn("axvspan", plotting_source)
        for name in ("geometry_change_MORIPF.csv", "geometry_change_RAFA16.csv", "atomic_force_MORIPF.csv", "atomic_force_RAFA16.csv"):
            columns = [column.lower() for column in pd.read_csv(CSV / name, nrows=1).columns]
            self.assertFalse(any("median" in column or "iqr" in column or "q25" in column or "q75" in column for column in columns))


class BaseCompatibilityTests(unittest.TestCase):
    def test_preserved_base_analysis_and_wrapper_exist(self):
        self.assertTrue((PROJECT / "base_analysis" / "src" / "run_analysis.py").exists())
        self.assertTrue((PROJECT / "base_analysis" / "outputs" / "reports" / "analysis_summary.md").exists())
        self.assertTrue((PROJECT / "src" / "run_analysis.py").exists())


if __name__ == "__main__":
    unittest.main()
