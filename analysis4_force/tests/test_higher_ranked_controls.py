"""Coverage and geometric/force invariants for individual control identification."""
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parents[1]/'outputs'

class HigherRankedControlsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ids=pd.read_csv(OUT/'csv/higher_rank_contact_identities.csv')
        cls.competitors=pd.read_csv(OUT/'csv/higher_rank_competitors.csv')

    def test_complete_original_rank_counts(self):
        ranks=pd.read_csv(OUT/'csv/contact_controls_within_ranks.csv')
        keys=['kind','system','window','target','pool','metric']
        counts=self.competitors.groupby(keys).size()
        for r in ranks.itertuples():
            self.assertEqual(counts.get(tuple(getattr(r,k) for k in keys),0),r.stronger_controls)
        c=self.competitors
        self.assertTrue((c.control_value.abs()>c.target_value.abs()+1e-12).all())
        self.assertFalse(c.duplicated(keys+['control_atom','control_pair']).any())
        self.assertTrue(c.loc[c.kind=='contact','interpretation_ja'].notna().all())
        self.assertTrue(c.loc[c.kind=='atom','role_ja'].notna().all())

    def test_identity_coverage_and_existing_metrics_unchanged(self):
        source=pd.read_csv(OUT/'csv/contact_controls_pairs.csv');source=source[source.eligible]
        keys=['system','window','pair']
        joined=source.merge(self.ids,on=keys,validate='one_to_one')
        self.assertEqual(len(joined),len(source))
        for metric in ['pair_atomic_work_kcal_mol','radial_work_kcal_mol','mean_abs_projection_mEh_bohr','mean_distance_A']:
            np.testing.assert_allclose(joined[metric+'_x'],joined[metric+'_y'],atol=1e-12)
        dictionary=pd.read_csv(OUT/'csv/higher_rank_atom_dictionary.csv')
        self.assertEqual(dictionary.groupby('system').size().tolist(),[57]*4)
        self.assertTrue(dictionary.role_ja.notna().all())

    def test_signed_force_components_are_consistent(self):
        d=self.ids
        self.assertTrue((d.mean_closing_projection_mEh_bohr>=-1e-12).all())
        self.assertTrue((d.mean_opening_projection_mEh_bohr>=-1e-12).all())
        np.testing.assert_allclose(d.mean_closing_projection_mEh_bohr+d.mean_opening_projection_mEh_bohr,d.mean_abs_projection_mEh_bohr,atol=1e-9)
        np.testing.assert_allclose(d.mean_opening_projection_mEh_bohr-d.mean_closing_projection_mEh_bohr,d.mean_signed_projection_mEh_bohr,atol=1e-9)
        self.assertTrue((d.mean_signed_projection_mEh_bohr.abs()<=d.mean_abs_projection_mEh_bohr+1e-9).all())

    def test_atomic_percentages_use_energy_differences_not_work_sum(self):
        values=pd.read_csv(OUT/'csv/atomic_work_barrier_percentages.csv')
        audit=pd.read_csv(OUT/'csv/atomic_work_barrier_percentage_audit.csv')
        source=pd.read_csv(OUT/'csv/contact_controls_atoms.csv')
        np.testing.assert_allclose(values.atomic_work_kcal_mol,source.atomic_work_kcal_mol)
        self.assertEqual(len(values),57*4*3)
        for system,g in values.groupby('system'):
            trajectory=pd.read_csv(OUT/f'csv/parsed_oriented_{system}.csv.gz')
            energy=trajectory.groupby('signed_irc_coordinate').energy_hartree.first().sort_index()
            barrier=(energy.loc[0.]-energy.iloc[0])*627.509474
            np.testing.assert_allclose(g.IRC_left_to_TS_kcal_mol,barrier,atol=1e-7)
            np.testing.assert_allclose(g.atomic_work_pct_IRC_endpoint_barrier,
                                       100*g.atomic_work_kcal_mol/barrier,atol=1e-7)
        windows=pd.read_csv(OUT/'csv/two_sided_window_metrics.csv')
        windows=windows[windows.H_choice=='primary'].set_index(['system','window'])
        for r in audit.itertuples():
            g=values[(values.system==r.system)&(values.window==r.window)]
            expected=windows.loc[(r.system,r.window)]
            self.assertEqual(len(g),57)
            self.assertAlmostEqual(r.window_energy_rise_kcal_mol,expected.window_energy_rise_kcal_mol)
            self.assertAlmostEqual(r.integration_residual_kcal_mol,expected.work_energy_residual_kcal_mol)
            self.assertAlmostEqual(g.atomic_work_pct_window_energy_rise.sum(),
                                   100*g.atomic_work_kcal_mol.sum()/r.window_energy_rise_kcal_mol)

    def test_target_geometry_matches_existing_independent_analysis(self):
        old=pd.read_csv(OUT/'csv/two_sided_window_metrics.csv')
        old=old[(old.window=='common_CN_2.50_to_2.10')&(old.H_choice=='primary')].set_index('system')
        ids=self.ids[self.ids.window=='common_CN_2.50_to_2.10']
        for system,g in ids.groupby('system'):
            near=g[g.nominated=='amine_contact'].iloc[0];opp=g[g.nominated=='opposite_contact'].iloc[0]
            self.assertAlmostEqual(near.TS_distance_A,old.loc[system,'near_TS_distance_A'],places=8)
            self.assertAlmostEqual(opp.TS_distance_A,old.loc[system,'opposite_TS_distance_A'],places=8)
            if opp.element_a=='H' or opp.element_b=='H':
                self.assertAlmostEqual(opp.TS_DHA_angle_deg,old.loc[system,'opposite_TS_CXO_angle_deg'],places=8)
            formed=g[g.contact_type=='形成中のC20–求核N結合']
            self.assertEqual(len(formed),1)
            self.assertEqual(formed.topology.iloc[0],'分子間')

if __name__=='__main__':unittest.main()
