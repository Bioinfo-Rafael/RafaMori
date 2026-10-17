"""Scientific checks for ortho-site assignment, ratio and work diagnostics."""
import sys,unittest,yaml
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
ROOT=Path(__file__).resolve().parents[2];CSV=ROOT/'analysis4_force/outputs/csv'
sys.path[:0]=[str(ROOT/'base_analysis/src'),str(ROOT/'analysis2_force/src'),str(ROOT/'analysis4_force')]
from two_sided_contacts import pair_projection

class TwoSidedTests(unittest.TestCase):
    def test_same_ortho_site_and_different_elements(self):
        a=pd.read_csv(CSV/'two_sided_atom_mapping.csv').set_index('system')
        self.assertEqual(a.opposite_atom.to_dict(),{'MoriPF':29,'RAFA14':56,'RAFA15':56,'RAFA16':29})
        self.assertEqual(a.opposite_element.to_dict(),{'MoriPF':'F','RAFA14':'H','RAFA15':'H','RAFA16':'F'})
        self.assertTrue((a.opposite_O==21).all());self.assertTrue((a.opposite_ring_C==25).all())
    def test_projection_sign_and_rigid_rotation(self):
        xyz=np.array([[[0.,0.,0.],[2.,0.,0.]]]);forces=np.array([[[1.,0.,0.],[-1.,0.,0.]]]);velocity=np.array([[[0.,0.,0.],[-.5,0.,0.]]])
        q,Q,A=pair_projection(xyz,forces,0,1,velocity)
        np.testing.assert_allclose(q,2);np.testing.assert_allclose(Q,-1);np.testing.assert_allclose(A,.5)
        rotation=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
        rotated=pair_projection(xyz@rotation+17,forces@rotation,0,1,velocity@rotation)
        for actual,expected in zip(rotated,(q,Q,A)):np.testing.assert_allclose(actual,expected,atol=1e-12)
    def test_fixed_window_coverage_and_ratios(self):
        p=pd.read_csv(CSV/'two_sided_force_profiles.csv');m=pd.read_csv(CSV/'two_sided_window_metrics.csv')
        self.assertEqual(len(m),4*3*2)
        for row in m.itertuples():
            g=p[(p.system==row.system)&(p.H_choice==row.H_choice)&(p.R>=row.R_start)&(p.R<=row.R_end)]
            self.assertLess(row.R_end,0);self.assertTrue(g.near_projected_force_hartree_bohr.notna().all());self.assertGreater(len(g),5)
            self.assertTrue(0<=row.opposite_atomic_force_share<=1);self.assertTrue(0<=row.opposite_projected_magnitude_share<=1)
            top=row.opposite_X_mean_force_mEh_bohr+row.opposite_O_mean_force_mEh_bohr
            bottom=top+row.near_X_mean_force_mEh_bohr+row.amine_H_mean_force_mEh_bohr
            self.assertAlmostEqual(row.opposite_atomic_force_share,top/bottom,places=11)
        self.assertFalse(any(c.endswith('_minus_MoriPF') for c in m.columns))
    def test_integrated_all_atom_work_matches_energy_change(self):
        m=pd.read_csv(CSV/'two_sided_window_metrics.csv')
        for system,g in m.groupby('system'):
            profile=pd.read_csv(CSV/f'reaction_force_{system}.csv')
            for row in g.itertuples():
                de=(np.interp(row.R_end,profile.R,profile.E_raw_hartree)-np.interp(row.R_start,profile.R,profile.E_raw_hartree))*627.509474
                self.assertLess(abs(row.all_atomic_work_kcal_mol-de)/abs(de),.002)
    def test_all_correlations_retained_and_uncentered(self):
        m=pd.read_csv(CSV/'two_sided_window_metrics.csv');c=pd.read_csv(CSV/'two_sided_correlations.csv')
        cfg=yaml.safe_load((ROOT/'analysis4_force/config/two_sided_contacts.yaml').read_text())
        self.assertEqual(len(c),3*2*len(cfg['correlation_descriptors'])*2)
        self.assertTrue(c.exploratory_only.all());self.assertFalse(c.independent_validation.any())
        for row in c.itertuples():
            g=m[(m.window==row.window)&(m.H_choice==row.H_choice)]
            r=pearsonr(g[row.descriptor],g[row.target]).statistic
            self.assertAlmostEqual(r,row.pearson_r,places=10)
    def test_atomic_partition_and_full_barrier_accounting(self):
        m=pd.read_csv(CSV/'two_sided_window_metrics.csv')
        np.testing.assert_allclose(m.near_atomic_work_kcal_mol+m.opposite_atomic_work_kcal_mol+m.remaining_atomic_work_kcal_mol,m.all_atomic_work_kcal_mol,atol=1e-12)
        np.testing.assert_allclose(m.near_X_atomic_work_kcal_mol+m.amine_H_atomic_work_kcal_mol,m.near_atomic_work_kcal_mol,atol=1e-12)
        np.testing.assert_allclose(m.opposite_X_atomic_work_kcal_mol+m.opposite_O_atomic_work_kcal_mol,m.opposite_atomic_work_kcal_mol,atol=1e-12)
        np.testing.assert_allclose(m.entry_energy_kcal_mol+m.window_energy_rise_kcal_mol+m.after_window_to_TS_kcal_mol,m.absolute_electronic_barrier_kcal_mol,atol=1e-8)
        self.assertLess(m.work_energy_residual_kcal_mol.abs().max(),.004)

    def test_primary_near_projection_matches_existing_force_calculation(self):
        p=pd.read_csv(CSV/'two_sided_force_profiles.csv');old=pd.read_csv(CSV/'force_motion_alignment.csv')
        for system in ['MoriPF','RAFA15','RAFA16']:
            new=p[(p.system==system)&(p.H_choice=='primary')].sort_values('R')
            original=old[(old.system==system)&(old.target=='hb_substituent')].sort_values('R')
            np.testing.assert_allclose(new.near_projected_force_hartree_bohr,original.projected_generalized_force_hartree_per_bohr,atol=1e-12,equal_nan=True)

if __name__=='__main__':unittest.main()
