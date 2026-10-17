"""Check all-contact controls independently against existing selected results."""
import sys,unittest,itertools
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2];CSV=ROOT/'analysis4_force/outputs/csv'
sys.path[:0]=[str(ROOT/'base_analysis/src'),str(ROOT/'analysis2_force/src'),str(ROOT/'analysis4_force')]
from contact_specificity import TARGETS,pairkey

class SpecificityTests(unittest.TestCase):
    def test_complete_atom_pair_accounting_and_site_mapping(self):
        a=pd.read_csv(CSV/'contact_controls_atoms.csv');p=pd.read_csv(CSV/'contact_controls_pairs.csv');m=pd.read_csv(CSV/'contact_controls_mapping.csv')
        self.assertTrue((a.groupby(['system','window']).size()==57).all())
        self.assertTrue((p.groupby(['system','window']).size()==57*56//2).all())
        self.assertTrue((m.groupby('system').cross_system_mapped.sum()==42).all())
        for _,g in m.groupby('system'):
            self.assertEqual(g.atom.nunique(),57);self.assertEqual(g[g.cross_system_mapped].site.nunique(),42)
        self.assertTrue((p[p.eligible].graph_separation>=3).all());self.assertTrue((p[p.eligible].min_distance_A<=4).all())
        self.assertTrue((p[p.nominated.isin(TARGETS)].eligible).all())
    def test_selected_pair_and_atom_values_match_existing_analysis(self):
        old=pd.read_csv(CSV/'two_sided_window_metrics.csv');old=old[old.H_choice=='primary'].set_index(['system','window'])
        p=pd.read_csv(CSV/'contact_controls_pairs.csv');a=pd.read_csv(CSV/'contact_controls_atoms.csv')
        for row in p[p.nominated.isin(TARGETS)].itertuples():
            prefix='near' if row.nominated=='amine_contact' else 'opposite';ref=old.loc[(row.system,row.window)]
            self.assertAlmostEqual(row.radial_work_kcal_mol,ref[prefix+'_pair_diagnostic_work_kcal_mol'],places=10)
            self.assertAlmostEqual(row.pair_atomic_work_kcal_mol,ref[prefix+'_atomic_work_kcal_mol'],places=10)
            self.assertAlmostEqual(row.mean_abs_projection_mEh_bohr,ref[prefix+'_mean_abs_projected_force_mEh_bohr'],places=10)
        for row in a[a.site.isin(['X24','NH_near','X25','core21'])].itertuples():
            prefix={'X24':'near_X','NH_near':'amine_H','X25':'opposite_X','core21':'opposite_O'}[row.site]
            self.assertAlmostEqual(row.atomic_work_kcal_mol,old.loc[(row.system,row.window),prefix+'_atomic_work_kcal_mol'],places=10)
        for key,g in a.groupby(['system','window']):
            self.assertAlmostEqual(g.atomic_work_kcal_mol.sum(),old.loc[key,'all_atomic_work_kcal_mol'],places=10)
    def test_common_correlations_use_four_eligible_corresponding_sites(self):
        p=pd.read_csv(CSV/'contact_controls_pairs.csv');c=pd.read_csv(CSV/'contact_controls_correlations.csv');e=pd.read_csv(CSV/'energy_comparison.csv').set_index('system')
        for row in c[c.kind=='contact'].itertuples():
            g=p[(p.window==row.window)&(p.pair==row.item)].set_index('system')
            self.assertEqual(len(g),4);self.assertTrue(g.eligible.all());self.assertTrue(g.cross_system_mapped.all())
            g=g.loc[e.index]
            r=np.corrcoef(g[row.metric],e[row.barrier])[0,1]
            self.assertAlmostEqual(r,row.pearson_r,places=10)
    def test_rankings_and_disjoint_counterexamples(self):
        a=pd.read_csv(CSV/'contact_controls_atoms.csv');r=pd.read_csv(CSV/'contact_controls_within_ranks.csv')
        for row in r[r.kind=='atom'].itertuples():
            g=a[(a.system==row.system)&(a.window==row.window)]
            expected=1+(abs(g[row.metric])>abs(row.value)+1e-12).sum()
            self.assertEqual(expected,row.magnitude_rank)
        ex=pd.read_csv(CSV/'contact_controls_counterexamples.csv')
        for row in ex[ex.pool=='disjoint'].itertuples():
            self.assertFalse(set(row.control_pair.split('|')) & set(TARGETS[row.target]))
    def test_scan_calibration_retains_all_label_permutations(self):
        x=pd.read_csv(CSV/'contact_controls_label_permutations.csv')
        self.assertTrue((x.n_label_permutations==24).all());self.assertFalse(x.confirmatory_p_value.any())
        self.assertTrue(x.permutations_with_scan_max_at_least_target.between(1,24).all())
        s=pd.read_csv(CSV/'contact_controls_specificity.csv')
        self.assertTrue((s.controls_at_least_as_correlated<=s.n_common_controls).all())
        self.assertTrue((s.stronger_and_at_least_as_correlated<=s.stronger_in_all_four).all())

if __name__=='__main__':unittest.main()
