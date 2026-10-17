"""Scientific invariants and regression checks on the generated four-system run."""
import sys,unittest,json,gzip,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'analysis4_force/outputs';CSV=OUT/'csv'
sys.path[:0]=[str(ROOT/'base_analysis/src'),str(ROOT/'analysis2_force/src'),str(ROOT/'analysis4_force')]
from parse_gaussian_irc import parse_log
from inputs import parse_stationary
from numerics4 import reaction_force_profile

class Analysis4Tests(unittest.TestCase):
    def test_original_hashes_and_compressed_roundtrip(self):
        for row in json.loads((ROOT/'data/completed_20261017/manifest.json').read_text()):
            path=ROOT/row['path']
            self.assertEqual(hashlib.sha256(gzip.decompress(Path(str(path)+'.gz').read_bytes())).hexdigest(),row['sha256'])
    def test_four_oriented_logs_and_ts(self):
        audit=pd.read_csv(CSV/'orientation_audit.csv')
        self.assertEqual(set(audit.system),{'MoriPF','RAFA14','RAFA15','RAFA16'})
        self.assertEqual(audit.coordinate_multiplier.tolist(),[-1,1,-1,-1])
        for row in audit.itertuples():
            self.assertGreater(row.reactant_side_CN_A,row.product_side_CN_A)
            f,meta=parse_log(ROOT/row.oriented_log,row.system,1,57,627.509474)
            self.assertFalse(f.irc_coordinate_estimated.any())
            self.assertEqual(len(f),57*row.point_count)
            self.assertEqual(f.groupby('signed_irc_coordinate').size().unique().tolist(),[57])
            ts=f[f.signed_irc_coordinate.abs()<1e-12]
            self.assertEqual(len(ts),57)
            self.assertLess(abs(row.TS_energy_difference_kcal_mol),.001)
            self.assertLess(row.TS_midpoint_alignment_RMSD_A,.001)
    def test_stationary_audit_catches_wrong_lys_level(self):
        a=json.loads((OUT/'reports/stationary_point_audit.json').read_text())
        lys=next(r for r in a if r['system']=='Lys')
        self.assertEqual(lys['method_basis'],'B3LYP/6-311+G(2D,P)')
        self.assertEqual(lys['frequency_count'],0)
        self.assertIsNone(lys['G_hartree'])
        for row in a:
            if row['role']=='ts':
                self.assertEqual(len(row['imaginary_frequencies']),1)
                self.assertEqual(row['frequency_count'],165)
                self.assertEqual(row['method_basis'],'M062X/6-311++G(D,P)')
    def test_barrier_cancellation_and_closure(self):
        e=pd.read_csv(CSV/'energy_comparison.csv').set_index('system');ref=e.loc['RAFA16']
        for _,row in e.iterrows():
            self.assertFalse(row.absolute_barrier_valid)
            self.assertAlmostEqual(row.IRC_left_to_alpha_kcal_mol+row.alpha_to_TS_kcal_mol,row.IRC_left_to_TS_kcal_mol,places=9)
            for common_lys in [-174.4,-174.6]:
                barrier=(row.TS_E_hartree-row.substrate_E_hartree-common_lys)-(ref.TS_E_hartree-ref.substrate_E_hartree-common_lys)
                self.assertAlmostEqual(barrier*627.509474,row.relative_electronic_barrier_kcal_mol,places=8)
            dW2=row.alpha_to_TS_kcal_mol-ref.alpha_to_TS_kcal_mol
            self.assertAlmostEqual(row.relative_preparation_component_kcal_mol+dW2,row.relative_electronic_barrier_kcal_mol,places=8)
    def test_analytic_derivatives_and_extrema(self):
        for path in CSV.glob('parsed_oriented_*.csv.gz'):
            f=pd.read_csv(path);p,m,s=reaction_force_profile(f,627.509474,4001,None)
            np.testing.assert_allclose(p.F_smooth_hartree_per_R,-s.derivative()(p.R),atol=1e-12)
            np.testing.assert_allclose(p.kappa_paper_hartree_per_R2,s.derivative(2)(p.R),atol=1e-12)
            np.testing.assert_allclose(p.kappa_paper_hartree_per_R2,-p.kappa_smooth_hartree_per_R2,atol=1e-12)
            self.assertLess(m.alpha_continuous,0);self.assertGreater(m.gamma_continuous,0)
            self.assertLess(abs(float(s.derivative(2)(m.alpha_continuous))),1e-7)
            self.assertLess(abs(float(s.derivative(2)(m.gamma_continuous))),1e-7)
    def test_force_motion_partition(self):
        v=pd.read_csv(CSV/'force_motion_validation_summary.csv')
        self.assertTrue((v.pearson_r>.99).all());self.assertTrue((v.normalized_RMSE<.03).all())
        frag=pd.read_csv(CSV/'fragment_contributions.csv');total=pd.read_csv(CSV/'force_motion_validation.csv').set_index(['system','R'])
        for key,g in frag[frag.is_disjoint_partition_fragment].groupby(['system','R']):
            self.assertEqual(g.atom_count.sum(),57)
            if g.C_fragment_hartree_per_R.notna().all():self.assertAlmostEqual(g.C_fragment_hartree_per_R.sum(),total.loc[key,'sum_atomic_C_hartree_per_R'],places=11)
    def test_all_contacts_and_fixed_selection(self):
        c=pd.read_csv(CSV/'all_NH_F_contacts.csv');s=pd.read_csv(CSV/'system_summary.csv').set_index('system')
        for system,g in c.groupby('system'):
            count=10 if system=='MoriPF' else 6
            self.assertEqual(g.groupby('R').size().unique().tolist(),[count])
            selected=g[g.selected_at_TS]
            self.assertEqual(len(selected[['F_atom','H_atom']].drop_duplicates()),1)
            self.assertEqual(selected.groupby('R').size().unique().tolist(),[1])
            self.assertAlmostEqual(g[g.R.abs()<1e-12].FH_A.min(),s.loc[system,'TS_nearest_FH_A'],places=10)
        self.assertEqual(s.ortho_F_count.to_dict(),{'MoriPF':2,'RAFA14':0,'RAFA15':1,'RAFA16':2})

if __name__=='__main__':unittest.main()
