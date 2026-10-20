"""Audit union coverage, correspondence and source values of the comparison."""
import json
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
from contact_specificity import pool_mask
from cross_system_work_tables import site_key,pair_key


class CrossSystemWorkTests(unittest.TestCase):
    def test_exact_union_and_source_values(self):
        out=HERE/'outputs'
        cells=pd.read_csv(out/'csv/cross_system_work_union.csv')
        audits=json.loads((out/'reports/cross_system_work_union_audit.json').read_text())
        cfg=yaml.safe_load((HERE/'config/contact_specificity.yaml').read_text())
        atoms=pd.read_csv(out/'csv/contact_controls_atoms.csv')
        pairs=pd.read_csv(out/'csv/contact_controls_pairs.csv').fillna({'nominated':''})
        atoms=atoms[atoms.window=='common_CN_2.50_to_2.10'].copy()
        pairs=pairs[pairs.window=='common_CN_2.50_to_2.10'].copy()
        atoms['cross_key']=[site_key(r.system,r.site) for r in atoms.itertuples()]
        pairs['cross_key']=[pair_key(r.system,r.site_a,r.site_b) for r in pairs.itertuples()]
        for audit in audits:
            kind,target,pool=[audit[k] for k in ['kind','target','pool']]
            table=cells[(cells.kind==kind)&(cells.target==target)&(cells.pool==pool)]
            expected_union=set()
            for system,chosen in audit['selected_keys'].items():
                src=(atoms if kind=='atom' else pairs)
                src=src[src.system==system]
                metric='atomic_work_kcal_mol' if kind=='atom' else 'radial_work_kcal_mol'
                if kind=='atom':
                    cand=src
                    targets=cand[cand.site.isin(['X24','NH_near','X25','core21'])]
                else:
                    t=src[src.nominated==target].iloc[0]
                    cand=src[pool_mask(src,t,pool,cfg)]
                    targets=cand[cand.nominated==target]
                ordered=cand.sort_values([metric,'atom' if kind=='atom' else 'pair'],
                    key=lambda x:x.abs() if x.name==metric else x,ascending=[False,True])
                stop=max(ordered.cross_key.tolist().index(k) for k in targets.cross_key)+3
                cutoff=abs(ordered.iloc[min(stop,len(ordered)-1)][metric])
                expected=set(ordered.loc[ordered[metric].abs()>=cutoff-1e-12,'cross_key'])
                self.assertEqual(set(chosen),expected)
                expected_union|=expected
                for r in table[table.system==system].itertuples():
                    found=src[src.cross_key==r.cross_key]
                    if not r.correspondence_available:
                        self.assertTrue(found.empty)
                        continue
                    self.assertEqual(len(found),1)
                    value=found.iloc[0][metric]
                    self.assertAlmostEqual(r.work_kcal_mol,value)
                    if r.in_candidate_pool:
                        self.assertEqual(r.rank,int((cand[metric].abs()>abs(value)+1e-12).sum())+1)
                    else:
                        self.assertTrue(np.isnan(r.rank))
            self.assertEqual(set(table.cross_key),expected_union)
            self.assertEqual(set(audit['union_keys']),expected_union)
            self.assertTrue((table.groupby('cross_key').size()==4).all())


if __name__=='__main__':unittest.main()
