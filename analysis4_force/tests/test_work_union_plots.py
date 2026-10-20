"""Validate complete subplot coverage and regression/dispersion statistics."""
import sys
import json
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
from work_union_plots import statistics


class WorkUnionPlotTests(unittest.TestCase):
    def test_degenerate_cases_and_signed_mean(self):
        self.assertTrue(np.isnan(statistics(np.array([1.]),np.array([10.]))['r_squared']))
        self.assertTrue(np.isnan(statistics(np.ones(4),np.arange(4.))['r_squared']))
        zero=statistics(np.array([-2.,-1.,1.,2.]),np.arange(4.))
        self.assertTrue(np.isnan(zero['std_over_signed_mean']))
        negative=statistics(np.array([-4.,-3.,-2.,-1.]),np.arange(4.))
        self.assertLess(negative['std_over_signed_mean'],0)
        self.assertAlmostEqual(negative['r_squared'],1)

    def test_all_union_rows_and_independent_statistics(self):
        out=HERE/'outputs'
        stats=pd.read_csv(out/'csv/work_union_subplot_statistics.csv')
        points=pd.read_csv(out/'csv/work_union_subplot_points.csv')
        union=pd.read_csv(out/'csv/cross_system_work_union.csv')
        energy=pd.read_csv(out/'csv/energy_comparison.csv').set_index('system')
        keys=['kind','target','pool','cross_key']
        expected=set(map(tuple,union[keys].drop_duplicates().to_numpy()))
        self.assertEqual(set(map(tuple,stats[keys].to_numpy())),expected)
        for row in stats.itertuples():
            p=points[points.panel_id==row.panel_id]
            self.assertEqual(len(p),row.n)
            source=union[(union.kind==row.kind)&(union.target==row.target)&
                         (union.pool==row.pool)&(union.cross_key==row.cross_key)]
            self.assertEqual(len(p),int(source.correspondence_available.sum()))
            for point in p.itertuples():
                original=source[source.system==point.system].iloc[0]
                self.assertAlmostEqual(original.work_kcal_mol,point.work_kcal_mol)
                self.assertAlmostEqual(energy.loc[point.system,'absolute_RRHO_G_barrier_kcal_mol'],
                                       point.forward_activation_G_kcal_mol)
            if row.n==4:
                x=p.work_kcal_mol.to_numpy();y=p.forward_activation_G_kcal_mol.to_numpy()
                fit=linregress(x,y)
                self.assertAlmostEqual(row.r_squared,fit.rvalue**2)
                self.assertAlmostEqual(row.slope,fit.slope,places=7)
                self.assertAlmostEqual(row.std_over_signed_mean,np.std(x,ddof=0)/np.mean(x))
            else:
                self.assertTrue(np.isnan(row.r_squared))
                self.assertTrue(np.isnan(row.std_over_signed_mean))
        manifest=json.loads((out/'reports/work_union_subplot_manifest.json').read_text())
        for kind in ['scatter','bars']:
            panels=[p for page in manifest if page['kind']==kind for p in page['panel_ids']]
            self.assertEqual(len(panels),len(stats))
            self.assertEqual(set(panels),set(stats.panel_id))
        for page in manifest:self.assertTrue((out/'figures'/page['path']).exists())


if __name__=='__main__':unittest.main()
