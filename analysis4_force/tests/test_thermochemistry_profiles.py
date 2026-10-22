"""Check thermodynamic identities, approximation provenance and C1 display joins."""
import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'analysis4_force'))
from energetics import entropy_terms, KCAL
from thermochemistry_profiles import connector, connector_width, OUT


class ThermochemistryTests(unittest.TestCase):
    def test_entropy_is_h_minus_g_and_not_e_minus_g(self):
        def state(e,h,g):return dict(energy_hartree=e,H_hartree=h,G_hartree=g,temperature_K=300.)
        result=entropy_terms(state(-100,-99,-99.5),state(-109,-107,-108),state(-110,-108,-109),state(-10,-9.5,-9.75))
        self.assertAlmostEqual(result['T_delta_S_activation_kcal_mol'],.25*KCAL)
        self.assertAlmostEqual(result['G_minus_E_activation_kcal_mol'],.25*KCAL)
        self.assertAlmostEqual(result['H_minus_E_activation_kcal_mol'],.5*KCAL)
        self.assertNotAlmostEqual(result['T_delta_S_activation_kcal_mol'],-result['G_minus_E_activation_kcal_mol'])

    def test_real_barrier_identity_and_independent_printed_entropy(self):
        e=pd.read_csv(OUT/'csv/energy_comparison.csv')
        np.testing.assert_allclose(e.absolute_RRHO_G_barrier_kcal_mol,
            e.absolute_electronic_barrier_kcal_mol+e.H_minus_E_activation_kcal_mol-e.T_delta_S_activation_kcal_mol,atol=1e-8)
        audit=pd.read_csv(OUT/'csv/thermochemistry_log_audit.csv')
        self.assertEqual(len(audit),17)
        self.assertTrue((audit[audit.role=='irc'].freq_thermochemistry_blocks==0).all())
        self.assertLess(audit.H_minus_G_minus_printed_TS_kcal_mol.abs().max(),.001)

    def test_reverse_barriers_use_current_product_and_close_thermodynamic_cycle(self):
        audit=json.loads((OUT/'reports/stationary_point_audit.json').read_text())
        barriers=pd.read_csv(OUT/'csv/directional_barriers.csv').set_index('system')
        energy=pd.read_csv(OUT/'csv/energy_comparison.csv').set_index('system')
        saved=json.loads((OUT/'reports/directional_barriers.json').read_text())
        for system,row in barriers.iterrows():
            records={r['role']:r for r in audit if r['system']==system}
            exported=next(r for r in saved if r['system']==system)
            self.assertEqual(exported['product_sha256'],records['product']['sha256'])
            for quantity,key,label in [('E','energy_hartree','electronic'),('H','H_hartree','enthalpy'),('G','G_hartree','RRHO_G')]:
                reverse=(records['ts'][key]-records['product'][key])*KCAL
                self.assertAlmostEqual(row[f'reverse_delta_{quantity}_kcal_mol'],reverse,places=8)
                self.assertAlmostEqual(energy.loc[system,f'reverse_{label}_barrier_kcal_mol'],reverse,places=8)
                self.assertAlmostEqual(exported[f'reverse_delta_{quantity}_kcal_mol'],reverse,places=8)
            for quantity in ['delta_E','delta_H','delta_G','T_delta_S']:
                self.assertAlmostEqual(row[f'forward_{quantity}_kcal_mol']-row[f'reverse_{quantity}_kcal_mol'],
                                       row[f'reaction_{quantity}_kcal_mol'],places=8)
        r=barriers.loc['RAFA15']
        self.assertEqual(Path(r.product_log).name,'RAFA15_LYS_OPTFREQ2.LOG')
        self.assertAlmostEqual(r.reverse_delta_E_kcal_mol,2.2629309403,places=8)
        self.assertAlmostEqual(r.reverse_delta_G_kcal_mol,.9964850448,places=8)
        self.assertTrue(r.forward_RRHO_G_provisional)
        self.assertFalse(r.reverse_RRHO_G_provisional)

    def test_irc_values_and_constant_correction_are_preserved(self):
        data=pd.read_csv(OUT/'csv/thermochemistry_display_profiles.csv')
        e=pd.read_csv(OUT/'csv/energy_comparison.csv').set_index('system')
        for system,row in e.iterrows():
            p=pd.read_csv(OUT/f'csv/reaction_force_{system}.csv')
            curves=data[(data.system==system)&(data.segment=='IRC')]
            ev=curves[curves.quantity=='E'].set_index('R_display').value_kcal_mol
            gv=curves[curves.quantity=='G'].set_index('R_display').value_kcal_mol
            sv=curves[curves.quantity=='T_delta_S'].set_index('R_display').value_kcal_mol
            np.testing.assert_allclose(ev.reindex(p.R).to_numpy(),p.E_relative_separated_reactants_kcal_mol,atol=1e-9)
            np.testing.assert_allclose(gv-ev,row.G_minus_E_activation_kcal_mol,atol=1e-9)
            np.testing.assert_allclose(sv,row.T_delta_S_activation_kcal_mol,atol=1e-9)
            self.assertAlmostEqual(ev.loc[0.],row.absolute_electronic_barrier_kcal_mol,places=8)
            self.assertAlmostEqual(gv.loc[0.],row.absolute_RRHO_G_barrier_kcal_mol,places=8)
        joins=data[data.segment!='IRC']
        self.assertTrue(joins.display_coordinate_only.all())
        self.assertEqual(set(joins.provenance),{'visual_connector_NOT_computed_IRC'})

    def test_joins_match_values_and_slopes_without_vertical_segments(self):
        audit=pd.read_csv(OUT/'csv/thermochemistry_connector_audit.csv')
        self.assertEqual(len(audit),24)
        for row in audit.itertuples():
            self.assertGreater(row.x_end,row.x_start)
            c=connector(row.x_start,row.y_start_kcal_mol,row.slope_start,row.x_end,row.y_end_kcal_mol,row.slope_end)
            np.testing.assert_allclose(c([row.x_start,row.x_end]),[row.y_start_kcal_mol,row.y_end_kcal_mol],atol=1e-10)
            np.testing.assert_allclose(c.derivative()([row.x_start,row.x_end]),[row.slope_start,row.slope_end],atol=1e-10)
        # Opposing level gap and slope require a turning point; do not silently
        # enforce monotonicity or claim this artificial minimum is physical.
        c=connector(-4,0,0,0,-2,1)
        roots=c.derivative().roots(extrapolate=False)
        self.assertTrue(np.any((roots>-4)&(roots<0)))
        self.assertTrue(np.isfinite(connector_width(8,10,0)))
        with self.assertRaises(ValueError):connector(0,0,0,0,1,0)
        for (_,side),group in audit.groupby(['system','side']):
            self.assertEqual(group.x_start.nunique(),1)
            self.assertEqual(group.x_end.nunique(),1)


if __name__=='__main__':unittest.main()
