"""Check complete scaffold correspondence and independent local-geometry audit."""
import unittest
from pathlib import Path
import pandas as pd

OUT=Path(__file__).resolve().parents[1]/'outputs'


class AtomCorrespondenceTests(unittest.TestCase):
    def test_complete_element_and_connectivity_correspondence(self):
        m=pd.read_csv(OUT/'csv/audited_atom_correspondence.csv')
        edges=pd.read_csv(OUT/'csv/contact_controls_covalent_edges.csv')
        ref=m[m.system=='MoriPF'].set_index('cross_site')
        expected=None
        for system,g in m.groupby('system'):
            self.assertEqual(len(g),57)
            self.assertEqual(g.cross_site.nunique(),57)
            self.assertEqual(set(g.cross_site),set(ref.index))
            self.assertEqual(g.newly_resolved.sum(),15)
            for r in g.itertuples():
                rr=ref.loc[r.cross_site]
                if r.cross_site.startswith('X'):
                    self.assertIn(r.element,['H','F'])
                else:self.assertEqual(r.element,rr.element)
                if not r.newly_resolved:self.assertEqual(r.cross_site,r.original_site)
            sites=g.set_index('atom').cross_site.to_dict()
            graph={tuple(sorted([sites[r.atom_a],sites[r.atom_b]]))
                   for r in edges[edges.system==system].itertuples()}
            if expected is None:expected=graph
            self.assertEqual(graph,expected)
            for refatom,atom in [(49,47),(50,48),(52,50),(53,51),(55,53),(56,54),(57,55)]:
                actual=g[g.reference_MoriPF_atom==refatom].atom.iloc[0]
                self.assertEqual(actual,refatom if system=='MoriPF' else atom)

    def test_geometry_assignment_stable_and_all_panels_have_four_points(self):
        audit=pd.read_csv(OUT/'csv/atom_correspondence_geometry_audit.csv')
        self.assertEqual(len(audit),36)
        self.assertTrue((audit.consistent_window_points==41).all())
        self.assertTrue(audit.TS_agrees.all())
        self.assertTrue((audit.second_best_gap_deg2>0).all())
        panels=pd.read_csv(OUT/'csv/work_union_subplot_statistics.csv')
        self.assertTrue((panels.n==4).all())
        self.assertTrue(panels.r_squared.notna().all())
        cells=pd.read_csv(OUT/'csv/cross_system_work_union.csv')
        self.assertTrue(cells.correspondence_available.all())


if __name__=='__main__':unittest.main()
