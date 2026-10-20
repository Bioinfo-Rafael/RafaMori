"""Draw 2D formulas with the original 1-based Gaussian atom numbers.

Connectivity is the same audited reactant-side graph used in the contact tables.
Bond orders are assigned by valence, and sanitized before depiction. The diagram
is a connectivity reference, not a drawing of the IRC or TS spatial geometry.
"""
from pathlib import Path
import hashlib
import json

import pandas as pd

HERE=Path(__file__).resolve().parent
OUT=HERE/'outputs'
SYSTEMS=('MoriPF','RAFA14','RAFA15','RAFA16')
COLORS={'core20':(1.,.82,.38),'core21':(.55,.75,1.),'X24':(1.,.55,.6),
        'NH_near':(1.,.76,.45),'X25':(.60,.85,.65),'amineN':(.77,.69,.94)}


def structure_section(system=None):
    systems=[system] if system else SYSTEMS
    text='## 原子番号付きの構造式\n\n'
    text+='数字は各系のGaussianログと表で使用した原子番号。たとえば20は反応C20、21はそのカルボニルO21を示す。炭素は通常の骨格式で表し、各頂点のそばに番号を付けた。Hも省略せず、全57原子の番号を表示している。\n\n'
    text+='色付きの丸は、黄＝反応C20、紫＝求核剤N、赤＝環のアミン側F/H、橙＝選択アミンH、緑＝反対側F/H、青＝O21。色は比較対象を見つける目印で、相互作用の強さを表さない。\n\n'
    text+='基質と求核剤の結合関係を示す2次元構造式であり、TSでの距離・配向を再現した図ではない。形成中のC20–Nや非共有接触は共有結合として結んでいない。F/Hやアミンの番号は系ごとに異なるため、該当する系の図で確認する。\n\n'
    for s in systems:
        if not system:text+=f'### {s}\n\n'
        text+=f'![{s}：全原子番号付き構造式](../figures/24_{s}_numbered_structure.png)\n\n'
        text+=f'[拡大可能なSVG](../figures/24_{s}_numbered_structure.svg) · [原子番号順を保持したMOLファイル](../structures/{s}_numbered.mol)\n\n'
    return text


def main():
    from rdkit import Chem,rdBase
    from rdkit.Chem import rdDetermineBonds,rdDepictor,rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Geometry import Point2D
    mapping=pd.read_csv(OUT/'csv/contact_controls_mapping.csv')
    edges=pd.read_csv(OUT/'csv/contact_controls_covalent_edges.csv')
    structs=OUT/'structures';structs.mkdir(exist_ok=True)
    molecules={};metadata={}
    for system in SYSTEMS:
        m=mapping[mapping.system==system].sort_values('atom');g=edges[edges.system==system]
        assert m.atom.tolist()==list(range(1,58))
        rw=Chem.RWMol()
        for row in m.itertuples():
            atom=Chem.Atom(row.element);atom.SetNoImplicit(True)
            atom.SetProp('atomNote',str(row.atom));atom.SetIntProp('gaussian_atom_number',row.atom)
            atom.SetProp('analysis_site',row.site)
            assert rw.AddAtom(atom)==row.atom-1
        for row in g.itertuples():rw.AddBond(row.atom_a-1,row.atom_b-1,Chem.BondType.SINGLE)
        mol=rw.GetMol()
        rdDetermineBonds.DetermineBondOrders(mol,charge=0,allowChargedFragments=False,embedChiral=False)
        Chem.SanitizeMol(mol)
        assert mol.GetNumAtoms()==57 and Chem.GetFormalCharge(mol)==0
        assert all(atom.GetNumRadicalElectrons()==0 and atom.GetNumImplicitHs()==0 for atom in mol.GetAtoms())
        assert len(Chem.GetMolFrags(mol))==2
        for a,b in [(12,19),(20,21)]:assert mol.GetBondBetweenAtoms(a-1,b-1).GetBondType()==Chem.BondType.DOUBLE
        original={tuple(sorted([r.atom_a,r.atom_b])) for r in g.itertuples()}
        drawn={tuple(sorted([b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1])) for b in mol.GetBonds()}
        assert drawn==original
        molecules[system]=mol
        metadata[system]={'formula':rdMolDescriptors.CalcMolFormula(mol),'smiles':Chem.MolToSmiles(Chem.RemoveHs(mol)),
                          'atom_count':57,'fragment_count':2,'atoms':[],
                          'bonds':[{'a':b.GetBeginAtomIdx()+1,'b':b.GetEndAtomIdx()+1,'type':str(b.GetBondType())} for b in mol.GetBonds()]}
    # Use the same common-scaffold layout in all four drawings. Alkyl H positions
    # are depiction choices; no new cross-system matching of these H is assumed.
    ref=molecules['RAFA15'];rdDepictor.Compute2DCoords(ref)
    refcoords={a.GetProp('analysis_site'):ref.GetConformer().GetAtomPosition(a.GetIdx()) for a in ref.GetAtoms()
               if not a.GetProp('analysis_site').startswith('unmapped_')}
    # Match only drawing slots for equivalent alkyl H attached to the same
    # heavy-atom site. This is not a stereochemical/IRC atom correspondence.
    hydrogen_slots={}
    for atom in ref.GetAtoms():
        if atom.GetProp('analysis_site').startswith('unmapped_'):
            parent=atom.GetNeighbors()[0].GetProp('analysis_site')
            hydrogen_slots.setdefault(parent,[]).append(ref.GetConformer().GetAtomPosition(atom.GetIdx()))
    for system,mol in molecules.items():
        if system!='RAFA15':
            coords={a.GetIdx():Point2D(refcoords[a.GetProp('analysis_site')].x,refcoords[a.GetProp('analysis_site')].y)
                    for a in mol.GetAtoms() if a.GetProp('analysis_site') in refcoords}
            rdDepictor.Compute2DCoords(mol,canonOrient=False,coordMap=coords)
            used={}
            for atom in mol.GetAtoms():
                site=atom.GetProp('analysis_site')
                if site in refcoords:
                    mol.GetConformer().SetAtomPosition(atom.GetIdx(),refcoords[site])
                else:
                    assert atom.GetSymbol()=='H' and atom.GetDegree()==1
                    parent=atom.GetNeighbors()[0].GetProp('analysis_site')
                    slot=used.get(parent,0);used[parent]=slot+1
                    mol.GetConformer().SetAtomPosition(atom.GetIdx(),hydrogen_slots[parent][slot])
        color={a.GetIdx():COLORS[a.GetProp('analysis_site')] for a in mol.GetAtoms() if a.GetProp('analysis_site') in COLORS}
        for ext,drawer in [('svg',rdMolDraw2D.MolDraw2DSVG(1800,1050)),('png',rdMolDraw2D.MolDraw2DCairo(1800,1050))]:
            opt=drawer.drawOptions();opt.annotationFontScale=.85;opt.baseFontSize=.75;opt.minFontSize=18;opt.maxFontSize=32
            opt.bondLineWidth=2.3;opt.padding=.07;opt.addStereoAnnotation=False;opt.legendFontSize=24
            # Highlight atoms only; do not imply highlighted bonds are contacts.
            drawer.DrawMolecule(mol,legend=f'{system} | Gaussian atom numbers | 2D connectivity, not TS geometry',
                                highlightAtoms=list(color),highlightBonds=[],highlightAtomColors=color)
            drawer.FinishDrawing();path=OUT/f'figures/24_{system}_numbered_structure.{ext}'
            if ext=='svg':path.write_text(drawer.GetDrawingText())
            else:drawer.WriteDrawingText(str(path))
        Chem.MolToMolFile(mol,str(structs/f'{system}_numbered.mol'))
        restored=Chem.MolFromMolFile(str(structs/f'{system}_numbered.mol'),removeHs=False)
        assert restored.GetNumAtoms()==57
        assert [a.GetSymbol() for a in restored.GetAtoms()]==[a.GetSymbol() for a in mol.GetAtoms()]
        for atom in mol.GetAtoms():
            pos=mol.GetConformer().GetAtomPosition(atom.GetIdx())
            metadata[system]['atoms'].append({'number':atom.GetIdx()+1,'element':atom.GetSymbol(),
                'site':atom.GetProp('analysis_site'),'annotation':atom.GetProp('atomNote'),'x2d':pos.x,'y2d':pos.y})
    report={'rdkit_version':rdBase.rdkitVersion,'source':'audited reactant-side covalent graph; original atom numbering',
            'bond_order_method':'RDKit DetermineBondOrders, neutral non-radical valence assignment; sanitized',
            'mapping_sha256':hashlib.sha256((OUT/'csv/contact_controls_mapping.csv').read_bytes()).hexdigest(),
            'edges_sha256':hashlib.sha256((OUT/'csv/contact_controls_covalent_edges.csv').read_bytes()).hexdigest(),'systems':metadata}
    (structs/'numbered_structure_audit.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print('Verified all 57 atom labels and covalent edges for MoriPF, RAFA14, RAFA15 and RAFA16.')


if __name__=='__main__':main()
