"""Audit all 57 positions using connectivity and local IRC geometry, not energies."""
from pathlib import Path
import itertools
import json
import numpy as np
import pandas as pd
from numbered_structures import SYSTEMS

OUT=Path(__file__).resolve().parent/'outputs'
ANCHORS={
    'core14':('core06','core05'), 'amineN':('amineC1','amineC2'),
    'core15':('core14','sub_alkylC1'), 'sub_alkylC1':('core15','core14'),
    'sub_alkylC2':('core14','sub_alkylC3'), 'sub_alkylC3':('sub_alkylC2','core14'),
    'amineC1':('amineN','amineC2'), 'amineC2':('amineC1','amineC3'),
    'amineC3':('amineC2','amineC1')}


def build_correspondence():
    mapping=pd.read_csv(OUT/'csv/contact_controls_mapping.csv')
    edges=pd.read_csv(OUT/'csv/contact_controls_covalent_edges.csv')
    descriptors={};members={};graphs={}
    for system in SYSTEMS:
        m=mapping[mapping.system==system];site=m.set_index('site').atom;atom=m.set_index('atom')
        graph={int(a):[] for a in m.atom}
        for r in edges[edges.system==system].itertuples():
            graph[r.atom_a].append(r.atom_b);graph[r.atom_b].append(r.atom_a)
        graphs[system]=graph
        f=pd.read_csv(OUT/f'csv/parsed_oriented_{system}.csv.gz').sort_values(['signed_irc_coordinate','atom_number'])
        rs=f.signed_irc_coordinate.unique();xyz=f[['x','y','z']].to_numpy().reshape(-1,57,3)
        cn=np.linalg.norm(xyz[:,19]-xyz[:,site['amineN']-1],axis=1)
        good=rs<0;grid=np.linspace(2.1,2.5,41)
        assert (np.diff(cn[good])<0).all()
        assert cn[good].min()<=2.1<2.5<=cn[good].max()
        xx=np.array([np.interp(grid,cn[good][::-1],col[good][::-1])
                     for col in xyz.reshape(len(rs),-1).T]).T.reshape(41,57,3)
        xx=np.concatenate([xx,xyz[np.argmin(abs(rs))][None]])
        for parent,(a,b) in ANCHORS.items():
            p=int(site[parent])
            hs=sorted(i for i in graph[p] if (atom.loc[i,'site'] in ['core15','sub_alkylC2']
                      if parent=='core14' else atom.loc[i,'element']=='H'))
            assert len(hs) in [2,3]
            members[system,parent]=hs
            u=xx[:,site[a]-1]-xx[:,p-1];u/=np.linalg.norm(u,axis=1)[:,None]
            v=xx[:,site[b]-1]-xx[:,p-1];v-=np.einsum('ij,ij->i',u,v)[:,None]*u
            assert np.linalg.norm(v,axis=1).min()>.1
            v/=np.linalg.norm(v,axis=1)[:,None];w=np.cross(u,v)
            vec=xx[:,np.array(hs)-1]-xx[:,p-1,None,:];vec/=np.linalg.norm(vec,axis=2)[:,:,None]
            descriptors[system,parent]=np.einsum('phi,pij->phj',vec,np.stack([u,v,w],axis=2))
    reference=mapping[mapping.system=='MoriPF'].set_index('atom')
    refsites={a:(f'H_ref{a}' if r.site.startswith('unmapped_') else r.site) for a,r in reference.iterrows()}
    assignments={};group_audit=[]
    for system in SYSTEMS:
        for parent in ANCHORS:
            ref=descriptors['MoriPF',parent];d=descriptors[system,parent]
            cost=np.arccos(np.clip(np.einsum('phi,pji->phj',ref,d),-1,1))**2
            permutations=list(itertools.permutations(range(d.shape[1])))
            scores=np.array([cost[:41,np.arange(len(pm)),pm].mean() for pm in permutations])
            order=np.argsort(scores);best=permutations[order[0]]
            perpoint=np.array([cost[:41,np.arange(len(pm)),pm].mean(axis=1) for pm in permutations])
            tsbest=int(np.argmin([cost[-1,np.arange(len(pm)),pm].mean() for pm in permutations]))
            stability=float(np.mean(perpoint.argmin(axis=0)==order[0]))
            assert stability==1. and tsbest==order[0],f'Unstable local match: {system} {parent}'
            assert scores[order[1]]-scores[order[0]]>1e-6
            assignment=dict(zip(members['MoriPF',parent],[members[system,parent][i] for i in best]))
            for refatom,current in assignment.items():assignments[system,current]=refsites[refatom]
            group_audit.append(dict(system=system,parent_site=parent,anchor1=ANCHORS[parent][0],
                anchor2=ANCHORS[parent][1],reference_to_current=json.dumps(assignment),
                rms_local_angle_deg=np.sqrt(scores[order[0]])*180/np.pi,
                second_best_gap_deg2=(scores[order[1]]-scores[order[0]])*(180/np.pi)**2,
                consistent_window_points=int(stability*41),total_window_points=41,TS_agrees=tsbest==order[0]))
    rows=[]
    for r in mapping.itertuples():
        cross=assignments.get((r.system,r.atom),r.site)
        if not r.site.startswith('unmapped_'):assert cross==r.site
        assert not cross.startswith('unmapped_')
        refatom=next(a for a,s in refsites.items() if s==cross)
        if not cross.startswith('X'):assert r.element==reference.loc[refatom,'element']
        else:assert r.element in ['H','F']
        rows.append(dict(system=r.system,atom=r.atom,element=r.element,original_site=r.site,
                         cross_site=cross,reference_MoriPF_atom=refatom,newly_resolved=r.site.startswith('unmapped_'),
                         method='local_geometry_and_connectivity' if (r.system,r.atom) in assignments
                         else 'fixed_scaffold_connectivity'))
    result=pd.DataFrame(rows)
    reference_edges=None
    for system,g in result.groupby('system',sort=False):
        assert len(g)==g.cross_site.nunique()==57 and set(g.cross_site)==set(refsites.values())
        sites=g.set_index('atom').cross_site.to_dict()
        mapped_edges={tuple(sorted([sites[a],sites[b]])) for a,ns in graphs[system].items() for b in ns}
        if reference_edges is None:reference_edges=mapped_edges
        assert mapped_edges==reference_edges,f'Connectivity differs: {system}'
    result.to_csv(OUT/'csv/audited_atom_correspondence.csv',index=False)
    pd.DataFrame(group_audit).to_csv(OUT/'csv/atom_correspondence_geometry_audit.csv',index=False)
    roles=pd.read_csv(OUT/'csv/higher_rank_atom_dictionary.csv').set_index(['system','atom'])
    report='''# 全57原子の対応を再確認した結果

以前の横並び表では、アルキル基の15個のHを一律に未対応としていた。そのため対応位置が存在するのに1点しか描かず、フィットをN/Aにしていた。今回の構造ではこの除外は不要だった。

**Lysモデルとクマリン骨格・そのアミノ置換基は、元素も含めて原子ごとに対応する。置換ベンゼン環では、環炭素とその結合位置を対応させ、位置によってF/Hが異なることを許す。全4系で57対57の一対一対応と、共有結合グラフの一致を確認した。**

MoriPFを基準とし、同じ親原子に結合するHを照合した。親原子と2つの固定骨格原子から右手系の局所座標を作り、各結合の単位方向ベクトルを比較した。CH₂は2通り、CH₃は6通りの対応をすべて評価し、C–N距離2.10〜2.50 Åの41点で角度差の二乗平均が最小になる対応を採用した。アミンN上の2つのH、C15上のH16/H17、基質N14の2つのエチル枝も同じ方法で再確認した。

**全グループで41点すべての最適対応が一致し、独立に調べたTSでも同じ対応になった。** 既存の対応は入れ替わらず、これまで未対応だった各系15個のHだけを追加できた。仕事・力・ΔG‡・R²は対応決定に使用していない。各系内部の仕事、順位、上位と直下3件の選択は変更せず、横並び表の行の統合と4点フィットを更新した。

下表の数字は元ログの原子番号である。「同じ原子位置」は「同じ番号」とは限らず、たとえばMoriPFのH49は他3系のH47に対応する。F/Hの置換により入力内の原子番号がずれるためである。

| 対応位置 | MoriPF | RAFA14 | RAFA15 | RAFA16 | MoriPFでの位置・役割 |
| --- | --- | --- | --- | --- | --- |
'''
    for a,site in refsites.items():
        names=[]
        for system in SYSTEMS:
            r=result[(result.system==system)&(result.cross_site==site)].iloc[0]
            names.append(f'{r.element}{r.atom}')
        report+='| '+' | '.join([site,*names,roles.loc[('MoriPF',a),'role_ja']])+' |\n'
    report+='\n[全原子の対応CSV](../csv/audited_atom_correspondence.csv) · [局所座標・候補差・41点とTSの一致の監査CSV](../csv/atom_correspondence_geometry_audit.csv)\n'
    report+='\n局所角度の差は分子間の配座の違いも反映する。たとえばRAFA15のN14のエチル枝は基準との角度差が大きいが、枝を交換する候補より小さく、41点とTSで対応は一貫している。構造が完全に重なるという意味ではない。\n'
    (OUT/'reports/atom_correspondence_audit_ja.md').write_text(report)
    return result


if __name__=='__main__':build_correspondence()
