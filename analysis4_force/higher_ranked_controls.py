"""Identify every control above nominated atoms/contacts, without barrier ranking.

Reuses all control values and selection rules. Adds atom identities, topology,
H-contact angles and signed projections from existing IRC coordinates/forces.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import yaml
from contact_specificity import pool_mask, POOLS
from numbered_structures import structure_section, SYSTEMS

HERE=Path(__file__).resolve().parent
OUT=HERE/'outputs'
PRIMARY='common_CN_2.50_to_2.10'
ATOM_TARGETS=['X24','NH_near','X25','core21']
PAIR_METRICS=['pair_atomic_work_kcal_mol','radial_work_kcal_mol','mean_abs_projection_mEh_bohr']
ATOM_METRICS=['atomic_work_kcal_mol','mean_force_mEh_bohr']
METRICS={'pair_atomic_work_kcal_mol':'2原子の仕事の合計（kcal/mol）','radial_work_kcal_mol':'接触方向の仕事（kcal/mol）',
         'mean_abs_projection_mEh_bohr':'接触方向の力の平均絶対値（mEh/Bohr）',
         'atomic_work_kcal_mol':'原子の仕事（kcal/mol）','mean_force_mEh_bohr':'原子に働く力の平均大きさ（mEh/Bohr）'}
NAMES={'amine_contact':'アミン側','opposite_contact':'反対側'}


def table(headers,rows):
    def clean(x):return str(x).replace('|','／').replace('\n',' ')
    return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(clean(x) for x in row)+' |\n' for row in rows)


def fmt(x):return '—' if pd.isna(x) else f'{x:.6g}'


def work_percentages(atoms):
    energy=pd.read_csv(OUT/'csv/energy_comparison.csv')
    windows=pd.read_csv(OUT/'csv/two_sided_window_metrics.csv')
    windows=windows[windows.H_choice=='primary'][['system','window','window_energy_rise_kcal_mol']]
    result=atoms.merge(windows,on=['system','window'],validate='many_to_one').merge(
        energy[['system','IRC_left_to_TS_kcal_mol','absolute_electronic_barrier_kcal_mol']],
        on='system',validate='many_to_one')
    assert (result.IRC_left_to_TS_kcal_mol>0).all()
    assert (result.window_energy_rise_kcal_mol>0).all()
    result['atomic_work_pct_IRC_endpoint_barrier']=100*result.atomic_work_kcal_mol/result.IRC_left_to_TS_kcal_mol
    result['atomic_work_pct_window_energy_rise']=100*result.atomic_work_kcal_mol/result.window_energy_rise_kcal_mol
    audit=result.groupby(['system','window'],sort=False).agg(
        atom_count=('atom','size'),sum_atomic_work_kcal_mol=('atomic_work_kcal_mol','sum'),
        window_energy_rise_kcal_mol=('window_energy_rise_kcal_mol','first'),
        IRC_endpoint_barrier_kcal_mol=('IRC_left_to_TS_kcal_mol','first'),
        separated_reactants_barrier_kcal_mol=('absolute_electronic_barrier_kcal_mol','first'),
        sum_pct_IRC_endpoint_barrier=('atomic_work_pct_IRC_endpoint_barrier','sum'),
        sum_pct_window_energy_rise=('atomic_work_pct_window_energy_rise','sum')).reset_index()
    assert (audit.atom_count==57).all()
    audit['integration_residual_kcal_mol']=audit.sum_atomic_work_kcal_mol-audit.window_energy_rise_kcal_mol
    result.to_csv(OUT/'csv/atomic_work_barrier_percentages.csv',index=False)
    audit.to_csv(OUT/'csv/atomic_work_barrier_percentage_audit.csv',index=False)
    return result


def percentage_summary():
    audit=pd.read_csv(OUT/'csv/atomic_work_barrier_percentage_audit.csv')
    audit=audit[audit.window==PRIMARY]
    text='## 原子仕事の合計と電子障壁に対する割合\n\n'
    text+='原子の仕事の表には、全57原子と合計行を載せた。分母を明示して、次の2種類の割合を併記する。どちらも符号付きの仕事を使用し、絶対値にはしていない。\n\n'
    text+=r'''$$
p_i^{\mathrm{IRC}}=100\,\frac{W_i^{[2.50\to2.10]}}{E_{\mathrm{TS}}-E_{\mathrm{IRC,start}}},
\qquad
p_i^{\mathrm{window}}=100\,\frac{W_i^{[2.50\to2.10]}}{E_{c=2.10}-E_{c=2.50}}
$$

前者は「IRC反応物側端点→TSの電子障壁の何％を、この部分区間でその原子が担ったか」、後者は「積分区間の電子エネルギー上昇の何％をその原子が担ったか」を表す。分母は以下の実際のエネルギー差であり、仕事の合計を分母にして100%へ強制的に合わせてはいない。

**全原子の仕事の和が対応するのは、積分区間の電子エネルギー差。現在の2.50→2.10 Åの仕事を足しても、端点→TSの障壁全体にはならない。** 端点からTSまで同じ積分を行えば、その区間の電子障壁に対応する。自由エネルギー障壁ΔG‡や熱補正込みの障壁には対応しない。

| 分子系 | IRC端点→TSのΔE‡ / kcal/mol | 積分区間ΔE / kcal/mol | 全57原子の仕事ΣW / kcal/mol | 積分残差ΣW−区間ΔE / kcal/mol | ΣW／IRC端点→TSのΔE‡ / % | ΣW／区間ΔE / % | 分離反応物＋Lys→TSのΔE‡ / kcal/mol |
| --- | --- | --- | --- | --- | --- | --- | --- |
'''
    for r in audit.itertuples():
        text+='| '+' | '.join([r.system,*[fmt(v) for v in [r.IRC_endpoint_barrier_kcal_mol,r.window_energy_rise_kcal_mol,
            r.sum_atomic_work_kcal_mol,r.integration_residual_kcal_mol,r.sum_pct_IRC_endpoint_barrier,
            r.sum_pct_window_energy_rise,r.separated_reactants_barrier_kcal_mol]]])+' |\n'
    text+='\nIRC端点は分離反応物の最適化構造とは異なる。したがって、端点→TSの障壁が分離反応物基準の障壁を必ず過小評価するとは限らない。今回、分離反応物基準のΔE‡はMoriPF・RAFA15・RAFA16で負であり、上の割合の分母には使用していない。負の原子仕事はエネルギー低下への寄与で、その割合も負になる。全原子の割合の和と100%の差は、区間の不足と積分残差を区別して読む。\n\n'
    text+='[全57原子×4系×3区間の割合CSV](../csv/atomic_work_barrier_percentages.csv) · [合計・分母・残差の監査CSV](../csv/atomic_work_barrier_percentage_audit.csv)\n\n'
    return text


def rank_context(frame,metric,key,targets,below=3):
    """All preceding rows, nominated rows, and at least three following rows.

    Ties at the final boundary are retained. Competition ranks use the same
    numerical tolerance as the source analysis.
    """
    ordered=frame.sort_values([metric,key],key=lambda x: x.abs() if x.name==metric else x,
                              ascending=[False,True]).copy()
    ordered['rank']=[int((frame[metric].abs()>abs(v)+1e-12).sum())+1 for v in ordered[metric]]
    ordered['is_target']=ordered[key].isin(targets)
    positions=np.flatnonzero(ordered.is_target.to_numpy())
    assert len(positions)==len(targets)
    last=min(len(ordered)-1,int(positions.max())+below)
    cutoff=abs(ordered.iloc[last][metric])
    return ordered[ordered[metric].abs()>=cutoff-1e-12]


def context_tables(system,atoms,pairs,dictionary,config):
    """Self-contained tables used by both the overview and individual reports."""
    a=atoms[(atoms.system==system)&(atoms.window==PRIMARY)].merge(
        dictionary[dictionary.system==system],on=['system','atom','element','site'])
    p=pairs[(pairs.system==system)&(pairs.window==PRIMARY)]
    text='### 原子の仕事・力：着目する4原子と、その前後\n\n'
    text+=f'この系の割合の分母：IRC端点→TSのΔE‡ = {a.IRC_left_to_TS_kcal_mol.iloc[0]:.6f} kcal/mol、積分区間ΔE = {a.window_energy_rise_kcal_mol.iloc[0]:.6f} kcal/mol。\n\n'
    text+='★が着目原子。仕事は全57原子と合計を示す。力は着目4原子の上位すべてと直下3件を含む範囲を示す。順位は絶対値の降順で、仕事と割合には符号を残す。割合の分母はIRC端点→TSの電子障壁と、C–N 2.50→2.10 Åの区間エネルギー差である。\n\n'
    for metric in ATOM_METRICS:
        iswork=metric=='atomic_work_kcal_mol'
        g=rank_context(a,metric,'atom',a.loc[a.target_atom,'atom'].tolist(),below=57 if iswork else 3)
        text+=f'<details>\n<summary>{METRICS[metric]}：上位・着目原子・下位を開く</summary>\n\n'
        headers=['順位','原子',METRICS[metric]]
        if iswork:headers+=['原子仕事／IRC端点→TSのΔE‡ / %','原子仕事／区間ΔE / %']
        headers+=['位置・役割','直接結合する原子']
        rows=[]
        for r in g.itertuples():
            row=[r.rank,'**★ '+f'{r.element}{r.atom}'+'**' if r.is_target else f'{r.element}{r.atom}',
                 f'**{fmt(getattr(r,metric))}**' if r.is_target else fmt(getattr(r,metric))]
            if iswork:row+=[fmt(r.atomic_work_pct_IRC_endpoint_barrier),fmt(r.atomic_work_pct_window_energy_rise)]
            rows.append(row+[r.role_ja,r.neighbors])
        if iswork:rows.append(['—','**全57原子の合計**',fmt(a.atomic_work_kcal_mol.sum()),
            fmt(a.atomic_work_pct_IRC_endpoint_barrier.sum()),fmt(a.atomic_work_pct_window_energy_rise.sum()),'C–N 2.50→2.10 Åの積分','—'])
        text+=table(headers,rows)+'\n</details>\n\n'
    text+='### 接触の仕事・力：着目する接触と、その前後\n\n'
    text+='表ごとに指定接触を★で示し、上位すべて・指定接触・直下3件を同じ表に並べた。同値が境界にある場合は同値の行も残す。下位候補が3件未満なら残りをすべて示す。順位は各表の比較条件の中で付けたもの。Q̄の負は距離を縮める向き、正は開く向きで、2原子間だけの引力・斥力を示すものではない。\n\n'
    atomvalues=a.set_index('atom').atomic_work_kcal_mol
    for target in ['amine_contact','opposite_contact']:
        t=p[p.nominated==target].iloc[0]
        text+=f'#### {NAMES[target]}：★ {t.pair_label}\n\n'
        for metric in PAIR_METRICS:
            text+=f'**{METRICS[metric]}**\n\n'
            for pool,label in [('distance_and_chemistry','種類・距離を絞った比較'),('same_category','分子内／分子間の区分を合わせた広い比較')]:
                candidates=p[pool_mask(p,t,pool,config)]
                g=rank_context(candidates,metric,'pair',[t.pair])
                rank=int(g.loc[g.is_target,'rank'].iloc[0])
                intro=f'{label}：指定接触は{rank}/{len(candidates)}位'
                if pool=='same_category':text+=f'<details>\n<summary>{intro}（上位全件・指定接触・直下3件）</summary>\n\n'
                else:text+=intro+'。\n\n'
                rows=[]
                for r in g.itertuples():
                    note=r.interpretation_ja
                    if metric=='pair_atomic_work_kcal_mol':
                        note=f'{r.element_a}{r.atom_a}の仕事={fmt(atomvalues[r.atom_a])} kcal/mol、{r.element_b}{r.atom_b}の仕事={fmt(atomvalues[r.atom_b])} kcal/mol。'+note
                    rows.append([r.rank,'**★ '+r.pair_label+'**' if r.is_target else r.pair_label,
                        f'**{fmt(getattr(r,metric))}**' if r.is_target else fmt(getattr(r,metric)),
                        f'{r.mean_distance_A:.3f}',fmt(r.mean_signed_projection_mEh_bohr),
                        f'{r.contact_type}／{r.topology}',f'{r.role_a} ＋ {r.role_b}',note])
                text+=table(['順位','原子対',METRICS[metric],'平均原子間距離 / Å','接触方向の力の符号付き平均Q / mEh/Bohr','種類・結合関係','原子の位置','個別の読み方'],rows)+'\n'
                if pool=='same_category':text+='</details>\n\n'
    return text


def graph_for(mapping,edges):
    graph={int(a):[] for a in mapping.atom}
    for r in edges.itertuples():graph[r.atom_a].append(r.atom_b);graph[r.atom_b].append(r.atom_a)
    return {a:sorted(v) for a,v in graph.items()}


def atom_roles(mapping,graph):
    m=mapping.set_index('atom');roles={}
    core={1:'クマリン縮合ベンゼン環C',2:'クマリン縮合ベンゼン環C',3:'クマリン縮合部C',4:'環内O18に結合するクマリン環C',
          5:'クマリン縮合ベンゼン環C',6:'置換アミノN14に結合する環C',10:'クマリン不飽和部C',12:'クマリンのラクトンカルボニルC',
          13:'反応C20に結合するクマリン骨格C',14:'基質のジアルキルアミノN',15:'基質アミノ置換基のN隣接CH2',
          18:'クマリンの環内O',19:'クマリンのラクトンカルボニルO',20:'求核攻撃を受ける反応カルボニルC',
          21:'反応カルボニルC20のO',22:'C20と置換ベンゼンC23をつなぐエステルO',23:'O22につながる置換ベンゼンipso C',
          24:'アミン側ortho C',25:'反対側ortho C',26:'置換ベンゼンmeta C',27:'置換ベンゼンmeta C',28:'置換ベンゼンpara C'}
    for a,r in m.iterrows():
        if r.site.startswith('core') and a in core:role=core[a]
        elif r.site.startswith('X'):
            c=int(r.site[1:]);pos={24:'アミン側ortho',25:'反対側ortho',26:'meta',27:'meta',28:'para'}[c]
            role=f'置換ベンゼンC{c}上の{pos} {r.element}'
        elif r.site=='amineN':role='求核剤アミンN'
        elif r.site=='NH_near':role='TSでX24に近い方のアミンN–H'
        elif r.site=='NH_other':role='もう一方のアミンN–H'
        elif r.site.startswith('amineC'):
            role={'amineC1':'求核剤のN隣接CH2','amineC2':'求核剤アルキル鎖の中央CH2','amineC3':'求核剤アルキル鎖の末端CH3'}[r.site]
        elif r.site.startswith('sub_alkylC'):
            role={'sub_alkylC1':'基質アミノ置換基の末端CH3（枝1）','sub_alkylC2':'基質アミノ置換基のN隣接CH2（枝2）','sub_alkylC3':'基質アミノ置換基の末端CH3（枝2）'}[r.site]
        else:role=''
        roles[a]=role
    for a,r in m.iterrows():
        if not roles[a]:
            assert r.element=='H' and len(graph[a])==1
            b=graph[a][0];roles[a]=f'{m.loc[b,"element"]}{b}（{roles[b]}）上のH'
    return roles


def atom_label(m,a):return f'{m.loc[a,"element"]}{a}'


def shortest_path(graph,start,end):
    queue=[(start,[start])];seen={start}
    for a,path in queue:
        if a==end:return path
        for b in graph[a]:
            if b not in seen:seen.add(b);queue.append((b,path+[b]))
    return []


def chemical_contact(row,m,graph):
    a,b=int(row.atom_a),int(row.atom_b);sites={row.site_a,row.site_b};ea,eb=row.element_a,row.element_b
    if sites=={'core20','amineN'}:kind='形成中のC20–求核N結合'
    elif 'H' in [ea,eb] and any(e in ['O','N','F'] for e in [ea,eb]):
        h,acc=(a,b) if ea=='H' else (b,a);donor=graph[h][0]
        kind=f'{m.loc[donor,"element"]}–H···{m.loc[acc,"element"]}接触'
    elif ea==eb=='H':kind='H···H近接'
    elif set([ea,eb])=={'F','O'}:kind='F···O近接（水素を含まない）'
    elif set([ea,eb])=={'F','C'}:kind='C···F近接'
    elif ea==eb=='C':kind='C···C近接'
    elif ea==eb=='O':kind='O···O近接'
    else:kind=f'{ea}···{eb}近接'
    path=shortest_path(graph,a,b)
    top='分子間' if not path else f'{len(path)-1}結合離れ（1–{len(path)}）'
    h=acc=donor=None
    if 'H' in [ea,eb] and any(e in ['O','N','F'] for e in [ea,eb]):
        h,acc=(a,b) if ea=='H' else (b,a);donor=graph[h][0]
    return kind,top,path,h,acc,donor


def geometry(frame,n):
    frame=frame.sort_values(['signed_irc_coordinate','atom_number'])
    r=frame.signed_irc_coordinate.unique();xyz=frame[['x','y','z']].to_numpy().reshape(len(r),57,3)
    forces=frame[['fx','fy','fz']].to_numpy().reshape(len(r),57,3)
    cn=np.linalg.norm(xyz[:,19]-xyz[:,n-1],axis=1);good=(r<0)&np.isfinite(forces).all(axis=(1,2))
    assert np.all(np.diff(cn[good])<0)
    return r,xyz,forces,cn,good


def window_average(values,cn,good,start,end):
    grid=np.linspace(end,start,501)
    return float(np.trapezoid(np.interp(grid,cn[good][::-1],values[good][::-1]),grid)/(start-end))


def enrich_pairs(pairs,mapping,edges,windows):
    enriched=[];dictionary=[]
    for system,g in pairs[pairs.eligible].groupby('system',sort=False):
        mm=mapping[mapping.system==system];m=mm.set_index('atom');graph=graph_for(mm,edges[edges.system==system]);roles=atom_roles(mm,graph)
        n=int(mm.loc[mm.site=='amineN','atom'].iloc[0])
        r,xyz,forces,cn,good=geometry(pd.read_csv(OUT/f'csv/parsed_oriented_{system}.csv.gz'),n);ts=np.argmin(abs(r))
        for a,mr in m.iterrows():
            dictionary.append(dict(system=system,atom=a,element=mr.element,site=mr.site,role_ja=roles[a],
                neighbors=', '.join(atom_label(m,b) for b in graph[a]),target_atom=mr.site in ATOM_TARGETS))
        for row in g.itertuples():
            a,b=row.atom_a,row.atom_b;start,end=windows[row.window]
            kind,top,path,h,acc,donor=chemical_contact(row,m,graph)
            delta=xyz[:,b-1]-xyz[:,a-1];dist=np.linalg.norm(delta,axis=1);unit=delta/dist[:,None]
            q=.5*np.einsum('ij,ij->i',forces[:,b-1]-forces[:,a-1],unit)*1000
            meanq=window_average(q,cn,good,start,end)
            # Absolute value is taken after interpolation, matching the original definition.
            grid=np.linspace(end,start,501)
            meanabs=float(np.trapezoid(abs(np.interp(grid,cn[good][::-1],q[good][::-1])),grid)/(start-end))
            if not np.isclose(meanabs,row.mean_abs_projection_mEh_bohr,atol=1e-9):raise ValueError('Projection does not reproduce the source value')
            angle=anglemean=np.nan;dh=''
            if h is not None:
                u=xyz[:,donor-1]-xyz[:,h-1];v=xyz[:,acc-1]-xyz[:,h-1]
                ang=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',u,v)/np.linalg.norm(u,axis=1)/np.linalg.norm(v,axis=1),-1,1)))
                angle=float(ang[ts]);anglemean=window_average(ang,cn,good,start,end)
                dh=f'{atom_label(m,donor)}–{atom_label(m,h)}···{atom_label(m,acc)}'
            notes=[]
            if kind=='形成中のC20–求核N結合':notes.append('主反応の結合形成を直接含むため、補助接触との大小を分けて読む')
            elif 'core20' in {row.site_a,row.site_b}:notes.append('反応C20の全原子力がこの対にも投影される')
            elif 'amineN' in {row.site_a,row.site_b}:notes.append('求核Nの全原子力がこの対にも投影される')
            if len(path)==4:notes.append('1–4の骨格内近接で、独立した非共有結合の強さとは扱わない')
            if h is not None:
                notes.append(f'{dh}の共通区間平均角は{anglemean:.1f}°。距離・角度だけで水素結合とは断定しない')
            elif kind=='H···H近接':notes.append('水素結合の供与体–受容体対ではない')
            elif 'F···O' in kind:notes.append('この分類や力の大きさだけではF–Oの斥力を判定できない')
            if not notes:notes.append('骨格・配向変化を含む全原子力の投影であり、この対だけの相互作用ではない')
            enriched.append(dict(row._asdict(),pair_label=f'{atom_label(m,a)}···{atom_label(m,b)}',role_a=roles[a],role_b=roles[b],
                neighbors_a=', '.join(atom_label(m,c) for c in graph[a]),neighbors_b=', '.join(atom_label(m,c) for c in graph[b]),
                contact_type=kind,topology=top,bond_path='–'.join(atom_label(m,c) for c in path) or '別分子',
                includes_reaction_C20='core20' in {row.site_a,row.site_b},includes_nucleophile_N='amineN' in {row.site_a,row.site_b},
                includes_O21='core21' in {row.site_a,row.site_b},mean_signed_projection_mEh_bohr=meanq,mean_closing_projection_mEh_bohr=.5*(meanabs-meanq),mean_opening_projection_mEh_bohr=.5*(meanabs+meanq),
                DH_A=dh,mean_DHA_angle_deg=anglemean,TS_DHA_angle_deg=angle,TS_distance_A=dist[ts],interpretation_ja='。'.join(notes)))
    return pd.DataFrame(enriched).drop(columns='Index'),pd.DataFrame(dictionary)


def enumerate_controls(atoms,pairs,config):
    rows=[]
    for (system,window),g in atoms.groupby(['system','window'],sort=False):
        for target in ATOM_TARGETS:
            t=g[g.site==target].iloc[0]
            for metric in ATOM_METRICS:
                stronger=g[(g.atom!=t.atom)&(g[metric].abs()>abs(t[metric])+1e-12)]
                for r in stronger.itertuples():
                    rows.append(dict(kind='atom',system=system,window=window,target=target,target_atom=t.atom,pool='all_57_atoms',metric=metric,
                        target_value=t[metric],control_atom=r.atom,control_pair='',control_value=getattr(r,metric),
                        control_rank=int((g[metric].abs()>abs(getattr(r,metric))+1e-12).sum())+1))
    for (system,window),g in pairs.groupby(['system','window'],sort=False):
        for _,t in g[g.nominated!=''].iterrows():
            for pool in POOLS:
                eligible=g[pool_mask(g,t,pool,config)&(g.pair!=t.pair)]
                for metric in PAIR_METRICS:
                    stronger=eligible[eligible[metric].abs()>abs(t[metric])+1e-12]
                    for r in stronger.itertuples():
                        rows.append(dict(kind='contact',system=system,window=window,target=t.nominated,target_atom=np.nan,pool=pool,metric=metric,
                            target_value=t[metric],control_atom=np.nan,control_pair=r.pair,control_value=getattr(r,metric),
                            control_rank=int((eligible[metric].abs()>abs(getattr(r,metric))+1e-12).sum())+1))
    result=pd.DataFrame(rows)
    # Audit all original ranks, including zero-competitor cases in every pool/window.
    ranks=pd.read_csv(OUT/'csv/contact_controls_within_ranks.csv')
    counts=result.groupby(['kind','system','window','target','pool','metric']).size()
    for r in ranks.itertuples():
        count=int(counts.get((r.kind,r.system,r.window,r.target,r.pool,r.metric),0))
        if count!=r.stronger_controls:raise ValueError(f'Incomplete competitor list: {r}')
    return result


def system_report(system,atoms,pairs,dictionary,controls,config):
    text=f'# {system}：着目する原子・接触と、その上位・直下3件の比較\n\n[全4系の考察と比較表](higher_ranked_controls_summary_ja.md)\n\n'
    text+='主区間はC–N 2.50→2.10 Å。順位は各量の絶対値が大きい順。「上位」は相互作用エネルギーが大きいという意味ではない。ΔG‡や障壁との相関では選別していない。\n\n'
    text+=structure_section(system)+context_tables(system,atoms,pairs,dictionary,config)
    text+='## データの追跡\n\n[全原子の対応・結合先](../csv/higher_rank_atom_dictionary.csv)、[各接触の距離・角度・符号付きQ](../csv/higher_rank_contact_identities.csv)、[全3区間・全7比較条件の上位一覧](../csv/higher_rank_competitors.csv)。\n'
    (OUT/f'reports/higher_ranked_controls_{system}_ja.md').write_text(text)


def main():
    cfg=yaml.safe_load((HERE/'config/contact_specificity.yaml').read_text())
    windows=yaml.safe_load((HERE/'config/two_sided_contacts.yaml').read_text())['windows']
    atoms=work_percentages(pd.read_csv(OUT/'csv/contact_controls_atoms.csv'))
    pairs=pd.read_csv(OUT/'csv/contact_controls_pairs.csv').fillna({'nominated':''})
    mapping=pd.read_csv(OUT/'csv/contact_controls_mapping.csv');edges=pd.read_csv(OUT/'csv/contact_controls_covalent_edges.csv')
    identified,dictionary=enrich_pairs(pairs,mapping,edges,windows)
    dictionary.to_csv(OUT/'csv/higher_rank_atom_dictionary.csv',index=False)
    identified.to_csv(OUT/'csv/higher_rank_contact_identities.csv',index=False)
    competitors=enumerate_controls(atoms,pairs,cfg)
    # Attach full chemical descriptions to each occurrence in a ranked list.
    lookup=identified[['system','window','pair','pair_label','role_a','role_b','contact_type','topology','bond_path','mean_distance_A','mean_DHA_angle_deg','mean_signed_projection_mEh_bohr','interpretation_ja']]
    competitors=competitors.merge(lookup,left_on=['system','window','control_pair'],right_on=['system','window','pair'],how='left').drop(columns='pair')
    competitors=competitors.merge(dictionary[['system','atom','role_ja','neighbors']],left_on=['system','control_atom'],right_on=['system','atom'],how='left').drop(columns='atom')
    competitors.to_csv(OUT/'csv/higher_rank_competitors.csv',index=False)
    from numbered_structures import main as draw_numbered_structures
    draw_numbered_structures()
    for system in SYSTEMS:system_report(system,atoms,identified,dictionary,competitors,cfg)
    write_overview(atoms,identified,competitors,cfg)
    print(f'{len(dictionary)} atom identities; {len(identified)} eligible contacts across 3 windows; {len(competitors)} higher-ranking entries (repeated across targets/metrics/pools).')


def write_overview(atoms,pairs,controls,config):
    primary=pairs[pairs.window==PRIMARY]
    a=atoms[atoms.window==PRIMARY]
    c=controls[(controls.window==PRIMARY)&(controls.kind=='contact')]
    text='''# 着目する原子・接触の順位と、比較対象の中身を調べた結果

作成日：2026-10-20（JST）。ΔG‡を説明指標にせず、原子の仕事・力と接触方向の仕事・力について、指定対象より上位のものを一つずつ同定した。新しいDFTは行っていない。

## 今回の確認方法

主区間はC–N 2.50→2.10 Å。指定した4原子（環のアミン側F/H、選択アミンH、反対側F/H、O21）は全57原子と比較した。指定した2接触は、従来と同じ分子内／分子間区分・距離・結合関係で選んだ比較対象と比べた。上位の原子・接触を都合よく除外せず、すべて残した。

比較候補は、元の分子内で直接結合する1–2対と、2本の結合でつながる1–3対を除き、区間内の最短距離が4.0 Å以下の対とした。1–4対は広い比較に含める。アミン側は分子間、反対側は基質分子内で比較する。種類・距離を絞る場合、アミン側はH···H/F/O/N、反対側はH/F···O/Nを残し、指定接触との平均距離差を±0.75 Å以内にする。形成中のC20–求核Nは元の分子間として広い比較に含まれるので、表では主反応の結合形成と明記した。

各原子の結合先を確認し、各原子対について、原子番号・元素・分子内の位置・最短結合経路・平均距離を記載した。H···O/N/Fについては、Hが結合する原子Dを特定してD–H···A角も調べた。さらに、これまで順位に使っていた力の絶対値に加え、符号付きの接触方向の力Qも既存の座標・力から計算した。

「上位」は定義した数値の絶対値が大きいことを意味する。原子の仕事、2原子の仕事の和、接触方向の仕事、力の平均絶対値は別の量として扱う。接触方向のQも、2原子間だけの力ではなく、両原子に働く全原子力の差から求める。

## 個別の全件一覧

'''
    for system in atoms.system.unique():
        text+=f'- [{system}：番号付き構造式・着目対象・上位全件・直下3件](higher_ranked_controls_{system}_ja.md)\n'
    text+='\n'+(HERE/'docs/force_work_definitions_ja.md').read_text()+'\n'+percentage_summary()
    text+='''
このまとめの後半に、4系すべての番号付き構造式と比較表を掲載した。各系の個別ページにも同じ表を置いた。★が着目対象で、上位だけでなく指定対象自身と直下3件も同じ表に載せた。広い比較の表は折りたたんである。

## 1. 原子の仕事で上位なのは、主反応と求核剤の動きを担う原子だった

全4系で原子仕事の1位は求核剤N、2位は反応C20だった。その後には求核剤のアルキルC、もう一方のアミンH、C20に結合するクマリン骨格C13などが並ぶ。したがって、指定F/HやO21より仕事が大きい原子の多くは、反応中心や求核剤全体の移動・変形に関係する原子である。補助接触がこれらすべてを上回る必要はない。

RAFA15のO21を例に、上位8原子・O21自身・直下3原子を示す。4系の全指定原子と力の表は後半に掲載した。

'''
    g=a[a.system=='RAFA15'].sort_values('atomic_work_kcal_mol',key=abs,ascending=False)
    d=pd.read_csv(OUT/'csv/higher_rank_atom_dictionary.csv');rd=d[d.system=='RAFA15'].set_index('atom')
    o=g[g.site=='core21'].iloc[0]
    h=rank_context(g,'atomic_work_kcal_mol','atom',[o.atom])
    text+=table(['順位','原子','原子の仕事 / kcal/mol','原子仕事／IRC端点→TSのΔE‡ / %','原子仕事／区間ΔE / %','位置・役割'],[[r.rank,('**★ O21**' if r.is_target else f'{r.element}{r.atom}'),fmt(r.atomic_work_kcal_mol),fmt(r.atomic_work_pct_IRC_endpoint_barrier),fmt(r.atomic_work_pct_window_energy_rise),rd.loc[r.atom,'role_ja']] for r in h.itertuples()])
    text+='''
2原子の仕事を足した指標では、仕事の大きいNやC20などと組み合わせた多くの対が上位になる。これは、その組合せの相互作用が強いことを直接示さない。個別表では、2原子それぞれの仕事も併記し、どちらが和を担うかを確認できるようにした。

## 2. 反対側の接触を4系で比較する

まず4系について、広い比較の順位と、その上位にある原子対の内訳を示す。原子対を一つずつ確認する表は後半にあり、指定接触と直下3件も含めている。

'''
    for metric in PAIR_METRICS[1:]:
        rows=[]
        for system in SYSTEMS:
            g=primary[primary.system==system];t=g[g.nominated=='opposite_contact'].iloc[0]
            broad=g[pool_mask(g,t,'same_category',config)]
            upper=broad[broad[metric].abs()>abs(t[metric])+1e-12]
            others=upper[~upper.includes_reaction_C20]
            rows.append([system,t.pair_label,fmt(t[metric]),
                         f'{len(upper)+1}/{len(broad)}',int(upper.includes_reaction_C20.sum()),
                         len(others),'、'.join(others.pair_label) or 'なし'])
        text+=table(['分子系','指定接触',METRICS[metric],'順位／候補数','上位のうちC20を含む対の数',
                     '上位のうちC20を含まない対の数','C20を含まない上位対'],rows)+'\n'
    text+='''
MoriPF・RAFA16では、仕事の上位にC20を含まない対も4つあり、どちらもC12···O22、C13···F29、O19···O21、O22···F29だった。RAFA14でもO22···H56などの骨格内近接が上位にある。**「指定接触より上位はすべてC20を含む」という特徴はRAFA15に限られる。** 反対側のH/F···O21を、4系すべてで一律に主要な接触と呼ぶことはできない。

### RAFA15の反対側H56···O21を詳しく見る

反対側の指定接触O21···H56は、接触方向の仕事が0.009085 kcal/mol、力の平均絶対値が0.513777 mEh/Bohr。広い基質分子内の比較では、それぞれ169候補中5位、3位だった。

### 接触方向の仕事：上位4対・指定接触・直下3対

'''
    g=primary[primary.system=='RAFA15'];t=g[g.nominated=='opposite_contact'].iloc[0]
    for metric,label in [('radial_work_kcal_mol',METRICS['radial_work_kcal_mol']),('mean_abs_projection_mEh_bohr',METRICS['mean_abs_projection_mEh_bohr'])]:
        if metric.startswith('mean'):text+='\n### 接触方向の力：上位2対・指定接触・直下3対\n\n'
        broad=g[pool_mask(g,t,'same_category',config)]
        h=rank_context(broad,metric,'pair',[t.pair])
        text+=table(['順位','原子対',label,'平均距離 / Å','結合関係','原子の位置'],[
            [r.rank,'**★ '+r.pair_label+'**' if r.is_target else r.pair_label,fmt(getattr(r,metric)),f'{r.mean_distance_A:.3f}',r.topology,f'{r.role_a} ＋ {r.role_b}'] for r in h.itertuples()])
    text+='''
上位にあるのは、C20とortho置換原子または環Cの組合せである。C20は求核攻撃を受ける反応点なので、その全原子力が複数の方向に投影されている可能性がある。これらを「H56···O21より強い独立した非共有結合」と呼ぶのは適切ではない。

この中身を踏まえると、**H56···O21が全体5位・3位であることは、その接触を軽視する理由にはならない**。従来の種類・距離を絞った13候補では、仕事・力ともに1位という結果も維持される。一方で、この数値だけでH···Oの安定化エネルギーを決めることはできない。平均C–H···O角は約90.4°であり、直線的な水素結合の配置とも異なる。

## 3. アミン側F···Hより上位の接触は、力の向きまで同じではなかった

Q<0は距離を縮める向き、Q>0は開く向きの全原子力の射影である。力の大きさだけの順位は両方向を区別しない。

まず、指定したアミン側接触の符号を示す。

'''
    chosen=primary[primary.nominated=='amine_contact']
    text+=table(['分子系','指定原子対','接触方向の力の平均絶対値 / mEh/Bohr','接触方向の力の符号付き平均Q / mEh/Bohr','平均D–H···A角 / °'],[
        [r.system,r.pair_label,fmt(r.mean_abs_projection_mEh_bohr),fmt(r.mean_signed_projection_mEh_bohr),fmt(r.mean_DHA_angle_deg)] for r in chosen.itertuples()])
    text+='''
Fを含むMoriPF・RAFA15・RAFA16では、指定接触の平均Qは負だった。RAFA14の指定接触はF···HではなくH···Hで、平均Qは正だった。この違いは、平均絶対値だけの順位では表せない。

以下は、種類・距離を絞った比較における上位全件・指定接触・直下3件である。★が指定接触で、順位はこの絞った候補内のもの。

'''
    for system in atoms.system.unique():
        g=primary[primary.system==system];t=g[g.nominated=='amine_contact'].iloc[0]
        h=rank_context(g[pool_mask(g,t,'distance_and_chemistry',config)],'mean_abs_projection_mEh_bohr','pair',[t.pair])
        text+=f'### {system}\n\n'
        text+=table(['順位','原子対','接触の種類','平均距離 / Å','平均D–H···A角 / °','接触方向の力の平均絶対値 / mEh/Bohr','接触方向の力の符号付き平均Q / mEh/Bohr'],[
            [r.rank,'**★ '+r.pair_label+'**' if r.is_target else r.pair_label,r.contact_type,f'{r.mean_distance_A:.3f}',fmt(r.mean_DHA_angle_deg),fmt(r.mean_abs_projection_mEh_bohr),fmt(r.mean_signed_projection_mEh_bohr)] for r in h.itertuples()])+'\n'
    text+='''
MoriPFとRAFA16では、この絞った比較で上位にある5対はすべて平均Qが正で、指定F···Hの負のQとは向きが異なる。N–H···OやC–H···Oという元素の組合せでも、この結果から「F···Hより強い水素結合」とは結論できない。

RAFA15では上位6対のうち、負のQを示すのは**同じF31ともう一方のアミンH45**の接触だった。指定F31···H44は平均距離2.603 Å、平均絶対値0.096767、平均角99.4°。F31···H45は距離2.655 Åと少し長いが、平均絶対値0.108108、平均角96.4°だった。いずれも直線的なN–H···F配置ではないが、Fの影響を追う際に片方のHだけに限定すると情報を落とす可能性がある。既存の「TSで近いHを固定する」定義は変更せず、もう一方も別の候補として示した。

RAFA15で値の大きいN–H···O候補の平均角は約66〜69°で、平均Qは正だった。これらは単純な水素結合強度の競合相手として扱うより、アミンと反応部位の配向・変形を含む量として検討する方がよい。

## 4. ここから指標をどう考え直すか

**以前の「アミン側F···Hは力の順位が低いので、強い接触という根拠が弱い」という評価は、比較対象の役割と力の向きを区別しておらず、粗すぎた。** 上位の中には主反応の結合形成、骨格内の1–4近接、距離を開く向きの大きな成分が含まれていた。

次の評価では、全候補を残したまま、少なくとも①反応C20・求核Nを直接含む対、②骨格内の1–4対、③それ以外の接触を区別する必要がある。また、Qの正負を平均絶対値にまとめず、縮める向き・開く向きの成分を別々に見る方が今回の仮説に対応する。

今回のCSVには、共通区間での平均 max(−Q,0) と平均 max(Q,0) も併記した。単位はmEh/Bohrで、両者の和が従来の平均絶対値、開く成分−縮める成分が平均Qになる。**これは今回の確認を受けて追加した記述用の分解であり、独立した検証やペア間の引力・斥力の分解ではない。** この新しい量で有利な候補だけを残した順位や、障壁に合う指標の選別は行っていない。

現在優先して追う対象は、アミン側では **Fと2つのN–Hそれぞれの接触**、反対側では **H/F···O21と、それより上位のC20を含む対、O22···H/Fなどの骨格内近接**である。距離・角度・符号付きQ・接触方向の仕事を同じC–N軸で追えば、Fがアミンの配置を保つ方向に働くという仮説と、単に反応中心の全体的な変形を拾っている可能性を、具体的な対照を使って比較できる。

## 原子番号の意味についての補足

旧説明で「クマリンO」とまとめていたO21は、結合関係では**反応カルボニルC20に結合するO**である。クマリン環内のOはO18、ラクトンのカルボニルOはO19、C20と置換ベンゼンC23をつなぐOはO22。この区別を各行の説明に反映した。結合グラフは反応物側IRC端点の距離に基づく既存の推定を使用し、形成中のC20–求核Nは元の分子間として別に識別している。原子の位置・結合先は確認できるが、各近接に独立した相互作用があることまで確認したわけではない。

'''
    text+='## 4系すべての構造式と比較表\n\n以下は共通C–N区間2.50→2.10 Åでの比較。原子は仕事と力の2指標、接触は2原子の仕事の合計・接触方向の仕事・接触方向の力の3指標を、各系について示す。\n\n'
    for system in SYSTEMS:
        text+=f'### {system}\n\n'
        text+=structure_section(system).replace('## 原子番号付きの構造式','#### 原子番号付きの構造式')
        body=context_tables(system,atoms,pairs,d,config)
        body=body.replace('#### ','##### ').replace('### 原子','#### 原子').replace('### 接触','#### 接触')
        text+=body
    text+='## 保存した範囲と検証\n\n'
    text+=f'全4系の{len(d)}原子の対応、全3区間の{len(pairs)}接触候補を同定した。上位一覧は{len(controls)}行で、同じ対が区間・指標・比較条件ごとに重複して現れる。独立した相互作用がこの数だけあるという意味ではない。主区間の広い接触比較だけでは{len(c[c.pool=="same_category"])}行である。\n\n'
    text+='''
- [原子の位置・結合先の辞書](../csv/higher_rank_atom_dictionary.csv)
- [接触ごとの位置・結合経路・距離・角度・Q](../csv/higher_rank_contact_identities.csv)
- [全3区間・全7条件について指定対象を上回る全件](../csv/higher_rank_competitors.csv)

従来の順位表にある全条件について、今回列挙した件数と「指定対象より大きい比較対象数」の一致を確認した。符号付きQの計算と同時に、既存の平均絶対値を全接触・全区間で再現することも確認した。力や仕事の既存値、比較対象の選び方、原子対応は変更していない。角度と力は同じC–N区間で平均し、TS角だけで代用していない。

再生成：`python analysis4_force/higher_ranked_controls.py`。既存の解析CSVとIRC座標・力を使い、量子化学計算は実行しない。
'''
    (OUT/'reports/higher_ranked_controls_summary_ja.md').write_text(text)



if __name__=='__main__':
    main()
    manifest_path=OUT/'reports/output_manifest.json'
    manifest=[dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
              for p in sorted(OUT.rglob('*')) if p.is_file() and p!=manifest_path]
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
