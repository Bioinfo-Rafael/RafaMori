"""Compare signed work using complete, geometry-audited atom correspondence."""
from pathlib import Path
from functools import lru_cache
import json
import pandas as pd
from contact_specificity import pool_mask
from numbered_structures import SYSTEMS

OUT=Path(__file__).resolve().parent/'outputs'


@lru_cache(maxsize=1)
def correspondence_lookup():
    d=pd.read_csv(OUT/'csv/audited_atom_correspondence.csv')
    return d.set_index(['system','original_site']).cross_site.to_dict()


def site_key(system,site):
    return correspondence_lookup()[system,site]


def pair_key(system,a,b):
    return '||'.join(sorted([site_key(system,a),site_key(system,b)]))


def build_cross_system_tables(atoms,identified,dictionary,config):
    from higher_ranked_controls import PRIMARY,ATOM_TARGETS,rank_context,table,fmt
    from atom_correspondence import build_correspondence
    build_correspondence()
    correspondence_lookup.cache_clear()
    atom_metric='atomic_work_kcal_mol';pair_metric='radial_work_kcal_mol'
    atoms=atoms[atoms.window==PRIMARY].copy()
    allpairs=pd.read_csv(OUT/'csv/contact_controls_pairs.csv').fillna({'nominated':''})
    allpairs=allpairs[allpairs.window==PRIMARY].copy()
    atoms['cross_key']=[site_key(r.system,r.site) for r in atoms.itertuples()]
    allpairs['cross_key']=[pair_key(r.system,r.site_a,r.site_b) for r in allpairs.itertuples()]
    atomlookup={s:g.set_index('cross_key') for s,g in atoms.groupby('system')}
    pairlookup={s:g.set_index('cross_key') for s,g in allpairs.groupby('system')}
    assert all(g.index.is_unique for g in [*atomlookup.values(),*pairlookup.values()])
    roles=dictionary.set_index(['system','atom'])
    labels={}
    special={'X24':'アミン側ortho F/H','X25':'反対側ortho F/H','X26':'meta位26のF/H',
             'X27':'meta位27のF/H','X28':'para位のF/H','amineN':'求核アミンN',
             'NH_near':'選択アミンH','NH_other':'もう一方のアミンH',
             'amineC1':'求核剤のN隣接CH₂','amineC2':'求核剤の中央CH₂','amineC3':'求核剤の末端CH₃',
             'sub_alkylC1':'基質アミノ置換基の末端CH₃（枝1）',
             'sub_alkylC2':'基質アミノ置換基のN隣接CH₂（枝2）',
             'sub_alkylC3':'基質アミノ置換基の末端CH₃（枝2）'}
    for r in atoms.itertuples():
        labels[r.cross_key]=special.get(r.site,
            f'{r.element}{r.atom}' if r.cross_system_mapped else
            f'対応H（MoriPF H{r.cross_key.removeprefix("H_ref")}）')
    records=[];audits=[]

    def emit(kind,target,pool,selections,candidates,lookup,metric):
        union=set().union(*(set(g.cross_key) for g in selections.values()))
        # Common ranking does not exist: sort by best rank among systems where
        # the row was selected. Signed work and original per-system ranks stay intact.
        best={k:min(int(g.loc[g.cross_key==k,'rank'].iloc[0]) for g in selections.values()
                    if k in set(g.cross_key)) for k in union}
        ordered=sorted(union,key=lambda k:(best[k],k))
        rows=[]
        for key in ordered:
            cells=[];details=[];target_any=False
            for system in SYSTEMS:
                selected=selections[system];cand=candidates[system];data=lookup[system]
                present=key in data.index;sel=key in set(selected.cross_key)
                if not present:
                    cells.append('—（原子の個体対応未確定）')
                    details.append(f'{system}：同じ原子と保証できないため比較しない')
                    records.append(dict(kind=kind,target=target,pool=pool,cross_key=key,system=system,
                                        selected=False,correspondence_available=False))
                    continue
                r=data.loc[key];value=float(r[metric])
                incandidate=key in set(cand.cross_key)
                rank=int((cand[metric].abs()>abs(value)+1e-12).sum())+1 if incandidate else None
                nominated=r.site in ATOM_TARGETS if kind=='atom' else r.nominated==target
                target_any|=nominated
                if kind=='atom':
                    name=f'{r.element}{int(r.atom)}';role=roles.loc[(system,int(r.atom)),'role_ja']
                else:
                    name=f'{r.element_a}{int(r.atom_a)}···{r.element_b}{int(r.atom_b)}'
                    role=''
                reason='着目対象' if nominated else ('選択範囲内' if sel else '他系で選ばれたため併記')
                if kind=='atom' and not nominated and sel:
                    targets=cand[cand.site.isin(ATOM_TARGETS)]
                    reason=('着目原子の少なくとも1つより上位' if
                            (targets[metric].abs()<abs(value)-1e-12).any()
                            else '最下位の着目原子の直下3件（境界同値を含む）')
                if kind=='contact' and incandidate and not nominated:
                    t=cand[cand.nominated==target].iloc[0]
                    reason=('着目接触より上位' if abs(value)>abs(t[metric])+1e-12 else
                            '着目接触と同順位' if abs(abs(value)-abs(t[metric]))<=1e-12 else
                            '着目接触の直下3件（境界同値を含む）' if sel else '直下3件より下位')
                ranklabel=f'{rank}/{len(cand)}位' if rank is not None else '条件外・順位なし'
                if not incandidate:
                    reason=('候補条件外（距離・結合関係）' if not bool(r.eligible)
                            else '種類・平均距離の条件外')
                cells.append(f'{"**★ " if nominated else ""}{fmt(value)}{"**" if nominated else ""}（{name}）')
                details.append(f'{system}：{ranklabel}；{reason}'+(f'。{role}' if kind=='atom' else ''))
                records.append(dict(kind=kind,target=target,pool=pool,cross_key=key,system=system,
                    selected=sel,correspondence_available=True,in_candidate_pool=incandidate,
                    nominated=nominated,atom_or_pair=name,work_kcal_mol=value,rank=rank,
                    pool_size=len(cand),selection_detail=reason))
            name=labels[key] if kind=='atom' else ' ··· '.join(labels[k] for k in key.split('||'))
            if target_any:name='**★ '+name+'**'
            rows.append([name,*cells,' ／ '.join(details)])
        audits.append(dict(kind=kind,target=target,pool=pool,union_count=len(union),
            selected_keys={s:sorted(g.cross_key) for s,g in selections.items()},union_keys=sorted(union)))
        return table(['対応する原子' if kind=='atom' else '対応する接触',
                      *[f'{s}：{"原子の仕事" if kind=="atom" else "接触方向の仕事"} / kcal/mol' for s in SYSTEMS],
                      '詳細説明'],rows)+'\n'

    text='## 分子系を横に並べた仕事の比較（和集合）\n\n'
    text+=(OUT.parent/'docs/work_comparison_groups_ja.md').read_text()+'\n\n'
    text+='共通のC–N距離2.50→2.10 Åについて、符号付きの仕事をMoriPF・RAFA14・RAFA15・RAFA16の列に並べた。MoriPFとの差には換算していない。★が着目対象。各セルには仕事と実際の原子番号を、詳細説明列には各系の候補内での絶対値順位／候補数と選択理由を示す。正の仕事はエネルギー上昇側、負は低下側の寄与であり、正負を区別して比較する。\n\n'
    text+='各系で「上位全件＋着目対象＋直下3件」を選び、その**和集合**を表の行にした。他系ではもっと下位になる行も削除せず、元の仕事と順位を示す。行の順は選択元の系での最も小さい順位を基準にした表示順で、4系共通の順位ではない。\n\n'
    text+='**原子対応を再確認し、全4系・全57原子の一対一対応を確認した。** Lysモデルとクマリン骨格は元素も含めて対応し、置換ベンゼンのF/Hは同じ結合位置として対応する。従来未対応としていた各系15個のアルキルHも、結合先と局所座標を共通IRC区間41点・TSで照合できた。以前の1点だけの行は対応する4系を統合し、フィットを再計算した。各系内部の仕事・順位・候補選択は変更していない。\n\n'
    text+='[全57位置の対応表と確認方法](atom_correspondence_audit_ja.md) · [原子対応CSV](../csv/audited_atom_correspondence.csv)。H_ref49などの位置ラベルは「MoriPFのH49に対応する位置」を意味し、各系の実際の原子番号はセル内に併記した。\n\n'
    text+='### 対応する原子の仕事\n\n4つの着目原子それぞれについて選択した行を、さらに4系で和集合にした。順位の母集団は各系の全57原子。各系での全57原子・合計・障壁に対する割合は後続の個別表にある。\n\n'
    selections={};candidates={}
    for s in SYSTEMS:
        g=atoms[atoms.system==s];candidates[s]=g
        selections[s]=rank_context(g,atom_metric,'atom',g.loc[g.site.isin(ATOM_TARGETS),'atom'].tolist())
    text+=emit('atom','four_target_atoms','all_57_atoms',selections,candidates,atomlookup,atom_metric)
    for target,title in [('amine_contact','アミン側'),('opposite_contact','反対側')]:
        text+=f'### {title}の接触方向の仕事\n\n'
        for pool,label in [('distance_and_chemistry','種類・平均距離を絞った比較'),('same_category','分子内／分子間の区分を合わせた広い比較')]:
            selections={};candidates={}
            for s in SYSTEMS:
                g=allpairs[allpairs.system==s];t=g[g.nominated==target].iloc[0]
                c=g[pool_mask(g,t,pool,config)];candidates[s]=c
                selections[s]=rank_context(c,pair_metric,'pair',[t.pair])
            text+=f'#### {label}\n\n'
            text+='他系で候補条件を外れる場合も、対応する原子対の計算値は掲載し、順位を「条件外」とする。原子の対応が確定しない場合とは区別する。\n\n'
            if pool=='same_category':text+='<details>\n<summary>広い比較の和集合を開く</summary>\n\n'
            text+=emit('contact',target,pool,selections,candidates,pairlookup,pair_metric)
            if pool=='same_category':text+='</details>\n\n'
    text+='[横並び表の全セルと選択理由CSV](../csv/cross_system_work_union.csv) · [各系の選択集合と和集合の監査JSON](cross_system_work_union_audit.json)\n\n'
    pd.DataFrame(records).to_csv(OUT/'csv/cross_system_work_union.csv',index=False)
    (OUT/'reports/cross_system_work_union_audit.json').write_text(json.dumps(audits,ensure_ascii=False,indent=2)+'\n')
    from work_union_plots import main as draw_work_union
    gallery=draw_work_union()
    (OUT/'reports/cross_system_work_comparison_ja.md').write_text('# 4系の原子仕事・接触方向の仕事の比較\n\n[全体の考察](higher_ranked_controls_summary_ja.md) · [散布図・棒グラフへ](#work-plots)\n\n'+text+gallery)
    return text+'\n[各原子・接触の散布図（ΔG‡と仕事）と4系の棒グラフ](cross_system_work_comparison_ja.md#work-plots)\n\n'
