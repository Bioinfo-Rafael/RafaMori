"""Report local prominence and specificity without selecting away controls."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting4 import SYSTEM_COLORS
from reporting import markdown as base_markdown
from contact_specificity import HERE,TARGETS,ENERGIES,pairkey,pool_mask

def markdown(df):
    # Small diagnostic works must not be rounded to apparent zeros.
    formatted=df.copy()
    for col in formatted.select_dtypes(include='floating'):
        formatted[col]=formatted[col].map(lambda value:f'{value:.6g}')
    return base_markdown(formatted)


def create_report(atoms,pairs,ranks,corr,spec,examples,perms,config,windows,energy):
    figdir=HERE/'outputs/figures';reportdir=HERE/'outputs/reports';primary=windows['primary_window'];systems=energy.index.tolist()
    a=atoms[atoms.window==primary];p=pairs[pairs.window==primary];r=ranks[ranks.window==primary];c=corr[corr.window==primary];s=spec[spec.window==primary]
    figs,axes=plt.subplots(2,2,figsize=(12,9))
    for ax,system in zip(axes.flat,systems):
        g=a[a.system==system].sort_values('atomic_work_kcal_mol',key=abs,ascending=False).reset_index(drop=True)
        ax.plot(np.arange(1,58),abs(g.atomic_work_kcal_mol),'o',color='lightgrey',ms=4,label='All 57 atoms')
        for site,color in zip(['X24','NH_near','X25','core21'],['#E45756','#F2A541','#54A24B','#4C78A8']):
            idx=g.index[g.site==site][0];v=abs(g.loc[idx,'atomic_work_kcal_mol'])
            ax.scatter(idx+1,v,color=color,s=55,zorder=4,label=f'{site}: atom {g.loc[idx,"atom"]}, rank {idx+1}')
        ax.set(yscale='log',xlabel='Rank by absolute atomic work (1 = largest)',ylabel='|Atomic work| (kcal/mol)',title=system);ax.legend(fontsize=8);ax.grid(alpha=.2)
    figs.suptitle('All atoms within each system; C-N 2.50 → 2.10 Å')
    figs.tight_layout();figs.savefig(figdir/'20_all_atom_work_controls.png',dpi=190);plt.close(figs)

    figs,axes=plt.subplots(2,2,figsize=(12,9))
    for row,target in enumerate(TARGETS):
        for col,(metric,label) in enumerate([('radial_work_kcal_mol','|Radial work diagnostic| (kcal/mol)'),('mean_abs_projection_mEh_bohr','Mean |Q| (mEh/Bohr)')]):
            ax=axes[row,col];vals=[];targetvals=[];ranklabels=[]
            for system in systems:
                local=p[p.system==system];t=local[local.nominated==target].iloc[0];control=local[pool_mask(local,t,'same_category',config)&(local.pair!=t.pair)]
                vals.append(abs(control[metric]).to_numpy());targetvals.append(abs(t[metric]));rr=r[(r.system==system)&(r.target==target)&(r.pool=='same_category')&(r.metric==metric)].iloc[0]
                ranklabels.append(f'{int(rr.magnitude_rank)}/{int(rr.n_controls)+1}')
            ax.boxplot(vals,tick_labels=systems,showfliers=True,flierprops={'marker':'.','markersize':3,'markerfacecolor':'grey','markeredgecolor':'grey'})
            for i,(v,system,ranklabel) in enumerate(zip(targetvals,systems,ranklabels)):
                ax.scatter(i+1,v,s=60,color=SYSTEM_COLORS[system],zorder=5);ax.annotate('rank '+ranklabel,(i+1,v),xytext=(8,5),textcoords='offset points',fontsize=8)
            ax.set(yscale='log',ylabel=label,title=target.replace('_',' ')+' / same intra- or intermolecular class');ax.grid(alpha=.2)
    figs.suptitle('Nominated contacts against all geometry-selected controls within each system\nBoxes: controls; colored points: nominated pair. Values are total-force diagnostics.')
    figs.tight_layout(rect=(0,0,1,.95));figs.savefig(figdir/'21_contact_strength_controls.png',dpi=190);plt.close(figs)

    figs,axes=plt.subplots(1,2,figsize=(12,5.5))
    for ax,target in zip(axes,TARGETS):
        key=pairkey(*TARGETS[target]);category=p[p.pair==key].category.iloc[0]
        common=c[(c.kind=='contact')&(c.metric=='radial_work_kcal_mol')&(c.barrier==ENERGIES[0])]
        allowed=set(p[(p.category==category)&p.eligible].pair);common=common[common.item.isin(allowed)]
        ax.scatter(common.median_magnitude,common.pearson_r,color='grey',alpha=.55,s=22,label='All mapped controls + nominated contact')
        t=common[common.item==key].iloc[0];ax.scatter(t.median_magnitude,t.pearson_r,color='#D62728',s=85,zorder=4,label=target.replace('_',' '))
        ax.annotate(f'{key}\nr={t.pearson_r:.3f}',(t.median_magnitude,t.pearson_r),xytext=(8,12),textcoords='offset points',fontsize=8)
        ax.set(xscale='log',xlabel='Median |radial work| across four systems (kcal/mol)',ylabel='Pearson r with electronic barrier ΔE‡',ylim=(-1.08,1.08),title=category.replace('_',' '));ax.legend(fontsize=7);ax.grid(alpha=.2)
    figs.suptitle('Strength and correlation together: every mapped contact qualifying in all four systems')
    figs.tight_layout();figs.savefig(figdir/'22_contact_correlation_specificity.png',dpi=190);plt.close(figs)

    atomic=r[(r.kind=='atom')&(r.metric=='atomic_work_kcal_mol')][['system','target','value','magnitude_rank','n_controls']]
    local=r[(r.kind=='contact')&(r.pool.isin(['same_category','distance_and_chemistry']))][['system','target','pool','metric','value','magnitude_rank','n_controls']]
    main=s[(s.metric=='radial_work_kcal_mol')&(s.barrier==ENERGIES[0])]
    original=s[(s.metric=='pair_atomic_work_kcal_mol')&(s.barrier==ENERGIES[0])]
    qmain=s[(s.metric=='mean_abs_projection_mEh_bohr')&(s.barrier==ENERGIES[0])]
    ex=examples[(examples.window==primary)&(examples.pool=='disjoint')&(examples.metric=='radial_work_kcal_mol')&(examples.barrier==ENERGIES[0])]
    perm=perms[(perms.window==primary)&(perms.pool=='same_category')&(perms.metric=='radial_work_kcal_mol')&(perms.barrier==ENERGIES[0])]
    sensitivity=spec[(spec.pool=='same_category')&(spec.metric=='radial_work_kcal_mol')][['window','target','barrier','target_pearson_r','n_common_controls','controls_at_least_as_correlated','stronger_and_at_least_as_correlated']]
    atomcorr=c[(c.kind=='atom')&(c.metric=='atomic_work_kcal_mol')&(c.barrier==ENERGIES[0])].copy()
    atomcorr=atomcorr.assign(absr=abs(atomcorr.pearson_r)).sort_values('absr',ascending=False)
    direct=s[(s.target=='opposite_contact')&(s.metric=='pair_atomic_work_kcal_mol')&(s.barrier==ENERGIES[0])].set_index('pool')
    radial=s[(s.target=='opposite_contact')&(s.metric=='radial_work_kcal_mol')&(s.barrier==ENERGIES[0])].set_index('pool')
    short=r[(r.kind=='contact')&(r.target=='opposite_contact')&(r.metric=='radial_work_kcal_mol')&(r.pool.isin(['same_category','distance_and_chemistry']))][['system','pool','value','magnitude_rank','n_controls']]
    text=f'''# 分子内順位と、他の原子・接触を対照にした特異性検証

追加DFTなし。2026-10-17。全57原子と全1596原子対を計算し、幾何・結合関係で対照を選んだ。障壁は各系のΔE‡・ΔG‡そのもの。MoriPFとの差やpercentage pointsは用いない。

## 結論

**反対側H/F···Oは、近い種類・距離の接触の中で応答が大きい。ただし、障壁に相関する接触がこれだけという結果ではない。アミン側は「大きさ」と「相関」を同時に満たす根拠が弱い。**

元の図15で示した反対側2原子の仕事和（電子障壁r≈−0.965）は、同じ分類で4系に共通する対照{int(direct.loc['same_category','n_common_controls'])}個中、{int(direct.loc['same_category','controls_at_least_as_correlated'])}個が同等以上の|r|を持つ。H/F···O/Nかつ平均距離差±0.75 Åに絞っても、対照{int(direct.loc['distance_and_chemistry','n_common_controls'])}個中{int(direct.loc['distance_and_chemistry','controls_at_least_as_correlated'])}個が上回る。前の高相関だけでは特異性を示せなかった。

接触方向の射影仕事（電子障壁r≈−0.967）では、同分類の対照{int(radial.loc['same_category','n_common_controls'])}個中{int(radial.loc['same_category','controls_at_least_as_correlated'])}個が同等以上の|r|を持つ。化学種・距離を絞っても対照{int(radial.loc['distance_and_chemistry','n_common_controls'])}個中{int(radial.loc['distance_and_chemistry','controls_at_least_as_correlated'])}個。したがってこの指標でも「ここだけ相関する」とはいえない。

一方、**RAFA15の反対側H56···O21は、化学種・距離を近づけた13候補の中で射影仕事・平均射影力とも1位**。同じ対照定義でMoriPFは2位、RAFA14は3位、RAFA16は2位であり、局所応答としては目立つ。2原子の仕事和では4系とも4/13位。指定接触を研究対象にする根拠は残るが、因果や唯一性は未確定。

{markdown(short)}

「強さと相関の両方で全4系の対象を上回る」対照は、反対側の2原子仕事和には存在しなかった。射影仕事には広い対照でC20···C24が1件あり（r≈−0.995）、H/F···O/Nに限定すると0件。これは指標・対照の選択により結論が変わることを示す。C20···C24は1–4のC···Cなので、同じ化学相互作用とは主張せず、1–4対を除いた感度比較も残した。

## 今回の問い

「指定した接触に相関がある」だけでなく、(1) 同じ系の他の原子・接触より応答が大きいか、(2) 他の接触には同様の障壁相関がないか、(3) 指定原子を共有しない対照でも同じ結果が出ないか、を確認する。

全原子や共有結合より最大であることは機構上の必要条件ではない。一方、同じ定義・近い距離・同じ分子内/分子間分類の対照でも平凡なら、「この接触だけが特別に強い」という説明は支持されない。本解析は全原子力からの間接量の比較であり、真のペア結合強度の順位ではない。

## 対照の固定ルール

- 反応C–N距離の3区間は従来と同じ。主区間2.50→2.10 Å、感度区間2.55→2.10、2.40→2.15 Å。
- 最も反応物側の構造で共有結合グラフを推定。距離が共有結合半径和×{config['covalent_radius_scale']}未満を結合とし、反応で形成中のC–Nは元の分子間の対として扱う。
- 直接結合(1–2)・結合角(1–3)の原子対は対照から除外。グラフ距離3以上（1–4対は含む）、区間内最短距離{config['contact_cutoff_A']:.1f} Å以下を非共有接触候補とする。これは操作的な接触定義で、相互作用の有無の判定ではない。
- アミン側指定接触は分子間、反対側は基質分子内。同じ分類内で比較する。全原子は別に57原子すべてを比較する。
- `same_category`：同じ分子内/分子間分類の全候補。`distance_matched`：さらに平均距離差±{config['matched_mean_distance_tolerance_A']:.2f} Å。`shared_anchor`：反対側はO21、アミン側は選択NHを共有。`disjoint`：指定接触の2原子をどちらも含まない対照。
- `chemically_related`：反対側はH/F···O/N、アミン側はH···H/F/O/Nに限定。`distance_and_chemistry`はさらに距離条件も課す。RAFA14の指定アミン側はH···Hなので、この対照からH···Hを除外しない。
- `exclude_1_4`：グラフ距離3の1–4対まで除外し、4結合以上離れた接触だけを比較する感度確認。
- 元素の種類・供与体角度・環内外は一律に同じではない。これは広い反証用対照で、精密な化学的マッチングではない。H/F置換のある同じ環位置を比較するため、元素一致を必須にして対照を空にすることは避けた。
- 同一分子内の強さ順位は対称Hも含む全原子対。系間の相関は番号だけでは対応させず、共通骨格、環上の置換位置、求核剤の対応位置を使用する。アルキル基上の複数Hは一意な化学的対応を仮定せず系間相関から除外した。原子対応と共有結合グラフをCSVで公開する。
- 系間比較では同じ対が**全4系で**接触条件を満たすものだけを採用。距離マッチも全4系で満たす。障壁値から対照を選ばない。

## 全原子内で指定原子はどこに位置するか

仕事の符号を保って保存し、大きさの順位は絶対値で比較。rank=1が最大、対照数56に対象1原子を加えて全57原子。

{markdown(atomic)}

![全原子の順位](../figures/20_all_atom_work_controls.png)

## 同じ系の他の接触に比べて強いか

`radial_work_kcal_mol`は接触方向の仕事診断、`mean_abs_projection_mEh_bohr`は接触方向の力の平均絶対値。強さのrankは大きい順、比較総数はn_controls+1。少数の順位差を精密な結合強度差とは読まない。

{markdown(local)}

![接触の強さ分布](../figures/21_contact_strength_controls.png)

## この接触だけが障壁に相関するか

主区間、電子障壁についての結果。相関は符号を保つ。対照との相関順位は|r|で比較し、正/負いずれも強い対応として数える。`controls_at_least_as_correlated`は対象以上の|r|を持つ対照数。`stronger_in_all_four`は全4系で対象以上の大きさを持つ対照数。最後の列は両条件を同時に満たす対照数。

{markdown(main)}

元の図15の「2原子の仕事和」でも同じ対照比較を行った。これは指定した2原子の仕事の和であり、2原子間の力ではない：

{markdown(original)}

アミン側の平均射影力は電子障壁との|r|順位が同分類で1位だが、各分子内での大きさは下位である。相関順位だけで強い接触と呼べない例になる。

同じ比較を力の平均絶対値で行った結果：

{markdown(qmain)}

![強さと相関の同時比較](../figures/22_contact_correlation_specificity.png)

### 指定した2原子をどちらも含まない対照の具体例

以下は相関の強い順に各5候補に加え、「全4系でより大きく、相関も同等以上」の候補を最大3つ表示した反例候補で、これらを新たな機構として採用したわけではない。全候補はCSVに残す。`coreNN`は共通骨格の原子番号NN、`X24`〜`X28`はその環C上のF/H、`NH_near`/`NH_other`はTSでアミン側Xに近い/遠いアミンH。各系の具体的な原子番号・元素はmappingを参照。

{markdown(ex)}

### 他の単一原子も相関するか

電子障壁と原子仕事の|r|が高い10原子。原子番号の対応が一意でないアルキルHは含めない。

{markdown(atomcorr.head(10).drop(columns='absr'))}

## 区間感度と探索の影響

{markdown(sensitivity)}

4系の障壁ラベルの全24通りの並べ替えについて、候補全体の最大|r|を求めた。対象の|r|以上の最大値が何通りで得られるかを記録する。下表は同分類の全候補を同時に走査した診断である。

{markdown(perm)}

これは候補を多数見たときに高い相関がどれほど出やすいかの記述で、事前登録された確証検定や独立したp値ではない。対照同士は原子を共有して互いに依存する。弱い相関の対照だけを選んで「ここだけ相関する」とは結論しない。全候補、全3区間、ΔE‡とΔG‡の双方を保存した。

## 次に既存logだけで詰めること

現在の結果からは、相関の良い指標を追加探索するより、RAFA15で大きい反対側応答がO21単独の移動・環の回転・C–N形成に伴う共通変形で説明できるかを調べる方が有用。累積仕事と二面角・接触角の変化を同じC–N軸で比較し、指定ペアとH11···O18などの高相関対照を同時に示す。複数の原子力射影は同じ構造運動を共有するので、相関が多数出ること自体も結果として扱う。

原子仕事の座標合わせや射影仕事の微小値への数値誤差の影響は未評価の部分が残る。今回の対照比較は3区間で行ったが、追加の物理計算による独立検証ではない。

## 再現と参照

`python analysis4_force/contact_specificity.py`、または通常の全解析を実行する。出力は `contact_controls_*.csv`、図20〜22。既存の指定接触・原子仕事と今回の全候補走査の値が一致することをテストする。

原子力と移動からの仕事配分には [Haas & Chu (2009), DOI:10.1063/1.3243080](https://doi.org/10.1063/1.3243080) の背景がある。今回の接触軸への射影、対照抽出、特異性評価は本プロジェクトのad hocな設計であり、論文によりこの仮説が証明されたという意味ではない。
'''
    (reportdir/'contact_specificity_summary_ja.md').write_text(text)
