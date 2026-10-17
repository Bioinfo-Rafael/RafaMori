"""Log-only work diagnostics with dimensional axes and uncentered barriers."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting4 import SYSTEM_COLORS
from reporting import markdown

TARGETS=['absolute_electronic_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol']

def correlation_figure(chosen,primary,descriptors,labels,titles,path,subtitle,targets=TARGETS):
    fig,axes=plt.subplots(len(targets),len(descriptors),figsize=(6.1*len(descriptors),4.6*len(targets)),squeeze=False)
    for i,target in enumerate(targets):
        for ax,metric,label,title in zip(axes[i],descriptors,labels,titles):
            rc=primary[(primary.descriptor==metric)&(primary.target==target)].iloc[0]
            for row in chosen.itertuples():
                x=getattr(row,metric);y=getattr(row,target)
                ax.scatter(x,y,color=SYSTEM_COLORS[row.system],s=65,zorder=3)
                ax.annotate(row.system,(x,y),textcoords='offset points',xytext=(6,7),fontsize=9)
            ax.set(title=f'{title}\nr={rc.pearson_r:.3f}; Spearman={rc.spearman_r:.2f}',xlabel=label,
                   ylabel=('Electronic barrier ΔE‡ (kcal/mol)' if i==0 else 'RRHO Gibbs barrier ΔG‡ (kcal/mol)'))
            ax.margins(x=.30,y=.27);ax.grid(alpha=.2)
    fig.suptitle(subtitle,fontsize=11)
    fig.tight_layout(rect=(0,0,1,.96));fig.savefig(path,dpi=190);plt.close(fig)


def create_report(profiles,metrics,corr,mapping,cfg,fig,reports):
    chosen=metrics[(metrics.window==cfg['primary_window'])&(metrics.H_choice=='primary')].copy()
    primary=corr[(corr.window==cfg['primary_window'])&(corr.H_choice=='primary')]
    systems=chosen.system.tolist()
    fig1,axes=plt.subplots(2,2,figsize=(12,9),sharex=True)
    for system in systems:
        p=profiles[(profiles.system==system)&(profiles.H_choice=='primary')&(profiles.R<0)]
        p=p[(p.CN_A<=2.55)&(p.CN_A>=2.10)].sort_values('CN_A',ascending=False)
        color=SYSTEM_COLORS[system]
        for ax,col,scale in zip(axes.flat,['near_distance_A','opposite_distance_A','near_projected_force_hartree_bohr','opposite_projected_force_hartree_bohr'],[1,1,1000,1000]):
            ax.plot(p.CN_A,p[col]*scale,color=color,label=system)
    for ax,title in zip(axes.flat,['Amine-facing ortho X...H (X=F or H)','Opposite ortho X...O21 (X=F or H)','Amine-side total-force projection','Opposite-side total-force projection']):
        ax.set_title(title);ax.grid(alpha=.2);ax.legend(fontsize=8)
    axes[0,0].invert_xaxis()
    for ax in axes[0]:ax.set_ylabel('Distance (Å)')
    for ax in axes[1]:ax.set(xlabel='Reaction C-N distance (Å); reaction progresses right',ylabel='Q = 0.5 (Fb-Fa)·eab (mEh/Bohr)');ax.axhline(0,color='black',lw=.7)
    fig1.suptitle('Fixed ortho sites: net atomic-force projection is not a pair-interaction force')
    fig1.tight_layout();fig1.savefig(fig/'14_two_sided_contacts.png',dpi=190);plt.close(fig1)

    correlation_figure(chosen,primary,['near_atomic_work_kcal_mol','opposite_atomic_work_kcal_mol'],
        ['Amine-side atomic work (kcal/mol)','Opposite-side atomic work (kcal/mol)'],
        ['Ortho X + amine H','Opposite ortho X + coumarin O21'],fig/'15_two_sided_barrier_correlations.png',
        'Log-only atomic work; C-N 2.50 → 2.10 Å; 4 systems\nBarriers from each separated substrate + Lys; RAFA15 RRHO provisional')
    correlation_figure(chosen,primary,['opposite_projected_magnitude_share','opposite_atomic_force_share'],
        ['Opposite / (amine + opposite) projected-force magnitude (dimensionless)',
         'Opposite / (amine + opposite) atomic-force magnitude (dimensionless)'],
        ['Original primary ratio (supplement)','Original secondary ratio (supplement)'],fig/'17_dimensionless_share_supplement.png',
        'Original exploratory ratios retained without centering or percentages; not interaction strengths',targets=TARGETS[:1])
    correlation_figure(chosen,primary,['near_pair_diagnostic_work_kcal_mol','opposite_pair_diagnostic_work_kcal_mol'],
        ['Amine-side radial work diagnostic (kcal/mol)','Opposite-side radial work diagnostic (kcal/mol)'],
        ['Contact-direction diagnostic','Contact-direction diagnostic'],fig/'18_radial_work_barrier_correlations.png',
        'Q from total atomic forces; radial work is not isolated pair-interaction energy; 4 systems')

    fig3,axes=plt.subplots(1,2,figsize=(12,4.8))
    for ax,target,title in zip(axes,TARGETS,['Electronic barrier ΔE‡','RRHO Gibbs barrier ΔG‡ (RAFA15 provisional)']):
        for descriptor,label,marker in [('near_atomic_work_kcal_mol','Amine-side work','o'),('opposite_atomic_work_kcal_mol','Opposite-side work','s'),('opposite_atomic_force_share','Original atomic-force ratio','^')]:
            subset=corr[(corr.descriptor==descriptor)&(corr.target==target)]
            labels=[w.replace('common_CN_','').replace('sensitivity_CN_','').replace('_to_',' → ')+'\n'+('nearest H' if h=='primary' else 'other H') for w,h in zip(subset.window,subset.H_choice)]
            ax.plot(np.arange(len(subset)),subset.pearson_r,marker=marker,label=label)
        ax.set_xticks(np.arange(len(subset)),labels,rotation=35,ha='right',fontsize=7)
        ax.set(title=title,ylabel='Pearson r',ylim=(-1.05,1.05));ax.axhline(0,color='black',lw=.7);ax.legend(fontsize=8)
    fig3.suptitle('All predefined windows and both amine H choices; exploratory sensitivity')
    fig3.tight_layout();fig3.savefig(fig/'16_two_sided_sensitivity.png',dpi=190);plt.close(fig3)

    fig4,axes=plt.subplots(1,2,figsize=(12,5))
    x=np.arange(len(chosen));width=.19
    for j,(key,label) in enumerate([('near_X','Amine-facing X'),('amine_H','Amine H'),('opposite_X','Opposite X'),('opposite_O','Coumarin O21')]):
        axes[0].bar(x+(j-1.5)*width,chosen[key+'_atomic_work_kcal_mol'],width,label=label)
    axes[0].set_xticks(x,systems);axes[0].set(ylabel='Atomic work (kcal/mol)',title='Selected atoms, shown separately');axes[0].legend(fontsize=8)
    positive=np.zeros(len(chosen));negative=np.zeros(len(chosen))
    for key,label in [('near','Amine-side atoms'),('opposite','Opposite-side atoms'),('remaining','All remaining atoms')]:
        v=chosen[key+'_atomic_work_kcal_mol'].to_numpy()
        bottom=np.where(v>=0,positive,negative)
        axes[1].bar(x,v,bottom=bottom,label=label)
        positive+=np.maximum(v,0);negative+=np.minimum(v,0)
    axes[1].scatter(x,chosen.all_atomic_work_kcal_mol,marker='_',s=180,color='black',label='Sum of all atoms')
    axes[1].set_xticks(x,systems);axes[1].set(ylabel='Work in common interval (kcal/mol)',title='Accounting for the complete interval');axes[1].legend(fontsize=8)
    fig4.suptitle('C-N 2.50 → 2.10 Å; atomic work allocation, not pair binding energy')
    fig4.tight_layout();fig4.savefig(fig/'19_atomic_work_partition.png',dpi=190);plt.close(fig4)

    disp=chosen[['system','near_atomic_work_kcal_mol','opposite_atomic_work_kcal_mol','absolute_electronic_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol']].copy()
    disp.columns=['系','アミン側原子仕事','反対側原子仕事','ΔE‡','ΔG‡（RRHO）']
    def selected_corr(names):
        return primary[primary.descriptor.isin(names)][['descriptor','target','pearson_r','spearman_r','LOO_prediction_RMSE_kcal_mol']]
    works=['near_atomic_work_kcal_mol','opposite_atomic_work_kcal_mol','near_pair_diagnostic_work_kcal_mol','opposite_pair_diagnostic_work_kcal_mol']
    raw=chosen[['system']+[n+'_atomic_work_kcal_mol' for n in ['near_X','amine_H','opposite_X','opposite_O','remaining']]+['all_atomic_work_kcal_mol']]
    accounting=chosen[['system','entry_energy_kcal_mol','window_energy_rise_kcal_mol','after_window_to_TS_kcal_mol','absolute_electronic_barrier_kcal_mol','work_energy_residual_kcal_mol']]
    geom=chosen[['system','near_TS_distance_A','opposite_TS_distance_A','opposite_TS_CXO_angle_deg']]
    tableforce=chosen[['system','near_X_mean_force_mEh_bohr','amine_H_mean_force_mEh_bohr','opposite_X_mean_force_mEh_bohr','opposite_O_mean_force_mEh_bohr','near_mean_signed_projected_force_mEh_bohr','opposite_mean_signed_projected_force_mEh_bohr']]
    allcor=primary[['descriptor','target','pearson_r','spearman_r','LOO_prediction_RMSE_kcal_mol']]
    sens=corr[corr.descriptor.isin(works+['opposite_atomic_force_share'])][['window','H_choice','descriptor','target','pearson_r','spearman_r']]
    text=f'''# 既存logによる両側接触・原子仕事と活性化障壁

## 結論と今回の変更

追加のDFT計算は行わず、保存済みの座標・Gaussian原子力・エネルギーだけを使用した。図15の横軸は **各側の原子仕事（kcal/mol）**、縦軸は **各系の分離反応物を基準とするΔE‡・ΔG‡** に変更した。MoriPFとの差やpercentage pointsは使わない。元の比率の相関は補足図17に残し、相互作用の強さとは解釈しない。

文献調査に基づく選択理由・今後のlog限定の計算案は [調査ノート](log_only_research_plan_ja.md) を参照。今回の原子仕事は「どの原子の移動にエネルギー上昇が配分されるか」という間接指標。F···HやH/F···Oの結合エネルギーではないが、単位・積分範囲・全体との加算関係を明示できる。

共通区間はC–N 2.50→2.10 Å。表の全数値はkcal/mol：

{markdown(disp)}

![原子仕事と障壁](../figures/15_two_sided_barrier_correlations.png)

{markdown(selected_corr(works))}

反対側原子仕事は電子障壁とr≈−0.965、Gibbs障壁とr≈−0.976。RAFA15の仕事が大きいほど障壁が低いという**逆向きの対応**であり、「反発コストが低いので速い」とは読めない。RAFA15では反対側仕事の約98%がO21由来で、H56単独の仕事は小さい。

相関の良さは指標を物理的な相互作用量に変えない。両側を別々に示し、相関が弱い側も残す。選択4原子の合計仕事、O21単独、残りの全原子も保存し、特定の接触ではなく反応中心全体の変化で説明できないかを確認する。

## 原子の対応

{markdown(mapping)}

C24上のX=F/Hと選択アミンHをアミン側、C25上のX=F/HとO21を反対側とする。RAFA15の反対側はベンゼン環H56···O21。RAFA14のアミン側はH57···H44でありF···Hではない。アミンHはTSでC24置換原子に近い方を固定し、もう一方も感度解析する。

{markdown(geom)}

RAFA15のC25–H56···O21角は約90.4°で、直線的なC–H···O配置ではない。近接・相互作用の仮説として扱い、水素結合と断定しない。

## 仕事の定義と解釈

```text
C_i(R) = F_i(R) · dr_i/dR
W_g = −∫[R(c=2.50) → R(c=2.10)] Σ(i∈g) C_i(R) dR
W_near + W_opposite + W_remaining = W_all ≈ E(end) − E(start)
```

座標はBohr、力はHartree/Bohrとして積分し、kcal/molに変換した。全57原子をTSへ剛体重ね合わせし、力ベクトルも同じ回転を施した。原子別配分はこの座標系・経路に依存し、独立した原子の状態関数ではない。全系の和と電子エネルギー差の一致を数値検証する。既存区間の仕事であり、分離反応物→TSの全障壁を分解したものではない。ΔG‡との比較も電子仕事の相関であって自由エネルギー分解ではない。

{markdown(raw)}

![原子別・全体の仕事](../figures/19_atomic_work_partition.png)

反応物側のIRCはエネルギー上昇方向なので、W>0を「その接触が反発している」と解釈しない。他の原子・結合・電子の応答が各原子力に含まれる。O21が主成分なら「反対側Oを含む構造変化のコスト」と表現する。全体との和が合うことは数値検証であり、H···Oの因果証明ではない。

## 局所仕事と全障壁を混同しないためのエネルギー収支

```text
ΔE‡ = [E(c=2.50) − E_separated]
      + [E(c=2.10) − E(c=2.50)] + [E_TS − E(c=2.10)]
```

{markdown(accounting)}

RAFA15は共通区間のエネルギー上昇が約4.031 kcal/molと大きい一方、その区間に入る時点のエネルギーが約−9.019 kcal/molと低い。したがって共通区間の仕事が大きくても全障壁は低くなれる。入口の値は会合・配座変化などを含む総エネルギーであり、特定のH···O安定化エネルギーではない。入口値と障壁の相関は上式を共有するため、独立した予測証拠に使わない。

## 接触方向に限定した補助指標

```text
Q_ab = 0.5 (F_b − F_a) · e_ab
W_radial,ab = −∫ Q_ab d r_ab
```

力を接触軸へ射影し、その距離変化に沿って積分する。単位はkcal/mol。共通C–N区間内での動きを評価するad hocな診断値で、単独ペアのポテンシャルではなく、全原子仕事に足してはいけない。Q>0は開く向き、Q<0は閉じる向き。全系の反対側平均Qは負であり、F···O固有の反発をこの符号から判定しない。

反対側の射影仕事は電子障壁とr≈−0.967（3区間で約−0.959〜−0.969）、アミン側は主区間でr≈−0.185。反対側の量は約−0.0016〜0.0091 kcal/molにすぎず、数kcal/molの全障壁差をこの値の大きさで説明しているわけではない。微小な構造・力の応答を拾う候補であり、原子仕事とは別定義である。反対側の定義はアミンHの選択に依存しないため、Hを替えて一致すること自体は独立の感度検証ではない。

![接触方向の仕事](../figures/18_radial_work_barrier_correlations.png)

![接触距離と射影力](../figures/14_two_sided_contacts.png)

平均力（mEh/Bohr）：

{markdown(tableforce)}

## 元の比率は補足に保持

```text
B_atomic = ∫ (|F_Xopp| + |F_O|) dc
           / ∫ (|F_Xnear| + |F_HN| + |F_Xopp| + |F_O|) dc
```

割合そのものを無次元値として表示する。電子障壁とのr≈0.958という相関はMoriPFを引かなくても変わらない。定数の差し引きは相関を変えないため、以前の相関がゼロ合わせで生まれたわけではない。ただし、分母・原子選択に依存し、物理的な解釈に乏しいという問題は残る。

RAFA15の反対側原子力の大きさは約97%がO21に由来する。この比率をH···O相互作用の強さと読むことはできない。

![元の無次元比率](../figures/17_dimensionless_share_supplement.png)

## 全指標と感度

元の12指標に、今回7つの原子仕事成分と両側の射影仕事2成分を追加した。新指標はユーザーのlog限定・有次元指標の要望を受けた探索的追加であり、事前登録や独立検証ではない。全{len(cfg['correlation_descriptors'])}指標×3区間×2通りのアミンH×2種類の障壁を保存する。区間を変えて最も高い相関だけを選ばない。

{markdown(allcor)}

![感度](../figures/16_two_sided_sensitivity.png)

{markdown(sens)}

4系しかないため、相関係数・順位と1系除外予測は記述的診断。複数指標を同じ4系で探索した交差検証は独立検証ではない。比率や線形結合の係数を障壁に合わせて最適化しない。RAFA15のRRHO Gibbs障壁は反応物の微小虚振動のため暫定。

## 再現・文献

`python analysis4_force/run_analysis.py` で全結果、`python analysis4_force/two_sided_contacts.py` でこの追加解析を再生成できる。TSの力欠測を補間せず、全系で実測力がある共通区間のみ使用した。

- Haas & Chu (2009), [Decomposition of energy and free energy changes by following the flow of work along reaction path](https://doi.org/10.1063/1.3243080)：原子力と移動による仕事配分の根拠。本解析で自由エネルギー経路を計算したという意味ではない。
- Vöhringer-Martinez & Toro-Labbé (2012), [Understanding the physics and chemistry of reaction mechanisms from atomic contributions: a reaction force perspective](https://doi.org/10.1021/jp303075k)：原子別の反応力・エネルギー寄与の考え方。
'''
    (reports/'two_sided_contact_summary_ja.md').write_text(text)
