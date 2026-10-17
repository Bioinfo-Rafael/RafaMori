"""Transparent exploratory report: show primary and secondary descriptors alike."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting4 import SYSTEM_COLORS
from reporting import markdown


def create_report(profiles,metrics,corr,mapping,cfg,fig,reports):
    chosen=metrics[(metrics.window==cfg['primary_window'])&(metrics.H_choice=='primary')].copy()
    primary=corr[(corr.window==cfg['primary_window'])&(corr.H_choice=='primary')]
    systems=chosen.system.tolist();ref=cfg['reference_system']
    fig1,axes=plt.subplots(2,2,figsize=(12,9),sharex=True)
    for system in systems:
        p=profiles[(profiles.system==system)&(profiles.H_choice=='primary')&(profiles.R<0)]
        p=p[(p.CN_A<=2.55)&(p.CN_A>=2.10)].sort_values('CN_A',ascending=False)
        color=SYSTEM_COLORS[system]
        axes[0,0].plot(p.CN_A,p.near_distance_A,color=color,label=system)
        axes[0,1].plot(p.CN_A,p.opposite_distance_A,color=color,label=system)
        axes[1,0].plot(p.CN_A,p.near_projected_force_hartree_bohr*1000,color=color,label=system)
        axes[1,1].plot(p.CN_A,p.opposite_projected_force_hartree_bohr*1000,color=color,label=system)
    for ax,title in zip(axes.flat,['Amine-facing ortho X...H (X=F or H)','Opposite ortho X...O21 (X=F or H)','Amine-side total-force projection','Opposite-side total-force projection']):
        ax.set_title(title);ax.grid(alpha=.2);ax.legend(fontsize=8)
    axes[0,0].invert_xaxis()
    for ax in axes[0]:ax.set_ylabel('Distance (Å)')
    for ax in axes[1]:ax.set(xlabel='Reaction C-N distance (Å); reaction progresses right',ylabel='Q = 0.5 (Fb-Fa)·eab (mEh/Bohr)');ax.axhline(0,color='black',lw=.7)
    fig1.suptitle('Fixed ortho sites: net atomic-force projection is not a pair-interaction force')
    fig1.tight_layout();fig1.savefig(fig/'14_two_sided_contacts.png',dpi=190);plt.close(fig1)
    fig2,axes=plt.subplots(1,2,figsize=(11,5.4))
    for ax,metric,title in zip(axes,['opposite_projected_magnitude_share','opposite_atomic_force_share'],['Primary: projected-force magnitude share','Secondary: atomic-force magnitude share']):
        rowcorr=primary[(primary.descriptor==metric)&(primary.target=='relative_electronic_barrier_kcal_mol')].iloc[0]
        for row in chosen.itertuples():
            x=getattr(row,metric+'_minus_MoriPF')*100;y=row.relative_electronic_barrier_kcal_mol
            ax.scatter(x,y,color=SYSTEM_COLORS[row.system],s=65);ax.annotate(row.system,(x,y),textcoords='offset points',xytext=(5,6),fontsize=9)
        ax.axhline(0,color='grey',lw=.7);ax.axvline(0,color='grey',lw=.7)
        ax.set(title=f'{title}\nr={rowcorr.pearson_r:.3f}, rank correlation={rowcorr.spearman_r:.1f}',xlabel=f'Opposite-side share minus {ref} (percentage points)',ylabel=f'Electronic barrier minus {ref} (kcal/mol)')
        ax.margins(.35);ax.grid(alpha=.2)
    fig2.suptitle('4 systems; C-N 2.50 → 2.10 Å; descriptive association, not independent validation',fontsize=11)
    fig2.tight_layout();fig2.savefig(fig/'15_two_sided_barrier_correlations.png',dpi=190);plt.close(fig2)
    fig3,axes=plt.subplots(1,2,figsize=(12,4.8))
    for ax,target,title in zip(axes,['relative_electronic_barrier_kcal_mol','relative_RRHO_G_barrier_kcal_mol'],['Electronic barrier','RRHO Gibbs barrier (RAFA15 provisional)']):
        subset=corr[(corr.descriptor=='opposite_atomic_force_share')&(corr.target==target)]
        labels=[w.replace('common_CN_','').replace('sensitivity_CN_','').replace('_to_',' → ')+'\n'+('nearest H' if h=='primary' else 'other H') for w,h in zip(subset.window,subset.H_choice)]
        ax.bar(np.arange(len(subset)),subset.pearson_r,color='#4C78A8');ax.set_xticks(np.arange(len(subset)),labels,rotation=35,ha='right',fontsize=7);ax.set(title=title,ylabel='Pearson r',ylim=(-1,1));ax.axhline(0,color='black',lw=.7)
    fig3.suptitle('Sensitivity of the secondary atomic-force-share descriptor (all tested windows shown)')
    fig3.tight_layout();fig3.savefig(fig/'16_two_sided_sensitivity.png',dpi=190);plt.close(fig3)
    disp=chosen[['system','opposite_atomic_force_share','opposite_atomic_force_share_minus_MoriPF','relative_electronic_barrier_kcal_mol','relative_RRHO_G_barrier_kcal_mol']].copy()
    disp.iloc[:,1:3]*=100
    disp.columns=['系','反対側原子力の割合 (%)','MoriPFとの差 (ポイント)','ΔΔE‡ (kcal/mol)','ΔΔG‡ (kcal/mol)']
    geom=chosen[['system','near_TS_distance_A','opposite_TS_distance_A','opposite_TS_CXO_angle_deg']]
    allcor=primary[['descriptor','target','pearson_r','spearman_r','LOO_prediction_RMSE_kcal_mol']]
    sens=corr[(corr.descriptor=='opposite_atomic_force_share')][['window','H_choice','target','pearson_r','spearman_r','LOO_prediction_RMSE_kcal_mol']]
    tableforce=chosen[['system','near_X_mean_force_mEh_bohr','amine_H_mean_force_mEh_bohr','opposite_X_mean_force_mEh_bohr','opposite_O_mean_force_mEh_bohr','near_mean_signed_projected_force_mEh_bohr','opposite_mean_signed_projected_force_mEh_bohr']]
    er=primary[(primary.descriptor=='opposite_atomic_force_share')&(primary.target=='relative_electronic_barrier_kcal_mol')].iloc[0]
    gr=primary[(primary.descriptor=='opposite_atomic_force_share')&(primary.target=='relative_RRHO_G_barrier_kcal_mol')].iloc[0]
    e_sensitive=corr[(corr.descriptor=='opposite_atomic_force_share')&(corr.target=='relative_electronic_barrier_kcal_mol')].pearson_r
    text=f'''# 両ortho位置の接触・力と活性化障壁の探索解析

## 結論

**障壁の順位と対応する「力の配分」の候補は得られた。ただし、接触方向の力の比率は弱い対応にとどまり、F···O反発／H···O引力を定量した結果ではない。**

ユーザーの仮説は、アミン側のF···Hに加え、反対側のortho H···クマリンOが有利に働く一方、同位置がFならOとの反発が不利に働く、という両側のバランス。前の解析は「最も近いF」を選んでいたため、RAFA15のH56···O21を含まなかった。今回は同じ環位置を対応させて追加解析した。

C–N結合長2.50→2.10 Åの共通区間で集計すると、**反対側のH/FとOに働く全原子力の大きさの割合**はRAFA15＜MoriPF＜RAFA16＜RAFA14で、電子障壁・Gibbs障壁の順位と一致する。

{markdown(disp)}

この二次的な候補のPearson相関は、電子障壁で **r={er.pearson_r:.3f}**、RRHO Gibbs障壁で **r={gr.pearson_r:.3f}**、順位相関はいずれも1.0。一方、最初に指定した主指標「接触方向に射影した力の大きさの比率」は電子障壁でr≈−0.393、順位相関−0.4で、十分な対応を示さなかった。良い指標だけを残さず、12種類の指標を全てCSVと下表に掲載した。

![両側比率と障壁](../figures/15_two_sided_barrier_correlations.png)

## 原子の対応と接触形状

{markdown(mapping)}

C24がアミン側、C25が反対側。RAFA15の反対側はH56–O21、MoriPF・RAFA16はF29–O21。RAFA14は両側のortho位置がHで、アミン側はH57···H44という対照になる。RAFA14のアミン側をF···Hと呼ぶことはできない。主解析のアミンHはTSでC24上の置換原子に近いHを固定し、もう一方のHでも感度を調べた。

{markdown(geom)}

RAFA15のH56···O21はTSで約2.764 Å、C25–H56···O21角は約90.4°。距離が近いことだけで水素結合・安定化を確定しない。ユーザーの「相互作用が有利」という仮説と、今回測れた幾何・力の指標は区別する。

![接触方向の力](../figures/14_two_sided_contacts.png)

## 比率の定義

`X_near` はC24上のF/H、`H_N` は選択したアミンH、`X_opp` はC25上のF/H、OはO21。Gaussianの原子力ベクトルを **F_i**、反応C–N距離をcとする。

原子力の大きさを共通のc区間で積分し、

```text
B_atomic = ∫ (|F_Xopp| + |F_O|) dc
           / ∫ (|F_Xnear| + |F_HN| + |F_Xopp| + |F_O|) dc
```

を計算した。表の割合はこの値の100倍。分子全体に対する割合ではなく、重複しない選択4原子内の配分である。分子全体に対する選択4原子の割合は `selected_atomic_force_share` として別途保存し、そちらの電子障壁との相関は約−0.632にとどまる。

各原子の平均力と射影力（mEh/Bohr）も保存している：

{markdown(tableforce)}

RAFA15の反対側2原子の力の大きさは、H56が約0.050、O21が約1.477 mEh/Bohrで、合計の約97%がO21に由来する。このため、指標が大きいことをH···Oの結合力が強いことに読み替えるのは不適切で、O21を含む反応中心の応答を捉えている可能性がある。

**このB_atomicは全原子力のノルムの比であり、H···FやF/H···Oのペア相互作用力の比ではない。** 原子力には共有結合、他の非共有結合、電子構造全体の効果が含まれる。FをHへ置き換えた原子種の違いやIRCの質量重み付けも影響し得る。新しい指標の名称は「反対側原子力割合」とし、「反発率」や「水素結合強度」とはしない。

## 射影力と力・運動の積

```text
Q_ab = 0.5 (F_b − F_a) · e_ab
A_ab = Q_ab × d(r_ab in Bohr)/dR
W_ab_diagnostic = −∫ A_ab dR
```

Q>0は距離を開く側、Q<0は距離を閉じる側への全原子力の射影。ただし、これをa–b間の独立した斥力／引力と同一視しない。実際、今回の反対側射影の平均はF···Oを含めて全系で負であり、この符号からF···O反発の仮説を確認することも否定することもできない。

`opposite_projected_magnitude_share = ⟨|Q_opp|⟩ / (⟨|Q_near|⟩ + ⟨|Q_opp|⟩)` を主指標として調べたが、原子力ノルムの比率ほど障壁に対応しなかった。TSでは全体の勾配が小さく、今回のTS力も未収録なので、TS一点のforce比は計算していない。

原子別 `C_i=F_i·dr_i/dR` を選択4原子で足し、その負積分も保存した。これは共通IRC区間の力・運動診断であり、分離反応物→TSの全障壁の相互作用エネルギー分解ではない。分離反応物からIRC左端までの会合・配座変化と、Gibbs障壁のエントロピー効果はこの積分に含まれない。

## 相関・感度・予測の限界

共通区間はC–N=2.50→2.10 Å。TSの力欠測に触れず全系に実データがある範囲で、相関を見る前に指定した。2.55→2.10 Å、2.40→2.15 Åと、アミンのもう一方のHを使う解析も全て実行した。MoriPFを引いた指標と障壁で比較し、回帰には切片を含める（MoriPFを通るよう強制しない）。

{markdown(allcor)}

高相関の候補 `B_atomic` について、電子障壁とのrは全感度条件で **{e_sensitive.min():.3f}〜{e_sensitive.max():.3f}**。順位は維持されるが、区間によって相関の強さは変わる。

{markdown(sens)}

![感度](../figures/16_two_sided_sensitivity.png)

主区間のB_atomicを使った1系除外の線形予測では、電子障壁のRMSEは{er.LOO_prediction_RMSE_kcal_mol:.2f} kcal/mol、Gibbs障壁では{gr.LOO_prediction_RMSE_kcal_mol:.2f} kcal/mol。4系しかなく、指標探索に同じ4系を使っているため、この交差検証も独立検証ではない。複数指標の中から高い相関が見つかった探索結果であり、「十分に予測できる」「機構を証明した」とはまだ判断しない。p値による有意性は主張しない。

## 次の検証

この候補が実用的かは、追加置換系や独立した配座で、今回固定した定義と区間のまま予測できるかを確認するのが重要。機構については、共通C–N距離での活性化ひずみ／相互作用エネルギー分解を追加すれば、構造変形と分子間相互作用の寄与を切り分けられる。ただし通常の分子間分解だけでは分子内の反対側H/F···Oを一意に分離できず、適切な配座対照や電子密度解析も必要になる。これらの追加量子化学計算は今回行っていない。

## 再現・参照

`python analysis4_force/two_sided_contacts.py` で再生成。設定は `config/two_sided_contacts.yaml`、全結果は `two_sided_*.csv`。通常の4系解析を実行した後の座標・力を使用する。

- [原子寄与からのReaction Force解析（2012）](https://pubmed.ncbi.nlm.nih.gov/22746292/)：原子の反応進行への寄与を考える背景。今回のペア射影比・原子力割合はこのプロジェクト独自の探索指標。
- [Distortion/Interaction–Activation Strainモデル（2017）](https://pubmed.ncbi.nlm.nih.gov/28447369/)：障壁の構造変形と相互作用への分解を考える背景。
'''
    (reports/'two_sided_contact_summary_ja.md').write_text(text)
