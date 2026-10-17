"""Figures and Japanese result report from computed CSVs."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting4 import SYSTEM_COLORS


def markdown(df, columns=None):
    if columns:df=df[columns]
    def fmt(v):
        if isinstance(v,(float,np.floating)):return f'{v:.3f}'
        return str(v).replace('|','/')
    return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))


def create_extra_figures(summary,energy,contacts,corrected,figdir,config):
    (figdir/'02_relative_barriers.png').unlink(missing_ok=True)
    fig,axes=plt.subplots(1,3,figsize=(14,4.7))
    for ax,col,title in zip(axes,['absolute_electronic_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol','alpha_to_TS_kcal_mol'],['Electronic barrier ΔE‡','RRHO Gibbs barrier ΔG‡','Work from alpha to TS']):
        bars=ax.bar(summary.system,summary[col],color=[SYSTEM_COLORS[s] for s in summary.system]);ax.bar_label(bars,fmt='%.2f',padding=3);ax.set(title=title,ylabel='kcal/mol');ax.axhline(0,color='black',lw=.7);ax.margins(y=.22)
    fig.suptitle('Barriers from each separated substrate + Lys; no inter-system centering\nRAFA15 RRHO values are provisional (small imaginary reactant mode)',fontsize=11)
    fig.tight_layout();fig.savefig(figdir/'02_barriers_and_work.png',dpi=190);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.8))
    for ax,col,title in zip(axes,['absolute_electronic_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol'],['Electronic barrier from separated reactants','RRHO Gibbs barrier from separated reactants']):
        bars=ax.bar(summary.system,summary[col],color=[SYSTEM_COLORS[s] for s in summary.system])
        ax.bar_label(bars,fmt='%.2f',padding=3);ax.axhline(0,color='black',lw=.7)
        ax.set(title=title,ylabel='kcal/mol');ax.margins(y=.25)
    fig.suptitle('M06-2X/6-311++G(d,p), SMD water; RRHO: 298.15 K, 1 atm\nRAFA15 RRHO values are provisional (small imaginary reactant mode)',fontsize=11)
    fig.tight_layout();fig.savefig(figdir/'03_absolute_barriers.png',dpi=190);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,9))
    for ax,system in zip(axes.flat,summary.system):
        c=contacts[(contacts.system==system)&contacts.selected_at_TS].sort_values('R')
        ax.plot(c.R,c.FH_A,color=SYSTEM_COLORS[system],label='Fixed TS-selected F...H distance')
        for (f,h),g in contacts[contacts.system==system].groupby(['F_atom','H_atom']):
            ax.plot(g.R,g.FH_A,color='grey',alpha=.20,lw=.8)
        ax2=ax.twinx();ax2.plot(c.R,c.NHF_angle_deg,color='black',ls='--',label='N-H...F angle');ax2.set_ylabel('N-H...F angle (degree)');ax2.set_ylim(0,180)
        ax.axvline(0,color='black',lw=.7);ax.set(title=system,xlabel='Oriented IRC coordinate R',ylabel='F...H distance (Å)');ax.legend(fontsize=7,loc='upper left');ax2.legend(fontsize=7,loc='upper right')
    fig.suptitle('All amine-H / F contacts; highlighted pair fixed at TS (no switching)')
    fig.tight_layout();fig.savefig(figdir/'11_NH_F_contacts.png',dpi=190);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.8))
    for s,g in contacts[contacts.selected_at_TS].groupby('system',sort=False):
        axes[0].plot(g.CN_A,g.FH_A,color=SYSTEM_COLORS[s],label=s)
        row=summary[summary.system==s].iloc[0]
        axes[1].scatter(row.TS_nearest_FH_A,row.absolute_electronic_barrier_kcal_mol,color=SYSTEM_COLORS[s]);axes[1].annotate(s,(row.TS_nearest_FH_A,row.absolute_electronic_barrier_kcal_mol),xytext=(5,5),textcoords='offset points')
    axes[0].set(xlabel='Reaction C-N distance (Å)',ylabel='Fixed F...H distance (Å)',title='Common chemical progress coordinate');axes[0].legend()
    axes[1].set(xlabel='TS nearest amine-H...F distance (Å)',ylabel='Electronic barrier ΔE‡ (kcal/mol)',title='Four systems: descriptive comparison only');axes[1].margins(.25)
    fig.tight_layout();fig.savefig(figdir/'12_contact_vs_progress_and_barrier.png',dpi=190);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
    for ax,(target,group) in zip(axes.flat,corrected.groupby('target',sort=False)):
        for s,g in group.groupby('system',sort=False):ax.plot(g.R,g.COM_corrected_A,color=SYSTEM_COLORS[s],label=s)
        ax.set(title=target,ylabel='TS-referenced distance minus COM change (Å)',xlabel='IRC coordinate R');ax.legend(fontsize=7);ax.axvline(0,color='black',lw=.7)
    fig.tight_layout();fig.savefig(figdir/'13_COM_corrected_distances.png',dpi=190);plt.close(fig)


def write_report(summary,energy,orientation,stationary,endpoints,validation,warnings,reports,config):
    display=summary[['system','ortho_F_count','absolute_electronic_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol','TS_nearest_FH_A','TS_NHF_angle_deg','alpha_to_TS_kcal_mol']].copy()
    display.columns=['系','ortho F数','ΔE‡ / kcal mol−1','ΔG‡（RRHO）','TS F···H / Å','TS N–H···F / °','α→TS / kcal mol−1']
    absolute=energy[['system','absolute_electronic_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol','electronic_reaction_energy_kcal_mol','RRHO_G_reaction_energy_kcal_mol']].copy()
    absolute.columns=['系','ΔE‡','ΔG‡（RRHO）','ΔE反応','ΔG反応（RRHO）']
    work=energy[['system','separated_reactants_to_alpha_kcal_mol','alpha_to_TS_kcal_mol','absolute_electronic_barrier_kcal_mol']].copy()
    lys=next(r for r in stationary if r['system']=='Lys')
    order=' < '.join(energy.sort_values('absolute_electronic_barrier_kcal_mol').system)
    values=summary.set_index('system')
    energies=energy.set_index('system')
    corr=pd.read_csv(reports.parent/'csv/descriptive_correlations.csv')
    transition_corr=corr[(corr.region=='transition')&(corr.x=='FH_A')&(corr.y=='CN_A')][['system','n','pearson_r','spearman_r']]
    text=f'''# 4分子系のIRC・Reaction Force解析 summary

実行日：2026-10-17（JST）。入力：完成.zip内の16ログと、同条件で再計算された `LYS (1).LOG`（保存名 `Lys.log`）。旧結果を上書き更新。既存の2系解析から独立した4系解析。

## 背景と問い

ユーザーの仮説は、フッ素置換フェニル基のortho Fが求核剤アミンのHと接触することで反応を促進し、ortho Fがない系では反応性が下がる、というもの。距離・角度、反応力、構造変化の同期性から間接的に検討する。ここで「Lys」は入力ファイルの名称であり、ログの実体はC3H9N（13原子）のアミンモデルである。

## 数値結果

電子障壁の低い順は **{order}**。各系の分離反応物（基質＋Lys）を基準とする障壁そのものを表示する。MoriPFを引く系間のゼロ合わせは行わない。生の全電子エネルギー（Hartree）も `energy_comparison.csv` に保持するが、組成の異なる系の全電子エネルギーの大小を反応性とは解釈しない。

{markdown(display)}

表の接触はTSで最短のF···Hを選び、その原子対をIRC全体で追跡した値。F位置はフェニル環C23をipsoとして判定した。

![障壁と仕事](../figures/02_barriers_and_work.png)

分離反応物基準の障壁・反応エネルギー（すべてkcal/mol）：

{markdown(absolute)}

`ΔE‡ = E_TS − E_substrate − E_Lys`、`ΔG‡ = G_TS − G_substrate − G_Lys`。生成物についても同様に計算した。Lysの電子エネルギーは{lys['energy_hartree']:.9f} Hartree、RRHO Gibbsエネルギーは{lys['G_hartree']:.6f} Hartree。

電子障壁が負でも、分離した2分子の電子エネルギーよりTSが低いという意味であり、反応障壁が存在しないという意味ではない。RRHO Gibbs障壁は全系で正で、IRCの反応物側からTSへもエネルギー上昇がある。生成物側は与えられたopt構造についての値であり、特にRAFA15のIRCとの対応には後述の制限がある。

![絶対障壁](../figures/03_absolute_barriers.png)

## 仮説に対する解釈

**最新の対照検証：[分子内の強さ順位と、他の接触の相関](contact_specificity_summary_ja.md)。指定接触は局所応答として目立つ場合があるが、障壁と相関する接触は他にもあり、特異的な機構の証拠としては不十分。**

追加DFTなしの計算方針は [文献調査とlog限定の評価案](log_only_research_plan_ja.md) を参照。反対側ortho位置のH/F···Oまで含む追加定量化は、[両側接触・力と障壁の探索summary](two_sided_contact_summary_ja.md)を参照。従来の `interference_F` は最も近いFを選ぶ指標で、RAFA15の反対側H56···O21を表していないため、新解析では同じ環位置のF/Hを対応させた。以下はH···F単独の指標についての解釈であり、両側を組み合わせた仮説を否定するものではない。

**今回の結果は、ortho Fによる近接配置の形成とは整合するが、F···H接触だけで4系の障壁差を説明することはできない。**

- RAFA14（ortho Fなし）の最短F···Hは{values.loc['RAFA14','TS_nearest_FH_A']:.2f} Å、電子障壁は{values.loc['RAFA14','absolute_electronic_barrier_kcal_mol']:.2f} kcal/mol。
- MoriPFとRAFA16は約2.35 Å、N–H···F角約120°の近接配置を持ち、電子障壁はそれぞれ{values.loc['MoriPF','absolute_electronic_barrier_kcal_mol']:.2f}、{values.loc['RAFA16','absolute_electronic_barrier_kcal_mol']:.2f} kcal/mol。
- RAFA15は最短F···Hが{values.loc['RAFA15','TS_nearest_FH_A']:.2f} Å、角度{values.loc['RAFA15','TS_NHF_angle_deg']:.1f}°でも、電子障壁は{values.loc['RAFA15','absolute_electronic_barrier_kcal_mol']:.2f} kcal/molと最も低い。接触の短さだけでは順位を説明できない。
- RAFA15のα→TS仕事は{energies.loc['RAFA15','alpha_to_TS_kcal_mol']:.2f} kcal/molと最大。分離反応物→αは{energies.loc['RAFA15','separated_reactants_to_alpha_kcal_mol']:.2f} kcal/mol。この区間には会合・配座変化が含まれ、特定の接触の安定化には分解できない。

ortho Fの有無、F···H距離、N–H···F角度と障壁を合わせて評価する。短い接触や大きな相関係数だけでは、水素結合・安定化エネルギー・因果関係の証明にならない。F置換に伴う電子効果、分子の配向、立体効果、他の接触も同時に変わる。また、MoriPFはF数自体がRAFA14–16と異なるため、位置だけを比較する対照ではない。

全F×アミンHの距離と角度を `all_NH_F_contacts.csv` に保存した。RAFA15ではアミンの配向も異なり、同じ化学的位置のFが存在しないRAFA14との単一ペア比較には限界がある。固定ペアの距離と反応C–N距離を同じ横軸で比較した図も作成した。

![全接触](../figures/11_NH_F_contacts.png)

![進行度と接触](../figures/12_contact_vs_progress_and_barrier.png)

## Reaction Forceと仕事の分解

`F(R) = −dE/dR`。αは反応物側でFが最小となる点、β=0はTS、γは生成物側でFが最大となる点。エネルギーを平滑化スプラインで近似し、解析的な微分を用いた。極値は全域の細かい探索と局所最適化で求めた。

**反応力定数は参照論文の `κ = d²E/dR² = −dF/dR` を図に使用する。旧解析のκは逆符号だった。** CSVには論文規約を `kappa_paper_*` として保存し、旧規約の列も比較用に保持する。

`α→TS = E_TS − E(α)` は各IRC内部の遷移前半の仕事であり、分離反応物からの全障壁ではない。`IRC_left_to_alpha` は計算済み左端からαまでの区間で、系によってIRCの到達範囲が異なる。これをそのまま全構造準備コストの順位として扱わない。分離反応物からαまでの値は `separated_reactants_to_alpha_kcal_mol` に保存する。

{markdown(work)}

比較図・個別系の図ともに `ΔE(R) = E(R) − E_substrate − E_Lys` を使い、TSの高さはその系の電子障壁を表す。

![反応力](../figures/01_four_systems_energy_force_kappa.png)

`descriptive_correlations.csv` は全経路・反応物側・遷移領域・TS前後のPearson/Spearman相関を収録する。IRC点は時系列状に連続しており独立標本ではないため、p値による有意性を主張しない。距離同士の相関には共通の反応進行が混入するため、符号付き微分同士の相関も併記した。4系のみから機構を回帰的に確定しない。

遷移領域のF···H距離とC–N距離の記述相関：

{markdown(transition_corr)}

遠い接触しかないRAFA14でも相関が約0.975に達する。これは「高相関そのものがH···F相互作用の証拠ではない」ことを、このデータ内で示す重要な対照である。RAFA15は負相関で、結合形成に伴って選択したF···Hがむしろ離れる方向に動く。

## IRC方向とTSの確認

全系でC20–求核N距離が長い側を負、短い側を正とした。MoriPF・RAFA15・RAFA16は元のparser座標を反転、RAFA14はそのままとした。RAFA16には最大0.0114 Åの局所的なC–N伸長があり、全点が厳密に単調減少するわけではないが、端点間の反応方向は明瞭である。元のGaussian出力は圧縮してバイト単位で保存し、解析用LOGを `outputs/oriented_logs/` に別途生成した。これは新規Gaussian計算ではなく、座標・エネルギー・力・反応座標を正規化した派生LOGである。解析はこのLOGを読み直して実施し、書き出し前後の値の一致を検証した。

{markdown(pd.DataFrame(orientation),['system','coordinate_multiplier','reactant_side_CN_A','product_side_CN_A','CN_monotonically_decreasing','point_count','TS_midpoint_alignment_RMSD_A'])}

TSは対応するfreqログの実座標・電子エネルギーを使用した。IRC内のTSエネルギーとの一致と、近傍構造への重ね合わせを確認した。全TSは虚振動1個。ただし虚振動固有ベクトルが反応C–N形成モードに対応するかの可視化確認はこの解析では実施していない。未収束末端点（反応座標の記録が欠ける点）は除外し、外挿しない。

## 計算条件・対応関係の制限

分子側の反応物・TS・生成物と新しいLysはすべてM06-2X/6-311++G(d,p)、SMD水。Lysはopt/freq正常終了、13原子・33振動モード・虚振動0個。組成、電荷、計算レベル、温度・圧力、熱化学の有無を確認してから絶対障壁を計算する。旧Lysに由来する異条件の算術値と警告は今回の出力から削除した。

ΔG‡はGaussian標準RRHO値による298.15 K・1 atmのモデル値であり、1 M標準状態補正や低振動数の準調和補正はしていない。1 atmから1 Mへの補正を行う場合、会合型反応の絶対値が変わるが、同一分子数変化の系間比較では共通項は相殺する。低振動数・配座差の誤差は残る。RAFA15反応物には−5.4036 cm⁻¹の小さな虚振動があり、この系のRRHO値は暫定としてフラグを残している。

`absolute_barrier_valid` は計算条件・組成が整合した算術値であることを示す。全反応物が厳密な極小であることや実験速度の再現を保証するフラグではない。別列 `reactant_minimum_verified` と `RRHO_G_provisional` を参照する。単一の与えられた配座・経路についての比較なので、速度定数や実験的反応性の順位は確定しない。

{markdown(endpoints)}

生成物optとIRC末端のRMSDは、原子番号をそのまま対応させた重ね合わせの診断値であり、対称等価原子の入れ替え探索はしていない。特にopt生成物がIRC末端より高エネルギーの場合、別配座／別経路の可能性がある。両者を滑らかな接続線で結ばない。

## 力・運動解析の検証

{markdown(validation)}

`Σ F_i · dr_i/dR` とエネルギー由来反応力の一致を確認する診断。原子・断片の寄与や2原子の力射影は、H···F相互作用エネルギーの一意な分解ではない。平滑化への依存は `smoothing_sensitivity.csv` に記録した。

## 次に確認したい点

- RAFA15反応物の微小虚振動を解消し、各系の配座探索とIRCからの生成物最適化で比較条件をそろえる。
- F···Hの安定化を直接論じるには、電子密度指標（NCI/QTAIM等）や適切なエネルギー分解・配座対照を追加する。これらは今回未計算。

## 実行時の警告

'''+'\n'.join('- '+w for w in warnings)+'''

## 参照文献

- Toro-Labbé et al. (2009), *The reaction force and the transition region of a reaction*, J. Mol. Model. 15, 707–710. [DOI](https://doi.org/10.1007/s00894-008-0431-8)、[PubMed要旨（κの定義を確認）](https://pubmed.ncbi.nlm.nih.gov/19085021/)。構造準備／遷移領域という解釈の理論的背景。
- Dunitz & Taylor (1997), *Organic Fluorine Hardly Ever Accepts Hydrogen Bonds*, Chem. Eur. J. 3, 89–98. [DOI](https://doi.org/10.1002/chem.19970030115)。有機Fの水素結合受容性を距離だけで判断しないための背景資料。

## 再現方法

ルートで `python -m pip install -r requirements.txt`、続いて `python analysis4_force/run_analysis.py`。圧縮原本は自動展開し、manifestのSHA-256と照合する。原子対応表・全数値・図・このsummaryはスクリプトから生成する。
'''
    (reports/'summary_ja.md').write_text(text)
