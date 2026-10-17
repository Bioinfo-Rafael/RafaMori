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


def create_extra_figures(summary,energy,contacts,corrected,figdir):
    fig,axes=plt.subplots(1,3,figsize=(14,4.7))
    for ax,col,title in zip(axes,['relative_electronic_barrier_kcal_mol','relative_RRHO_G_barrier_kcal_mol','alpha_to_TS_kcal_mol'],['Relative electronic barrier (RAFA16=0)','Relative RRHO Gibbs barrier (provisional)','IRC alpha-to-TS work W2']):
        bars=ax.bar(summary.system,summary[col],color=[SYSTEM_COLORS[s] for s in summary.system]);ax.bar_label(bars,fmt='%.2f',padding=3);ax.set(title=title,ylabel='kcal/mol');ax.axhline(0,color='black',lw=.7);ax.margins(y=.22)
    fig.suptitle('Common Lys term cancels in relative barriers; absolute barriers unavailable')
    fig.tight_layout();fig.savefig(figdir/'02_relative_barriers.png',dpi=190);plt.close(fig)
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
        axes[1].scatter(row.TS_nearest_FH_A,row.relative_electronic_barrier_kcal_mol,color=SYSTEM_COLORS[s]);axes[1].annotate(s,(row.TS_nearest_FH_A,row.relative_electronic_barrier_kcal_mol),xytext=(5,5),textcoords='offset points')
    axes[0].set(xlabel='Reaction C-N distance (Å)',ylabel='Fixed F...H distance (Å)',title='Common chemical progress coordinate');axes[0].legend()
    axes[1].set(xlabel='TS nearest amine-H...F distance (Å)',ylabel='Relative electronic barrier (kcal/mol)',title='Four systems: descriptive comparison only');axes[1].margins(.25)
    fig.tight_layout();fig.savefig(figdir/'12_contact_vs_progress_and_barrier.png',dpi=190);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
    for ax,(target,group) in zip(axes.flat,corrected.groupby('target',sort=False)):
        for s,g in group.groupby('system',sort=False):ax.plot(g.R,g.COM_corrected_A,color=SYSTEM_COLORS[s],label=s)
        ax.set(title=target,ylabel='TS-referenced distance minus COM change (Å)',xlabel='IRC coordinate R');ax.legend(fontsize=7);ax.axvline(0,color='black',lw=.7)
    fig.tight_layout();fig.savefig(figdir/'13_COM_corrected_distances.png',dpi=190);plt.close(fig)


def write_report(summary,energy,orientation,stationary,endpoints,validation,warnings,reports,config):
    display=summary[['system','ortho_F_count','relative_electronic_barrier_kcal_mol','relative_RRHO_G_barrier_kcal_mol','TS_nearest_FH_A','TS_NHF_angle_deg','alpha_to_TS_kcal_mol']].copy()
    display.columns=['系','ortho F数','相対ΔE‡ / kcal mol−1','相対ΔG‡（暫定）','TS F···H / Å','TS N–H···F / °','α→TS / kcal mol−1']
    order=' < '.join(energy.sort_values('relative_electronic_barrier_kcal_mol').system)
    values=summary.set_index('system')
    energies=energy.set_index('system')
    r15_substrate_gap=(energies.loc['RAFA15','substrate_E_hartree']-energies.loc['RAFA16','substrate_E_hartree'])*627.509474
    r15_ts_gap=(energies.loc['RAFA15','TS_E_hartree']-energies.loc['RAFA16','TS_E_hartree'])*627.509474
    corr=pd.read_csv(reports.parent/'csv/descriptive_correlations.csv')
    transition_corr=corr[(corr.region=='transition')&(corr.x=='FH_A')&(corr.y=='CN_A')][['system','n','pearson_r','spearman_r']]
    text=f'''# 4分子系のIRC・Reaction Force解析 summary

実行日：2026-10-17（JST）。入力：完成.zip内の16ログとLys.log。既存の2系解析から独立した4系解析。

## 背景と問い

ユーザーの仮説は、フッ素置換フェニル基のortho Fが求核剤アミンのHと接触することで反応を促進し、ortho Fがない系では反応性が下がる、というもの。距離・角度、反応力、構造変化の同期性から間接的に検討する。ここで「Lys」は入力ファイルの名称であり、ログの実体はC3H9N（13原子）のアミンモデルである。

## 数値結果

電子エネルギーによる相対障壁の低い順は **{order}**。RAFA16を0とする。同じLys反応物を仮定すると、系間の障壁差ではLysのエネルギーが消えるため、Lysの異なる計算条件を混ぜずに比較できる。

{markdown(display)}

相対障壁の定義：`ΔΔE‡(i) = [E_TS(i) − E_substrate(i)] − [E_TS(RAFA16) − E_substrate(RAFA16)]`。表のTS接触は各系でTSにおける最短のF···H（アミンNに結合した2つのHを対象）を選び、その同じ原子対をIRC全体で追跡した値。Fの位置はフェニル環C23をipsoとして結合グラフから判定した。

![相対障壁](../figures/02_relative_barriers.png)

## 仮説に対する解釈

**今回の結果は、ortho Fによる近接配置の形成とは整合するが、F···H接触だけで4系の障壁差を説明することはできない。**

- RAFA14（ortho Fなし）は最短F···Hが{values.loc['RAFA14','TS_nearest_FH_A']:.2f} Åと遠く、RAFA16より相対電子障壁が{values.loc['RAFA14','relative_electronic_barrier_kcal_mol']:.2f} kcal/mol高い。この2系の対比は提案された仮説と整合する。
- MoriPFとRAFA16はともに約2.35 Å、N–H···F角約120°の近接配置を持つ。MoriPFの相対障壁はRAFA16より約0.76 kcal/mol低いが、他のF置換や電子効果との切り分けはできない。
- **RAFA15が単純な説明への反例になる。** 最短F···Hは{values.loc['RAFA15','TS_nearest_FH_A']:.2f} Å、角度は{values.loc['RAFA15','TS_NHF_angle_deg']:.1f}°で、MoriPF/RAFA16より遠く、直線からも外れる。それでも相対電子障壁はRAFA16より{abs(values.loc['RAFA15','relative_electronic_barrier_kcal_mol']):.2f} kcal/mol低い。最も短い接触やortho Fの個数だけでは順位を説明できない。
- RAFA15は反応物基準では低障壁でも、α→TSの仕事は約2.58 kcal/molで4系中最大。この区間の仕事を全活性化障壁と同一視すると逆の結論になる。反応物基準・配座・経路全体を合わせて見る必要がある。

同組成のRAFA15とRAFA16を比較すると、RAFA15のTSそのものはRAFA16より{r15_ts_gap:.2f} kcal/mol高い。一方、与えられたRAFA15反応物は{r15_substrate_gap:.2f} kcal/mol高い。その差し引きにより相対障壁が低くなる。したがって今回のRAFA15の低障壁を、TSでのH···Fによる強い安定化の証拠とは読めない。

ortho Fの有無、F···H距離、N–H···F角度と障壁を合わせて評価する。短い接触や大きな相関係数だけでは、水素結合・安定化エネルギー・因果関係の証明にならない。F置換に伴う電子効果、分子の配向、立体効果、他の接触も同時に変わる。また、MoriPFはF数自体がRAFA14–16と異なるため、位置だけを比較する対照ではない。

全F×アミンHの距離と角度を `all_NH_F_contacts.csv` に保存した。RAFA15ではアミンの配向も異なり、同じ化学的位置のFが存在しないRAFA14との単一ペア比較には限界がある。固定ペアの距離と反応C–N距離を同じ横軸で比較した図も作成した。

![全接触](../figures/11_NH_F_contacts.png)

![進行度と接触](../figures/12_contact_vs_progress_and_barrier.png)

## Reaction Forceと仕事の分解

`F(R) = −dE/dR`。αは反応物側でFが最小となる点、β=0はTS、γは生成物側でFが最大となる点。エネルギーを平滑化スプラインで近似し、解析的な微分を用いた。極値は全域の細かい探索と局所最適化で求めた。

**反応力定数は参照論文の `κ = d²E/dR² = −dF/dR` を図に使用する。旧解析のκは逆符号だった。** CSVには論文規約を `kappa_paper_*` として保存し、旧規約の列も比較用に保持する。

`α→TS = E_TS − E(α)` は各IRC内部の遷移前半の仕事であり、分離反応物からの全障壁ではない。`IRC_left_to_alpha` は計算済み左端からαまでの区間で、系によってIRCの到達範囲が異なる。これをそのまま全構造準備コストの順位として扱わない。分離反応物からαまでの系間差は `relative_preparation_component_kcal_mol` に別途記録した。

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

分子側の反応物・TS・生成物はM06-2X/6-311++G(d,p)、SMD水。`Lys.log` はB3LYP/6-311+G(2d,p)、水PCMのoptであり、freq・熱化学補正はない。したがって絶対ΔE‡・ΔG‡は有効な結果として報告しない。異条件の単純和は監査用CSVの `mixed_level_barrier_INVALID_kcal_mol` に明示的に隔離した。

相対ΔG‡は各分子のGaussian標準RRHO値の差から共通Lys項を相殺した暫定値。298.15 K・1 atmの出力値であり、1 M標準状態補正や低振動数の準調和補正はしていない。同一分子数変化の系間差では共通の標準状態項は相殺するが、低振動数・配座差の誤差は残る。RAFA15反応物に−5.4036 cm⁻¹の小さな虚振動があり、厳密な極小構造としての再確認が必要。電子障壁の比較も単一の与えられた配座・経路についての比較であり、速度定数や実験的反応性の順位を確定しない。

{markdown(endpoints)}

生成物optとIRC末端のRMSDは、原子番号をそのまま対応させた重ね合わせの診断値であり、対称等価原子の入れ替え探索はしていない。特にopt生成物がIRC末端より高エネルギーの場合、別配座／別経路の可能性がある。両者を滑らかな接続線で結ばない。

## 力・運動解析の検証

{markdown(validation)}

`Σ F_i · dr_i/dR` とエネルギー由来反応力の一致を確認する診断。原子・断片の寄与や2原子の力射影は、H···F相互作用エネルギーの一意な分解ではない。平滑化への依存は `smoothing_sensitivity.csv` に記録した。

## 次に確認したい点

- Lysを同じM06-2X/6-311++G(d,p)、SMD水でopt/freqし、絶対障壁を整合した条件で評価する。
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
