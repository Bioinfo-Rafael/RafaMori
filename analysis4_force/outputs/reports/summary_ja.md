# 4分子系のIRC・Reaction Force解析 summary

実行日：2026-10-17（JST）。入力：完成.zip内の16ログと、同条件で再計算された `LYS (1).LOG`（保存名 `Lys.log`）。旧結果を上書き更新。既存の2系解析から独立した4系解析。

## 背景と問い

ユーザーの仮説は、フッ素置換フェニル基のortho Fが求核剤アミンのHと接触することで反応を促進し、ortho Fがない系では反応性が下がる、というもの。距離・角度、反応力、構造変化の同期性から間接的に検討する。ここで「Lys」は入力ファイルの名称であり、ログの実体はC3H9N（13原子）のアミンモデルである。

## 数値結果

電子エネルギーによる相対障壁の低い順は **RAFA15 < MoriPF < RAFA16 < RAFA14**。系間比較はすべて **MoriPF＝0** とする。新しいLysは分子側と計算条件が一致するため、分離反応物（分子＋Lys）基準の絶対障壁も算出した。以下のΔΔ値はMoriPFとの差、Δ値は分離反応物との差を表す。

| 系 | ortho F数 | 相対ΔE‡ / kcal mol−1 | 相対ΔG‡（RRHO） | TS F···H / Å | TS N–H···F / ° | α→TS / kcal mol−1 |
| --- | --- | --- | --- | --- | --- | --- |
| MoriPF | 2 | 0.000 | 0.000 | 2.358 | 120.526 | 1.468 |
| RAFA14 | 0 | 1.818 | 1.845 | 4.795 | 118.059 | 2.022 |
| RAFA15 | 1 | -3.245 | -5.141 | 2.622 | 101.335 | 2.578 |
| RAFA16 | 2 | 0.763 | 0.940 | 2.345 | 120.620 | 1.603 |

相対障壁の定義：`ΔΔE‡(i) = [E_TS(i) − E_substrate(i)] − [E_TS(MoriPF) − E_substrate(MoriPF)]`。表のTS接触は各系でTSにおける最短のF···H（アミンNに結合した2つのHを対象）を選び、その同じ原子対をIRC全体で追跡した値。Fの位置はフェニル環C23をipsoとして結合グラフから判定した。

![相対障壁](../figures/02_relative_barriers.png)

分離反応物基準の障壁・反応エネルギー（すべてkcal/mol）：

| 系 | ΔE‡ | ΔG‡（RRHO） | ΔE反応 | ΔG反応（RRHO） |
| --- | --- | --- | --- | --- |
| MoriPF | -1.342 | 14.083 | -5.386 | 11.347 |
| RAFA14 | 0.476 | 15.928 | -1.461 | 15.725 |
| RAFA15 | -4.587 | 8.943 | -2.837 | 12.233 |
| RAFA16 | -0.579 | 15.023 | -3.653 | 12.700 |

`ΔE‡ = E_TS − E_substrate − E_Lys`、`ΔG‡ = G_TS − G_substrate − G_Lys`。生成物についても同様に計算した。Lysの電子エネルギーは-174.451556589 Hartree、RRHO Gibbsエネルギーは-174.357693 Hartree。

電子障壁が負でも、分離した2分子の電子エネルギーよりTSが低いという意味であり、反応障壁が存在しないという意味ではない。RRHO Gibbs障壁は全系で正で、IRCの反応物側からTSへもエネルギー上昇がある。生成物側は与えられたopt構造についての値であり、特にRAFA15のIRCとの対応には後述の制限がある。

![絶対障壁](../figures/03_absolute_barriers.png)

## 仮説に対する解釈

**今回の結果は、ortho Fによる近接配置の形成とは整合するが、F···H接触だけで4系の障壁差を説明することはできない。**

- RAFA14（ortho Fなし）は最短F···Hが4.79 Åと遠く、MoriPFより相対電子障壁が1.82 kcal/mol高い。この2系の対比は提案された仮説と整合する。
- MoriPFとRAFA16はともに約2.35 Å、N–H···F角約120°の近接配置を持つ。RAFA16の相対障壁はMoriPFより約0.76 kcal/mol高いが、他のF置換や電子効果との切り分けはできない。
- **RAFA15が単純な説明への反例になる。** 最短F···Hは2.62 Å、角度は101.3°で、MoriPF/RAFA16より遠く、直線からも外れる。それでも相対電子障壁はMoriPFより3.24 kcal/mol低い。最も短い接触やortho Fの個数だけでは順位を説明できない。
- RAFA15は反応物基準では低障壁でも、α→TSの仕事は約2.58 kcal/molで4系中最大。この区間の仕事を全活性化障壁と同一視すると逆の結論になる。反応物基準・配座・経路全体を合わせて見る必要がある。

MoriPF基準でRAFA15の全電子障壁は3.24 kcal/mol低い一方、α→TSの仕事は1.11 kcal/mol大きい。分離反応物からαまでのエネルギー差はMoriPFより4.35 kcal/mol低く、この差が全障壁に寄与している。ただしこの区間には分子間会合や配座変化が含まれ、H···F安定化だけには分解できない。MoriPFとRAFA14–16は組成が異なるので、生のTSエネルギーや分子単体エネルギーをそのまま差し引いた比較は行わない。

ortho Fの有無、F···H距離、N–H···F角度と障壁を合わせて評価する。短い接触や大きな相関係数だけでは、水素結合・安定化エネルギー・因果関係の証明にならない。F置換に伴う電子効果、分子の配向、立体効果、他の接触も同時に変わる。また、MoriPFはF数自体がRAFA14–16と異なるため、位置だけを比較する対照ではない。

全F×アミンHの距離と角度を `all_NH_F_contacts.csv` に保存した。RAFA15ではアミンの配向も異なり、同じ化学的位置のFが存在しないRAFA14との単一ペア比較には限界がある。固定ペアの距離と反応C–N距離を同じ横軸で比較した図も作成した。

![全接触](../figures/11_NH_F_contacts.png)

![進行度と接触](../figures/12_contact_vs_progress_and_barrier.png)

## Reaction Forceと仕事の分解

`F(R) = −dE/dR`。αは反応物側でFが最小となる点、β=0はTS、γは生成物側でFが最大となる点。エネルギーを平滑化スプラインで近似し、解析的な微分を用いた。極値は全域の細かい探索と局所最適化で求めた。

**反応力定数は参照論文の `κ = d²E/dR² = −dF/dR` を図に使用する。旧解析のκは逆符号だった。** CSVには論文規約を `kappa_paper_*` として保存し、旧規約の列も比較用に保持する。

`α→TS = E_TS − E(α)` は各IRC内部の遷移前半の仕事であり、分離反応物からの全障壁ではない。`IRC_left_to_alpha` は計算済み左端からαまでの区間で、系によってIRCの到達範囲が異なる。これをそのまま全構造準備コストの順位として扱わない。分離反応物からαまでの系間差は `relative_preparation_component_kcal_mol` に別途記録した。

| system | separated_reactants_to_alpha_kcal_mol | alpha_to_TS_kcal_mol | relative_preparation_component_kcal_mol | relative_alpha_to_TS_kcal_mol |
| --- | --- | --- | --- | --- |
| MoriPF | -2.810 | 1.468 | 0.000 | 0.000 |
| RAFA14 | -1.547 | 2.022 | 1.264 | 0.554 |
| RAFA15 | -7.165 | 2.578 | -4.355 | 1.110 |
| RAFA16 | -2.182 | 1.603 | 0.629 | 0.134 |

上の `relative_*` 列はMoriPFとの差。比較図のエネルギー軸は各系で分離反応物を差し引いた後、MoriPFの電子障壁を引き、MoriPFのTSを0とする。各系単独の形状診断図だけはその系のTS＝0を使う。

![反応力](../figures/01_four_systems_energy_force_kappa.png)

`descriptive_correlations.csv` は全経路・反応物側・遷移領域・TS前後のPearson/Spearman相関を収録する。IRC点は時系列状に連続しており独立標本ではないため、p値による有意性を主張しない。距離同士の相関には共通の反応進行が混入するため、符号付き微分同士の相関も併記した。4系のみから機構を回帰的に確定しない。

遷移領域のF···H距離とC–N距離の記述相関：

| system | n | pearson_r | spearman_r |
| --- | --- | --- | --- |
| MoriPF | 30 | 0.981 | 1.000 |
| RAFA14 | 30 | 0.975 | 0.982 |
| RAFA15 | 29 | -0.978 | -1.000 |
| RAFA16 | 30 | 0.979 | 0.998 |

遠い接触しかないRAFA14でも相関が約0.975に達する。これは「高相関そのものがH···F相互作用の証拠ではない」ことを、このデータ内で示す重要な対照である。RAFA15は負相関で、結合形成に伴って選択したF···Hがむしろ離れる方向に動く。

## IRC方向とTSの確認

全系でC20–求核N距離が長い側を負、短い側を正とした。MoriPF・RAFA15・RAFA16は元のparser座標を反転、RAFA14はそのままとした。RAFA16には最大0.0114 Åの局所的なC–N伸長があり、全点が厳密に単調減少するわけではないが、端点間の反応方向は明瞭である。元のGaussian出力は圧縮してバイト単位で保存し、解析用LOGを `outputs/oriented_logs/` に別途生成した。これは新規Gaussian計算ではなく、座標・エネルギー・力・反応座標を正規化した派生LOGである。解析はこのLOGを読み直して実施し、書き出し前後の値の一致を検証した。

| system | coordinate_multiplier | reactant_side_CN_A | product_side_CN_A | CN_monotonically_decreasing | point_count | TS_midpoint_alignment_RMSD_A |
| --- | --- | --- | --- | --- | --- | --- |
| MoriPF | -1 | 2.617 | 1.639 | True | 64 | 0.000 |
| RAFA14 | 1 | 2.581 | 1.652 | True | 60 | 0.000 |
| RAFA15 | -1 | 2.590 | 1.624 | True | 57 | 0.000 |
| RAFA16 | -1 | 2.603 | 1.641 | False | 71 | 0.000 |

TSは対応するfreqログの実座標・電子エネルギーを使用した。IRC内のTSエネルギーとの一致と、近傍構造への重ね合わせを確認した。全TSは虚振動1個。ただし虚振動固有ベクトルが反応C–N形成モードに対応するかの可視化確認はこの解析では実施していない。未収束末端点（反応座標の記録が欠ける点）は除外し、外挿しない。

## 計算条件・対応関係の制限

分子側の反応物・TS・生成物と新しいLysはすべてM06-2X/6-311++G(d,p)、SMD水。Lysはopt/freq正常終了、13原子・33振動モード・虚振動0個。組成、電荷、計算レベル、温度・圧力、熱化学の有無を確認してから絶対障壁を計算する。旧Lysに由来する異条件の算術値と警告は今回の出力から削除した。

ΔG‡はGaussian標準RRHO値による298.15 K・1 atmのモデル値であり、1 M標準状態補正や低振動数の準調和補正はしていない。1 atmから1 Mへの補正を行う場合、会合型反応の絶対値が変わるが、同一分子数変化のMoriPF基準の系間差では共通項は相殺する。低振動数・配座差の誤差は残る。RAFA15反応物には−5.4036 cm⁻¹の小さな虚振動があり、この系のRRHO値は暫定としてフラグを残している。

`absolute_barrier_valid` は計算条件・組成が整合した算術値であることを示す。全反応物が厳密な極小であることや実験速度の再現を保証するフラグではない。別列 `reactant_minimum_verified` と `RRHO_G_provisional` を参照する。単一の与えられた配座・経路についての比較なので、速度定数や実験的反応性の順位は確定しない。

| system | product_opt_vs_IRC_RMSD_A | product_opt_minus_IRC_endpoint_kcal_mol | product_opt_CN_A | IRC_product_CN_A | composition_balanced_with_Lys |
| --- | --- | --- | --- | --- | --- |
| MoriPF | 0.774 | -1.804 | 1.582 | 1.639 | True |
| RAFA14 | 0.393 | -0.489 | 1.622 | 1.652 | True |
| RAFA15 | 2.111 | 3.708 | 1.584 | 1.624 | True |
| RAFA16 | 1.071 | -1.160 | 1.588 | 1.641 | True |

生成物optとIRC末端のRMSDは、原子番号をそのまま対応させた重ね合わせの診断値であり、対称等価原子の入れ替え探索はしていない。特にopt生成物がIRC末端より高エネルギーの場合、別配座／別経路の可能性がある。両者を滑らかな接続線で結ばない。

## 力・運動解析の検証

| system | n | pearson_r | normalized_RMSE |
| --- | --- | --- | --- |
| MoriPF | 63 | 1.000 | 0.003 |
| RAFA14 | 59 | 1.000 | 0.002 |
| RAFA15 | 56 | 1.000 | 0.002 |
| RAFA16 | 70 | 1.000 | 0.023 |

`Σ F_i · dr_i/dR` とエネルギー由来反応力の一致を確認する診断。原子・断片の寄与や2原子の力射影は、H···F相互作用エネルギーの一意な分解ではない。平滑化への依存は `smoothing_sensitivity.csv` に記録した。

## 次に確認したい点

- RAFA15反応物の微小虚振動を解消し、各系の配座探索とIRCからの生成物最適化で比較条件をそろえる。
- F···Hの安定化を直接論じるには、電子密度指標（NCI/QTAIM等）や適切なエネルギー分解・配座対照を追加する。これらは今回未計算。

## 実行時の警告

- RAFA15 reactant: imaginary frequencies [-5.4036] cm^-1; minimum not fully verified
- MoriPF: omitted 1 unconverged terminal point(s); no extrapolated coordinates used
- RAFA15: omitted 1 unconverged terminal point(s); no extrapolated coordinates used
- RAFA15: supplied optimized product is higher than the IRC product-side endpoint; conformer/path correspondence needs checking
- RAFA16: local C-N increases up to 0.01138 A; overall orientation follows long-to-short endpoints

## 参照文献

- Toro-Labbé et al. (2009), *The reaction force and the transition region of a reaction*, J. Mol. Model. 15, 707–710. [DOI](https://doi.org/10.1007/s00894-008-0431-8)、[PubMed要旨（κの定義を確認）](https://pubmed.ncbi.nlm.nih.gov/19085021/)。構造準備／遷移領域という解釈の理論的背景。
- Dunitz & Taylor (1997), *Organic Fluorine Hardly Ever Accepts Hydrogen Bonds*, Chem. Eur. J. 3, 89–98. [DOI](https://doi.org/10.1002/chem.19970030115)。有機Fの水素結合受容性を距離だけで判断しないための背景資料。

## 再現方法

ルートで `python -m pip install -r requirements.txt`、続いて `python analysis4_force/run_analysis.py`。圧縮原本は自動展開し、manifestのSHA-256と照合する。原子対応表・全数値・図・このsummaryはスクリプトから生成する。
