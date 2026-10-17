# 両ortho位置の接触・力と活性化障壁の探索解析

## 結論

**障壁の順位と対応する「力の配分」の候補は得られた。ただし、接触方向の力の比率は弱い対応にとどまり、F···O反発／H···O引力を定量した結果ではない。**

ユーザーの仮説は、アミン側のF···Hに加え、反対側のortho H···クマリンOが有利に働く一方、同位置がFならOとの反発が不利に働く、という両側のバランス。前の解析は「最も近いF」を選んでいたため、RAFA15のH56···O21を含まなかった。今回は同じ環位置を対応させて追加解析した。

C–N結合長2.50→2.10 Åの共通区間で集計すると、**反対側のH/FとOに働く全原子力の大きさの割合**はRAFA15＜MoriPF＜RAFA16＜RAFA14で、電子障壁・Gibbs障壁の順位と一致する。

| 系 | 反対側原子力の割合 (%) | MoriPFとの差 (ポイント) | ΔΔE‡ (kcal/mol) | ΔΔG‡ (kcal/mol) |
| --- | --- | --- | --- | --- |
| MoriPF | 58.224 | 0.000 | 0.000 | 0.000 |
| RAFA14 | 59.999 | 1.775 | 1.818 | 1.845 |
| RAFA15 | 57.113 | -1.111 | -3.245 | -5.141 |
| RAFA16 | 58.987 | 0.763 | 0.763 | 0.940 |

この二次的な候補のPearson相関は、電子障壁で **r=0.958**、RRHO Gibbs障壁で **r=0.922**、順位相関はいずれも1.0。一方、最初に指定した主指標「接触方向に射影した力の大きさの比率」は電子障壁でr≈−0.393、順位相関−0.4で、十分な対応を示さなかった。良い指標だけを残さず、12種類の指標を全てCSVと下表に掲載した。

![両側比率と障壁](../figures/15_two_sided_barrier_correlations.png)

## 原子の対応と接触形状

| system | near_ring_C | near_atom | near_element | primary_amine_H | alternative_amine_H | opposite_ring_C | opposite_atom | opposite_element | opposite_O |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MoriPF | 24 | 33 | F | 46 | 47 | 25 | 29 | F | 21 |
| RAFA14 | 24 | 57 | H | 44 | 45 | 25 | 56 | H | 21 |
| RAFA15 | 24 | 31 | F | 44 | 45 | 25 | 56 | H | 21 |
| RAFA16 | 24 | 31 | F | 44 | 45 | 25 | 29 | F | 21 |

C24がアミン側、C25が反対側。RAFA15の反対側はH56–O21、MoriPF・RAFA16はF29–O21。RAFA14は両側のortho位置がHで、アミン側はH57···H44という対照になる。RAFA14のアミン側をF···Hと呼ぶことはできない。主解析のアミンHはTSでC24上の置換原子に近いHを固定し、もう一方のHでも感度を調べた。

| system | near_TS_distance_A | opposite_TS_distance_A | opposite_TS_CXO_angle_deg |
| --- | --- | --- | --- |
| MoriPF | 2.358 | 2.973 | 77.568 |
| RAFA14 | 2.358 | 2.807 | 87.127 |
| RAFA15 | 2.622 | 2.764 | 90.436 |
| RAFA16 | 2.345 | 2.991 | 78.380 |

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

| system | near_X_mean_force_mEh_bohr | amine_H_mean_force_mEh_bohr | opposite_X_mean_force_mEh_bohr | opposite_O_mean_force_mEh_bohr | near_mean_signed_projected_force_mEh_bohr | opposite_mean_signed_projected_force_mEh_bohr |
| --- | --- | --- | --- | --- | --- | --- |
| MoriPF | 0.370 | 0.239 | 0.331 | 0.518 | -0.046 | -0.111 |
| RAFA14 | 0.058 | 0.315 | 0.048 | 0.510 | 0.039 | -0.182 |
| RAFA15 | 0.768 | 0.380 | 0.050 | 1.477 | -0.097 | -0.514 |
| RAFA16 | 0.460 | 0.262 | 0.438 | 0.601 | -0.060 | -0.120 |

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

| descriptor | target | pearson_r | spearman_r | LOO_prediction_RMSE_kcal_mol |
| --- | --- | --- | --- | --- |
| opposite_projected_magnitude_share | relative_electronic_barrier_kcal_mol | -0.393 | -0.400 | 3.427 |
| opposite_projected_magnitude_share | relative_RRHO_G_barrier_kcal_mol | -0.485 | -0.400 | 4.679 |
| opposite_atomic_force_share | relative_electronic_barrier_kcal_mol | 0.958 | 1.000 | 1.410 |
| opposite_atomic_force_share | relative_RRHO_G_barrier_kcal_mol | 0.922 | 1.000 | 2.761 |
| selected_atomic_force_share | relative_electronic_barrier_kcal_mol | -0.632 | -0.400 | 3.158 |
| selected_atomic_force_share | relative_RRHO_G_barrier_kcal_mol | -0.553 | -0.400 | 5.625 |
| selected_atomic_work_kcal_mol | relative_electronic_barrier_kcal_mol | -0.937 | -0.800 | 1.648 |
| selected_atomic_work_kcal_mol | relative_RRHO_G_barrier_kcal_mol | -0.961 | -0.800 | 2.990 |
| pair_diagnostic_work_kcal_mol | relative_electronic_barrier_kcal_mol | -0.807 | -0.400 | 3.981 |
| pair_diagnostic_work_kcal_mol | relative_RRHO_G_barrier_kcal_mol | -0.866 | -0.400 | 4.882 |
| near_mean_abs_projected_force_mEh_bohr | relative_electronic_barrier_kcal_mol | -0.932 | -0.800 | 1.559 |
| near_mean_abs_projected_force_mEh_bohr | relative_RRHO_G_barrier_kcal_mol | -0.939 | -0.800 | 2.693 |
| opposite_mean_abs_projected_force_mEh_bohr | relative_electronic_barrier_kcal_mol | -0.881 | -0.200 | 6.981 |
| opposite_mean_abs_projected_force_mEh_bohr | relative_RRHO_G_barrier_kcal_mol | -0.927 | -0.200 | 7.866 |
| near_mean_signed_projected_force_mEh_bohr | relative_electronic_barrier_kcal_mol | 0.831 | 0.800 | 3.136 |
| near_mean_signed_projected_force_mEh_bohr | relative_RRHO_G_barrier_kcal_mol | 0.776 | 0.800 | 5.347 |
| opposite_mean_signed_projected_force_mEh_bohr | relative_electronic_barrier_kcal_mol | 0.873 | 0.200 | 6.377 |
| opposite_mean_signed_projected_force_mEh_bohr | relative_RRHO_G_barrier_kcal_mol | 0.920 | 0.200 | 7.261 |
| opposite_TS_distance_A | relative_electronic_barrier_kcal_mol | 0.439 | 0.400 | 3.493 |
| opposite_TS_distance_A | relative_RRHO_G_barrier_kcal_mol | 0.533 | 0.400 | 4.715 |
| opposite_TS_CXO_angle_deg | relative_electronic_barrier_kcal_mol | -0.475 | -0.200 | 3.571 |
| opposite_TS_CXO_angle_deg | relative_RRHO_G_barrier_kcal_mol | -0.568 | -0.200 | 4.782 |
| near_TS_distance_A | relative_electronic_barrier_kcal_mol | -0.937 | -0.400 | 3.717 |
| near_TS_distance_A | relative_RRHO_G_barrier_kcal_mol | -0.970 | -0.400 | 3.108 |

高相関の候補 `B_atomic` について、電子障壁とのrは全感度条件で **0.760〜0.993**。順位は維持されるが、区間によって相関の強さは変わる。

| window | H_choice | target | pearson_r | spearman_r | LOO_prediction_RMSE_kcal_mol |
| --- | --- | --- | --- | --- | --- |
| common_CN_2.50_to_2.10 | primary | relative_electronic_barrier_kcal_mol | 0.958 | 1.000 | 1.410 |
| common_CN_2.50_to_2.10 | primary | relative_RRHO_G_barrier_kcal_mol | 0.922 | 1.000 | 2.761 |
| sensitivity_CN_2.55_to_2.10 | primary | relative_electronic_barrier_kcal_mol | 0.985 | 1.000 | 0.842 |
| sensitivity_CN_2.55_to_2.10 | primary | relative_RRHO_G_barrier_kcal_mol | 0.966 | 1.000 | 1.935 |
| sensitivity_CN_2.40_to_2.15 | primary | relative_electronic_barrier_kcal_mol | 0.760 | 1.000 | 2.857 |
| sensitivity_CN_2.40_to_2.15 | primary | relative_RRHO_G_barrier_kcal_mol | 0.685 | 1.000 | 4.639 |
| common_CN_2.50_to_2.10 | alternative_amine_H | relative_electronic_barrier_kcal_mol | 0.993 | 1.000 | 0.528 |
| common_CN_2.50_to_2.10 | alternative_amine_H | relative_RRHO_G_barrier_kcal_mol | 0.983 | 1.000 | 1.523 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | relative_electronic_barrier_kcal_mol | 0.982 | 1.000 | 0.632 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | relative_RRHO_G_barrier_kcal_mol | 0.989 | 1.000 | 1.191 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | relative_electronic_barrier_kcal_mol | 0.866 | 1.000 | 2.351 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | relative_RRHO_G_barrier_kcal_mol | 0.807 | 1.000 | 4.023 |

![感度](../figures/16_two_sided_sensitivity.png)

主区間のB_atomicを使った1系除外の線形予測では、電子障壁のRMSEは1.41 kcal/mol、Gibbs障壁では2.76 kcal/mol。4系しかなく、指標探索に同じ4系を使っているため、この交差検証も独立検証ではない。複数指標の中から高い相関が見つかった探索結果であり、「十分に予測できる」「機構を証明した」とはまだ判断しない。p値による有意性は主張しない。

## 次の検証

この候補が実用的かは、追加置換系や独立した配座で、今回固定した定義と区間のまま予測できるかを確認するのが重要。機構については、共通C–N距離での活性化ひずみ／相互作用エネルギー分解を追加すれば、構造変形と分子間相互作用の寄与を切り分けられる。ただし通常の分子間分解だけでは分子内の反対側H/F···Oを一意に分離できず、適切な配座対照や電子密度解析も必要になる。これらの追加量子化学計算は今回行っていない。

## 再現・参照

`python analysis4_force/two_sided_contacts.py` で再生成。設定は `config/two_sided_contacts.yaml`、全結果は `two_sided_*.csv`。通常の4系解析を実行した後の座標・力を使用する。

- [原子寄与からのReaction Force解析（2012）](https://pubmed.ncbi.nlm.nih.gov/22746292/)：原子の反応進行への寄与を考える背景。今回のペア射影比・原子力割合はこのプロジェクト独自の探索指標。
- [Distortion/Interaction–Activation Strainモデル（2017）](https://pubmed.ncbi.nlm.nih.gov/28447369/)：障壁の構造変形と相互作用への分解を考える背景。
