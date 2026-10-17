# 既存logによる両側接触・原子仕事と活性化障壁

## 結論と今回の変更

追加のDFT計算は行わず、保存済みの座標・Gaussian原子力・エネルギーだけを使用した。図15の横軸は **各側の原子仕事（kcal/mol）**、縦軸は **各系の分離反応物を基準とするΔE‡・ΔG‡** に変更した。MoriPFとの差やpercentage pointsは使わない。元の比率の相関は補足図17に残し、相互作用の強さとは解釈しない。

文献調査に基づく選択理由・今後のlog限定の計算案は [調査ノート](log_only_research_plan_ja.md) を参照。今回の原子仕事は「どの原子の移動にエネルギー上昇が配分されるか」という間接指標。F···HやH/F···Oの結合エネルギーではないが、単位・積分範囲・全体との加算関係を明示できる。

共通区間はC–N 2.50→2.10 Å。表の全数値はkcal/mol：

| 系 | アミン側原子仕事 | 反対側原子仕事 | ΔE‡ | ΔG‡（RRHO） |
| --- | --- | --- | --- | --- |
| MoriPF | 0.085 | 0.043 | -1.342 | 14.083 |
| RAFA14 | 0.098 | 0.027 | 0.476 | 15.928 |
| RAFA15 | 0.141 | 0.140 | -4.587 | 8.943 |
| RAFA16 | 0.096 | 0.054 | -0.579 | 15.023 |

![原子仕事と障壁](../figures/15_two_sided_barrier_correlations.png)

| descriptor | target | pearson_r | spearman_r | LOO_prediction_RMSE_kcal_mol |
| --- | --- | --- | --- | --- |
| near_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.837 | -0.200 | 5.044 |
| near_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.886 | -0.200 | 6.230 |
| opposite_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965 | -0.800 | 0.746 |
| opposite_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.976 | -0.800 | 1.437 |
| near_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.185 | 0.200 | 3.795 |
| near_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.286 | 0.200 | 5.284 |
| opposite_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.967 | -0.800 | 1.220 |
| opposite_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.943 | -0.800 | 2.559 |

反対側原子仕事は電子障壁とr≈−0.965、Gibbs障壁とr≈−0.976。RAFA15の仕事が大きいほど障壁が低いという**逆向きの対応**であり、「反発コストが低いので速い」とは読めない。RAFA15では反対側仕事の約98%がO21由来で、H56単独の仕事は小さい。

相関の良さは指標を物理的な相互作用量に変えない。両側を別々に示し、相関が弱い側も残す。選択4原子の合計仕事、O21単独、残りの全原子も保存し、特定の接触ではなく反応中心全体の変化で説明できないかを確認する。

## 原子の対応

| system | near_ring_C | near_atom | near_element | primary_amine_H | alternative_amine_H | opposite_ring_C | opposite_atom | opposite_element | opposite_O |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MoriPF | 24 | 33 | F | 46 | 47 | 25 | 29 | F | 21 |
| RAFA14 | 24 | 57 | H | 44 | 45 | 25 | 56 | H | 21 |
| RAFA15 | 24 | 31 | F | 44 | 45 | 25 | 56 | H | 21 |
| RAFA16 | 24 | 31 | F | 44 | 45 | 25 | 29 | F | 21 |

C24上のX=F/Hと選択アミンHをアミン側、C25上のX=F/HとO21を反対側とする。RAFA15の反対側はベンゼン環H56···O21。RAFA14のアミン側はH57···H44でありF···Hではない。アミンHはTSでC24置換原子に近い方を固定し、もう一方も感度解析する。

| system | near_TS_distance_A | opposite_TS_distance_A | opposite_TS_CXO_angle_deg |
| --- | --- | --- | --- |
| MoriPF | 2.358 | 2.973 | 77.568 |
| RAFA14 | 2.358 | 2.807 | 87.127 |
| RAFA15 | 2.622 | 2.764 | 90.436 |
| RAFA16 | 2.345 | 2.991 | 78.380 |

RAFA15のC25–H56···O21角は約90.4°で、直線的なC–H···O配置ではない。近接・相互作用の仮説として扱い、水素結合と断定しない。

## 仕事の定義と解釈

```text
C_i(R) = F_i(R) · dr_i/dR
W_g = −∫[R(c=2.50) → R(c=2.10)] Σ(i∈g) C_i(R) dR
W_near + W_opposite + W_remaining = W_all ≈ E(end) − E(start)
```

座標はBohr、力はHartree/Bohrとして積分し、kcal/molに変換した。全57原子をTSへ剛体重ね合わせし、力ベクトルも同じ回転を施した。原子別配分はこの座標系・経路に依存し、独立した原子の状態関数ではない。全系の和と電子エネルギー差の一致を数値検証する。既存区間の仕事であり、分離反応物→TSの全障壁を分解したものではない。ΔG‡との比較も電子仕事の相関であって自由エネルギー分解ではない。

| system | near_X_atomic_work_kcal_mol | amine_H_atomic_work_kcal_mol | opposite_X_atomic_work_kcal_mol | opposite_O_atomic_work_kcal_mol | remaining_atomic_work_kcal_mol | all_atomic_work_kcal_mol |
| --- | --- | --- | --- | --- | --- | --- |
| MoriPF | 0.016 | 0.069 | 0.013 | 0.030 | 2.238 | 2.366 |
| RAFA14 | 0.004 | 0.095 | 0.003 | 0.024 | 3.020 | 3.145 |
| RAFA15 | 0.043 | 0.098 | 0.003 | 0.137 | 3.749 | 4.029 |
| RAFA16 | 0.021 | 0.075 | 0.019 | 0.036 | 2.477 | 2.627 |

![原子別・全体の仕事](../figures/19_atomic_work_partition.png)

反応物側のIRCはエネルギー上昇方向なので、W>0を「その接触が反発している」と解釈しない。他の原子・結合・電子の応答が各原子力に含まれる。O21が主成分なら「反対側Oを含む構造変化のコスト」と表現する。全体との和が合うことは数値検証であり、H···Oの因果証明ではない。

## 局所仕事と全障壁を混同しないためのエネルギー収支

```text
ΔE‡ = [E(c=2.50) − E_separated]
      + [E(c=2.10) − E(c=2.50)] + [E_TS − E(c=2.10)]
```

| system | entry_energy_kcal_mol | window_energy_rise_kcal_mol | after_window_to_TS_kcal_mol | absolute_electronic_barrier_kcal_mol | work_energy_residual_kcal_mol |
| --- | --- | --- | --- | --- | --- |
| MoriPF | -3.791 | 2.367 | 0.082 | -1.342 | -0.001 |
| RAFA14 | -2.913 | 3.148 | 0.242 | 0.476 | -0.003 |
| RAFA15 | -9.019 | 4.031 | 0.401 | -4.587 | -0.002 |
| RAFA16 | -3.326 | 2.628 | 0.120 | -0.579 | -0.001 |

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

| system | near_X_mean_force_mEh_bohr | amine_H_mean_force_mEh_bohr | opposite_X_mean_force_mEh_bohr | opposite_O_mean_force_mEh_bohr | near_mean_signed_projected_force_mEh_bohr | opposite_mean_signed_projected_force_mEh_bohr |
| --- | --- | --- | --- | --- | --- | --- |
| MoriPF | 0.370 | 0.239 | 0.331 | 0.518 | -0.046 | -0.111 |
| RAFA14 | 0.058 | 0.315 | 0.048 | 0.510 | 0.039 | -0.182 |
| RAFA15 | 0.768 | 0.380 | 0.050 | 1.477 | -0.097 | -0.514 |
| RAFA16 | 0.460 | 0.262 | 0.438 | 0.601 | -0.060 | -0.120 |

## 元の比率は補足に保持

```text
B_atomic = ∫ (|F_Xopp| + |F_O|) dc
           / ∫ (|F_Xnear| + |F_HN| + |F_Xopp| + |F_O|) dc
```

割合そのものを無次元値として表示する。電子障壁とのr≈0.958という相関はMoriPFを引かなくても変わらない。定数の差し引きは相関を変えないため、以前の相関がゼロ合わせで生まれたわけではない。ただし、分母・原子選択に依存し、物理的な解釈に乏しいという問題は残る。

RAFA15の反対側原子力の大きさは約97%がO21に由来する。この比率をH···O相互作用の強さと読むことはできない。

![元の無次元比率](../figures/17_dimensionless_share_supplement.png)

## 全指標と感度

元の12指標に、今回7つの原子仕事成分と両側の射影仕事2成分を追加した。新指標はユーザーのlog限定・有次元指標の要望を受けた探索的追加であり、事前登録や独立検証ではない。全21指標×3区間×2通りのアミンH×2種類の障壁を保存する。区間を変えて最も高い相関だけを選ばない。

| descriptor | target | pearson_r | spearman_r | LOO_prediction_RMSE_kcal_mol |
| --- | --- | --- | --- | --- |
| opposite_projected_magnitude_share | absolute_electronic_barrier_kcal_mol | -0.393 | -0.400 | 3.427 |
| opposite_projected_magnitude_share | absolute_RRHO_G_barrier_kcal_mol | -0.485 | -0.400 | 4.679 |
| opposite_atomic_force_share | absolute_electronic_barrier_kcal_mol | 0.958 | 1.000 | 1.410 |
| opposite_atomic_force_share | absolute_RRHO_G_barrier_kcal_mol | 0.922 | 1.000 | 2.761 |
| selected_atomic_force_share | absolute_electronic_barrier_kcal_mol | -0.632 | -0.400 | 3.158 |
| selected_atomic_force_share | absolute_RRHO_G_barrier_kcal_mol | -0.553 | -0.400 | 5.625 |
| selected_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.937 | -0.800 | 1.648 |
| selected_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.961 | -0.800 | 2.990 |
| pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.807 | -0.400 | 3.981 |
| pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.866 | -0.400 | 4.882 |
| near_mean_abs_projected_force_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.932 | -0.800 | 1.559 |
| near_mean_abs_projected_force_mEh_bohr | absolute_RRHO_G_barrier_kcal_mol | -0.939 | -0.800 | 2.693 |
| opposite_mean_abs_projected_force_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.881 | -0.200 | 6.981 |
| opposite_mean_abs_projected_force_mEh_bohr | absolute_RRHO_G_barrier_kcal_mol | -0.927 | -0.200 | 7.866 |
| near_mean_signed_projected_force_mEh_bohr | absolute_electronic_barrier_kcal_mol | 0.831 | 0.800 | 3.136 |
| near_mean_signed_projected_force_mEh_bohr | absolute_RRHO_G_barrier_kcal_mol | 0.776 | 0.800 | 5.347 |
| opposite_mean_signed_projected_force_mEh_bohr | absolute_electronic_barrier_kcal_mol | 0.873 | 0.200 | 6.377 |
| opposite_mean_signed_projected_force_mEh_bohr | absolute_RRHO_G_barrier_kcal_mol | 0.920 | 0.200 | 7.261 |
| opposite_TS_distance_A | absolute_electronic_barrier_kcal_mol | 0.439 | 0.400 | 3.493 |
| opposite_TS_distance_A | absolute_RRHO_G_barrier_kcal_mol | 0.533 | 0.400 | 4.715 |
| opposite_TS_CXO_angle_deg | absolute_electronic_barrier_kcal_mol | -0.475 | -0.200 | 3.571 |
| opposite_TS_CXO_angle_deg | absolute_RRHO_G_barrier_kcal_mol | -0.568 | -0.200 | 4.782 |
| near_TS_distance_A | absolute_electronic_barrier_kcal_mol | -0.937 | -0.400 | 3.717 |
| near_TS_distance_A | absolute_RRHO_G_barrier_kcal_mol | -0.970 | -0.400 | 3.108 |
| near_X_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.957 | -0.800 | 1.167 |
| near_X_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.945 | -0.800 | 2.380 |
| amine_H_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.349 | -0.200 | 3.595 |
| amine_H_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.448 | -0.200 | 4.912 |
| opposite_X_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.332 | 0.400 | 3.400 |
| opposite_X_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | 0.421 | 0.400 | 4.683 |
| opposite_O_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.954 | -0.800 | 2.980 |
| opposite_O_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.977 | -0.800 | 1.469 |
| near_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.837 | -0.200 | 5.044 |
| near_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.886 | -0.200 | 6.230 |
| opposite_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965 | -0.800 | 0.746 |
| opposite_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.976 | -0.800 | 1.437 |
| remaining_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.656 | -0.200 | 3.946 |
| remaining_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.733 | -0.200 | 5.094 |
| near_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.185 | 0.200 | 3.795 |
| near_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.286 | 0.200 | 5.284 |
| opposite_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.967 | -0.800 | 1.220 |
| opposite_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.943 | -0.800 | 2.559 |

![感度](../figures/16_two_sided_sensitivity.png)

| window | H_choice | descriptor | target | pearson_r | spearman_r |
| --- | --- | --- | --- | --- | --- |
| common_CN_2.50_to_2.10 | primary | opposite_atomic_force_share | absolute_electronic_barrier_kcal_mol | 0.958 | 1.000 |
| common_CN_2.50_to_2.10 | primary | opposite_atomic_force_share | absolute_RRHO_G_barrier_kcal_mol | 0.922 | 1.000 |
| common_CN_2.50_to_2.10 | primary | near_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.837 | -0.200 |
| common_CN_2.50_to_2.10 | primary | near_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.886 | -0.200 |
| common_CN_2.50_to_2.10 | primary | opposite_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965 | -0.800 |
| common_CN_2.50_to_2.10 | primary | opposite_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.976 | -0.800 |
| common_CN_2.50_to_2.10 | primary | near_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.185 | 0.200 |
| common_CN_2.50_to_2.10 | primary | near_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.286 | 0.200 |
| common_CN_2.50_to_2.10 | primary | opposite_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.967 | -0.800 |
| common_CN_2.50_to_2.10 | primary | opposite_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.943 | -0.800 |
| sensitivity_CN_2.55_to_2.10 | primary | opposite_atomic_force_share | absolute_electronic_barrier_kcal_mol | 0.985 | 1.000 |
| sensitivity_CN_2.55_to_2.10 | primary | opposite_atomic_force_share | absolute_RRHO_G_barrier_kcal_mol | 0.966 | 1.000 |
| sensitivity_CN_2.55_to_2.10 | primary | near_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.865 | -0.400 |
| sensitivity_CN_2.55_to_2.10 | primary | near_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.908 | -0.400 |
| sensitivity_CN_2.55_to_2.10 | primary | opposite_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.967 | -0.800 |
| sensitivity_CN_2.55_to_2.10 | primary | opposite_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.974 | -0.800 |
| sensitivity_CN_2.55_to_2.10 | primary | near_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.182 | 0.200 |
| sensitivity_CN_2.55_to_2.10 | primary | near_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.283 | 0.200 |
| sensitivity_CN_2.55_to_2.10 | primary | opposite_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.960 | -0.800 |
| sensitivity_CN_2.55_to_2.10 | primary | opposite_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.933 | -0.800 |
| sensitivity_CN_2.40_to_2.15 | primary | opposite_atomic_force_share | absolute_electronic_barrier_kcal_mol | 0.760 | 1.000 |
| sensitivity_CN_2.40_to_2.15 | primary | opposite_atomic_force_share | absolute_RRHO_G_barrier_kcal_mol | 0.685 | 1.000 |
| sensitivity_CN_2.40_to_2.15 | primary | near_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.805 | -0.200 |
| sensitivity_CN_2.40_to_2.15 | primary | near_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.862 | -0.200 |
| sensitivity_CN_2.40_to_2.15 | primary | opposite_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.959 | -0.800 |
| sensitivity_CN_2.40_to_2.15 | primary | opposite_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.976 | -0.800 |
| sensitivity_CN_2.40_to_2.15 | primary | near_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.166 | 0.200 |
| sensitivity_CN_2.40_to_2.15 | primary | near_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.268 | 0.200 |
| sensitivity_CN_2.40_to_2.15 | primary | opposite_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.969 | -0.800 |
| sensitivity_CN_2.40_to_2.15 | primary | opposite_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.947 | -0.800 |
| common_CN_2.50_to_2.10 | alternative_amine_H | opposite_atomic_force_share | absolute_electronic_barrier_kcal_mol | 0.993 | 1.000 |
| common_CN_2.50_to_2.10 | alternative_amine_H | opposite_atomic_force_share | absolute_RRHO_G_barrier_kcal_mol | 0.983 | 1.000 |
| common_CN_2.50_to_2.10 | alternative_amine_H | near_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.915 | -0.400 |
| common_CN_2.50_to_2.10 | alternative_amine_H | near_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.946 | -0.400 |
| common_CN_2.50_to_2.10 | alternative_amine_H | opposite_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965 | -0.800 |
| common_CN_2.50_to_2.10 | alternative_amine_H | opposite_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.976 | -0.800 |
| common_CN_2.50_to_2.10 | alternative_amine_H | near_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.214 | 0.200 |
| common_CN_2.50_to_2.10 | alternative_amine_H | near_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | 0.115 | 0.200 |
| common_CN_2.50_to_2.10 | alternative_amine_H | opposite_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.967 | -0.800 |
| common_CN_2.50_to_2.10 | alternative_amine_H | opposite_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.943 | -0.800 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | opposite_atomic_force_share | absolute_electronic_barrier_kcal_mol | 0.982 | 1.000 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | opposite_atomic_force_share | absolute_RRHO_G_barrier_kcal_mol | 0.989 | 1.000 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | near_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.925 | -0.400 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | near_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.953 | -0.400 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | opposite_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.967 | -0.800 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | opposite_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.974 | -0.800 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | near_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.122 | 0.200 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | near_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | 0.020 | 0.200 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | opposite_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.960 | -0.800 |
| sensitivity_CN_2.55_to_2.10 | alternative_amine_H | opposite_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.933 | -0.800 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | opposite_atomic_force_share | absolute_electronic_barrier_kcal_mol | 0.866 | 1.000 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | opposite_atomic_force_share | absolute_RRHO_G_barrier_kcal_mol | 0.807 | 1.000 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | near_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.896 | -0.400 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | near_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.932 | -0.400 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | opposite_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.959 | -0.800 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | opposite_atomic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.976 | -0.800 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | near_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.387 | 0.200 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | near_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | 0.294 | 0.200 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | opposite_pair_diagnostic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.969 | -0.800 |
| sensitivity_CN_2.40_to_2.15 | alternative_amine_H | opposite_pair_diagnostic_work_kcal_mol | absolute_RRHO_G_barrier_kcal_mol | -0.947 | -0.800 |

4系しかないため、相関係数・順位と1系除外予測は記述的診断。複数指標を同じ4系で探索した交差検証は独立検証ではない。比率や線形結合の係数を障壁に合わせて最適化しない。RAFA15のRRHO Gibbs障壁は反応物の微小虚振動のため暫定。

## 再現・文献

`python analysis4_force/run_analysis.py` で全結果、`python analysis4_force/two_sided_contacts.py` でこの追加解析を再生成できる。TSの力欠測を補間せず、全系で実測力がある共通区間のみ使用した。

- Haas & Chu (2009), [Decomposition of energy and free energy changes by following the flow of work along reaction path](https://doi.org/10.1063/1.3243080)：原子力と移動による仕事配分の根拠。本解析で自由エネルギー経路を計算したという意味ではない。
- Vöhringer-Martinez & Toro-Labbé (2012), [Understanding the physics and chemistry of reaction mechanisms from atomic contributions: a reaction force perspective](https://doi.org/10.1021/jp303075k)：原子別の反応力・エネルギー寄与の考え方。
