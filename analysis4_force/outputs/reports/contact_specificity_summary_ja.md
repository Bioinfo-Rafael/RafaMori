# 分子内順位と、他の原子・接触を対照にした特異性検証

追加DFTなし。2026-10-17。全57原子と全1596原子対を計算し、幾何・結合関係で対照を選んだ。障壁は各系のΔE‡・ΔG‡そのもの。MoriPFとの差やpercentage pointsは用いない。

## 結論

**反対側H/F···Oは、近い種類・距離の接触の中で応答が大きい。ただし、障壁に相関する接触がこれだけという結果ではない。アミン側は「大きさ」と「相関」を同時に満たす根拠が弱い。**

元の図15で示した反対側2原子の仕事和（電子障壁r≈−0.965）は、同じ分類で4系に共通する対照96個中、18個が同等以上の|r|を持つ。H/F···O/Nかつ平均距離差±0.75 Åに絞っても、対照6個中1個が上回る。前の高相関だけでは特異性を示せなかった。

接触方向の射影仕事（電子障壁r≈−0.967）では、同分類の対照96個中34個が同等以上の|r|を持つ。化学種・距離を絞っても対照6個中4個。したがってこの指標でも「ここだけ相関する」とはいえない。

一方、**RAFA15の反対側H56···O21は、化学種・距離を近づけた13候補の中で射影仕事・平均射影力とも1位**。同じ対照定義でMoriPFは2位、RAFA14は3位、RAFA16は2位であり、局所応答としては目立つ。2原子の仕事和では4系とも4/13位。指定接触を研究対象にする根拠は残るが、因果や唯一性は未確定。

| system | pool | value | magnitude_rank | n_controls |
| --- | --- | --- | --- | --- |
| MoriPF | same_category | 0.00283647 | 12 | 151 |
| MoriPF | distance_and_chemistry | 0.00283647 | 2 | 12 |
| RAFA14 | same_category | -0.00162221 | 17 | 151 |
| RAFA14 | distance_and_chemistry | -0.00162221 | 3 | 12 |
| RAFA15 | same_category | 0.00908524 | 5 | 168 |
| RAFA15 | distance_and_chemistry | 0.00908524 | 1 | 12 |
| RAFA16 | same_category | 0.00311402 | 12 | 152 |
| RAFA16 | distance_and_chemistry | 0.00311402 | 2 | 12 |

「強さと相関の両方で全4系の対象を上回る」対照は、反対側の2原子仕事和には存在しなかった。射影仕事には広い対照でC20···C24が1件あり（r≈−0.995）、H/F···O/Nに限定すると0件。これは指標・対照の選択により結論が変わることを示す。C20···C24は1–4のC···Cなので、同じ化学相互作用とは主張せず、1–4対を除いた感度比較も残した。

## 今回の問い

「指定した接触に相関がある」だけでなく、(1) 同じ系の他の原子・接触より応答が大きいか、(2) 他の接触には同様の障壁相関がないか、(3) 指定原子を共有しない対照でも同じ結果が出ないか、を確認する。

全原子や共有結合より最大であることは機構上の必要条件ではない。一方、同じ定義・近い距離・同じ分子内/分子間分類の対照でも平凡なら、「この接触だけが特別に強い」という説明は支持されない。本解析は全原子力からの間接量の比較であり、真のペア結合強度の順位ではない。

## 対照の固定ルール

- 反応C–N距離の3区間は従来と同じ。主区間2.50→2.10 Å、感度区間2.55→2.10、2.40→2.15 Å。
- 最も反応物側の構造で共有結合グラフを推定。距離が共有結合半径和×1.2未満を結合とし、反応で形成中のC–Nは元の分子間の対として扱う。
- 直接結合(1–2)・結合角(1–3)の原子対は対照から除外。グラフ距離3以上（1–4対は含む）、区間内最短距離4.0 Å以下を非共有接触候補とする。これは操作的な接触定義で、相互作用の有無の判定ではない。
- アミン側指定接触は分子間、反対側は基質分子内。同じ分類内で比較する。全原子は別に57原子すべてを比較する。
- `same_category`：同じ分子内/分子間分類の全候補。`distance_matched`：さらに平均距離差±0.75 Å。`shared_anchor`：反対側はO21、アミン側は選択NHを共有。`disjoint`：指定接触の2原子をどちらも含まない対照。
- `chemically_related`：反対側はH/F···O/N、アミン側はH···H/F/O/Nに限定。`distance_and_chemistry`はさらに距離条件も課す。RAFA14の指定アミン側はH···Hなので、この対照からH···Hを除外しない。
- `exclude_1_4`：グラフ距離3の1–4対まで除外し、4結合以上離れた接触だけを比較する感度確認。
- 元素の種類・供与体角度・環内外は一律に同じではない。これは広い反証用対照で、精密な化学的マッチングではない。H/F置換のある同じ環位置を比較するため、元素一致を必須にして対照を空にすることは避けた。
- 同一分子内の強さ順位は対称Hも含む全原子対。系間の相関は番号だけでは対応させず、共通骨格、環上の置換位置、求核剤の対応位置を使用する。アルキル基上の複数Hは一意な化学的対応を仮定せず系間相関から除外した。原子対応と共有結合グラフをCSVで公開する。
- 系間比較では同じ対が**全4系で**接触条件を満たすものだけを採用。距離マッチも全4系で満たす。障壁値から対照を選ばない。

## 全原子内で指定原子はどこに位置するか

仕事の符号を保って保存し、大きさの順位は絶対値で比較。rank=1が最大、対照数56に対象1原子を加えて全57原子。

| system | target | value | magnitude_rank | n_controls |
| --- | --- | --- | --- | --- |
| MoriPF | X24 | 0.0156109 | 18 | 56 |
| MoriPF | NH_near | 0.0690206 | 8 | 56 |
| MoriPF | X25 | 0.012798 | 20 | 56 |
| MoriPF | core21 | 0.0298387 | 12 | 56 |
| RAFA14 | X24 | 0.00381192 | 33 | 56 |
| RAFA14 | NH_near | 0.0945811 | 8 | 56 |
| RAFA14 | X25 | 0.00314581 | 34 | 56 |
| RAFA14 | core21 | 0.0239114 | 17 | 56 |
| RAFA15 | X24 | 0.0427993 | 12 | 56 |
| RAFA15 | NH_near | 0.0978649 | 10 | 56 |
| RAFA15 | X25 | 0.00290557 | 30 | 56 |
| RAFA15 | core21 | 0.137244 | 9 | 56 |
| RAFA16 | X24 | 0.0207225 | 16 | 56 |
| RAFA16 | NH_near | 0.0749908 | 8 | 56 |
| RAFA16 | X25 | 0.0187642 | 18 | 56 |
| RAFA16 | core21 | 0.0356022 | 12 | 56 |

![全原子の順位](../figures/20_all_atom_work_controls.png)

## 同じ系の他の接触に比べて強いか

`radial_work_kcal_mol`は接触方向の仕事診断、`mean_abs_projection_mEh_bohr`は接触方向の力の平均絶対値。強さのrankは大きい順、比較総数はn_controls+1。少数の順位差を精密な結合強度差とは読まない。

| system | target | pool | metric | value | magnitude_rank | n_controls |
| --- | --- | --- | --- | --- | --- | --- |
| MoriPF | opposite_contact | same_category | pair_atomic_work_kcal_mol | 0.0426367 | 26 | 151 |
| MoriPF | opposite_contact | same_category | radial_work_kcal_mol | 0.00283647 | 12 | 151 |
| MoriPF | opposite_contact | same_category | mean_abs_projection_mEh_bohr | 0.119227 | 14 | 151 |
| MoriPF | opposite_contact | distance_and_chemistry | pair_atomic_work_kcal_mol | 0.0426367 | 4 | 12 |
| MoriPF | opposite_contact | distance_and_chemistry | radial_work_kcal_mol | 0.00283647 | 2 | 12 |
| MoriPF | opposite_contact | distance_and_chemistry | mean_abs_projection_mEh_bohr | 0.119227 | 2 | 12 |
| MoriPF | amine_contact | same_category | pair_atomic_work_kcal_mol | 0.0846315 | 39 | 53 |
| MoriPF | amine_contact | same_category | radial_work_kcal_mol | -0.00490808 | 44 | 53 |
| MoriPF | amine_contact | same_category | mean_abs_projection_mEh_bohr | 0.0456488 | 48 | 53 |
| MoriPF | amine_contact | distance_and_chemistry | pair_atomic_work_kcal_mol | 0.0846315 | 3 | 7 |
| MoriPF | amine_contact | distance_and_chemistry | radial_work_kcal_mol | -0.00490808 | 6 | 7 |
| MoriPF | amine_contact | distance_and_chemistry | mean_abs_projection_mEh_bohr | 0.0456488 | 6 | 7 |
| RAFA14 | opposite_contact | same_category | pair_atomic_work_kcal_mol | 0.0270572 | 42 | 151 |
| RAFA14 | opposite_contact | same_category | radial_work_kcal_mol | -0.00162221 | 17 | 151 |
| RAFA14 | opposite_contact | same_category | mean_abs_projection_mEh_bohr | 0.182389 | 14 | 151 |
| RAFA14 | opposite_contact | distance_and_chemistry | pair_atomic_work_kcal_mol | 0.0270572 | 4 | 12 |
| RAFA14 | opposite_contact | distance_and_chemistry | radial_work_kcal_mol | -0.00162221 | 3 | 12 |
| RAFA14 | opposite_contact | distance_and_chemistry | mean_abs_projection_mEh_bohr | 0.182389 | 3 | 12 |
| RAFA14 | amine_contact | same_category | pair_atomic_work_kcal_mol | 0.098393 | 42 | 51 |
| RAFA14 | amine_contact | same_category | radial_work_kcal_mol | 0.00299599 | 44 | 51 |
| RAFA14 | amine_contact | same_category | mean_abs_projection_mEh_bohr | 0.0385756 | 48 | 51 |
| RAFA14 | amine_contact | distance_and_chemistry | pair_atomic_work_kcal_mol | 0.098393 | 4 | 7 |
| RAFA14 | amine_contact | distance_and_chemistry | radial_work_kcal_mol | 0.00299599 | 5 | 7 |
| RAFA14 | amine_contact | distance_and_chemistry | mean_abs_projection_mEh_bohr | 0.0385756 | 7 | 7 |
| RAFA15 | opposite_contact | same_category | pair_atomic_work_kcal_mol | 0.14015 | 26 | 168 |
| RAFA15 | opposite_contact | same_category | radial_work_kcal_mol | 0.00908524 | 5 | 168 |
| RAFA15 | opposite_contact | same_category | mean_abs_projection_mEh_bohr | 0.513777 | 3 | 168 |
| RAFA15 | opposite_contact | distance_and_chemistry | pair_atomic_work_kcal_mol | 0.14015 | 4 | 12 |
| RAFA15 | opposite_contact | distance_and_chemistry | radial_work_kcal_mol | 0.00908524 | 1 | 12 |
| RAFA15 | opposite_contact | distance_and_chemistry | mean_abs_projection_mEh_bohr | 0.513777 | 1 | 12 |
| RAFA15 | amine_contact | same_category | pair_atomic_work_kcal_mol | 0.140664 | 39 | 48 |
| RAFA15 | amine_contact | same_category | radial_work_kcal_mol | 0.00152599 | 47 | 48 |
| RAFA15 | amine_contact | same_category | mean_abs_projection_mEh_bohr | 0.0967671 | 43 | 48 |
| RAFA15 | amine_contact | distance_and_chemistry | pair_atomic_work_kcal_mol | 0.140664 | 7 | 8 |
| RAFA15 | amine_contact | distance_and_chemistry | radial_work_kcal_mol | 0.00152599 | 9 | 8 |
| RAFA15 | amine_contact | distance_and_chemistry | mean_abs_projection_mEh_bohr | 0.0967671 | 7 | 8 |
| RAFA16 | opposite_contact | same_category | pair_atomic_work_kcal_mol | 0.0543665 | 25 | 152 |
| RAFA16 | opposite_contact | same_category | radial_work_kcal_mol | 0.00311402 | 12 | 152 |
| RAFA16 | opposite_contact | same_category | mean_abs_projection_mEh_bohr | 0.129816 | 16 | 152 |
| RAFA16 | opposite_contact | distance_and_chemistry | pair_atomic_work_kcal_mol | 0.0543665 | 4 | 12 |
| RAFA16 | opposite_contact | distance_and_chemistry | radial_work_kcal_mol | 0.00311402 | 2 | 12 |
| RAFA16 | opposite_contact | distance_and_chemistry | mean_abs_projection_mEh_bohr | 0.129816 | 2 | 12 |
| RAFA16 | amine_contact | same_category | pair_atomic_work_kcal_mol | 0.0957133 | 39 | 53 |
| RAFA16 | amine_contact | same_category | radial_work_kcal_mol | -0.00664926 | 43 | 53 |
| RAFA16 | amine_contact | same_category | mean_abs_projection_mEh_bohr | 0.0597148 | 48 | 53 |
| RAFA16 | amine_contact | distance_and_chemistry | pair_atomic_work_kcal_mol | 0.0957133 | 3 | 7 |
| RAFA16 | amine_contact | distance_and_chemistry | radial_work_kcal_mol | -0.00664926 | 6 | 7 |
| RAFA16 | amine_contact | distance_and_chemistry | mean_abs_projection_mEh_bohr | 0.0597148 | 6 | 7 |

![接触の強さ分布](../figures/21_contact_strength_controls.png)

## この接触だけが障壁に相関するか

主区間、電子障壁についての結果。相関は符号を保つ。対照との相関順位は|r|で比較し、正/負いずれも強い対応として数える。`controls_at_least_as_correlated`は対象以上の|r|を持つ対照数。`stronger_in_all_four`は全4系で対象以上の大きさを持つ対照数。最後の列は両条件を同時に満たす対照数。

| window | target | pool | metric | barrier | target_pearson_r | target_spearman_r | n_common_controls | controls_at_least_as_correlated | abs_correlation_rank | stronger_in_all_four | stronger_and_at_least_as_correlated |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| common_CN_2.50_to_2.10 | amine_contact | same_category | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.185312 | 0.2 | 28 | 28 | 29 | 23 | 23 |
| common_CN_2.50_to_2.10 | amine_contact | distance_matched | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.185312 | 0.2 | 11 | 11 | 12 | 10 | 10 |
| common_CN_2.50_to_2.10 | amine_contact | shared_anchor | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.185312 | 0.2 | 6 | 6 | 7 | 5 | 5 |
| common_CN_2.50_to_2.10 | amine_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.185312 | 0.2 | 20 | 20 | 21 | 18 | 18 |
| common_CN_2.50_to_2.10 | amine_contact | chemically_related | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.185312 | 0.2 | 5 | 5 | 6 | 4 | 4 |
| common_CN_2.50_to_2.10 | amine_contact | distance_and_chemistry | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.185312 | 0.2 | 1 | 1 | 2 | 1 | 1 |
| common_CN_2.50_to_2.10 | amine_contact | exclude_1_4 | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.185312 | 0.2 | 28 | 28 | 29 | 23 | 23 |
| common_CN_2.50_to_2.10 | opposite_contact | same_category | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.966567 | -0.8 | 96 | 34 | 35 | 4 | 1 |
| common_CN_2.50_to_2.10 | opposite_contact | distance_matched | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.966567 | -0.8 | 56 | 26 | 27 | 3 | 1 |
| common_CN_2.50_to_2.10 | opposite_contact | shared_anchor | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.966567 | -0.8 | 7 | 1 | 2 | 0 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.966567 | -0.8 | 84 | 33 | 34 | 3 | 1 |
| common_CN_2.50_to_2.10 | opposite_contact | chemically_related | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.966567 | -0.8 | 8 | 5 | 6 | 0 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | distance_and_chemistry | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.966567 | -0.8 | 6 | 4 | 5 | 0 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | exclude_1_4 | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.966567 | -0.8 | 21 | 9 | 10 | 2 | 0 |

元の図15の「2原子の仕事和」でも同じ対照比較を行った。これは指定した2原子の仕事の和であり、2原子間の力ではない：

| window | target | pool | metric | barrier | target_pearson_r | target_spearman_r | n_common_controls | controls_at_least_as_correlated | abs_correlation_rank | stronger_in_all_four | stronger_and_at_least_as_correlated |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| common_CN_2.50_to_2.10 | amine_contact | same_category | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.836934 | -0.2 | 28 | 3 | 4 | 26 | 3 |
| common_CN_2.50_to_2.10 | amine_contact | distance_matched | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.836934 | -0.2 | 11 | 0 | 1 | 11 | 0 |
| common_CN_2.50_to_2.10 | amine_contact | shared_anchor | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.836934 | -0.2 | 6 | 1 | 2 | 4 | 1 |
| common_CN_2.50_to_2.10 | amine_contact | disjoint | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.836934 | -0.2 | 20 | 1 | 2 | 20 | 1 |
| common_CN_2.50_to_2.10 | amine_contact | chemically_related | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.836934 | -0.2 | 5 | 3 | 4 | 5 | 3 |
| common_CN_2.50_to_2.10 | amine_contact | distance_and_chemistry | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.836934 | -0.2 | 1 | 0 | 1 | 1 | 0 |
| common_CN_2.50_to_2.10 | amine_contact | exclude_1_4 | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.836934 | -0.2 | 28 | 3 | 4 | 26 | 3 |
| common_CN_2.50_to_2.10 | opposite_contact | same_category | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965135 | -0.8 | 96 | 18 | 19 | 21 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | distance_matched | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965135 | -0.8 | 56 | 12 | 13 | 12 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | shared_anchor | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965135 | -0.8 | 7 | 0 | 1 | 3 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | disjoint | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965135 | -0.8 | 84 | 18 | 19 | 16 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | chemically_related | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965135 | -0.8 | 8 | 1 | 2 | 3 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | distance_and_chemistry | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965135 | -0.8 | 6 | 1 | 2 | 3 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | exclude_1_4 | pair_atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.965135 | -0.8 | 21 | 3 | 4 | 4 | 0 |

アミン側の平均射影力は電子障壁との|r|順位が同分類で1位だが、各分子内での大きさは下位である。相関順位だけで強い接触と呼べない例になる。

同じ比較を力の平均絶対値で行った結果：

| window | target | pool | metric | barrier | target_pearson_r | target_spearman_r | n_common_controls | controls_at_least_as_correlated | abs_correlation_rank | stronger_in_all_four | stronger_and_at_least_as_correlated |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| common_CN_2.50_to_2.10 | amine_contact | same_category | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.932259 | -0.8 | 28 | 0 | 1 | 25 | 0 |
| common_CN_2.50_to_2.10 | amine_contact | distance_matched | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.932259 | -0.8 | 11 | 0 | 1 | 11 | 0 |
| common_CN_2.50_to_2.10 | amine_contact | shared_anchor | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.932259 | -0.8 | 6 | 0 | 1 | 5 | 0 |
| common_CN_2.50_to_2.10 | amine_contact | disjoint | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.932259 | -0.8 | 20 | 0 | 1 | 19 | 0 |
| common_CN_2.50_to_2.10 | amine_contact | chemically_related | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.932259 | -0.8 | 5 | 0 | 1 | 4 | 0 |
| common_CN_2.50_to_2.10 | amine_contact | distance_and_chemistry | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.932259 | -0.8 | 1 | 0 | 1 | 1 | 0 |
| common_CN_2.50_to_2.10 | amine_contact | exclude_1_4 | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.932259 | -0.8 | 28 | 0 | 1 | 25 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | same_category | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.881125 | -0.2 | 96 | 44 | 45 | 2 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | distance_matched | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.881125 | -0.2 | 56 | 29 | 30 | 1 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | shared_anchor | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.881125 | -0.2 | 7 | 4 | 5 | 0 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | disjoint | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.881125 | -0.2 | 84 | 40 | 41 | 1 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | chemically_related | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.881125 | -0.2 | 8 | 4 | 5 | 0 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | distance_and_chemistry | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.881125 | -0.2 | 6 | 2 | 3 | 0 | 0 |
| common_CN_2.50_to_2.10 | opposite_contact | exclude_1_4 | mean_abs_projection_mEh_bohr | absolute_electronic_barrier_kcal_mol | -0.881125 | -0.2 | 21 | 12 | 13 | 2 | 0 |

![強さと相関の同時比較](../figures/22_contact_correlation_specificity.png)

### 指定した2原子をどちらも含まない対照の具体例

以下は相関の強い順に各5候補に加え、「全4系でより大きく、相関も同等以上」の候補を最大3つ表示した反例候補で、これらを新たな機構として採用したわけではない。全候補はCSVに残す。`coreNN`は共通骨格の原子番号NN、`X24`〜`X28`はその環C上のF/H、`NH_near`/`NH_other`はTSでアミン側Xに近い/遠いアミンH。各系の具体的な原子番号・元素はmappingを参照。

| window | target | pool | metric | barrier | control_pair | pearson_r | spearman_r | stronger_in_all_four | MoriPF_value | RAFA14_value | RAFA15_value | RAFA16_value | MoriPF_atoms | RAFA14_atoms | RAFA15_atoms | RAFA16_atoms |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| common_CN_2.50_to_2.10 | amine_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | amineN/core23 | 0.969268 | 0.8 | True | 0.34706 | 0.383787 | 0.178856 | 0.39148 | 23-45 | 23-43 | 23-43 | 23-43 |
| common_CN_2.50_to_2.10 | amine_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | amineC1/core19 | -0.937481 | -0.4 | False | 0.00430654 | 0.00529835 | 0.128989 | 0.00389499 | 19-48 | 19-46 | 19-46 | 19-46 |
| common_CN_2.50_to_2.10 | amine_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | amineN/core19 | -0.915912 | -0.4 | True | 0.0464877 | 0.0737347 | 0.422328 | 0.0459466 | 19-45 | 19-43 | 19-43 | 19-43 |
| common_CN_2.50_to_2.10 | amine_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | amineC1/core13 | -0.905195 | -0.2 | True | 0.0529354 | 0.0791875 | 0.307793 | 0.0576549 | 13-48 | 13-46 | 13-46 | 13-46 |
| common_CN_2.50_to_2.10 | amine_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | amineN/core12 | -0.901864 | -0.2 | True | 0.0685257 | 0.115021 | 0.491738 | 0.071698 | 12-45 | 12-43 | 12-43 | 12-43 |
| common_CN_2.50_to_2.10 | opposite_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | core19/core20 | 0.998345 | 1 | False | 0.00563499 | 0.0085164 | 0.00145398 | 0.00678192 | 19-20 | 19-20 | 19-20 | 19-20 |
| common_CN_2.50_to_2.10 | opposite_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | core11/core18 | 0.99833 | 1 | False | 0.000442432 | 0.00067749 | 4.29776e-05 | 0.000575353 | 11-18 | 11-18 | 11-18 | 11-18 |
| common_CN_2.50_to_2.10 | opposite_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | core08/core14 | -0.998204 | -1 | False | -0.000124524 | -0.000234582 | 6.11572e-05 | -0.000187279 | 8-14 | 8-14 | 8-14 | 8-14 |
| common_CN_2.50_to_2.10 | opposite_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | core10/core18 | 0.997924 | 1 | False | 0.00011401 | 0.000174294 | 2.65881e-05 | 0.000145723 | 10-18 | 10-18 | 10-18 | 10-18 |
| common_CN_2.50_to_2.10 | opposite_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | core07/core22 | -0.99778 | -1 | False | -0.00210884 | -0.00356314 | 8.49079e-05 | -0.00258647 | 7-22 | 7-22 | 7-22 | 7-22 |
| common_CN_2.50_to_2.10 | opposite_contact | disjoint | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | core20/core24 | -0.994598 | -1 | True | 0.00541022 | 0.00368614 | 0.00987585 | 0.00510686 | 20-24 | 20-24 | 20-24 | 20-24 |

### 他の単一原子も相関するか

電子障壁と原子仕事の|r|が高い10原子。原子番号の対応が一意でないアルキルHは含めない。

| kind | window | item | metric | barrier | pearson_r | spearman_r | leave_one_out_r_min | leave_one_out_r_max | median_magnitude |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| atom | common_CN_2.50_to_2.10 | core02 | atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.996859 | 1 | 0.985659 | 0.999517 | 0.00533346 |
| atom | common_CN_2.50_to_2.10 | core01 | atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.996258 | 1 | 0.997085 | 0.999777 | 0.00256298 |
| atom | common_CN_2.50_to_2.10 | core18 | atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.994978 | 1 | 0.993562 | 0.999997 | 0.0076048 |
| atom | common_CN_2.50_to_2.10 | core04 | atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.99461 | 1 | 0.984258 | 0.999298 | 0.00147852 |
| atom | common_CN_2.50_to_2.10 | core03 | atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.994135 | 1 | 0.980769 | 0.998863 | 0.00538645 |
| atom | common_CN_2.50_to_2.10 | core19 | atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.994096 | 1 | 0.980064 | 0.999991 | 0.0570321 |
| atom | common_CN_2.50_to_2.10 | core08 | atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.994091 | 1 | 0.995373 | 0.999597 | 0.000368373 |
| atom | common_CN_2.50_to_2.10 | core09 | atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.990826 | 1 | 0.975431 | 0.998649 | 0.000687417 |
| atom | common_CN_2.50_to_2.10 | sub_alkylC1 | atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | -0.977178 | -0.8 | -0.996501 | -0.78308 | 0.000314287 |
| atom | common_CN_2.50_to_2.10 | X28 | atomic_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 0.974304 | 1 | 0.944516 | 0.996199 | 0.00118799 |

## 区間感度と探索の影響

| window | target | barrier | target_pearson_r | n_common_controls | controls_at_least_as_correlated | stronger_and_at_least_as_correlated |
| --- | --- | --- | --- | --- | --- | --- |
| common_CN_2.50_to_2.10 | amine_contact | absolute_electronic_barrier_kcal_mol | -0.185312 | 28 | 28 | 23 |
| common_CN_2.50_to_2.10 | amine_contact | absolute_RRHO_G_barrier_kcal_mol | -0.286262 | 28 | 27 | 23 |
| common_CN_2.50_to_2.10 | opposite_contact | absolute_electronic_barrier_kcal_mol | -0.966567 | 96 | 34 | 1 |
| common_CN_2.50_to_2.10 | opposite_contact | absolute_RRHO_G_barrier_kcal_mol | -0.943007 | 96 | 45 | 2 |
| sensitivity_CN_2.55_to_2.10 | amine_contact | absolute_electronic_barrier_kcal_mol | -0.182013 | 28 | 26 | 22 |
| sensitivity_CN_2.55_to_2.10 | amine_contact | absolute_RRHO_G_barrier_kcal_mol | -0.283316 | 28 | 26 | 22 |
| sensitivity_CN_2.55_to_2.10 | opposite_contact | absolute_electronic_barrier_kcal_mol | -0.959718 | 96 | 36 | 2 |
| sensitivity_CN_2.55_to_2.10 | opposite_contact | absolute_RRHO_G_barrier_kcal_mol | -0.932885 | 96 | 48 | 2 |
| sensitivity_CN_2.40_to_2.15 | amine_contact | absolute_electronic_barrier_kcal_mol | -0.166379 | 28 | 28 | 23 |
| sensitivity_CN_2.40_to_2.15 | amine_contact | absolute_RRHO_G_barrier_kcal_mol | -0.267729 | 28 | 28 | 23 |
| sensitivity_CN_2.40_to_2.15 | opposite_contact | absolute_electronic_barrier_kcal_mol | -0.969323 | 96 | 30 | 2 |
| sensitivity_CN_2.40_to_2.15 | opposite_contact | absolute_RRHO_G_barrier_kcal_mol | -0.946931 | 96 | 42 | 2 |

4系の障壁ラベルの全24通りの並べ替えについて、候補全体の最大|r|を求めた。対象の|r|以上の最大値が何通りで得られるかを記録する。下表は同分類の全候補を同時に走査した診断である。

| window | target | pool | metric | barrier | n_label_permutations | observed_target_abs_r | permutations_with_scan_max_at_least_target | confirmatory_p_value |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| common_CN_2.50_to_2.10 | amine_contact | same_category | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 24 | 0.185312 | 24 | False |
| common_CN_2.50_to_2.10 | opposite_contact | same_category | radial_work_kcal_mol | absolute_electronic_barrier_kcal_mol | 24 | 0.966567 | 17 | False |

これは候補を多数見たときに高い相関がどれほど出やすいかの記述で、事前登録された確証検定や独立したp値ではない。対照同士は原子を共有して互いに依存する。弱い相関の対照だけを選んで「ここだけ相関する」とは結論しない。全候補、全3区間、ΔE‡とΔG‡の双方を保存した。

## 次に既存logだけで詰めること

現在の結果からは、相関の良い指標を追加探索するより、RAFA15で大きい反対側応答がO21単独の移動・環の回転・C–N形成に伴う共通変形で説明できるかを調べる方が有用。累積仕事と二面角・接触角の変化を同じC–N軸で比較し、指定ペアとH11···O18などの高相関対照を同時に示す。複数の原子力射影は同じ構造運動を共有するので、相関が多数出ること自体も結果として扱う。

原子仕事の座標合わせや射影仕事の微小値への数値誤差の影響は未評価の部分が残る。今回の対照比較は3区間で行ったが、追加の物理計算による独立検証ではない。

## 再現と参照

`python analysis4_force/contact_specificity.py`、または通常の全解析を実行する。出力は `contact_controls_*.csv`、図20〜22。既存の指定接触・原子仕事と今回の全候補走査の値が一致することをテストする。

原子力と移動からの仕事配分には [Haas & Chu (2009), DOI:10.1063/1.3243080](https://doi.org/10.1063/1.3243080) の背景がある。今回の接触軸への射影、対照抽出、特異性評価は本プロジェクトのad hocな設計であり、論文によりこの仮説が証明されたという意味ではない。
