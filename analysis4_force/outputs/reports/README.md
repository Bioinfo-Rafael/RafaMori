# 解析資料の案内：どのMarkdownを読めばよいか

更新日：2026-10-25（JST）。現在の対象はMoriPF・RAFA14・RAFA15・RAFA16の4系。ここを文書の入口とし、目的別に必要な資料を開く。似た内容の報告は解析を追加した順に残っているため、古い総括と最新の比較結果を区別する。

## 最初に読む3つ

1. [仕事の大きさ・ばらつき・相関を合わせた考察](work_candidate_assessment_ja.md)：現在の結論と、Fによる反応促進についてどこまで言えるか。
2. [4系を横に並べた仕事の比較](cross_system_work_comparison_ja.md)：原子・接触の仕事、Pearson r、R²、STD/mean、散布図・棒グラフ。数値表と詳細説明表は分離してある。
3. [正方向・逆方向の活性化障壁](directional_barriers_discussion_ja.md)：電子エネルギー・自由エネルギー・エントロピー項の障壁と、RAFA15の反応物・生成物更新による両方向障壁の再計算。

原子番号が分からなくなったら[全57位置の対応表](atom_correspondence_audit_ja.md)、図の読み方が分からなければ[図の説明](figure_guide_ja.md)へ進む。

## 現在、Fの効果をどこまで定量化できているか

**「何も定量化できていない」わけではない。F置換位置に伴う障壁、原子の仕事、接触方向の仕事、4系間のばらつきと相関は数値化できている。一方、その違いをF···H相互作用による安定化に帰属し、『この相互作用を強めれば障壁が下がる』と因果的に結論するところまではできていない。**

今回比べているのはFの有無だけを変えた単一条件ではなく、F/Hの置換位置が異なる4系である。RAFA14にもFは存在する。「Fが存在する効果」と「ortho位置にFがある効果」「他の位置との組合せ」を一緒にしない。

以下の仕事は、反応点C20–求核Nの距離が2.50 Åから2.10 Åになる共通区間の値である。相関の縦軸は分離反応物基準の正方向ΔG‡、点数は4系。区間内の仕事と反応全体の障壁は区別する。数値の出典は[候補評価CSV](../csv/work_candidate_assessment.csv)。

| 数値化できたこと | 現在の読み方 |
| --- | --- |
| 4系の電子障壁ΔE‡と自由エネルギー障壁ΔG‡ | 計算結果として反応性の差を比較できる。分離反応物＋Lysを共通の定義とし、各系の障壁値そのものを表示する。 |
| アミン側置換F/H原子の仕事とΔG‡ | r=+0.679、R²=0.461。新しい反応物基準で再計算した記述的相関。局所的な相互作用エネルギーや因果効果ではない。 |
| 反応カルボニルO21の原子仕事とΔG‡ | r=+0.874、R²=0.763。新しい反応物基準で再計算した記述的相関。局所的な相互作用エネルギーや因果効果ではない。 |
| ラクトンカルボニルO19の原子仕事とΔG‡ | r=-0.755、R²=0.570。新しい反応物基準で再計算した記述的相関。局所的な相互作用エネルギーや因果効果ではない。 |
| 指定アミン側F/H···アミンHの接触方向の仕事 | r=+0.714、R²=0.509。新しい反応物基準で再計算した記述的相関。局所的な相互作用エネルギーや因果効果ではない。 |
| 反対側F/H···O21の接触方向の仕事 | r=+0.598、R²=0.357。新しい反応物基準で再計算した記述的相関。局所的な相互作用エネルギーや因果効果ではない。 |

新基準ではO21・アミン側F/Hの相関が正、O19は負となった。旧基準での相関の符号と候補評価は引き継がない。O21は原子仕事全体の最大成分ではない。原子仕事の平均は求核Nが約0.990、C20が約0.636、O21が約0.0566 kcal/mol。C20–求核Nの接触方向の仕事は平均約1.460 kcal/molで、これは原子仕事とは別の量である。大きさと相関を一緒に見る場合も、原子仕事同士、接触方向の仕事同士で比較する。

4系の相関は探索的な関係であり、局所的な仕事の大きさをそのまま相互作用エネルギーと呼ばない。正の原子仕事は共通区間内の電子エネルギー上昇への寄与で、「安定化エネルギー」ではない。

## 障壁差は「1 kcal/mol程度」なのか

現在の計算値は次のとおり。単位はkcal/mol。ΔE‡とΔG‡は、それぞれTSから分離反応物（基質＋Lys）を引いた差である。

| 分子系 | 正方向ΔE‡ / kcal/mol | 正方向ΔG‡ / kcal/mol |
| --- | ---: | ---: |
| MoriPF | −1.342 | 14.083 |
| RAFA14 | 0.476 | 15.928 |
| RAFA15 | 4.235 | 18.063 |
| RAFA16 | −0.579 | 15.023 |

2026-10-22更新：RAFA15生成物を `RAFA15_LYS_OPTFREQ2.LOG` に置換。逆向きΔE‡ = 2.263、ΔG‡ = 0.996 kcal/molとなり、旧ログの負値は解消した。その後2026-10-25に反応物も更新し、正方向はΔE‡ = 4.235、ΔG‡ = 18.063へ上昇した。今回の反応物置換では逆向きの値は変わらない。

出典：[energy_comparison.csv](../csv/energy_comparison.csv)。RRHO自由エネルギーは298.15 K・1 atm。RAFA15の新反応物は虚振動なしで、従来の暫定フラグを解除した。電子障壁が負でも、自由エネルギー障壁がないという意味ではない。

| 比較する2系 | 正方向ΔG‡の差の大きさ / kcal/mol |
| --- | ---: |
| MoriPFとRAFA16 | 0.940 |
| RAFA16とRAFA14 | 0.905 |
| MoriPFとRAFA14 | 1.845 |
| RAFA15とMoriPF | 3.980 |
| RAFA15とRAFA14 | 2.135 |

したがって、**RAFA15以外の隣り合う順位同士は約1 kcal/mol差だが、4系全体の障壁差が1 kcal/molしかないわけではない。** また、同じ温度・標準状態・透過係数を仮定する遷移状態理論では、低障壁側と高障壁側の速度定数の比は

$$
\frac{k_{\mathrm{low}}}{k_{\mathrm{high}}}
=\exp\left(\frac{\Delta G^\ddagger_{\mathrm{high}}-\Delta G^\ddagger_{\mathrm{low}}}{RT}\right)
$$

となる。298.15 Kでは、**ΔG‡が1 kcal/mol低いことは約5.4倍の速度定数に対応**する。1 kcal/molは、反応速度への影響として必ずしも小さくない。ただし、これは理論式による換算であり、今回の計算値がその精度で実測速度を予測できると確認したものではない。[IUPAC：transition state theory](https://goldbook.iupac.org/terms/view/T06470)

ΔE‡の差をそのままこの速度比に代入せず、ΔG‡の差を使う。今の議論の限界は、単に障壁差が小さいことよりも、4系という点数、配座・熱補正への依存、仕事が特定の原子対の相互作用エネルギーではないことにある。

## 現在の4系解析：各報告の役割

| Markdown | 何が書かれているか | どんなときに読むか・位置付け |
| --- | --- | --- |
| [README.md](README.md) | この案内。読む順番、各文書の役割、現在の結論、障壁差の大きさ。 | 最初に開く入口。 |
| [work_candidate_assessment_ja.md](work_candidate_assessment_ja.md) | 仕事の大きさ・絶対／相対変動・相関を同時に比較。優先候補、対照候補、1系除外時の相関、分子設計への解釈。 | 現在の考察を読む。 |
| [cross_system_work_comparison_ja.md](cross_system_work_comparison_ja.md) | 対応位置を行、4系を列にした仕事の和集合表。r・R²・STD/mean。各系の順位・選択理由は別表。図25・26も収録。 | 系間差を数値・図で確認する主なページ。 |
| [higher_ranked_controls_summary_ja.md](higher_ranked_controls_summary_ja.md) | 番号付き構造式、仕事・Qの数式、障壁に対する割合、横並び表、各系内の上位・下位比較をまとめた大きな資料。 | 定義から個別根拠まで一つの文書で追う。長いので通常は上の2文書から読む。 |
| [higher_ranked_controls_MoriPF_ja.md](higher_ranked_controls_MoriPF_ja.md) | MoriPFの番号付き構造式、全57原子の仕事、着目接触の上位全件と直下3件。 | MoriPF内で何が大きいかを調べる。 |
| [higher_ranked_controls_RAFA14_ja.md](higher_ranked_controls_RAFA14_ja.md) | RAFA14について同じ構成。 | ortho位置がHの系の比較を調べる。 |
| [higher_ranked_controls_RAFA15_ja.md](higher_ranked_controls_RAFA15_ja.md) | RAFA15について同じ構成。 | アミン側F／反対側Hの系を調べる。 |
| [higher_ranked_controls_RAFA16_ja.md](higher_ranked_controls_RAFA16_ja.md) | RAFA16について同じ構成。 | 両側ortho位置がFの系を調べる。 |
| [atom_correspondence_audit_ja.md](atom_correspondence_audit_ja.md) | 全57位置の4系間対応。結合関係、局所配置、共通IRC区間41点・TSによるH対応の確認。 | 同じ行が本当に同じ位置か、原子番号がどうずれるかを確認する。最新の横並び表・図25/26はこの対応を使用。 |
| [directional_barriers_discussion_ja.md](directional_barriers_discussion_ja.md) | 正方向・逆方向のΔE‡・ΔG‡・TΔSの表と考察。RAFA15新生成物とIRC末端の対応改善も説明。 | 障壁・生成物の安定性・逆反応を読む。 |
| [thermochemistry_profiles_ja.md](thermochemistry_profiles_ja.md) | E・G・TΔSを3×3で表示。IRCのみ、両端の反応物／生成物準位、滑らかな接続を比較。 | エネルギープロファイルを見る。IRC上のG・TΔSはTS熱補正固定の参考曲線で、各点のfreq計算ではない。 |
| [figure_guide_ja.md](figure_guide_ja.md) | 図01〜26の目的、軸・単位、色、パネル、解釈、amine/opposite・matched/broadの意味。 | 図の読み方を確認する。 |
| [summary_ja.md](summary_ja.md) | 4系解析の基礎的な総括。障壁、IRC方向、反応力、距離、計算条件と制約。 | 全体像と入力計算を知る資料。後から追加した和集合・全H対応・候補評価までは網羅しない。 |
| [two_sided_contact_summary_ja.md](two_sided_contact_summary_ja.md) | アミン側と反対側の仕事・距離・力、Hの選び方、集計区間への感度。 | 両側の指標を導入した段階の定義と結果を確認する。現在の全候補比較は横並び表を優先。 |
| [contact_specificity_summary_ja.md](contact_specificity_summary_ja.md) | 他の原子・接触との順位比較、相関の特異性、比較候補の条件。 | 対照比較の方法を確認する。初期の対応範囲での結果を含み、最新の全H対応による和集合の値は横並び表・図25/26を参照。 |
| [log_only_research_plan_ja.md](log_only_research_plan_ja.md) | 文献調査、追加DFTなしで計算可能な指標、原子仕事・接触方向の力の位置付け。 | 指標を選んだ理由と文献を確認する。計画の一部は既に実装済み。 |

## 数式・区分の説明原稿

以下は報告書に埋め込む説明文の元ファイルであり、別の解析結果ではない。定義を変更するときは、この原稿を更新して報告書を再生成する。

| Markdown | 内容 |
| --- | --- |
| [force_work_definitions_ja.md](../../docs/force_work_definitions_ja.md) | 原子仕事、2原子の仕事の和、接触方向の仕事、Q、平均絶対値、符号付き平均、積分・単位。 |
| [work_comparison_groups_ja.md](../../docs/work_comparison_groups_ja.md) | amine/opposite、matched/broad、距離・原子種・結合関係の比較条件。 |

## プロジェクト説明・旧解析

| Markdown | 内容・位置付け |
| --- | --- |
| [ルートREADME](../../../README.md) | プロジェクト概要、ディレクトリ、実行方法、理論の参照元。この文書への入口。 |
| [analysis4_force/README.md](../../README.md) | 現在の4系解析の実装、入力、再現性、IRC処理、数値計算、検証方法。 |
| [図フォルダのREADME](../figures/README.md) | 図ファイルの一覧。各図の詳しい説明はfigure_guideにある。 |
| [base_analysis/README.md](../../../base_analysis/README.md) | 旧2系の基礎解析の実行方法・設定。parserなど一部は4系解析でも再利用。 |
| [base_analysisのanalysis_summary.md](../../../base_analysis/outputs/reports/analysis_summary.md) | 旧2系の基礎解析結果。現在の4系比較の結論には使わない。 |
| [analysis2_force/README.md](../../../analysis2_force/README.md) | 旧2系の反応力解析の目的・理論・数値処理・注意点。 |
| [analysis2_force/README_ANALYSIS.md](../../../analysis2_force/README_ANALYSIS.md) | 旧解析の実装記録。原子対応、図の設計、平均・標準偏差の扱いなど。 |

旧`analysis2_force/outputs/`はGit管理外のため、この案内では対象外とする。

## 重複と更新の扱い

数値の追跡には各報告からリンクされたCSVを使う。現在の系間比較は`cross_system_work_comparison_ja.md`、その解釈は`work_candidate_assessment_ja.md`、原子対応は`atom_correspondence_audit_ja.md`を優先する。初期報告の相関や候補数は、対象・区間・比較条件が異なることがあるため混ぜない。

今回の整理では、入口となるREADMEを追加し、既存報告は参照可能な形で残した。図・仕事・障壁の計算値は変更していない。
