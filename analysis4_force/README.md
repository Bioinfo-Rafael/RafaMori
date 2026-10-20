# 4系IRC・Reaction Force解析

**上位候補の個別確認：[原子・接触の位置、結合関係、角度、力の向きまで調べた結果](outputs/reports/higher_ranked_controls_summary_ja.md)**。4系それぞれについて、指定対象より上位の全件を掲載しています。

**正方向・逆方向の障壁：[方向別の表とRAFA15生成物についての考察](outputs/reports/directional_barriers_discussion_ja.md)**

**E・G・TΔSと反応経路：[3×3の図と熱力学成分の説明](outputs/reports/thermochemistry_profiles_ja.md)**。TSの熱補正を固定した参考曲線と、両端準位への滑らかな表示用接続を比較できます。

**図の説明：[全24種類・41枚の図の見方](outputs/reports/figure_guide_ja.md)**。各図の目的、軸の単位、点・線・色の意味、読み取れることを説明しています。

[解析結果](outputs/reports/summary_ja.md)が結論と解釈の入口です。ルートで `python analysis4_force/run_analysis.py` を実行すると、入力照合、解析用LOGの生成・再読込、CSV・図・summaryの生成まで実施します。

## 入力と再現性

- `../data/completed_20261017/manifest.json` に原本のファイル名、サイズ、SHA-256を保存。
- 新しい `LYS (1).LOG` は `Lys.log` として旧ファイルを置換。manifestに元ファイル名を記録。17ログを個別にgzip圧縮し、Gitに保存。展開後のバイト列をハッシュで照合する。
- `config/systems.yaml` は4系のIRC・反応物・生成物・TSと求核Nの対応。
- `config/atom_map.csv` と `analysis_targets.yaml` はTS構造から再生成。反応C20、O19、O21、フェニル環C23–C28という共有骨格の番号は入力に基づく前提で、元素と環の結合性を検査する。
- アミンNの共有結合Hを2個抽出。全Fとの組合せを保存し、TSで最短のF···Hペアを図示対象として固定する。途中で最短ペアを切り替えない。
- O19接触はTSで最も近いアミンH、O21接触はTSで最も近いFを選ぶ。各系で異なる化学的位置になる場合があるため、等価な相互作用を自動的に仮定しない。
- phenyl fragmentはC23–C28とその置換F・H。H置換原子も含め、アミン・残りのcoumarin側と57原子の非重複分割を構成する。旧解析の重原子中心のfragmentとは数値が異なり得る。

## IRC処理

元の共通parserを使用し、重複したpath/pointを整理する。反応座標がなく外挿された未収束末端点を除外。TS座標と電子エネルギーをfreqログの実値に置換する。原子元素列、TSエネルギー（許容差0.01 kcal/mol）、構造重ね合わせ（RMSD 0.10 Å以内）を確認する。

負側・正側のC–N距離を比較して符号を決める。解析用 `.log` は元の点番号を維持しつつpath 1を正、path 2を負として書き出す。力のベクトルは反応座標の符号反転では変わらない。TSの力は取得せず欠損として保持する。書き出し後の再読込で座標・エネルギー・力が10⁻¹⁰以内で一致することを検証する。

このLOGは抽出・整形された解析用の派生データであり、Gaussian計算の完了を意味しない。原本の終了状態や未収束情報は原本と監査JSONに残す。

## 数値計算と指標

- エネルギーはTS値を引いてからSciPy `make_smoothing_spline` に渡す。平滑化係数はGCVで決定し、オフセットを戻す。有限差分の反応力も保持。
- α/γは各側の4001点探索で全候補を拾い、局所最適化した極値と端点候補を比較して決定。
- κの図は論文規約 `+d²E/dR²`。旧規約 `−d²E/dR²` も既存列に保存し、論文規約には `kappa_paper_*` を付ける。
- 平滑化0（補間）、GCV、λ=10⁻⁴、λ=10⁻²でα/γ・α→TS仕事を比較。
- 原子の力と運動はTSへの剛体重ね合わせ後に `C_i = F_i · dr_i/dR` を計算。ΣC_iと反応力を比較する。原子・断片の相互作用エネルギーとは解釈しない。
- COM補正距離は `[d(R)−d(TS)]−[D_COM(R)−D_COM(TS)]`。
- 相関は同一IRC内の記述値。点の独立性を仮定する有意差検定はしない。全経路だけでなく領域別・TS前後に分けて保存する。

## エネルギーの基準

反応物分子＋Lysの組成・電荷がTS・生成物に一致し、計算条件・温度・圧力が一致することを検証する。新しいLysはM06-2X/6-311++G(d,p)・SMD水のopt/freq正常終了、33モード・虚振動なし。`energetics.py` で参照条件の不整合や熱化学の欠落を拒否する。

分離反応物基準の `ΔE‡ = E_TS − E_substrate − E_Lys`、`ΔG‡ = G_TS − G_substrate − G_Lys` と、生成物の反応エネルギーを計算する。各系の値をそのまま表示する。MoriPFを引く系間のゼロ合わせとそのCSVは廃止した。全電子エネルギー（Hartree）も保存する。

Gibbs値は298.15 K・1 atmのRRHOモデルで、1 Mや準調和への補正は行わない。RAFA15の微小虚振動は `RRHO_G_provisional=true` として残す。`absolute_barrier_valid` は条件・組成が整合した算術値のフラグであり、反応物の厳密な極小確認は `reactant_minimum_verified` に分ける。

4系比較・個別系のエネルギー図ともに各系の分離反応物基準のΔE(R)を使い、TSの高さはその系のΔE‡を表す。反応力は定数オフセットに依存しない。数値平滑化でTS値を一時的に引く処理は丸め誤差対策で、表示のゼロ合わせではない。

IRC左端からTSへの障壁・αまでの仕事は計算範囲に依存する。分離反応物基準の相対障壁やα→TSの仕事とは別々に保存する。生成物optはIRC末端と直接接続せず、エネルギー差と原子番号を固定したRMSDを監査する。

## 出力の入口

| ファイル／ディレクトリ | 内容 |
|---|---|
| `outputs/reports/summary_ja.md` | 日本語の結果報告、図、仮説評価、制限、参照文献 |
| `outputs/reports/stationary_point_audit.json` | 全opt/freqの条件、エネルギー、周波数、抽出行番号 |
| `outputs/reports/orientation_audit.json` | 符号、C–N距離、除外点、TS対応の検証 |
| `outputs/csv/system_summary.csv` | 4系の比較表 |
| `outputs/csv/energy_comparison.csv` | 各系の障壁、反応エネルギー、区間の仕事、生のHartree値 |
| `outputs/csv/all_NH_F_contacts.csv` | 全アミンH×Fの距離・角度 |
| `outputs/csv/descriptive_correlations.csv` | 距離・符号付き微分・反応力の領域別相関 |
| `outputs/csv/fluorine_positions.csv` | Fの原子番号とortho/meta/paraの対応 |
| `outputs/csv/endpoint_correspondence.csv` | 生成物optとIRC末端の対応診断 |
| `outputs/csv/smoothing_sensitivity.csv` | 平滑化に対する極値と仕事の感度 |
| `outputs/csv/force_motion_validation_summary.csv` | 原子寄与和と反応力の整合性 |
| `outputs/csv/parsed_oriented_*.csv.gz` | 実際に解析した全点の座標・力・エネルギー |
| `outputs/figures/` | 4系比較と各系のPNG図 |
| `outputs/oriented_logs/` | 方向を統一しTSを補った解析用LOG |

## 検証

```sh
python -m unittest discover -s analysis4_force/tests -v
```

圧縮原本の一致、4系の方向・TS、新Lysの条件一致と不整合入力の拒否、絶対障壁・MoriPF差・共通項相殺、仕事の加法性、微分と極値、断片和、全接触と固定ペアを検証する。

## 両側ortho接触の追加探索

`two_sided_contacts.py` はC24上のF/H–アミンHと、C25上のF/H–O21を対応づける。旧 `interference_F` 指標はF限定なので、RAFA15のH56–O21はこの追加解析で扱う。主区間は全系共通のC–N=2.50→2.10 Å。3区間×アミンHの2選択について、21種類の指標と電子／RRHO Gibbs障壁の相関を全て保存する。各系の実値を使い、MoriPFとの差へ変換しない。

`two_sided_atom_mapping.csv`、`two_sided_force_profiles.csv`、`two_sided_window_metrics.csv`、`two_sided_correlations.csv` と図14–16、[追加summary](outputs/reports/two_sided_contact_summary_ja.md)を生成する。主指標の射影力比と、探索上高相関だった各原子に働く合計の力ノルム比を区別する。力は全相互作用の合計であり、F/H···Oの独立した引力・斥力ではない。

通常の `run_analysis.py` から自動実行される。既存の4系結果から追加部分だけ再計算する場合は `python analysis4_force/two_sided_contacts.py`。

## 既存log限定の指標と文献調査

[調査・計算方針](outputs/reports/log_only_research_plan_ja.md) と [両側仕事の結果](outputs/reports/two_sided_contact_summary_ja.md) を参照。追加DFTなしで、共通C–N区間の原子仕事 `−∫F_i·dr_i` をアミン側・反対側・残りに分割し、接触方向の仕事も別に保存する。21指標×3区間×2通りのアミンH×2障壁を記録する。比率は補足、主図15はkcal/molの仕事と各系の障壁。全原子仕事と区間エネルギー上昇、区間前後を含む全障壁の収支をテストする。

## 全原子・全接触の対照比較

`python analysis4_force/contact_specificity.py` で、全57原子・全1596原子対について同じ区間の原子仕事、2原子の仕事の合計、接触方向の仕事、接触方向の力の平均絶対値を計算する。対照は共有結合グラフ・4 Å接触条件から抽出し、分子内/分子間を分け、距離・化学種・原子共有・1–4対除外の感度比較を行う。42個の対応が一意なサイトのみ系間相関に使い、アルキル基上の対称Hは分子内順位には含めるが系間相関から除く。条件は `config/contact_specificity.yaml`。通常の全解析にも組み込んだ。

[結果と制限](outputs/reports/contact_specificity_summary_ja.md)、`contact_controls_*.csv`、図20〜22を参照。

## 日本語の図解だけを更新する

`python analysis4_force/figure_explanations.py` で、既存の結果報告の図説明、全図の案内、出力ファイルのハッシュ一覧を更新できます。数値解析や図の再計算は行いません。通常の解析を再実行した場合も、結果報告に図の説明が付きます。

## 指定対象より上位の原子・接触を同定する

`python analysis4_force/higher_ranked_controls.py` で、既存の全原子・原子対の比較結果に、各原子の役割、結合経路、共通C–N区間での接触角、符号付きの力を付けます。指定対象より上位の全件を3区間・7比較条件について保存し、各系の日本語一覧と考察を生成します。通常の全解析にも組み込んであります。
