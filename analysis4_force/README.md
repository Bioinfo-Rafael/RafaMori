# 4系IRC・Reaction Force解析

[結果summary](outputs/reports/summary_ja.md)が結論と解釈の入口です。ルートで `python analysis4_force/run_analysis.py` を実行すると、入力照合、解析用LOGの生成・再読込、CSV・図・summaryの生成まで実施します。

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

分離反応物基準の `ΔE‡ = E_TS − E_substrate − E_Lys`、`ΔG‡ = G_TS − G_substrate − G_Lys` と、生成物の反応エネルギーを計算する。系間差はすべて `config/systems.yaml` の **MoriPF** を基準にする。`MoriPF_baseline_differences.csv` にTSの距離・角度・α→TS仕事のMoriPFとの差も保存する。

Gibbs値は298.15 K・1 atmのRRHOモデルで、1 Mや準調和への補正は行わない。RAFA15の微小虚振動は `RRHO_G_provisional=true` として残す。`absolute_barrier_valid` は条件・組成が整合した算術値のフラグであり、反応物の厳密な極小確認は `reactant_minimum_verified` に分ける。

4系比較エネルギー図では、各系の分離反応物基準の電子エネルギーからMoriPFの電子障壁を引き、MoriPFのTSを0とする。個別系の形状診断図では引き続き各TS＝0を使用し、系間の相対障壁とは区別する。反応力は定数オフセットに依存しない。

IRC左端からTSへの障壁・αまでの仕事は計算範囲に依存する。分離反応物基準の相対障壁やα→TSの仕事とは別々に保存する。生成物optはIRC末端と直接接続せず、エネルギー差と原子番号を固定したRMSDを監査する。

## 出力の入口

| ファイル／ディレクトリ | 内容 |
|---|---|
| `outputs/reports/summary_ja.md` | 日本語summary、図、仮説評価、制限、参照文献 |
| `outputs/reports/stationary_point_audit.json` | 全opt/freqの条件、エネルギー、周波数、抽出行番号 |
| `outputs/reports/orientation_audit.json` | 符号、C–N距離、除外点、TS対応の検証 |
| `outputs/csv/system_summary.csv` | 4系の比較表 |
| `outputs/csv/energy_comparison.csv` | 絶対・MoriPF基準の相対障壁、反応エネルギー、区間の仕事 |
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
