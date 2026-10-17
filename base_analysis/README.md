# Focused Gaussian IRC Analysis

MoriPFとRAFA16について、設定した4距離・2角度・3種類のfragment COM距離、およびenergy/reaction force/Cartesian force projectionを比較するパイプラインです。exhaustive search、ranking、Top5、union、heatmapは行いません。

## 実行方法

Windows PowerShell:

```powershell
cd C:\Users\mugi\Desktop\FRIS_SATO\20260619_compare4\TS\IRC_ANALYSIS
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python src\run_analysis.py
```

LOGを必ず読み直す場合:

```powershell
python src\run_analysis.py --force-reparse
```

設定されたLOGが `data/` に無い場合だけ、同名ファイルをDownloadsから移動します。両方に存在する場合は上書きしません。parse cacheはLOGのSHA-256、反応方向符号、parser versionが一致する場合に `outputs/parsed/` から再利用します。

## 設定ファイル

### `config/series_config.yaml`

系列名、LOGパス、reaction directionを定義します。内部規約はnegative side、TS=0、positive sideです。化学的なreactant/productが逆の場合は、その系列の `reaction_direction_sign` を `-1` にしてください。

現在の設定では、保存済み第2解析と同じreactant-left規約を維持するためMoriPFに `-1`、今回明示指定されたRAFA16にも `reaction_direction_sign: -1` を適用します。したがって両系列のsorting、微分、COM補正、全CSV、全図は同じ負側reactant／正側product規約を使用します。

### `config/atom_map.csv`

解析対象positionのaliasと、各系列におけるGaussian atom numberを対応付けます。`expected_element` はmapping validationに使われます。コードはMoriPF/RAFA16の番号を直接保持せず、この表から解決します。

### `config/analysis_targets.yaml`

以下を定義します。

- `distances`: 解析する2原子aliasと使用するCOM baseline
- `angles`: 解析する3原子alias。中央のaliasが角度の頂点
- `com_baselines`: 2つのfragment名
- `fragments`: 系列ごとのGaussian atom number list

fragment listはコード内に埋め込まれていません。実行時に、原子番号がLOGに存在するか、重複がないか、baselineを構成するfragment同士が重複していないかを検査し、`outputs/tables/configuration_validation.csv` とreportへ記録します。

## COM距離とbaseline補正

本パイプラインのCOMは、各原子の質量を用いたcenter of mass（質量重心）です。

各対象距離 `d(ξ)` とbaseline COM距離 `D_COM(ξ)` について:

```text
Δd(ξ)     = d(ξ) - d(TS)
ΔD_COM(ξ) = D_COM(ξ) - D_COM(TS)

d_residual(ξ) = Δd(ξ) - ΔD_COM(ξ)
```

`d_residual` を「重心ベースライン補正距離」として図示します。

- reaction bond: `com_amine_coumarin_core` (`amine_side` vs `coumarin_core`)
- substituent–amine H-bond: `com_amine_substituent` (`amine_side` vs `substituent_side`)
- coumarin–amine H-bond: `com_amine_coumarin_core` (`amine_side` vs `coumarin_core`)
- interference距離: `com_substituent_coumarin_core` (`substituent_side` vs `coumarin_core`)

角度は長さと次元が異なるためCOM補正しません。角度にはraw値とTS基準deltaだけを出力します。

## CSVと図の対応

| 内容 | CSV | 図 |
|---|---|---|
| 4距離 raw | `distances_raw_long/wide.csv` | `summary/distances_all_raw.png`, `distances/raw/` |
| 4距離 TS delta | `distances_delta_from_ts_long.csv` | `distances/delta_from_ts/` |
| COM距離 | `com_distances_long/wide.csv` | `summary/com_distances.png`, `com/` |
| COM補正距離 | `distances_com_corrected_long/wide.csv` | `summary/distances_all_com_corrected.png`, `distances/com_corrected/` |
| 2角度 raw | `angles_raw_long/wide.csv` | `summary/angles_all_raw.png`, `angles/raw/` |
| 2角度 TS delta | `angles_delta_from_ts_long.csv` | `angles/delta_from_ts/` |
| Relative energy | `energy_profile.csv` | `energy/relative_energy.png` |
| Reaction force | `reaction_force.csv` | `energy/reaction_force.png` |
| Reaction force constant | `reaction_force_constant.csv` | `energy/reaction_force_constant.png` |
| Cartesian force projection | `projected_forces.csv` | `forces/projected_force_*.png` |

wide CSVは `series + signed_irc_coordinate + path_number + point_number` を1行とし、targetを列に展開しています。系列間でIRC sampling coordinateが厳密には一致しないため、系列同士を無理に同じcoordinate行へouter joinしていません。

## Energy derivativeとCartesian force projection

Energy解析では次を保存します。

```text
ΔE(ξ) = E(ξ) - E(TS)                       [kcal/mol]
F(ξ)  = -dE/dξ                              [kcal mol-1 IRC-1]
κ(ξ)  = dF/dξ = -d²E/dξ²                   [kcal mol-1 IRC-2]
```

raw derivativeは実際のsigned IRC coordinate間隔を使ったNumPy gradientです。3点以上では端点に2次one-sided differenceを使います。smoothed derivativeはenergy trajectoryをSavitzky–Golay filter（既定 `window_length=7`, `polyorder=3`, `mode=interp`）で平滑化してから同じ数値微分を行います。parametersは `config/series_config.yaml` で変更でき、実際に使用したwindow/polyorderもCSVに保存されます。

2原子間距離 `r_ab` に対するCartesian force projectionは、両端原子の相対forceを使って次のように定義します。

```text
e_ab = (r_b - r_a) / |r_b - r_a|
F_projected = 0.5 * (F_b - F_a) · e_ab
```

単一原子forceのdot productではありません。単位はHartree/Bohrです。TSなどCartesian forceが無いpointはNAとして保持し、energyや距離・角度解析は停止しません。

## Parser上の注意

実LOGにある `Point Number / Path Number`、forward/reverse marker、`NET REACTION COORDINATE`、orientation、SCF energy、`Forces (Hartrees/Bohr)` を使ってpointを再構成します。merged/restartによる同一 `(path, point)` は完全性とsource lineで重複除去します。

今回の両LOGではTS energyは `Energy From Chk` にありますが、TS Cartesian geometryはLOG本文にありません。このため:

- TS energy: checkpoint記録の値
- TS geometry: 最初のforward/reverse構造のCartesian midpointによる推定

を使用します。またMoriPFの未収束reverse terminal pointはNET reaction coordinateが無いため、直前5区間のmedian stepで外挿し、parsed CSVの `irc_coordinate_estimated` をtrueにします。これらはreportにもwarningとして残ります。

## 別系列の追加

1. LOGを `data/` に置く。
2. `series_config.yaml` に系列を追加する。
3. `atom_map.csv` に同名の列を追加し、8 aliasのGaussian atom numberを記入する。
4. `analysis_targets.yaml` の `fragments` に同名系列の4 fragment listを追加する。
5. `python src\run_analysis.py --force-reparse` を実行する。

系列名に依存する原子番号の分岐はありません。新系列でも同じ4距離、2角度、3 COM baseline、energy/force図、CSVが自動生成されます。
