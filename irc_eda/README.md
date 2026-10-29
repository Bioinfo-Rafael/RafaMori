# Gaussian IRC → IQA / sobEDA Adapter

別のWindows PCへの移行は[Windows手順](WINDOWS.md)を参照してください。共有用ZIPは`python irc_eda/scripts/package_transfer.py`で作成できます。

Gaussianの既存計算を読み、公式APOST3D `ENPART` / easyEDA同梱`EDA_run.sh`に接続するCLIです。IQA/EDAの積分・電子状態計算は実装していません。APOST3Dの`EDAIQA`は使用せず、sobEDAとの成分の合算もしません。

**今回の実行結果:** 4系252点のIRC抽出・XYZ・固定座標SP入力・エネルギープロファイルを生成。IQA/EDAの新規数値計算は0点です。詳細は[実行報告](reports/execution_report_ja.md)。M06-2XのAPOST3D非対応、SMDの分解未検証、実行ソフト不足を明示して停止しています。

## 構成

| パス | 内容 |
|---|---|
| `parsing.py` | cclibとGaussian完了マーカーによるIRC抽出、fchk情報、原子順序固定の照合 |
| `upstream.py` | 公式入力・出力仕様、対応手法判定、環境調査 |
| `pipeline.py` | inventory、準備、実行、検証、resume、集計、3図 |
| `samples.py` | 公式サンプルと保存済み公式出力の検証を区別 |
| `vendor/APOST3D/` | clone済みGPL-3.0実装。ローカル保持、Git追跡対象外 |
| `vendor/easyEDA/` | clone済みMIT実装。ローカル保持、Git追跡対象外 |
| `upstream-lock.json` | 再取得できるURLと固定commit |
| `config/project.yaml` | このリポジトリの4系の設定 |
| `outputs/project/` | 実際の生成物。元データと分離、Git追跡対象外 |
| `reports/` | 実行報告・環境/検証記録 |

既存の`analysis4_force/`等は変更していません。`01_Vault/project.md` / `docs/project.md`は存在せず、MTD-RL Vault向けの保存規則は適用対象外です。

## セットアップとCLI

リポジトリルートから実行します。Python 3.10以降が必要です。この環境では専用venvを作成し、Python 3.14で検証しました。

```sh
bash irc_eda/scripts/bootstrap.sh
bash irc_eda/scripts/run.sh inspect
bash irc_eda/scripts/run.sh prepare
bash irc_eda/scripts/run.sh check-samples
bash irc_eda/scripts/run.sh run-iqa
bash irc_eda/scripts/run.sh run-eda
bash irc_eda/scripts/run.sh analyze
bash irc_eda/scripts/run.sh run-all
irc_eda/.venv/bin/python -m pytest irc_eda/tests -q
```

すべてのコマンドは`-c /path/to/config.yaml`を受け付けます。直接`python -m irc_eda ...`も使用可能です。`run-all`はinspect → prepare →公式サンプル確認 → IQA → EDA → analyzeを実行します。依存不足はJSONに`BLOCKED`として記録し、残りの処理は続けます。外部計算失敗は`FAILED`、入力不備はCLI終了コード2です。コマンド終了コード0だけを量子計算成功の判定に使わず、`execution_summary.json`と各`latest.json`を確認してください。

### 設定

パスはYAMLのあるディレクトリから解決します。`systems`を変更して他のIRCに再利用できます。

```yaml
input_root: /path/to/originals
output_dir: /path/to/new_adapter_output
scope: pilot
parallel_jobs: 1
threads: 4
memory: 4GB
executables:
  # PATHにある場合は省略可能。環境は先にsourceしておく。
  g16: /licensed/gaussian/g16/g16
  formchk: /licensed/gaussian/g16/formchk
  unfchk: /licensed/gaussian/g16/unfchk
  Multiwfn: /path/to/Multiwfn
  apost3d: /path/to/APOST3D/apost3d
  get_energy_g16: /path/to/APOST3D/utils/get_energy_g16
systems:
  molecule:
    irc: IRC.log
    # forward/reverseが別ファイルの場合はircの代わりにirc_filesを指定。
    # irc_files: [forward.log, reverse.log]
    ts: TS.log
    coordinate_sign: 1
    # 省略時は座標最小点を参照。化学的な反応物側は利用者が指定する。
    reference_point: molecule_p2_0035
    fragments: null
```

`inventory_exclude`はinput_rootからの相対ディレクトリ指定です。付属設定では自身の生成物`irc_eda`を除外します。`.log/.out/.chk/.fchk/.fch/.gjf/.com/.xyz`とgzip圧縮版を、大文字拡張子も含めて調査します。圧縮版と展開版は別ファイルのinventory行として残します。

### fragment定義

未指定時には`outputs/project/fragments.template.yaml`を生成します。化学的な分割は推定しません。編集した定義を設定の`systems.<name>.fragments`へ転記します。

```yaml
fragments:
  - atoms: [1, 2, 3]   # 例示。実系の分割ではない。
    charge: 0
    multiplicity: 1
  - atoms: [4, 5, 6]
    charge: 0
    multiplicity: 1
```

全原子を1回ずつ割り当て、全電荷、電子数と多重度のパリティ、全スピン射影を検査します。原子番号はGaussian元配列の1始まりです。sobEDAで多重度を負にするとフラグメントのスピンを反転するため、その符号も指定できます。設定の化学的妥当性までは自動認定しません。Multiwfn出力でフラグメント単位に原子が並び替わった場合も、明示した分割による順列だけを検証し、対応表を保存します。公開CSV/XYZのindexは元順序を維持します。

## 外部プログラム

- [APOST3D README](https://github.com/mgimferrer/APOST3D/blob/3ba2d95b58831d7a6adaaedc40cc4ed950eb16c5/README.md)、[DOCUMENTATION](https://github.com/mgimferrer/APOST3D/blob/3ba2d95b58831d7a6adaaedc40cc4ed950eb16c5/DOCUMENTATION.md)に従い、Intel oneAPI 2024の`ifort`・`icx`と互換Libxc 4.2.3を使用します。
- `bash irc_eda/scripts/build_apost3d.sh`は公式`compile_libxc.sh`・`make_compile.sh`を呼びます。元ソースは変更しません。固定commitにはREADMEが言及する`libxc-4.2.3.tar.gz`が含まれていないため、未配置の場合は明示停止します。新しいLibxcやgfortranへ黙って差し替えません。
- [easyEDA](https://github.com/kangmg/easyEDA)はGaussian ≥16 A.03、Multiwfn ≥2023-06-23に加え、実コード上は`formchk`、`unfchk`、`bc`、`dos2unix`、GNU sed、bashが必要です。実行対象はLinux/WSLです。GUI用の絶対パスを書き換える`first_install.py`は使用せず、公式`source/EDA_run.sh`を点ごとの作業領域でそのまま呼びます。
- Gaussianライセンス環境は利用者が通常どおり読み込んでおきます。独自の代替量子化学ソフトには切り替えません。
- `environment.json`に実行ファイル位置・検出したバージョン・共有ライブラリ情報・バイナリSHA256・OS・upstream commitを保存します。起動せずには確認できない製品バージョンは`unknown`と記録します。`software_versions`へ確認済みバージョンを指定できます。Gaussianログの生成元バージョンはinventoryに別途記録します。

### 手法の制約

IQAはTFVCを既定とし、HF、B3LYP、BP86の公式キーワードを使用できます。その他のLIBRARY指定は`iqa.libxc`にGaussian手法名・family（LDA/GGA/HYB_GGA）・ID・互換性確認の根拠`verified_equivalence_source`が必要です。M06-2X等のmeta-GGAは、固定commitの評価処理が未実装なので禁止しています。単にLibxcのIDが存在するだけでは対応としません。

SMD/PCMの場合は、両手法ともこのAdapterで妥当性確認済みの分解として実行しません。溶媒反応場・非静電寄与の扱いを明確にせず、ガス相の総エネルギーと比較しないためです。ガス相への条件変更は明示的に別研究条件として用意する必要があります。この制約を外すだけのスイッチはありません。

`ChkBasis` / `Gen`のfchkでは基底名の一致を自動認定しません。原子順序、整列RMSD、エネルギー、基底シェル/指数/係数の指紋を保存します。別途確認できた場合のみ、`wavefunction_conditions`でinput_rootからの相対ファイル名をキーに`basis`と`evidence`を指定できます。現行M06-2X/SMDのBLOCKEDは基底名の確認後も解除されません。

## 抽出と検証の規約

- cclibを用いてSCFエネルギーと構造を独立照合し、不足する完了点番号・方向・反応座標をログのマーカーで補います。
- 完了`Point Number:`に続く`NET REACTION COORDINATE`のある点だけを採用します。SCF反復数はIRC点数ではありません。
- 座標とエネルギーは、採用点の直前の**同じSCF評価**から取得します。HPC補正後の補間端点そのものが別途印字されていない場合は、最終SCF評価構造であることを`geometry_kind`に明記し、補正端点エネルギーとの差も保存します。
- forward正・reverse負を標準に、設定の`coordinate_sign`を掛けます。単位はGaussianの質量重み付き座標`√amu·bohr`、CartesianはÅです。今回の4系の方向設定は既存解析と一致しています。非質量重み付きIRCの独自形式は検証対象外です。
- TSを補間しません。明示的TSログを指定した場合は手法・基底・電荷・多重度・溶媒・原子順序・IRCチェックポイントエネルギー・隣接構造RMSDを検査します。指定なしでは、IRC内に実際の初期構造が印字されている場合だけ採用します。
- 重複点の不一致はエラー、欠落点は欠落のままです。合成ログの最後の正常終了だけで全経路正常とは判断せず、正常/異常終了イベント数を別々に保存します。
- `.chk`は`formchk`で出力領域に変換し、変換後の内容を照合します。basenameだけの対応付けはしません。
- 追加Gaussian計算は元の電子条件を保持した固定座標SPです。最適化・IRC・チェックポイント座標読み出しを外します。未解決の汎用基底・ECPカード等は自動生成しません。
- APOST3Dの公式`get_energy_g16`による追記は、対応が検証された単一ジョブ正常終了ログから、作業用fchkコピーにだけ行います。
- IQAは対角＋各非対角ペア1回の和を総エネルギーと比較します。sobEDAは`E_els + E_x + E_rep + E_orb + E_DFTc + E_dc = E_tot`を検査します。`E_xrep`、`E_c`は重複する合成項です。
- IQAのHartree→kcal/molは627.5094740631、sobEDA印字値の逆換算は公式コードの627.51です。sobEDAの通常出力は0.01 kcal/mol丸めなので、closure許容値は0.05 kcal/molです。
- 不明項は空欄/JSON null、未計算点は状態表に記録します。実計算CSVに架空の行は挿入しません。

## 状態・resume

`readiness.csv/json`は全252点について`READY` / `CONVERSION_REQUIRED` / `RECALC_REQUIRED` / `BLOCKED`と理由を記録します。必要データの状態`*_data_status`と、ソフト・手法を含めた実行可能性`*_status`を分けています。sobEDAはXYZがあっても内部で量子計算を要するので、準備が整った状態は`RECALC_REQUIRED`です。

既定の`scope: pilot`は各系の両端とTS付近、検証済み既存波動関数の点を選びます。`scope: all`は同じ条件・出力ハッシュでpilotがすべて成功していない限り停止します。成功点は入力、コード、設定、外部バイナリ、commit、出力ハッシュが一致すると再実行しません。scratchは点・手法・試行ごとに独立し、失敗試行も保存します。同一点を同時起動するとロックします。プロセスが強制終了した場合は、生存確認後に該当`running.lock`を手動解除してください。

## 生成物

[実際の出力ディレクトリ](outputs/project/)に以下を保存します。

| ファイル | 内容 |
|---|---|
| `inventory.csv/json` | ファイル単位の条件、終了状態、原子数、hash、点対応 |
| `irc_profile.csv/json`、`points.json`、`xyz/` | 点の番号・座標・エネルギー・出典・参照差分 |
| `readiness.csv/json` | 全点の両手法の実行可能性と理由 |
| `preparation.csv/json`、`prepared/` | 固定座標SP等の入力と生成できなかった理由 |
| `iqa_atomic.csv`、`iqa_pairwise.csv` | 検証成功したIQA原子/ペア成分、Hartree/kcal/mol、参照差分 |
| `eda_components.csv`、`eda_stages.csv` | sobEDA全9項、任意sobEDAw印字項、各段階の印字されたHartree成分 |
| `validation_report.csv/json` | IRCの照合、計算状態、closure、元IRCエネルギーとの差 |
| `01_irc_profile.png` | 4系の原IRCのSCFエネルギー。指定IRC端点基準で、分離反応物基準の障壁ではない |
| `02_iqa_interactions.png` | 大きく変化するペアまたは`plot_pairs`。数値なしの場合はその旨を表示 |
| `03_eda_components.png` | 同一図で成分差分を比較。数値なしの場合はその旨を表示 |
| `official_samples/` | 保存済み出力の検証と、新規サンプル実行の状態を分離 |
| `runs/*/*/latest.json` | 点別実行状態。成功出力のhashと試行先 |
| `execution_summary.json` | 実際の計算成功数、SP入力数、Gaussian実行数 |

各点の条件・元ファイルと、外部バージョン/commitを含む`provenance.json`への参照をCSV/JSONに保持します。参照点の分解値が未計算なら差分も空欄です。IQAとsobEDAの結果は別々のスキームとして保存します。

原著: [sobEDA, J. Phys. Chem. A (2023)](https://doi.org/10.1021/acs.jpca.3c04374)、[APOST3D, J. Chem. Phys. 160, 172502 (2024)](https://doi.org/10.1063/5.0206187)。各upstreamのLICENSEと引用案内をvendor内に保持しています。
