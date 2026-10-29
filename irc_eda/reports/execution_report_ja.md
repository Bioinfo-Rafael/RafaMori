# Gaussian IRC → IQA / sobEDA 実装・実行報告

実施日：2026-10-29 JST。作業ディレクトリ：`/Users/cls-lab/Git/RafaMori`。

**Adapter、CLI、入力生成、inventory、集計・可視化、テストを実装し、既存IRC 4系252点を実際に処理しました。新しいIQA・EDA数値計算の成功数はともに0点です。** 不足ソフトに加え、採用したAPOST3DのM06-2X非対応と、SMD分解の未検証があるため、計算条件を変更して実行することはしていません。

## 到達範囲

| 区分 | 結果 | 根拠 |
|---|---|---|
| 外部実装取得 | 完了、2リポジトリ、tracked source変更なし | [固定commit](../upstream-lock.json)、[環境](../outputs/project/environment.json) |
| Adapterコード・CLI | 実装済み。外部実機接続は依存不足で未検証 | [README](../README.md)、[テスト結果](tests.log) |
| Python環境 | 専用venvを構築、cclib 1.8.1等を導入 | [依存バージョン](../requirements-lock.txt) |
| APOST3Dビルド | 未完了。実行スクリプトを起動し、ifort不足で停止 | [ビルド記録](build_attempt.log) |
| 公式APOST3Dサンプルの実行 | 未実施、実行バイナリなし | [サンプル状態](../outputs/project/official_samples/sample_status.json) |
| 公式easyEDAサンプルの新規計算 | 未実施、Gaussian/Multiwfn等なし | 同上 |
| 公式easyEDA保存済み出力の読み取り | 成功。9項、closure確認 | [成分](../outputs/project/official_samples/easyEDA_recorded_components.csv)、[検証](../outputs/project/official_samples/easyEDA_recorded_validation.json) |
| 実IRCの抽出・XYZ・プロファイル | 252点で完了 | [irc_profile.csv](../outputs/project/irc_profile.csv)、[図](../outputs/project/01_irc_profile.png) |
| 実データIQA成功 | **0 / 252点** | [実行集計](../outputs/project/execution_summary.json) |
| 実データEDA成功 | **0 / 252点** | 同上 |
| Gaussian追加SP | 入力252件、実行0件 | [入力一覧](../outputs/project/preparation.csv) |
| 原データ保全 | 選択した8ログとfchk/圧縮版4ファイル、計12ファイルのhash不変 | [保全検証](source_integrity.json) |

[run-allの実行記録](run_all.log)を保存しました。各系の反応物側・TS・生成物側の計12点をpilot対象として実行制御に渡し、IQA・EDAそれぞれのBLOCKEDを記録しました。これはエンジン実行失敗ではなく、事前検査での停止です。エンジン失敗試行数も0です。

## 調査したデータと抽出結果

[最終inventory](../outputs/project/inventory.csv)はワークスペースの164ファイルを記録しています。内訳はログ122、Gaussian入力22、XYZ 16、fchk 4です。拡張子は大文字小文字を区別せず、gzipを除いた種類で集計しています。圧縮版・展開版・派生ファイルも含むため、164を独立計算数とは解釈しません。`.chk`、`.fch`、`.out`、`.com`はこの探索では検出しませんでした。自分の生成物`irc_eda/`、vendor、仮想環境、Git内部は対象外です。

IRCの本解析対象は設定で指定した`data/completed_20261017/IRC/`の4ログと、対応する`TS/`の4ログです。全系の電子条件はM06-2X/6-311++G(d,p)、全電荷0、一重項、SMD水です。

| 系 | forward完了点 | reverse完了点 | 明示TS | 合計 | 元ログの終了情報 |
|---|---:|---:|---:|---:|---|
| MoriPF | 35 | 28 | 1 | 64 | 異常終了1件。reverse 29は完了記録なし |
| RAFA14 | 24 | 35 | 1 | 60 | 合成ログに終了行なし。冒頭にforward 25の失敗除外の注記あり |
| RAFA15 | 35 | 21 | 1 | 57 | 異常終了1件。reverse 22は完了記録なし |
| RAFA16 | 35 | 35 | 1 | 71 | 正常終了行1件。ただし端点極小の証明ではない |
| 合計 | 129 | 119 | 4 | **252** | |

収束した部分経路は、後続点の失敗と区別して保存しました。欠落点・未収束点を補間していません。TSの座標は各TSログから取得し、IRCログの`Energy From Chk`、電子条件、原子順序、隣接点との距離を確認しています。

非TS 248点の座標とSCFエネルギーはcclibと照合しました。生成XYZの原子順序・座標を再読み取りし、252点の出力座標差は0 Åでした。既存`analysis4_force/outputs/csv/reaction_force_*.csv`に対して、反応座標差・電子エネルギー差の最大値は全4系で0でした。[数値比較表](project_validation.csv)。これは同じ元データの抽出・出力整合性の確認であり、別の量子化学計算による独立検証ではありません。

HPC補正過程では、最終SCFで評価された構造と、補間された補正端点が概念上区別されます。Adapterは印字された同一SCFの座標とエネルギーを保持しています。`Corrected End Point Energy`と採用SCFエネルギーの最大差はMoriPFで1.0800002×10⁻⁶ Hartree（約0.000678 kcal/mol）、RAFA14で5.1000006×10⁻⁷、RAFA15で4.9999994×10⁻⁷でした。RAFA16ではこの差を確認できる印字値がなく、未確認を0と解釈してはいけません。各点の差は[irc_profile.csv](../outputs/project/irc_profile.csv)で、未出力値は空欄です。印字されていない補正端点のCartesian座標は再構成していません。

反応座標は`√amu·bohr`、Cartesian座標はÅです。今回の符号設定はMoriPF/RAFA15/RAFA16が元ログのforwardを負側、RAFA14がreverseを負側に置きます。参照点はそれぞれ設定に明記した反応物側IRC端点です。[図1](../outputs/project/01_irc_profile.png)の端点基準エネルギー差は、分離反応物基準の活性化障壁や自由エネルギー障壁ではありません。

## 既存fchkの照合

4つのfchkファイル記録は、実質2つのTS波動関数とそれぞれのgzip版です。

| fchk | 幾何対応するIRC点 | 原子数 | 原子順序固定の整列RMSD / Å | fchk − TSログのエネルギー / Hartree |
|---|---|---:|---:|---:|
| `RAFA14_TS2_NBO.fchk` | `RAFA14_p0_0000` | 57 | 4.76999×10⁻⁷ | +1.80194×10⁻⁹ |
| `RAFA15_TS6_NBO.fchk` | `RAFA15_p0_0000` | 57 | 4.64282×10⁻⁷ | −7.46923×10⁻¹⁰ |

この数値は原子順序、幾何、エネルギーが非常によく一致する観測です。一方、fchkのrouteは`m062x/chkbasis`、ヘッダは`Gen`であり、基底名を単独で確定する根拠は不足しています。したがって`GEOMETRY_MATCH_CONDITIONS_UNVERIFIED`とし、名前やエネルギー一致だけで`READY`にしていません。基底関数879個とシェル・指数・係数の指紋も保存しましたが、関数数の一致だけは基底の完全一致の証明になりません。

残り250点には、照合できる波動関数がありません。波動関数の必要性だけを見ると250点が`RECALC_REQUIRED`、上記2点が条件未確定による`BLOCKED`です。ソフトと手法を含めた全体判定はIQA・EDAとも252点すべて`BLOCKED`です。[状態一覧](../outputs/project/readiness.csv)。

固定座標SP入力は252件生成しました。250件は波動関数欠如のため、2件は既存fchkの基底同一性を解決できない場合の代替準備です。既存2件については、基底情報の確認で再計算を避けられる可能性があります。入力は[prepared/](../outputs/project/prepared/)の`single_point.gjf`、一覧は[preparation.csv](../outputs/project/preparation.csv)です。元の条件を保持し、座標最適化は要求していません。

**これらのSPを実行しても、次のM06-2X非対応と溶媒分解問題は解決しません。今回は追加Gaussian計算を1件も実行していません。**

## 実行を止めた根拠

1. **APOST3Dのmeta-GGA未実装。** 固定commit `3ba2d95b58831d7a6adaaedc40cc4ed950eb16c5`の[enpart_dft.f](../vendor/APOST3D/sources/enpart_dft.f)では、`XC_FAMILY_MGGA` / `XC_FAMILY_HYB_MGGA`に対し1113–1123行付近で開発中と出力し、`stop`します。評価ルーチン呼び出しもコメント化されています。M06-2XのLibxc ID指定だけでは解決せず、B3LYPへの置換もしていません。
2. **溶媒条件の分解が未検証。** IQAにSMD反応場・溶媒和寄与をどのように含めるか、今回の公式実装での妥当性確認ができていません。easyEDAの[EDA_run.sh](../vendor/easyEDA/source/EDA_run.sh)はGaussianのSCF・L608項からsobEDA成分を算出し、独立したSMD/PCM項を扱う処理がありません。したがって「SMD計算もこのまま正しく分解できる」とは認定していません。これは未検証に基づく停止で、あらゆる溶媒EDAが理論的に不可能という主張ではありません。
3. **実行環境。** 現在はmacOS/Darwin arm64です。PATH、明示設定、ローカルvendorビルド、代表的な設置場所を確認しましたが、g16、formchk、unfchk、Multiwfn、ifort、ifx、icx、gfortran、apost3d、get_energy_g16、dos2unixは見つかりませんでした。リモート計算ホストは設定されていません。GNU sedでない点もsobEDA実行の制約です。[環境調査](../outputs/project/environment.json)。
4. **APOST3Dのビルド資材。** READMEの手順はoneAPI 2024のifortと対応Libxcを必要とします。さらにcloneしたcommitには`compile_libxc.sh`が要求する`libxc-4.2.3.tar.gz`がありません。今回のビルド試行はifort不足の事前検査で終了し、コンパイル完了とは報告していません。
5. **fragment未指定。** 各系の分割、電荷、多重度は推定していません。[編集用テンプレート](../outputs/project/fragments.template.yaml)を生成しました。指定後も1–4の該当問題は残ります。

## 出力仕様とテスト

公式easyEDAの保存済み`vendor/easyEDA/scr/1.log`から以下を読めました。これは今回の4系の結果ではありません。

| 項 | kcal/mol |
|---|---:|
| E_tot | −127.49 |
| E_els | −141.30 |
| E_x | −76.83 |
| E_rep | 287.00 |
| E_xrep = E_x + E_rep | 210.17 |
| E_orb | −172.95 |
| E_DFTc | −20.38 |
| E_dc | −3.03 |
| E_c = E_DFTc + E_dc | −23.41 |

非重複6項の和とE_totの差は印字精度で0.00 kcal/molです。E_xrepやE_cをさらに加算すると二重計上になるため、別項として保持してもclosureには重複させません。公式の`EDA_SUMMARY.log`で`E_Pauli`と呼ばれる列は、コード上はE_xrepです。AdapterはE_repと区別します。各段階のHartree値も抽出します。一般的な別流派のEDA項へ改名していません。

20テストが通過しました。[tests.log](tests.log)。主な検証は、途中SCF構造を点としないこと、未収束末尾の除外、重複不一致の拒否、別ファイルreverseの方向、元順序の保持、SP条件の保存、meta-GGA/SMD拒否、fragmentの被覆・電荷・スピン、公式sobEDA出力9項、APOST3Dのブロック行列と一意ペア和、失敗記録、成功済みresumeと出力改変検出、4実系のプロファイル一致、TS fchkの内容照合です。

APOST3Dの行列パーサーと実行制御の一部は**明示的な合成テストfixture**で検証しています。これをAPOST3D実機実行、公式サンプル成功、実データIQA成功と同一視していません。公式APOST3Dの実機出力形式との最終確認は未完了です。

`iqa_atomic.csv`、`iqa_pairwise.csv`、`eda_components.csv`は列定義のみで数値行なしです。未計算を0と記録していません。[図2](../outputs/project/02_iqa_interactions.png)と[図3](../outputs/project/03_eda_components.png)は「検証済み数値なし」と明示した表示です。実数値のプロファイルは図1のみです。

## 続行に必要なもの

現在のM06-2X/SMD条件を維持するには、まずその条件を正しく扱える公式実装・方法上の根拠が必要です。単に実行ファイルを追加するだけではIQA/EDAは開始しません。互換性が確認できた公式実装を採用する場合には、固定commit・対応判定・出力パーサーを再検証します。

HF/B3LYP等の対応する気相IRCへ本CLIを再利用する場合は、Linux/WSLのGaussianライセンス環境・Multiwfn・APOST3Dを設定し、fragmentを明示して公式サンプルとpilotを実行できます。成功済みpilotの条件・コード・出力hashを確認できた場合にのみ`scope: all`を許可します。現行4系を別手法や気相へ変更した結果は、元のM06-2X/SMDと同じ計算結果として比較してはいけません。
