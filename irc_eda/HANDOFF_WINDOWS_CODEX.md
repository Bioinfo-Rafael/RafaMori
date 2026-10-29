# Windows PCのCodexへの引継ぎ指示

このWindows PCにはWindows版Gaussian 16がインストールされています。GitHubの `https://github.com/Bioinfo-Rafael/RafaMori` にあるGaussian IRC → IQA / sobEDA Adapterの作業を引き継ぎ、実行可能な処理を実際に進めてください。計画だけの回答では終わらせないでください。

## 最初に読むもの

1. リポジトリのAGENTS.mdなど、このPCで適用される作業規則。
2. `irc_eda/README.md`
3. `irc_eda/WINDOWS.md`
4. `irc_eda/reports/execution_report_ja.md`
5. `irc_eda/config/project.yaml`、`irc_eda/upstream-lock.json`

未取得ならcloneしてください。既存checkoutなら未コミット変更を保護し、fetchして対象ブランチを確認してください。共有履歴を書き換えたり、他の研究作業の変更を破棄したりしないでください。

## 前PCで確認できた事実

- Adapter・CLI・入力生成・結果パーサー・resumeを実装済み。20テスト通過。
- 実IRCはMoriPF 64、RAFA14 60、RAFA15 57、RAFA16 71、計252点（各TSを含む）。非TS 248点はcclibと座標・SCFエネルギーを照合済み。
- 原データは `data/completed_20261017/` のgzip圧縮版。必要なIRC/TSログ8件とNBOログ/fchk4件はリポジトリに含まれる。`irc_eda/source-data-manifest.json`に相対パスとSHA-256を記録。
- RAFA14/RAFA15のNBO fchkはそれぞれTSと原子順序・座標・エネルギーが一致する。ただしChkBasis/Genなので基底名の同一性は未確定。
- 他の250点には対応する波動関数が見つからなかった。SP入力は、上記2点の代替準備も含め252件生成可能。
- 前PCはmacOSでGaussian・Multiwfn・ifort・APOST3D実行ファイルなし。新規IQA/EDA計算成功は各0点。公式サンプルの新規計算も未実施。
- easyEDA同梱 `scr/1.log` の9成分とclosureを検証済み。ただし保存済み出力の読み取りであり、こちらで計算した結果ではない。

## このPCで実施すること

1. OS、Python、Git、Gaussian 16の実体・版・ライセンス環境、formchk/unfchk、Multiwfn、WSL/Linux環境、APOST3Dとコンパイラの有無を確認して記録する。Windows版Gaussianが存在することを、WSL内でLinux版Gaussianが利用可能という意味に読み替えない。
2. 必要データのSHA-256をmanifestと照合する。元ログ・チェックポイント・圧縮原本は上書き・削除しない。
3. `irc_eda/scripts/setup_windows.ps1`で専用Python環境と固定commitの外部実装を取得する。スクリプトはWindows実機未検証なので、実際に確認し、必要な互換性修正を行う。Gaussian等の有償ソフトを無断で別ソフトへ置換しない。
4. `config/project.yaml`を同じディレクトリの`config/local.yaml`へコピーし、確認した実行ファイルパスを設定する。local.yamlはGit対象外。inspect → prepare → analyzeを実行して、252点の抽出・元順序・エネルギー・XYZ出力を再確認する。実入力テストが未展開LOG不在を理由にskipした場合も、それだけで実データ検証済みとせず、gzip原本を読むCLIの結果を確認する。
5. 以下の手法上の障害を、公式文書・ソースコードを根拠に検討する。前PCのBLOCKED判定を機械的に繰り返すだけでなく、このPCの実環境で解決できることを進める。ただし根拠なく判定を解除しない。
6. 元条件で実行可能な経路が確認できたら、公式サンプル → 両端・TSの少数点 → 数値closureと元IRCエネルギー照合の順に実行する。全点の大量計算はpilot成功後。既存TS波動関数を優先し、必要性の確定前に252件のSPを一括投入しない。
7. Windows版Gaussianの起動インターフェースは実機で確認する。stdin方式やLinux版と同じ動作を仮定しない。必要な固定座標SPは、独立作業ディレクトリで1点ずつ試験し、入力条件・原子順序・エネルギーを確認してから拡大する。
8. エンジンを実行できなくても、環境調査、データ照合、入力準備、互換性修正、テストは完了させる。新たに判明した未解決事項と次に必要な判断を報告する。

## 解除してはいけない前提条件

- 現行データは **M06-2X/6-311++G(d,p)、電荷0、一重項、SMD水**。固定commitのAPOST3Dはmeta-GGA/hybrid-meta-GGA選択時に開発中として停止する。該当根拠は `sources/enpart_dft.f`。Gaussianが利用可能になっても、この非対応は解消しない。
- APOST3DのENPART/TFVCをIQAに使う。EDAIQAとsobEDAは別手法。IQA/EDA数値計算アルゴリズムは自作しない。
- M06-2XをB3LYPへ、SMDを気相へ黙って置き換えない。公式の新しい対応実装等を調査することは可能だが、手法・基底・溶媒・参照状態が変わる計算の採用は、違いを説明してユーザーと合意する。異なる条件の結果を元の結果と同一扱いしない。
- sobEDAのSMD寄与の扱いは未検証。独立した溶媒項を持たない元実装へSMD総エネルギーを単に混ぜない。正しい扱いを確認できなければ明示停止する。
- 現在のsobEDA実行制御はLinux/WSL＋GNUツール群を前提とする。Windowsネイティブ対応やWindows GaussianとWSLの橋渡しは未実装。利用可能な環境を調べ、実装するなら公式計算本体を変えずAdapterとして検証する。
- fragmentの原子集合・電荷・スピン多重度はユーザー未指定。テンプレートを用意し、化学的根拠なく推測して実行しない。必要な定義を尋ねても、独立したIQA/環境確認作業は続ける。
- SPでは元の固定座標を使い、最適化しない。各入力のchk名が同じため、点ごとに作業ディレクトリを分ける。
- 未実行・エラー・未出力を0と記録しない。成功済みresume、入力/出力hash、元ファイル、条件、外部ソフト版、Git commit、closureを記録する。

## 最終報告

「環境構築」「Adapter検証」「公式出力読み取り」「公式サンプルの新規計算」「実データIQA成功」「実データEDA成功」を分けて、実行点数と保存先を報告する。追加SPについて必要性・対象点・生成件数・実行件数を明示する。科学的な観測、算術的な整合性、未検証の仮説を区別する。

作業記録とWindows対応の変更は `irc_eda/` 内にまとめ、既存の他の研究解析には不要な変更を加えない。
