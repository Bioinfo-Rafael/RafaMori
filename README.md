# MoriPF・RAFA14・RAFA15・RAFA16のIRC／反応力解析

GaussianのIRC（固有反応座標）計算ログから、フッ素置換位置が異なる4分子系の反応経路を比較するPythonプロジェクトです。エネルギー、反応力、原子間距離・角度、分子断片の質量重心（COM）、原子の力と運動を解析します。

**最新の結果：[4系解析の日本語summary](analysis4_force/outputs/reports/summary_ja.md)**

## 今回の比較でわかったこと

2026年10月17日に受領した `完成.zip` の16ログ（4系×IRC・反応物・生成物・TS）と `Lys.log` を解析しました。RAFA16を基準とする電子エネルギーの相対障壁は次の通りです。

| 分子系 | ortho F数 | 相対電子障壁（kcal/mol） | TSでの最短アミンH···F距離（Å） |
|---|---:|---:|---:|
| MoriPF | 2 | −0.763 | 2.358 |
| RAFA14 | 0 | +1.055 | 4.795 |
| RAFA15 | 1 | −4.008 | 2.622 |
| RAFA16 | 2 | 0.000 | 2.345 |

RAFA14とRAFA16の比較は「ortho FによるアミンHとの近接配置が反応に関係する」という仮説と整合します。一方、RAFA15は接触がより長いのに相対障壁が低く、H···F距離だけでは4系の障壁差を説明できません。接触の相関は水素結合や安定化エネルギーの直接的な証明ではありません。

`Lys.log` は他のログと計算条件が異なり、freqもありません。絶対障壁は報告せず、共通のLys項が消える系間差を比較しています。RAFA15反応物には微小虚振動があり、配座と生成物対応にも確認事項があります。詳細はsummaryを参照してください。

![相対障壁](analysis4_force/outputs/figures/02_relative_barriers.png)

## ディレクトリ構成

| ディレクトリ | 内容 |
|---|---|
| `analysis4_force/` | 最新の4系解析、設定、検証コード |
| `analysis4_force/outputs/` | 4系のCSV、グラフ、方向を統一した解析用LOG、日本語summary |
| `data/completed_20261017/` | 今回の17入力ログの圧縮原本（`.gz`）とSHA-256 manifest |
| `base_analysis/` | 従来の2系の基礎解析。共通のIRC parserもここに保存 |
| `analysis2_force/` | 従来の2系の反応力解析。共通の構造・力の解析関数を再利用 |
| `src/` | 従来の基礎解析を起動する入口 |

従来の `analysis2_force/outputs/` はローカルに残し、Gitの追跡対象から除外しています。

## 実行方法

Pythonの仮想環境を有効化し、プロジェクトのルートで実行します。

```sh
python -m pip install -r requirements.txt
python analysis4_force/run_analysis.py
```

圧縮原本は自動展開され、SHA-256で照合されます。展開後のLOGはGitで追跡せず、元の内容を保持した `.gz` を保存しています。再実行時にDownloadsフォルダは不要です。

検証：

```sh
python -m unittest discover -s analysis4_force/tests -v
```

従来の2系解析の実行方法：

```sh
python src/run_analysis.py
python analysis2_force/run_analysis.py
```

後者は別途Downloads内の反応物・生成物ログを必要とします。[旧解析README](analysis2_force/README.md)を参照してください。

## IRCの方向とTS

C20–求核Nの距離が長い側から短い側へ進むように、各系の方向を検証しました。MoriPF・RAFA15・RAFA16は元のparser座標を反転し、RAFA14はそのまま使っています。RAFA16には小さな局所的伸長がありますが、端点間の方向は明瞭です。

`outputs/oriented_logs/` のLOGは、原本から抽出した座標・力・エネルギーの**解析用派生ファイル**です。新規Gaussian計算の出力ではありません。解析はこのファイルを読み直して実施します。実際のTS構造・エネルギーを対応するfreqログから読み込み、未収束で反応座標が欠けた末端点は除外します。

## 理論の参照元

Toro-Labbé, Gutiérrez-Oliva, Murray & Politzer (2009), “The reaction force and the transition region of a reaction,” *Journal of Molecular Modeling* **15**, 707–710. [DOI](https://doi.org/10.1007/s00894-008-0431-8)、[PubMed](https://pubmed.ncbi.nlm.nih.gov/19085021/)。

反応力は `F = −dE/dR`、論文の反応力定数は `κ = d²E/dR² = −dF/dR` です。旧解析のκは逆符号だったため、4系解析では論文の規約を図と `kappa_paper_*` 列に使用します。

構造変化の統計、原子・分子断片の力と運動の積、H···F接触の相関は本プロジェクトの追加解析です。

## 詳しい資料

- [4系の結果と解釈・制限](analysis4_force/outputs/reports/summary_ja.md)
- [4系解析の実装・出力ガイド](analysis4_force/README.md)
- [基礎解析の説明](base_analysis/README.md)
- [旧2系の実装記録](analysis2_force/README_ANALYSIS.md)
