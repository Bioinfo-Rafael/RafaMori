# 別のWindows PCへの移行

このコードはGitHubリポジトリ [Bioinfo-Rafael/RafaMori](https://github.com/Bioinfo-Rafael/RafaMori) の `irc_eda/` に置いています。**ローカルでファイルを作っただけではGitHubへ反映されません。commit/pushした後にclone/pullするか、以下の持ち運び用ZIPを利用します。**

## Windows版Gaussianがある場合の範囲

| 処理 | 現時点での扱い |
|---|---|
| inventory、IRC抽出、XYZ・SP入力生成、集計・図 | Pythonで実行可能。今回の移行確認はmacOS上の別ディレクトリで行い、Windows実機は未検証 |
| Windows版GaussianによるSP | 生成したgjfをGaussianの通常のGUI/バッチ操作で実行可能。Gaussianのインストール先や起動方式はPCごとに確認する |
| APOST3Dの実行 | 実行バイナリ／公式ビルド環境が別途必要。取得した版のM06-2X非対応はPCを変えても解消しない |
| easyEDA/sobEDAの一括実行 | 現在のAdapterはLinux/WSL上のツール群を前提とし、Windowsネイティブ実行は未対応 |

Windows版Gaussianがインストールされていることと、WSL内でLinux版Gaussianが使えることは別です。WSLを入れただけでsobEDAの実行環境が揃ったとは扱いません。Windows版GaussianとLinux向けシェルの橋渡しも未実装です。

今回のM06-2X/SMD水のIQA/EDAは、ソフト不足以外に手法対応と溶媒分解の問題があります。まず環境を検出して結果を確認し、252件のSPを一括投入する前に必要性を判断してください。

## 持ち込み方法

GitHubに必要ファイルをpushした後なら、Git for Windowsで次のように取得します。

```powershell
git clone https://github.com/Bioinfo-Rafael/RafaMori.git
cd RafaMori
```

GitHubへ送る前でも、`irc_eda/outputs/transfer/RafaMori-irc-eda-transfer.zip`をコピーして展開すれば、独立した作業フォルダとして使えます。ZIPにはAdapter・設定・手順・元データ12圧縮ファイルを同梱します。Gaussian本体、仮想環境、マシン固有設定、過去の実行結果は含めません。APOST3D/easyEDAはセットアップ時に固定commitから取得します。GitHub取得用のネット接続が必要です。

ZIPの展開先で検証：

```powershell
py -3 .\irc_eda\scripts\verify_transfer.py
```

## セットアップと実行

Python 3.10以降とGit for Windowsを用意します。PowerShellでリポジトリまたはZIP展開先へ移動し、次を実行します。

```powershell
powershell -ExecutionPolicy Bypass -File .\irc_eda\scripts\setup_windows.ps1
.\irc_eda\scripts\run_windows.ps1 inspect
.\irc_eda\scripts\run_windows.ps1 prepare
.\irc_eda\scripts\run_windows.ps1 analyze
```

圧縮原本はそのまま読めるので、手動展開は不要です。`requirements-lock.txt`は今回のMac環境の記録です。Windowsでは`requirements.txt`から環境を構築します。`run_windows.ps1`はPythonのUTF-8モードを指定します。直接起動する場合も `python -X utf8 -m irc_eda ...` としてください。Windowsの既定コードページではUTF-8設定や出力の読み書きに失敗することがあります。

確認するファイル：

- `irc_eda/outputs/project/environment.json`：実行ソフトの検出結果。
- `irc_eda/outputs/project/readiness.csv`：点ごとの不足条件。
- `irc_eda/outputs/project/irc_profile.csv`：4系252点の抽出結果。
- `irc_eda/outputs/project/prepared/*/single_point.gjf`：固定座標SP入力。
- `irc_eda/outputs/project/fragments.template.yaml`：編集用fragment定義。

## マシン固有の設定

共有設定を直接書き換える代わりに、同じディレクトリへコピーします。`config/local.yaml`はGitと持ち運び用ZIPの対象外です。

```powershell
Copy-Item .\irc_eda\config\project.yaml .\irc_eda\config\local.yaml
```

`local.yaml`の`executables`に、実際に存在するファイルを指定できます。パスはYAMLの単一引用符で囲みます。以下は実在を確認していない例です。

```yaml
executables:
  g16: 'C:\G16W\g16.exe'
  formchk: 'C:\G16W\formchk.exe'
  unfchk: 'C:\G16W\unfchk.exe'
```

```powershell
.\irc_eda\scripts\run_windows.ps1 inspect -c .\irc_eda\config\local.yaml
```

この指定でWindows版Gaussianの無人起動方法まで検証済みになるわけではありません。現在のM06-2X/SMD対象では、実行前チェックが停止します。

## Gaussianで必要な点だけ実行する場合

実行が必要と判断した`prepared/<point_id>/single_point.gjf`を、点ごとに独立したフォルダへコピーして、Gaussianの通常の操作で実行します。各入力のチェックポイント名は`wavefunction.chk`なので、複数点を同じ作業フォルダに混在させないでください。最初は1点で正常終了、元IRCとの座標・原子順序・条件・エネルギーの照合を行います。

計算後はログとchk/fchkを点ごとのフォルダごと保存します。`formchk`が利用できれば、Gaussian付属の操作方法に従ってfchk化します。入力や元ログを上書きしないでください。生成した波動関数は有用ですが、M06-2XのAPOST3D非対応を解消するものではありません。

## 共有対象

GitHubへ反映すべきものは`irc_eda/`のコード・設定・スクリプト・説明と、未追跡だった次の4ファイルです。IRC/TSの原ログ圧縮版8件はすでにGit管理されています。

- `data/completed_20261017/TS/RAFA14_TS2_NBO.LOG.gz`
- `data/completed_20261017/TS/RAFA14_TS2_NBO.fchk.gz`
- `data/completed_20261017/TS/RAFA15_TS6_NBO.LOG.gz`
- `data/completed_20261017/TS/RAFA15_TS6_NBO.fchk.gz`

外部リポジトリ全体やソフト本体をこのリポジトリへ重ねて登録する必要はありません。固定commitと取得スクリプトで再現します。
