# 固定断片SMD変形・相互作用差分pilot

2026-10-30：RAFA14・MoriPF・RAFA15・RAFA16の各3点pilotが完了。成功した断片計算24件を検証し、[4系比較報告](reports/smd_distortion_interaction/all_report_ja.md)、CSV、図を保存した。IQA/sobEDAとは独立した解析で、[方法の定義](reports/smd_distortion_interaction_definition.md)に従う。

## 別環境で結果を取得・再解析

成功24件の `.gjf` と `.out.gz` は [data/smd_pilots](data/smd_pilots/README.md) に追跡済み。従来の `outputs/` は引き続きローカル生成物用であり、Gitには含めない。入力・gzip・展開出力のhashをindex.jsonに記録し、入力の改行変換を無効化している。

Python環境に `requirements.txt` の依存関係を導入し、リポジトリルートから：

```text
python -X utf8 irc_eda/scripts/archive_smd_pilots.py verify
python -X utf8 irc_eda/scripts/summarize_smd_all_pilots.py
```

Gaussianや元outputsがない環境でも、アーカイブから検証・CSV・図・報告を再生成できる。追加Gaussian計算は起動しない。全系・参照エネルギーと出典は各manifest.jsonに保持し、元全系・参照ログは既存data配下にある。

## Windowsでの計算手順

このPCで作成済みのmanifestがある場合、PowerShellから手動起動する：

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force
.\irc_eda\scripts\run_smd_manual.ps1 -System RAFA15
```

各計算4コア・4 GB。最大2画面で異なる系を実行する。同一系を重複起動しない。Windows Gaussian `C:\g16w\g16.exe` をファイル名指定で起動する。独立作業フォルダ内で固定座標SPを逐次実行し、正常終了・SCF・電子数・条件・座標・原子順序を検証する。成功済みはhash照合して再利用し、未検証出力は上書きしない。別環境の後処理には上のアーカイブ用コマンドを使う。

初回の入力作成は `smd_distortion_pilot.py prepare <system>` で行い、既存 `outputs/project/points.json` が必要。元ログ・参照極小・原子対応を検査する。RAFA14には失敗試行からの移行履歴があるため、このPCのmanifestをprepareで上書きしない。

## 継続実行スクリプトと履歴

説明で言及していた `scripts/continue_smd_pilots.py` と `scripts/continue_smd_pilots_visible.ps1` を今回収録した。前者は既存RAFA14 workerを1枠として最大2系を処理し、後者は表示されたWindows PowerShellから起動するためのラッパー。Windows専用で、ローカルmanifest・Gaussianを必要とする。

これらは過去の継続実行手順を保存するためのもので、現在は起動していない。ユーザーによる手動実行へ切り替え済みで、完了結果の再解析には不要。過去のPIDやworkflow.lockを別環境で再利用しない。初期の失敗・中断は実行記録に保持し、成功24件とは区別する。

## 解釈の範囲

M06-2X/6-311++G(d,p)、SMD水、固定中性一重項断片。差分にはSMD空洞と溶媒寄与の変化を含み、BSSE未補正。ZPE・熱補正なし。SCF安定性や全配座探索は未実施。恒等式の一致は算術検査であり、物理的妥当性の証明ではない。端点は利用可能な受理済みIRC点で、最適化反応物・生成物を意味しない。全252点への展開はしていない。
