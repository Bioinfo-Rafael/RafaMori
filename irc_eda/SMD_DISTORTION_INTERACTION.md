# 固定断片SMD変形・相互作用差分pilot

2026-10-29 追加指示による更新：RAFA14完了を待たずMoriPFを開始する。既存RAFA14 workerを1枠として数え、もう1枠でMoriPFを実行する。空いた枠へRAFA15、RAFA16を順に投入し、常に最大2分子とする。表示されたWindows PowerShellから `scripts/continue_smd_pilots_visible.ps1 -ExistingWorkerPid <既存RAFA14 worker PID>` を実行する方式へ移行した。以下の「RAFA14報告後に他系開始」という旧記述はこの追加指示で置き換える。PIDは稼働中プロセスを確認して指定し、制御を重複起動しない。

2026-10-29 並列数の更新：ユーザー指定により、1計算4コア・4 GB、最大2分子を並列実行する制御に変更した（合計最大8コア・8 GB、Gaussian以外のメモリを除く）。RAFA14の既存workerは中断せず、その完了・検証・報告後にMoriPFとRAFA15を開始し、空いた枠でRAFA16を実行する。各分子内の6断片SPは逐次実行。以下の旧「全系逐次」の記述は、この更新により置き換える。新制御ログは `workflow_parallel_stdout.log` / `workflow_parallel_stderr.log`。失敗時には新たな計算を投入せず、既に動いている他方の計算の終了を待つ。252点への展開は行わない。

この経路はIQA/sobEDA CLIとは独立している。定義は [ノート](reports/smd_distortion_interaction_definition.md)。対象はRAFA14の3点を先行し、その報告後にMoriPF・RAFA15・RAFA16の各3点。全252点を実行する機能はない。

```powershell
# リポジトリルートから。既存 outputs/project/points.json が必要。
.\irc_eda\.venv\Scripts\python.exe -X utf8 irc_eda/scripts/smd_distortion_pilot.py prepare RAFA14
.\irc_eda\.venv\Scripts\python.exe -X utf8 irc_eda/scripts/smd_distortion_pilot.py run RAFA14
.\irc_eda\.venv\Scripts\python.exe -X utf8 irc_eda/scripts/smd_distortion_pilot.py analyze RAFA14
```

prepareは元ログhash、SCF値と対応座標、手法・基底・溶媒、電荷/多重度、参照の正の振動数、結合グラフと原子対応を確認する。2個の参照反応物と全系3点を再利用し、断片6件の固定座標SPのみを用意する。Gaussianの設定はWindows G16Wの実機パスを使用するこのPC用runnerである。

runはファイル名指定でGaussianを**逐次**起動する。1件4スレッド・4GB。正常終了、単一SCF、電子数、原子順と座標を検証し、validation.jsonを保存する。成功済みは入力/出力hashを照合して再利用する。未検証出力がある場合は上書きせず停止する。プロセス終了コードだけでは成功とみなさない。元ログを変更せず、Gaussian生成物はoutputs以下に分離する。

初回RAFA14で同時起動した3件が正常終了前にexit 1となった。失敗フォルダは `outputs/smd_distortion_interaction/RAFA14/` に保持し、manifestのprevious_attempt_folderに記録した。短いパス `outputs/smd_di_jobs/RAFA14/` で逐次再試行した。最初のTS Lys成功出力は再利用している。この特別な移行を行ったため、RAFA14では既存manifestを使ってrun/analyzeを再開し、prepareで再作成しない。

analyzeは成功した6件を再検証し、原子対分解を行わず断片エネルギーの差と和のみを計算する。`reports/smd_distortion_interaction/<system>/` にCSV/JSON/PNGを保存する。入力とGaussianログ本体はGit対象外。正常終了欠落の結合IRCについては受理済み点を利用したことをsource_auditに明示し、端点極小とは呼ばない。

SMDモデルによる差分、BSSE未補正、固定中性一重項の仮想断片という限界がある。SCF安定性解析・全配座探索は実施しない。恒等式のclosureは算術検査であり、これらの限界を解消しない。

## 継続実行

`scripts/continue_smd_pilots.py` を非表示の独立プロセスとして起動済み。先行RAFA14 workerの終了を待ち、正常出力を再検証してRAFA14報告を作成し、その後MoriPF→RAFA15→RAFA16を逐次実行する。再起動する場合は、既存の制御プロセス・Gaussianが動いていないことと `outputs/smd_distortion_interaction/workflow.lock` を確認する。稼働中に同じ制御スクリプトを重複起動しない。

進捗は `outputs/smd_distortion_interaction/workflow_status.json`、実行記録は同ディレクトリの `workflow_stdout.log` / `workflow_stderr.log`。各計算の `validation.json` にはRUNNING/SUCCESS/FAILEDを区別して保存する。全4系が完了すれば `reports/smd_distortion_interaction/all_pilots.csv` に12点を集約する。計算完了前にはこの集約CSVは存在しない。
