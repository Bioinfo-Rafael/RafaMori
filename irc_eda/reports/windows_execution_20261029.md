# Windows 実行記録（2026-10-29 JST）

対象 main 74fbe82。既存の別 checkout は変更せず、新規 clone。原本12件の SHA-256 は windows_source_integrity.json で全一致。作業用出力は ../outputs/project/ に保存。

- Python 3.12.6 専用 venv と固定 APOST3D 3ba2d95、easyEDA 08d7a7d を setup_windows.ps1 で取得。Windows のコードページで UTF-8 設定読み込みが失敗したため run_windows.ps1 に UTF-8 モードを指定。pytest は15 passed、5 skipped（未展開実データに依存）。実データは gzip 経由で別途検証。
- Gaussian 16W C.01: C:\g16w\g16.exe、formchk.exe、unfchk.exe。インストール記録と実行ログで版確認。専用フォルダの H2 HF/STO-3G 単点はファイル名指定で正常終了（../outputs/windows_gaussian_file_test2/smoke.out、SCF -1.11675930751 Hartree）。最初の整数座標入力は Gaussian が拒否した（../outputs/windows_gaussian_file_test/smoke.out）。stdin 試験は exit 1 で完了せず、起動方式の同等性は未確認。これはライセンス全般の保証ではない。
- RAFA14 既存 TS fchk のコピーを unfchk と formchk で往復し、fchk_info の座標・原子番号・エネルギー・基底指紋・条件が一致（../outputs/windows_wavefunction_check/）。変換結果は基底名の同一性を証明しない。
- Ubuntu および Ubuntu-20.04 の WSL2 で PATH 上の g16、formchk、unfchk、Multiwfn、apost3d、ifort、ifx、icx、gfortran、bc、dos2unix は見つからず。探索は PATH に限る。Windows 側 PATH には gfortran 13.2、GNU make/sed、dos2unix があり、Multiwfn、ifort、ifx、icx、apost3d は検出なし。公式 APOST3D ビルドの Intel コンパイラ環境は未確認。
- inspect/prepare: 64+60+57+71=252点。SP入力252件、実IRC SP実行0件。別途 windows_project_validation.csv は既存プロファイルとの座標・エネルギー差最大0、XYZ原子順序・座標差最大0、非TS 248点のcclib照合を記録。同一原本に対する抽出照合であり、独立計算ではない。
- 候補フラグメントは既存 analysis4_force/config/analysis_targets.yaml の amine_side と analysis4_force/outputs/reports/atom_correspondence_audit_ja.md の求核剤割当、基質+Lys基準の解析に基づく。MoriPF Lys 45–57、基質1–44。他3系 Lys 43–55、基質1–42と56–57。各0、一重項。Lys.log.gz と各 Reactant ログの元素組成・電荷・多重度はそれぞれ対応する部分集合と一致。局所設定として記録し、SMD未検証のため EDA実行には使っていない。
- 固定版および2026-10-29確認時の公式 APOST3D enpart_dft.f は MGGA/HYB_MGGA で stop。easyEDA EDA_run.sh に独立 SMD/PCM 項の処理を確認できず、正しい SMD 分解は未検証。元条件の IQA/EDA は実行不可と判定。
- 公式 easyEDA 保存済み出力読み取り1件、closure 0.00 kcal/mol。公式サンプル新規計算0件（APOST3D 実行ファイル、Multiwfn等が不足）。実データ IQA 0/252、EDA 0/252。追加 SP の必要点は手法対応の確定後に再判定。既存 TS 波動関数2点を優先し、基底同一性は未確定。
