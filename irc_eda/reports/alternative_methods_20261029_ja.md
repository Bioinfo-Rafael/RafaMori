# 元条件を保った IQA / EDA の代替実装調査

調査日: 2026-10-29 JST。対象: M06-2X/6-311++G(d,p)、SMD 水、電荷0・一重項。Windows Gaussian 16 C.01 を使用。Windows の実行記録は commit `d100aedd51931971ba789cf15a4fd22d726b2dfe` として origin/main に push 済み。本記録は、その後の調査・基底照合を追加する。

## 結論と確認の段階

**既存 TS fchk を最小の追加計算で解析する候補は AIMAll、溶媒を含む相互作用エネルギーを分解する候補は Q-Chem ALMO-EDA(solv)。** 前者は SMD 分極密度上の溶質内部エネルギーに適用範囲を限定する。後者は別エンジンでの SCF が必要であり、Gaussian の元エネルギーとの一致を未検証のまま仮定しない。

文書上の対応を確認したことと、この TS で計算成功したことは異なる。AIMAll・Q-Chem の利用可能な実行環境は探索範囲内で未検出。今回の公式サンプル新規計算、TS IQA、TS EDA の成功数はいずれも **0**。入力準備は行ったが、実行成功とは数えない。252件のSPは投入していない。汎関数・溶媒の置換も行っていない。

## 候補の比較

| 既存実装 | M06-2X | SMD 波動関数・溶媒寄与 | 実行環境 | 既存 TS fchk | 追加 Gaussian | 分解量と IRC 比較の範囲 |
|---|---|---|---|---|---|---|
| AIMAll 19.10.12 IQA | 公式に対応。ただし原子間 XC の配分に規約あり | SMD で分極した密度を保持。明示的な反応場/CDS 分解は確認できない | Windows native、Linux、macOS | 入力可能。モデル認識と積分精度の検証が必要 | RAFA14 では原則不要。入力不足時だけ再検討 | 原子別の溶質内部エネルギーとその成分。Gaussian SMD 全エネルギーの分解ではない |
| Q-Chem 6.3 ALMO-EDA(solv) | EDA2 に M06-2X 公式例。DFT+SMD 対応。今回の組合せは実行未検証 | 各中間状態を溶媒中で評価し、静電・非静電溶媒効果を扱う | 有効なライセンスと配布版が必要。公式インストール文書は Windows/Linux/macOS を挙げる | 座標・比較基準として利用。既存密度だけでは EDA 完結不可 | 不要。ただし Q-Chem 内で全系・断片・拘束 SCF が必要 | FRZ/POL/CT、凍結項の細分。溶液中の断片相互作用エネルギー。IRC 差には断片変形分も必要 |
| BDF sobEDA / sobEDAw | 実装の対応表には含まれるが、物理解釈・パラメータ条件を満たすとは限らない | 明示的な implicit solvation energy を実装。SMD 非静電項にも言及 | 公式手順は Linux/WSL。BDF 本体の入手・導入が必要 | 完了 EDA の入力として直接流用できる根拠なし。座標を使って再計算 | 原則不要、BDF 内の計算は必要 | 静電・交換・反発・軌道・DFT相関・分散補正・溶媒項。現条件での採用は保留 |
| PROMOLDEN の COSMO-IQA 拡張 | 確認論文の検証は HF。M06-2X の同等実装は未確認 | COSMO 表面電荷による静電溶媒寄与。SMD CDS までの対応根拠なし | 著者実装・配布条件の確認が必要。Windows 可否未確認 | fchk 単独で成立する根拠なし。表面電荷などが必要 | 必要性未確定。少なくとも現 fchk のみでは不足 | 原子/群の静電溶媒寄与。現状では元 SMD 全エネルギーへの適用不可 |

### AIMAll: 既存の密度を利用する限定的な解析

[AIMQB 公式マニュアル](https://aim.tkgristmill.com/manual/aimqb/aimqb.html)は M062X の XC 汎関数を原子領域で積分することを明記している。ただし DFT の原子間 XC 配分は一意ではなく、実装は原子間部分に HF 交換式、原子内部分に残差を用いる。従って M06-2X の一意な「原子対相関エネルギー」とは解釈しない。領域は QTAIM であり、旧 Adapter の TFVC とは分割法が異なる。

記載された IQA 原子エネルギーの和には、明示的な PCM 反応場/CDS 項がない。このため採用範囲は `E_solute[rho_SMD]` の原子分解とする。密度を気相で再最適化する操作ではない。Gaussian SMD 全エネルギーとの差を、そのまま積分誤差または完全な溶媒分解と呼ばない。非静電項を引くだけで閉じるとも仮定しない。

fchk から WFX への変換を含め、モデルが M062X と認識されることを確認する。公式同梱サンプルの再計算を先に行い、その後 RAFA14 を `aimqb -nogui -nproc=1 -naat=1 -encomp=2 -wsp=false RAFA14_TS2_NBO.fchk` で原子エネルギーから検証する。これは未実行の予定コマンドであり、原子対計算はさらに精度・費用を確認してから行う。

[ライセンス](https://aim.tkgristmill.com/license.html)と[FAQ](https://aim.tkgristmill.com/faqs.html)によれば Standard は無料で、原子数等に関する制限は並列処理・可視化に関係する。57原子の逐次数値計算に Professional が必須とは判断しない。ただし[ダウンロード](https://aim.tkgristmill.com/download/download.html)は HTTP 401、[登録案内](https://aim.tkgristmill.com/register.html)は配布認証の取得を要求する。登録メール送信や資格情報の推測は行っていない。実行ファイルまたは正規取得済みインストーラが必要。

### Q-Chem: 溶媒を含む EDA の優先候補

[公式 SMD 対応節と Example 12.18](https://manual.q-chem.com/6.3/sect_eda2_sol.html)は、溶媒中の孤立断片・凍結・分極・全系状態を用いる ALMO-EDA(solv) を示す。[EDA2 の入力説明](https://manual.q-chem.com/6.3/sect_eda2_job.html)には M06-2X の例もある。これらは対応の根拠だが、M06-2X/6-311++G(d,p)/SMD の本 TS の成功実績ではない。

相互作用エネルギーは `E_full^s − Σ E_fragment^s = FRZ + POL + CT`。凍結項の静電相互作用には溶媒遮蔽を含み、非静電溶媒寄与も扱う。分散の溶媒遮蔽などには近似がある。詳細は [Mao らの原著、Chemical Science (2021), DOI 10.1039/D0SC05327A](https://doi.org/10.1039/D0SC05327A)。これは気相 EDA に溶媒補正を後付けしただけの方法とは区別する。

公式サンプルは water–Mg2+、wB97M-V/6-31+G(d)、SMD 水であり、サンプルの仕様どおり保存した。これは対象系の条件変更ではない。RAFA14 用入力は元条件を維持し、`EDA2=2`、`EDA_CLS_DISP=true`、BSSE 補正なしとした。`EDA2=2` は AO 部分空間を用いる分極で、FERF 分極とは定義が異なる。大きな基底での分極項の安定性も確認が必要。

6-311++G** と球面基底の指定は[基底一覧](https://manual.q-chem.com/6.3/A2.S3.SS1.html)・[PURECART 説明](https://manual.q-chem.com/6.1/Ch8.S4.SS2.html)に基づく。Gaussian 元ログは5D7F、879基底関数、独立873（6方向を除外）。Q-Chem の全系計算で、基底、線形従属の処理、積分グリッド、SMD 空洞/半径、SCF 収束を照合する。同じ理論レベル名でも同じ数値になるとは限らない。[導入条件](https://manual.q-chem.com/6.1/Ch2.S1.SS2.html)を確認したが、利用可能なライセンス・配布版は未確認。

EDA の相互作用差だけを IRC 障壁と比較してはいけない。固定した同じ断片定義・参照で、`ΔE_full = ΔE_interaction + Σ ΔE_fragment` として断片変形を加える。元の比較対象が SCF/SMD エネルギーである場合、ZPE・RRHO 熱補正は混ぜない。

### BDF: 溶媒項は存在するが、元条件の解釈に制約

[現行中国語公式ガイドの sobEDA 節](https://bdf-manual.readthedocs.io/zh-cn/latest/User%20Guide.html)には、溶媒項を持つ既存実装がある。対応汎関数表上 M062X は含まれる一方、同じ文書は sobEDA の意味ある利用を「分散を表現しない汎関数＋分散補正」の組合せに限定する。sobEDAw は理論レベル依存の3パラメータを必要とし、非共有結合解析向けである。M06-2X の本基底用パラメータや結合形成 TS への妥当性は確認できなかった。B3LYP 等のパラメータを転用しない。

従って「SMD を含む実装が存在しない」という従来の一般化は訂正できるが、「本 TS にそのまま妥当」とは言えない。[入力仕様](https://bdf-manual.readthedocs.io/zh-cn/latest/Input%20Library.html)、[導入手順](https://bdf-manual.readthedocs.io/zh-cn/latest/Installation.html)も参照。手法・汎関数を変えてこの経路を選ぶ場合は、変更の意味を説明して別途判断が必要。

### PROMOLDEN: 溶媒 IQA の先行実装としての位置付け

[Jiménez-Grávalos ら (2018), DOI 10.1002/cphc.201800733](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cphc.201800733)と[著者公開原稿](https://digibuo.uniovi.es/dspace/bitstream/handle/10651/50641/Interacting.pdf?sequence=1)は COSMO の静電溶媒エネルギーを原子/群に分割する。検証条件 COSMO-HF/aug-cc-pVTZ を M06-2X/SMD 対応の証拠に読み替えない。非静電 SMD 項と配布実装の適用条件が未解決なので、今回の実行候補にはしなかった。

## 実施した TS 検証・保存物

RAFA14 の元 fchk は Gen/ChkBasis 表記だった。Windows Gaussian で同じ原子順・座標から `M062X/6-311++G(d,p) Guess=(Only,Save,Huckel) NoSymm SCRF=(SMD,Solvent=Water)` を実行し、formchk で基底を比較した。**879基底関数、殻の種類・原始関数数・原子対応・指数・収縮係数の指紋が完全一致**。記録は [RAFA14_basis_audit_20261029.json](RAFA14_basis_audit_20261029.json)。この証拠は RAFA14 に限定し、RAFA15 に外挿しない。

最初の `Guess=Only` も正常終了したが必要な checkpoint 情報が保存されず、Save を付けて再実行した。両方とも SCF Done は0件。Hückel 初期推測の密度は解析に用いず、IQA には元の収束済み SMD fchk を使う。この2ジョブは基底構築検証であり、対象 SP 成功数には含めない。

準備物は Git 対象外の `irc_eda/outputs/alternatives/` に保存した。

- `aimall/RAFA14_TS2_NBO.fchk`: 元 gzip からの展開コピー。
- `qchem/official_example_12_18.in`: 公式 SMD サンプル入力。
- `qchem/RAFA14_TS_candidate.in`: 元 TS 座標の候補入力。未実行。
- `qchem/atom_order.json`: Gaussian 原子番号への対応。
- `RAFA14_basis_check_saved/`: 基底構築検証の入力・出力・chk/fchk。
- `sources/`: 閲覧文書のローカルキャッシュ。本文一式は commit しない。

断片は既存 `analysis4_force/config/analysis_targets.yaml` と原子対応監査、Lys/Reactant 元ログに基づく。RAFA14 は基質1–42,56,57、Lys43–55、それぞれ電荷0・一重項。Q-Chem 入力では基質を先に並べ、対応表を保存した。推定をユーザーへ丸投げせず、既存ファイルから根拠を確認した選択である。反応途中での断片参照の意味と SCF 状態は、実計算時にも点検する。

## 実行数と残る作業

| 区分 | 結果 |
|---|---|
| 元データの整合性、IRC 抽出・入力準備 | 12原本ハッシュ一致、252点準備済み（先行 Windows 記録） |
| Gaussian 起動試験 | H2 HF/STO-3G 1件成功。対象系の成功数とは別 |
| 今回の基底構築 | 2ジョブ正常終了、比較可能な fchk 取得1件、RAFA14 基底一致 |
| 対象 IRC の追加収束 SP | 0件 |
| 代替手法の公式サンプル新規実行 | 0件。Q-Chem 入力取得、AIMAll サンプルは配布物未取得 |
| 対象系 IQA / EDA 成功 | 各0件。数値結果・closure は未取得（0のエネルギーではない） |

AIMAll または Q-Chem の正規に利用可能な実行環境が必要。導入先の確認を依頼済み。AIMAll なら公式サンプル→既存 RAFA14 のモデル認識・積分・電子数を確認し、溶質内部エネルギーとして報告する。Q-Chem なら公式 SMD サンプル→RAFA14 全系エネルギー照合→EDA closure と断片参照を検証する。これらが成立するまで252件の追加計算を開始しない。

ソース取得の URL・ハッシュは [alternative_sources_20261029.json](alternative_sources_20261029.json)、準備物のハッシュは [alternative_pilot_manifest_20261029.json](alternative_pilot_manifest_20261029.json)。HTTP 403 と記録された BDF の直接取得は失敗したが、同 URL の内容は Web 閲覧ツールで確認した。商用実装の非公開ソースを確認したとは主張せず、対応判断は上記公式文書・論文に基づく。
