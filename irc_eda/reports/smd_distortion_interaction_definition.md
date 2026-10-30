# SMD変形・相互作用差分エネルギー解析の定義

2026-10-29。IQA・sobEDAとは独立した固定断片の差分解析。M06-2X/6-311++G(d,p)、SMD水、各系の元IRC固定座標を維持する。

基質AとLysモデルBの原子集合・中性一重項を固定する。参照は既存の各分離反応物のSMD最適化構造（局所極小、全配座中の最小を保証しない）。各点の孤立断片は切り出した座標で再SCFするが幾何最適化しない。

E_int(R) = E_full(R) − E_A(R_A) − E_B(R_B)

E_dist(R) = E_A(R_A) − E_A(ref) + E_B(R_B) − E_B(ref)

E_full(R) − E_A(ref) − E_B(ref) = E_int(R) + E_dist(R)

各EにはGaussianのSMD SCF全エネルギー（反応場とCDSを含む）を使う。ログ中のCDSについて「included in total energy above」を確認し、二重加算しない。ZPE、RRHO、並進/回転/振動の熱補正・標準状態補正は加えない。溶媒モデル由来の自由エネルギー寄与を含むが、反応の完全なGibbs自由エネルギーではない。

全系と孤立断片では空洞・分極・非静電溶媒寄与が変わる。従ってE_intは純粋な溶質間相互作用でも静電項でもなく、**それぞれ独立に溶媒和した系の差分**。E_distにも変形に伴う溶媒和変化が含まれる。通常の溶液ASMで溶質と溶媒寄与を分離する方式とは異なる定義を明示する。

BSSEは未補正。断片自身の基底のみを使い、ghost原子を置かない。有限基底に由来するBSSEがE_intに残る。counterpoiseを後付けすると基底と空洞の参照定義が変わるため今回の3点pilotでは混在させない。将来のCP値は別列・別定義で検証が必要。

原子の所属は反応に沿って変更しない。C–N形成時の断片は仮想的な変形反応物であり、全系の原子電荷を断片の整数電荷に置き換える操作ではない。プロトン移動で親原子が変わる、断片内結合が大きく切れる、電子状態が不安定となる場合には、同じ分解の物理解釈を再検討する。SCF収束だけで電子状態の安定性を保証しない。

再構成恒等式の一致は算術と出力処理の検査であり、断片選択の一意性・BSSEの小ささ・反応機構・SCF状態の妥当性の証明ではない。各点の有限差分値をpilotとして評価し、全経路の単調性や速度定数まで推論しない。

参照資料：

- [van Zeist & Bickelhaupt, OBC 2010, DOI 10.1039/B926828F](https://doi.org/10.1039/B926828F)：反応座標に沿う歪みと相互作用の分解。
- [Bickelhaupt & Houk, Angew. Chem. Int. Ed. 2017, DOI 10.1002/anie.201701486](https://pmc.ncbi.nlm.nih.gov/articles/PMC5601271/)：参照反応物、歪み/相互作用、溶液中の取り扱い。
- [Marenich, Cramer & Truhlar, JPCB 2009, DOI 10.1021/jp810292n](https://doi.org/10.1021/jp810292n)：SMDの静電および空洞・分散・溶媒構造項。
- [Gaussian SCRF](https://gaussian.com/scrf/)、[Counterpoise](https://gaussian.com/counterpoise/)：今回Web取得は502。インストール付属公式ヘルプ `C:/g16w/help/solvation.htm` と実際のGaussian SMD出力を併用し、取得できなかったWeb本文を読めたとは扱わない。
