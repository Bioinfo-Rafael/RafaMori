"""Postprocessing only: validate and compare RAFA14/MoriPF; never launch Gaussian."""
import csv
import json
import re
from decimal import Decimal
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from irc_eda.scripts.smd_distortion_pilot import OUT,REPORT,analyze,validate,save
from irc_eda.parsing import HARTREE_KCAL,read_text

def main():
    systems=['RAFA14','MoriPF'];allrows=[];checks=[];barriers=[]
    for system in systems:
        manifest=json.loads((OUT/system/'manifest.json').read_text())
        for job in manifest['jobs']:
            result=validate(job)
            folder=ROOT/job['folder'];record=json.loads((folder/'validation.json').read_text())
            assert record['status']=='SUCCESS'
            assert record['output_sha256']==result['output_sha256']
            text=read_text(folder/'fragment.out')
            energy=Decimal(re.findall(r'SCF Done:.*?=\s*([-\d.]+)',text)[-1])
            assert abs(energy-Decimal(str(result['energy_hartree'])))<Decimal('1e-9')
            checks.append({'system':system,'role':job['role'],'fragment':job['fragment'],**result})
        analyze(system)
        rows=json.loads((REPORT/system/'pilot.json').read_text())
        for row in rows:
            f=Decimal(str(row['full_energy_hartree']))
            fragments=Decimal(str(row['substrate_hartree']))+Decimal(str(row['lys_hartree']))
            refs=sum(Decimal(str(v['energy_hartree'])) for v in manifest['references'].values())
            assert (f-refs)-((fragments-refs)+(f-fragments))==0
            assert abs(row['distortion_kcal_mol']-row['substrate_distortion_kcal_mol']-row['lys_distortion_kcal_mol'])<1e-6
            allrows.append({'system':system,**row})
        barriers.append({'system':system,**{col:rows[1][col]-rows[0][col] for col in ['distortion_kcal_mol','interaction_difference_kcal_mol','total_relative_kcal_mol']}})
    dest=REPORT/'RAFA14_MoriPF';dest.mkdir(exist_ok=True)
    for name,rows in [('pilot_comparison',allrows),('endpoint_to_TS',barriers)]:
        with (dest/(name+'.csv')).open('w',newline='',encoding='utf8') as f:
            w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    save(dest/'validation.json',{'successful_fragment_outputs':len(checks),'checks':checks,'decimal_identity':'exactly zero for all 6 points using decimal representations','fragment_distortion_sum_tolerance_kcal_mol':1e-6,'gaussian_calculations_launched_by_this_script':0})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),sharey=True)
    for ax,system in zip(axes,systems):
        rows=[r for r in allrows if r['system']==system]
        for col,label in [('distortion_kcal_mol','Distortion'),('interaction_difference_kcal_mol','Interaction difference'),('total_relative_kcal_mol','Total vs separated references')]:
            ax.plot([r['C20_N_distance_A'] for r in rows],[r[col] for r in rows],'o--',label=label)
        for r in rows:ax.annotate(r['role'].replace('_side',''),(r['C20_N_distance_A'],r['total_relative_kcal_mol']),xytext=(0,7),textcoords='offset points',fontsize=8)
        ax.invert_xaxis();ax.axhline(0,color='grey',lw=.5);ax.set(title=system,xlabel='Forming C-N distance / Angstrom')
    axes[0].set_ylabel('Energy / kcal mol$^{-1}$');axes[1].legend(fontsize=8)
    fig.suptitle('SMD fixed-fragment differences: 3 sampled points per system')
    fig.tight_layout();fig.savefig(dest/'comparison.png',dpi=200);plt.close(fig)
    text='# RAFA14・MoriPF：SMD固定断片の変形／相互作用差分解析\n\n'
    text+='対象は各3点（反応物側の利用可能IRC端点・TS・生成物側の利用可能IRC端点）。12断片出力の正常終了、SCF収束、手法・基底・SMD水、電荷0・一重項、電子数、原子順序、座標、入出力hashを再検証した。RAFA14は新規5件と既存TS Lys成功1件、MoriPFは新規6件。過去の失敗・中断試行は成功件数に含めない。今回の後処理によるGaussian追加起動は0件。\n\n'
    text+='M06-2X/6-311++G(d,p)、SMD水、元固定座標を維持。全系6値と各系の分離参照値は既存出力を再利用した。独立な参照構造は基質2種と共通Lysモデルの3構造。原本・参照情報は各系のsource_audit.json、今回の検証証跡はvalidation.jsonを参照。\n\n'
    text+='## 数値（kcal/mol）\n\n| 系 | 点 | 変形 | 相互作用差分 | 全系−分離参照 |\n|---|---|---:|---:|---:|\n'
    for r in allrows:text+=f"| {r['system']} | {r['role']} | {r['distortion_kcal_mol']:.3f} | {r['interaction_difference_kcal_mol']:.3f} | {r['total_relative_kcal_mol']:.3f} |\n"
    text+='\n## 反応物側IRC端点からTSへの変化\n\n| 系 | Δ変形 | Δ相互作用差分 | Δ全系 |\n|---|---:|---:|---:|\n'
    for r in barriers:text+=f"| {r['system']} | {r['distortion_kcal_mol']:.3f} | {r['interaction_difference_kcal_mol']:.3f} | {r['total_relative_kcal_mol']:.3f} |\n"
    text+='\n両系でC–N形成に伴う変形項増加を、負方向へ変化する相互作用差分が部分的に相殺する。TSの変形項はRAFA14 14.105、MoriPF 12.035 kcal/molで、基質寄与はそれぞれ13.734、11.657 kcal/mol、Lys寄与は0.371、0.378 kcal/mol。TSでの変形差約2.070 kcal/molは主に基質側にある。\n\n'
    text+='端点→TSの全系差はRAFA14 4.053、MoriPF 3.144 kcal/molで、MoriPFが約0.909 kcal/mol小さい。これはそれぞれの異なるIRC端点を起点とする有限差分であり、同じ反応座標における厳密な比較や、最適化反応物からの活性化自由エネルギーではない。TSのC–N距離も1.978 Å対2.031 Åと異なるため、MoriPFの変形コストが経路全体で一様に小さいとは結論しない。\n\n'
    text+='生成物側の変形は34.459／35.235、相互作用差分は−35.432／−38.816 kcal/mol。強く変形した仮想的な中性一重項断片という定義への依存が大きい。負の相互作用差分を純粋な静電・共有結合成分と解釈しない。\n\n'
    text+='## 整合性と解釈の限界\n\nE_int=E_full−ΣE_frag(R)、E_dist=Σ[E_frag(R)−E_frag(ref)]。全6点で E_full−ΣE_ref=E_int+E_dist の残差は保存精度で0 Eh。独立なDecimal演算でも一致し、断片別変形和も1e−6 kcal/mol以内で一致した。これは算術整合性であり、物理的妥当性や積分精度の証明ではない。エネルギー値の桁数はGaussianログの出力精度に制限される。\n\n'
    text+='SMDの反応場とCDSを含むSCFエネルギーを使用し、CDSを二重加算していない。全系と断片で空洞が変化するため溶媒寄与変化を含む差分である。BSSE未補正。ZPE・熱・標準状態補正なし。SCF安定性検査、別スピン状態、参照配座の網羅探索は未実施。原子対応と参照の局所極小確認、選択3点におけるN–H所属の保持は既存監査を参照。全IRCにおけるプロトン移動を否定する結果ではない。\n\n'
    text+='元の結合IRCログに正常終了がない場合、受理済み点の収束SCFを利用している。利用可能端点を最適化反応物／生成物とは呼ばない。図の破線は点を結ぶ目安で、中間点の計算結果ではない。系間の絶対エネルギーを直接差し引かず、それぞれの参照に対する差分だけを比較した。IQA・sobEDA解析ではない。\n\n'
    text+='この2系の3点pilotについて後処理を完了した。全252点への展開と他2系のGaussian起動は実施していない。定義と原著資料は[方法ノート](../smd_distortion_interaction_definition.md)ではなく、レポート親階層から見て `../../smd_distortion_interaction_definition.md` を参照。\n'
    text=text.replace('[方法ノート](../smd_distortion_interaction_definition.md)ではなく、レポート親階層から見て `../../smd_distortion_interaction_definition.md` を参照','[方法ノート](../../smd_distortion_interaction_definition.md)を参照')
    (dest/'report_ja.md').write_text(text,encoding='utf8')
    print('Validated 12 fragment outputs; 6-point comparison saved to',dest)

if __name__=='__main__':main()
