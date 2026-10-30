"""Reanalyze tracked Gaussian archive and summarize all four 3-point pilots."""
import csv
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from irc_eda.scripts.archive_smd_pilots import verify,SYSTEMS,REPORT

def main():
    verify(reanalyze=True)
    rows=[];barriers=[]
    for system in SYSTEMS:
        points=json.loads((REPORT/system/'pilot.json').read_text())
        assert [p['role'] for p in points]==['reactant_side','TS','product_side']
        for p in points:
            assert abs(p['identity_residual_hartree'])<1e-10
            assert abs(p['distortion_kcal_mol']-p['substrate_distortion_kcal_mol']-p['lys_distortion_kcal_mol'])<1e-6
            rows.append({'system':system,**p})
        barriers.append({'system':system,**{k:points[1][k]-points[0][k] for k in ['distortion_kcal_mol','interaction_difference_kcal_mol','total_relative_kcal_mol']}})
    for name,data in [('all_pilots',rows),('all_endpoint_to_TS',barriers)]:
        with (REPORT/(name+'.csv')).open('w',newline='',encoding='utf8') as f:
            writer=csv.DictWriter(f,fieldnames=data[0]);writer.writeheader();writer.writerows(data)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8),sharey=True)
    for ax,system in zip(axes.flat,SYSTEMS):
        data=[r for r in rows if r['system']==system]
        for k,label in [('distortion_kcal_mol','Distortion'),('interaction_difference_kcal_mol','Interaction difference'),('total_relative_kcal_mol','Total vs separated references')]:
            ax.plot([r['C20_N_distance_A'] for r in data],[r[k] for r in data],'o--',label=label)
        ax.invert_xaxis();ax.axhline(0,color='grey',lw=.5)
        ax.set(title=system,xlabel='Forming C-N distance / Angstrom',ylabel='Energy / kcal mol$^{-1}$')
        ax.legend(fontsize=7)
    fig.suptitle('SMD fixed-fragment differences: three points per system')
    fig.tight_layout();fig.savefig(REPORT/'all_comparison.png',dpi=180);plt.close(fig)
    text='# 4分子のSMD固定断片・変形／相互作用差分pilot\n\n2026-10-30。RAFA15・RAFA16の各6件を追加検証し、RAFA14・MoriPFを含め計24件の成功断片出力を確認した。4系×3点、全系12値と各分離参照は既存値を再利用。今回Gaussianの追加起動は0件。\n\n'
    text+='## 反応物側IRC端点→TS（kcal/mol）\n\n| 系 | Δ変形 | Δ相互作用差分 | Δ全系 |\n|---|---:|---:|---:|\n'
    for r in barriers:text+=f"| {r['system']} | {r['distortion_kcal_mol']:.3f} | {r['interaction_difference_kcal_mol']:.3f} | {r['total_relative_kcal_mol']:.3f} |\n"
    text+='\n## RAFA15・RAFA16の3点（kcal/mol）\n\n| 系 | 点 | 変形 | 相互作用差分 | 全系−分離参照 |\n|---|---|---:|---:|---:|\n'
    for r in rows:
        if r['system'] in ['RAFA15','RAFA16']:text+=f"| {r['system']} | {r['role']} | {r['distortion_kcal_mol']:.3f} | {r['interaction_difference_kcal_mol']:.3f} | {r['total_relative_kcal_mol']:.3f} |\n"
    text+='\nRAFA15の端点→TSの変形増加はRAFA16より大きく、相互作用差分の負方向への変化は両者で近い。この定義でRAFA15の全系差が大きくなる主因は、相対的に大きい変形増加である。TSでの基質変形はRAFA15 14.747、RAFA16 12.336 kcal/mol、Lysモデルは0.482、0.402 kcal/molで、変形差は主に基質側にある。\n\n'
    text+='端点→TSの全系差はMoriPF < RAFA16 < RAFA14 < RAFA15の順。ただし起点は各経路の利用可能端点であり最適化反応物ではない。TSのC–N距離も異なり、この順位を反応速度や活性化自由エネルギーの順位と同一視しない。3点の間を結ぶ線は視覚的補助であり、未計算点の傾向を保証しない。\n\n'
    text+='## 検証と保存\n\n24件とも正常終了、単一SCF収束、M06-2X/6-311++G(d,p)、SMD水、中性一重項、電子数、原子順序・固定座標、入出力SHA-256を確認。全12点でE_full−ΣE_ref=E_int+E_distの残差は保存精度で0 Eh。これは算術検証であり、分解の一意性・物理的妥当性の証明ではない。\n\n'
    text+='[入力・圧縮出力アーカイブ](../../data/smd_pilots/README.md)に各24ファイルを保存。gzip展開後のバイト列と元出力のSHA-256が一致することを検証した。元outputsフォルダを使わず、アーカイブから一時フォルダに再構成して全系の解析を再生成した。検証記録はarchive_validation.json。入力の改行コードはGit属性で保持する。失敗・中断試行やchk/rwfは成功アーカイブに含めない。\n\n'
    text+='全系・孤立断片それぞれのSMD空洞と溶媒寄与を含む差分であり、純粋な溶質間静電相互作用ではない。BSSE未補正、ZPE・熱補正なし。SCF安定性、別スピン状態、参照配座の網羅探索は未検証。参照の原子対応・局所極小と元IRC受理済み点の根拠は各source_audit.jsonを参照。元IRC全体の正常終了とは区別する。IQA/sobEDAではなく、全252点への展開も実施していない。\n\n'
    text+='再現：`python -X utf8 irc_eda/scripts/summarize_smd_all_pilots.py`。Gaussian本体・元作業フォルダは不要。Python依存関係はirc_eda/requirements.txtを利用する。\n'
    (REPORT/'all_report_ja.md').write_text(text,encoding='utf8')
    print(json.dumps(barriers,indent=2))

if __name__=='__main__':main()
