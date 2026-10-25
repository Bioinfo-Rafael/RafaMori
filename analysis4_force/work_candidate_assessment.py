"""Assess displayed union candidates without selecting thresholds to fit targets."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from numbered_structures import SYSTEMS

OUT=Path(__file__).resolve().parent/'outputs'


def main():
    from cross_system_work_tables import pair_key
    from higher_ranked_controls import table
    union=pd.read_csv(OUT/'csv/cross_system_work_union.csv').drop_duplicates(['kind','cross_key','system'])
    energy=pd.read_csv(OUT/'csv/energy_comparison.csv').set_index('system')
    pairs=pd.read_csv(OUT/'csv/contact_controls_pairs.csv')
    pairs=pairs[pairs.window=='common_CN_2.50_to_2.10']
    categories={pair_key(r.system,r.site_a,r.site_b):r.category for r in pairs.itertuples()}
    rows=[]
    for (kind,key),g in union.groupby(['kind','cross_key'],sort=False):
        g=g.set_index('system').loc[list(SYSTEMS)]
        x=g.work_kcal_mol.to_numpy();y=energy.loc[list(SYSTEMS),'absolute_RRHO_G_barrier_kcal_mol'].to_numpy()
        r=np.corrcoef(x,y)[0,1];loo=[np.corrcoef(np.delete(x,i),np.delete(y,i))[0,1] for i in range(4)]
        rows.append(dict(kind=kind,cross_key=key,domain='all_atoms' if kind=='atom' else categories[key],
            mean_work_kcal_mol=np.mean(x),mean_abs_work_kcal_mol=np.mean(abs(x)),
            std_work_kcal_mol=np.std(x),relative_variation=np.std(x)/abs(np.mean(x)),
            pearson_r=r,abs_r=abs(r),r_squared=r*r,loo_r_min=min(loo),loo_r_max=max(loo),
            r_without_RAFA15=loo[2],MoriPF_atoms=g.loc['MoriPF','atom_or_pair'],
            RAFA15_atoms=g.loc['RAFA15','atom_or_pair']))
    data=pd.DataFrame(rows)
    for i,t in data.iterrows():
        controls=data[(data.domain==t.domain)&(data.cross_key!=t.cross_key)]
        for metric,name in [('relative_variation','relative'),('std_work_kcal_mol','absolute')]:
            stronger=controls[(controls.mean_abs_work_kcal_mol>=t.mean_abs_work_kcal_mol)&
                              (controls[metric]>=t[metric])&(controls.abs_r>=t.abs_r)]
            data.loc[i,f'dominators_{name}_count']=len(stronger)
            data.loc[i,f'dominators_{name}_keys']=json.dumps(stronger.cross_key.tolist())
    data.to_csv(OUT/'csv/work_candidate_assessment.csv',index=False)
    chosen=[('contact','amineN||core20'),('atom','amineN'),('atom','core20'),
            ('atom','X24'),('atom','NH_near'),('atom','X25'),('atom','core21'),
            ('atom','core19'),('atom','core22'),('contact','NH_near||X24'),('contact','X25||core21')]
    names=['C20–求核N','求核N','反応C20','アミン側の置換F/H','選択アミンH',
           '反対側の置換F/H','反応カルボニルO21','ラクトンカルボニルO19','エステルO22',
           'アミン側F/H···選択アミンH','反対側F/H···O21']
    selected=[data[(data.kind==k)&(data.cross_key==key)].iloc[0] for k,key in chosen]
    by_key={(r.kind,r.cross_key):r for r in selected}
    def value(kind,key,field):return by_key[kind,key][field]
    text=f'''# 仕事の大きさ・系間変動・ΔG‡との相関から、どの候補が残るか

## 結論

**2026-10-25のRAFA15反応物更新により、仕事と正方向ΔG‡の相関が大きく変わった。O21とアミン側F/H原子は正相関、O19は負相関となる。旧解析の負相関に基づく「仕事が大きいほど低障壁」という説明は維持できない。** 新反応物の安定化で障壁が上がる一方、IRCの区間仕事は同じであり、基準配座への依存を明示して比較する。

今回の期待は「反応の中心となるC20–求核Nは大きな共通成分を担い、補助接触は系ごとの差を担う」というもの。前半は相対変動について部分的に合うが、C20–Nの変動そのものが小さいわけではない。以下はその違いを数字で確認した結果である。

## 比較する量

- 仕事の大きさM：4系の平均絶対値 mean(|W|)。単位はkcal/mol。符号の相殺で小さくならない尺度。
- 絶対的な変動σ：4系の母標準偏差。単位はkcal/mol。系間で仕事がどの程度変わったかを表す。
- 相対変動C：σ/|mean(W)|。無次元。既存表のSTD/meanの絶対値である。平均が小さいと大きくなるため、M・σも必ず一緒に見る。
- 相関：符号付き仕事Wと正方向ΔG‡のPearson r。相関の強さの比較には|r|を使うが、解釈では正負を残す。

原子仕事と接触方向の仕事は別の量である。全原子の仕事は区間ΔEへ足し合わせられるが、各原子対の接触方向の仕事を足してエネルギー分解にはできない。またMが大きいことを、その対だけの相互作用エネルギーが大きいこととは解釈しない。

'''
    text+=table(['対象','仕事の種類','M＝平均絶対値 / kcal/mol','σ / kcal/mol',
                 'C＝σ/絶対平均（無次元）','ΔG‡とのr（無次元）','R²（無次元）'],[
        [name,'原子仕事' if r.kind=='atom' else '接触方向の仕事',f'{r.mean_abs_work_kcal_mol:.6f}',
         f'{r.std_work_kcal_mol:.6f}',f'{r.relative_variation:.3f}',f'{r.pearson_r:.3f}',f'{r.r_squared:.3f}']
        for name,r in zip(names,selected)])+'\n'
    text+=f'''C20–求核Nの仕事は平均1.460 kcal/mol、σ=0.318 kcal/mol、相対変動0.218。対してO21の原子仕事はM=0.056649、σ=0.046715、相対変動0.825である。**O21は相対的には大きく変わるが、絶対的な変動量はC20–Nの方が大きい。** したがって「主反応は一定で、差をほぼ補助接触だけが担う」とは言えない。

アミン側F/H原子の相関はr={value('atom','X24','pearson_r'):+.3f}、求核Nはr={value('atom','amineN','pearson_r'):+.3f}、選択アミンHはr={value('atom','NH_near','pearson_r'):+.3f}。指定アミン側接触の仕事はr={value('contact','NH_near||X24','pearson_r'):+.3f}。**F側原子の仕事とF···アミンH接触の仕事は別の指標であり、F原子だけが障壁との強い対応を示すわけではない。**

反対側F/H···O21はC={value('contact','X25||core21','relative_variation'):.3f}、r={value('contact','X25||core21','pearson_r'):+.3f}。M={value('contact','X25||core21','mean_abs_work_kcal_mol'):.6f}、σ={value('contact','X25||core21','std_work_kcal_mol'):.6f} kcal/molと小さく、数kcal/molの障壁差を直接説明するエネルギー項ではない。

## 他の候補と同時に比較した場合

表示済みの和集合から重複を除いた原子40候補・分子間接触65候補・基質内接触26候補を使用する。全ての可能な接触を網羅した母集団ではない。原子は原子同士、分子間接触は分子間同士、基質内接触は基質内同士で比較する。

恣意的な閾値や重み付きの総合点を作らず、M・変動・|r|の3量がすべて指定対象以上の候補を数えた。相対変動Cで比べる場合と、絶対変動σで比べる場合の両方を示す。候補を残すために片方だけを採用しない。

'''
    text+=table(['対象','M・C・|r|がすべて以上の他候補数','M・σ・|r|がすべて以上の他候補数'],[
        [name,int(r.dominators_relative_count),int(r.dominators_absolute_count)] for name,r in zip(names,selected)])+'\n'
    text+=f'''O21のrは{value('atom','core21','pearson_r'):+.3f}、O19は{value('atom','core19','pearson_r'):+.3f}。新基準でO21をM・C・|r|のすべてで上回る原子は{int(value('atom','core21','dominators_relative_count'))}件、M・σ・|r|では{int(value('atom','core21','dominators_absolute_count'))}件となった。候補の優劣は新しい表の件数と具体的な対照名で評価し、旧障壁での順位を引き継がない。

全候補の値と、各対象を上回る候補名は[比較CSV](../csv/work_candidate_assessment.csv)に保存した。

## 4点の相関が特定の系に依存していないか

相関の安定性を記述するため、1系ずつ除いた3点のrも計算した。独立な検証ではない。

'''
    text+=table(['対象','1系除外時のr最小値（無次元）','1系除外時のr最大値（無次元）','RAFA15を除いたr（無次元）'],[
        [name,f'{r.loo_r_min:.3f}',f'{r.loo_r_max:.3f}',f'{r.r_without_RAFA15:.3f}'] for name,r in zip(names,selected)])+'\n'
    text+=f'''1系除外の表は、基準配座を更新しても4点相関が独立検証にはならないことを確認するための診断である。RAFA15を除いた3点のrは今回の更新では変わらない一方、RAFA15を含む相関は符号・大きさとも変わる。新反応物の虚振動に由来する留保は解消したが、配座集団・TSへの接続は未検証である。

## 「仕事が大きいほど反応しやすい」を安定化と呼べるか

この解析の正の仕事は、共通区間内の電子エネルギー上昇への寄与である。仕事と全障壁の相関が正でも負でも、そのまま特定原子対の安定化エネルギーとは読めない。今回、反応物基準だけで相関の符号が変わったことも、この区別が必要な理由である。

実際、全体の障壁には共通区間の外側の電子エネルギー差と熱補正も入る。次の3項の和がΔG‡になる。

'''
    windows=pd.read_csv(OUT/'csv/two_sided_window_metrics.csv')
    windows=windows[(windows.window=='common_CN_2.50_to_2.10')&(windows.H_choice=='primary')].set_index('system')
    erows=[]
    for system in SYSTEMS:
        e=energy.loc[system];w=windows.loc[system,'window_energy_rise_kcal_mol']
        erows.append([system,f'{w:.3f}',f'{e.absolute_electronic_barrier_kcal_mol-w:.3f}',
                      f'{e.absolute_RRHO_G_barrier_kcal_mol-e.absolute_electronic_barrier_kcal_mol:.3f}',
                      f'{e.absolute_RRHO_G_barrier_kcal_mol:.3f}'])
    text+=table(['分子系','共通区間ΔE / kcal/mol','ΔE‡−共通区間ΔE / kcal/mol',
                 'ΔG‡−ΔE‡ / kcal/mol','正方向ΔG‡ / kcal/mol'],erows)+'\n'
    text+=f'''ΔG‡−ΔE‡には熱補正全体が含まれ、TΔSそのものではない。区間外の項が異なるため、共通区間の仕事だけを障壁全体の因果的な分解とは扱えない。

## 現段階での議論の書き方

「F/H置換位置に伴う原子仕事・接触方向の仕事には系間差がある。RAFA15の反応物をより低い配座へ更新すると、同じIRCの仕事に対する障壁相関が変化した。したがって、局所仕事だけからF···H接触の強化による障壁低下を主張することはできず、反応物配座の安定性とTSへの接続を含めて評価する必要がある。」

分子設計の仮説としては、**単にF···Hを近づけることではなく、アミン側F/Hと反対側F/H···O21の配置が、O21/O19を含む反応中心全体の再編とどう結び付くか**を追う方が、今回の結果に合っている。既存ログで次に比較する対象は、O21・O19・C20、アミン側置換原子、2つのアミンH、および上記のC20を含む対照接触である。仕事が変化するC–N区間と、距離・角度・符号付き力が変化する区間の対応を見ることで、4点相関を補う構造的な説明を検討できる。

再生成：`python analysis4_force/work_candidate_assessment.py`。今回の比較には新しいDFT計算を使っていない。
'''
    (OUT/'reports/work_candidate_assessment_ja.md').write_text(text)


if __name__=='__main__':main()
