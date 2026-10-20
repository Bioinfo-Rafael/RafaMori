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
    text='''# 仕事の大きさ・系間変動・ΔG‡との相関から、どの候補が残るか

## 結論

**原子の仕事ではO21が優先候補として残る。アミン側の置換F/H原子も相対的な変動と相関が大きい。一方、指定したF···アミンHの接触方向の仕事は相関が弱く、その接触を強めれば反応が速くなる、という結論にはまだつながらない。**

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
    text+='''C20–求核Nの仕事は平均1.460 kcal/mol、σ=0.318 kcal/mol、相対変動0.218。対してO21の原子仕事はM=0.056649、σ=0.046715、相対変動0.825である。**O21は相対的には大きく変わるが、絶対的な変動量はC20–Nの方が大きい。** したがって「主反応は一定で、差をほぼ補助接触だけが担う」とは言えない。

アミン側F/H原子はC=0.682、r=−0.945で、求核N（C=0.130、r=−0.663）より期待に近い。ただし、選択アミンHはC=0.147、r=−0.448、指定アミン側接触はr=−0.286である。**F側原子の仕事が相関することと、F···アミンH接触の仕事が相関することは同じではない。**

反対側F/H···O21はC=1.135、r=−0.943と期待に近いが、M=0.004164、σ=0.003806 kcal/molと小さい。これだけで数kcal/molの障壁差を直接説明するエネルギー項とは扱えない。

## 他の候補と同時に比較した場合

表示済みの和集合から重複を除いた原子40候補・分子間接触65候補・基質内接触26候補を使用する。全ての可能な接触を網羅した母集団ではない。原子は原子同士、分子間接触は分子間同士、基質内接触は基質内同士で比較する。

恣意的な閾値や重み付きの総合点を作らず、M・変動・|r|の3量がすべて指定対象以上の候補を数えた。相対変動Cで比べる場合と、絶対変動σで比べる場合の両方を示す。候補を残すために片方だけを採用しない。

'''
    text+=table(['対象','M・C・|r|がすべて以上の他候補数','M・σ・|r|がすべて以上の他候補数'],[
        [name,int(r.dominators_relative_count),int(r.dominators_absolute_count)] for name,r in zip(names,selected)])+'\n'
    text+='''- **O21**：どちらの比較でも3量すべてで上回る原子はなかった。優先候補として残るが、唯一の候補という意味ではない。
- **アミン側F/H原子**：相対変動を使うとO21が3量とも上回る。絶対変動を使うとO19も上回る。F原子だけに特異的な説明とはまだ言えない。
- **O19**：r=+0.999で、O21よりさらに強い相関を示す。符号が逆で、障壁が低い系ほどO19の仕事は小さい。O21とO19を一緒に追い、カルボニル周辺の変形・運動の違いを比較する必要がある。
- **反対側F/H···O21**：相対変動を使うと3量とも上回る基質内接触はなかった。しかし絶対変動では、C20···アミン側F/HとC20···H7の2対が上回る。どちらも反応C20を含む。反対側接触の特異性と、反応中心の力が複数方向へ投影された結果を区別する対照になる。
- **指定アミン側F/H···H**：相関が弱く、絶対変動を使う比較では多数の分子間接触が上回る。今の接触方向の仕事を中心に促進機構を主張する根拠は弱い。

全候補の値と、各対象を上回る候補名は[比較CSV](../csv/work_candidate_assessment.csv)に保存した。

## 4点の相関が特定の系に依存していないか

相関の安定性を記述するため、1系ずつ除いた3点のrも計算した。独立な検証ではない。

'''
    text+=table(['対象','1系除外時のr最小値（無次元）','1系除外時のr最大値（無次元）','RAFA15を除いたr（無次元）'],[
        [name,f'{r.loo_r_min:.3f}',f'{r.loo_r_max:.3f}',f'{r.r_without_RAFA15:.3f}'] for name,r in zip(names,selected)])+'\n'
    text+='''O21の相関はRAFA15を除くと弱くなる一方、O19は強い相関が残る。4系全体のrだけで「この相互作用だけが障壁を決める」とは言えない。RAFA15反応物の微小虚振動に由来するRRHO自由エネルギーの留保も残る。

## 「仕事が大きいほど反応しやすい」を安定化と呼べるか

この解析の正の仕事は、共通区間内の電子エネルギー上昇への寄与である。正の仕事が大きい系ほど全体のΔG‡が低いという負の相関は、そのまま「この原子がTSを安定化した」という意味にはならない。

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
    text+='''ΔG‡−ΔE‡には熱補正全体が含まれ、TΔSそのものではない。区間外の項が異なるため、共通区間の仕事だけを障壁全体の因果的な分解とは扱えない。

## 現段階での議論の書き方

「F置換位置の違いに伴い、アミン側置換原子と反応カルボニルO21の原子仕事は、求核Nより大きな相対変動を示した。これらの仕事は4系の活性化自由エネルギーと負に相関し、F置換による反応部位周辺の運動・変形の変化が反応性に関与する可能性を示す。一方、指定F···アミンH接触の仕事自体の相関は弱く、O19などにも強い相関が認められた。したがって、特定のF···H相互作用の増強が障壁を低下させるという因果関係は、この解析だけでは確定しない。」

分子設計の仮説としては、**単にF···Hを近づけることではなく、アミン側F/Hと反対側F/H···O21の配置が、O21/O19を含む反応中心全体の再編とどう結び付くか**を追う方が、今回の結果に合っている。既存ログで次に比較する対象は、O21・O19・C20、アミン側置換原子、2つのアミンH、および上記のC20を含む対照接触である。仕事が変化するC–N区間と、距離・角度・符号付き力が変化する区間の対応を見ることで、4点相関を補う構造的な説明を検討できる。

再生成：`python analysis4_force/work_candidate_assessment.py`。今回の比較には新しいDFT計算を使っていない。
'''
    (OUT/'reports/work_candidate_assessment_ja.md').write_text(text)


if __name__=='__main__':main()
