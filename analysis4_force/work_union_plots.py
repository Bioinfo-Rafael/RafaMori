"""Small multiples for every row of the five work-union comparison tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from plotting4 import SYSTEM_COLORS
from numbered_structures import SYSTEMS

OUT=Path(__file__).resolve().parent/'outputs'
GROUPS=[
    ('atom','four_target_atoms','all_57_atoms','atoms','原子の仕事'),
    ('contact','amine_contact','distance_and_chemistry','amine_matched','アミン側・種類と距離を限定'),
    ('contact','amine_contact','same_category','amine_broad','アミン側・広い比較'),
    ('contact','opposite_contact','distance_and_chemistry','opposite_matched','反対側・種類と距離を限定'),
    ('contact','opposite_contact','same_category','opposite_broad','反対側・広い比較')]


def statistics(x,y):
    """Only four genuinely corresponding systems receive a fit or CV."""
    n=len(x)
    result=dict(n=n,mean_work_kcal_mol=float(np.mean(x)),
                population_std_work_kcal_mol=float(np.std(x,ddof=0)),
                std_over_signed_mean=np.nan,slope=np.nan,intercept_kcal_mol=np.nan,
                pearson_r=np.nan,r_squared=np.nan)
    if n!=4:return dict(result,status='incomplete_correspondence')
    tol=1e-12*max(1.,float(np.max(np.abs(x))))
    if abs(result['mean_work_kcal_mol'])>tol:
        result['std_over_signed_mean']=result['population_std_work_kcal_mol']/result['mean_work_kcal_mol']
    if np.ptp(x)<=tol:return dict(result,status='constant_work_no_fit')
    slope,intercept=np.polyfit(x,y,1)
    predicted=slope*x+intercept
    result.update(slope=float(slope),intercept_kcal_mol=float(intercept),
                  pearson_r=float(np.corrcoef(x,y)[0,1]),
                  r_squared=float(1-np.sum((y-predicted)**2)/np.sum((y-y.mean())**2)))
    return dict(result,status='ok' if np.isfinite(result['std_over_signed_mean']) else 'mean_near_zero')


def main():
    manifest_path=OUT/'reports/work_union_subplot_manifest.json'
    previous=json.loads(manifest_path.read_text()) if manifest_path.exists() else []
    data=pd.read_csv(OUT/'csv/cross_system_work_union.csv')
    energy=pd.read_csv(OUT/'csv/energy_comparison.csv').set_index('system')
    barriers=energy.absolute_RRHO_G_barrier_kcal_mol
    assert set(SYSTEMS)<=set(barriers.index)
    # The energy table contains freq-based, separated-reactant forward barriers.
    audits=json.loads((OUT/'reports/cross_system_work_union_audit.json').read_text())
    fitrows=[];plotpoints=[];figures=[]
    gallery='<a id="work-plots"></a>\n\n## 仕事と活性化自由エネルギーの散布図・4系の棒グラフ\n\n'
    gallery+=(OUT.parent/'docs/work_comparison_groups_ja.md').read_text()+'\n\n'
    gallery+='[全散布図のPDF](../figures/25_work_union_scatter_all.pdf) · [全棒グラフのPDF](../figures/26_work_union_bars_all.pdf)\n\n'
    gallery+='全5表の和集合を、表ごとに最大4列×4行のパネルへ分けた。1パネルが1つの対応する原子または原子対。表の間で重複する対象は同じ計算値を再掲する。各パネル左上のIDと位置ラベル、下側の原子番号対応を表・CSVと照合できる。青=MoriPF、緑=RAFA14、紫=RAFA15、橙=RAFA16。\n\n'
    gallery+='**散布図**：横軸はC–N 2.50→2.10 Å区間の符号付き仕事（kcal/mol）、縦軸はfreq計算からの正方向のΔG‡ = G(TS)−G(基質)−G(Lys)（kcal/mol）。4系の対応が確定したパネルだけ、切片付き最小二乗直線を当てはめ、タイトルにR²を小数点以下3桁で示す。赤い枠は着目対象。比較条件外の系も、対応する仕事が得られる場合は4点の比較に含める。条件内の点を丸、条件外を×で区別した。順位は元の候補条件のままであり、フィットのために選び直していない。\n\n'
    gallery+=r'''$$
\widehat{\Delta G^\ddagger}=aW+b,\qquad
R^2=1-\frac{\sum_s(\Delta G_s^\ddagger-\widehat{\Delta G_s^\ddagger})^2}
{\sum_s(\Delta G_s^\ddagger-\overline{\Delta G^\ddagger})^2}
$$

**棒グラフ**：横軸は4分子系、縦軸は符号付き仕事（kcal/mol）。タイトルのSTD/meanは、4系の母標準偏差（ddof=0）を符号付き算術平均で割った無次元量で、小数点以下3桁を表示する。

$$
\overline W=\frac14\sum_sW_s,\qquad
\sigma_W=\sqrt{\frac14\sum_s(W_s-\overline W)^2},\qquad
\mathrm{STD/mean}=\frac{\sigma_W}{\overline W}
$$

平均が負ならSTD/meanも負になる。平均が0に近い場合は比が大きくなり、正負の相殺にも影響されるため、値だけを「効果の強さ」とは読まない。平均の絶対値が `1e-12 × max(1, max|W|)` 以下ならN/Aとする。今回の対応再確認により、全行に4系の対応する仕事が入り、各散布図はn=4となった。以前の対応未確定による1点だけの図は修正済み。

全パネルでΔG‡の縦軸範囲は共通。仕事の軸は小さな値も読めるようにパネルごとに調整しているため、傾きの見た目や棒の高さだけでは異なるパネルを比較できない。各軸の数値を確認する。**4系だけの探索的なフィットであり、高いR²だけで因果関係や相互作用の特異性を示すものではない。** 2026-10-25のRAFA15反応物更新を反映し、障壁と回帰を再計算した。新反応物は虚振動なしだが、配座・反応経路の不確かさは残る。

'''
    gallery+='[全パネルのR²・回帰係数・平均・STD/mean](../csv/work_union_subplot_statistics.csv) · [各点の仕事・ΔG‡・原子番号・条件内外](../csv/work_union_subplot_points.csv)\n\n'
    legends=[Patch(color=SYSTEM_COLORS[s],label=s) for s in SYSTEMS]
    pdf_scatter=OUT/'figures/25_work_union_scatter_all.pdf'
    pdf_bars=OUT/'figures/26_work_union_bars_all.pdf'
    with PdfPages(pdf_scatter) as pdfs,PdfPages(pdf_bars) as pdfb:
        for kind,target,pool,slug,japanese in GROUPS:
            group=data[(data.kind==kind)&(data.target==target)&(data.pool==pool)]
            keys=list(dict.fromkeys(group.cross_key))
            audit=next(a for a in audits if a['kind']==kind and a['target']==target and a['pool']==pool)
            assert set(keys)==set(audit['union_keys'])
            panels=[]
            for k,key in enumerate(keys,1):
                rows=group[group.cross_key==key].set_index('system').reindex(SYSTEMS)
                valid=rows.correspondence_available.eq(True)&rows.work_kcal_mol.notna()
                systems=list(rows.index[valid]);x=rows.loc[valid,'work_kcal_mol'].to_numpy(float)
                y=barriers.loc[systems].to_numpy(float)
                assert len(x)>0
                stat=statistics(x,y);panel_id=f'{slug}-{k:03d}'
                fitrows.append(dict(panel_id=panel_id,kind=kind,target=target,pool=pool,cross_key=key,**stat))
                for s in systems:
                    plotpoints.append(dict(panel_id=panel_id,system=s,cross_key=key,
                        atom_or_pair=rows.loc[s,'atom_or_pair'],work_kcal_mol=rows.loc[s,'work_kcal_mol'],
                        forward_activation_G_kcal_mol=barriers[s],in_candidate_pool=rows.loc[s,'in_candidate_pool']))
                panels.append((panel_id,key,rows,systems,x,y,stat))
            gallery+=f'### {japanese}：{len(panels)}パネル\n\n'
            for start in range(0,len(panels),16):
                chunk=panels[start:start+16];page=start//16+1;nrows=(len(chunk)+3)//4
                fs,axs=plt.subplots(nrows,4,figsize=(17,4.0*nrows),squeeze=False)
                fb,axb=plt.subplots(nrows,4,figsize=(17,4.0*nrows),squeeze=False)
                for index,(panel_id,key,rows,systems,x,y,stat) in enumerate(chunk):
                    ax=axs.flat[index];bx=axb.flat[index]
                    nominated=rows.nominated.eq(True).any()
                    label=key.replace('||',' ~ ').replace('unmapped_','H?')
                    for s,xx,yy in zip(systems,x,y):
                        inside=bool(rows.loc[s,'in_candidate_pool'])
                        ax.scatter(xx,yy,s=38,c=SYSTEM_COLORS[s],marker='o' if inside else 'x',zorder=3)
                    if np.isfinite(stat['r_squared']):
                        xx=np.linspace(x.min(),x.max(),100)
                        ax.plot(xx,stat['slope']*xx+stat['intercept_kcal_mol'],color='#444444',lw=1)
                    ax.set_title(f"R² = {stat['r_squared']:.3f}" if np.isfinite(stat['r_squared']) else 'R² = N/A',fontsize=12)
                    ax.set_ylim(barriers.min()-.8,barriers.max()+3.3)
                    ax.margins(x=.15)
                    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
                    ax.set_xlabel('Atomic work / kcal mol⁻¹' if kind=='atom' else 'Radial work / kcal mol⁻¹',fontsize=9)
                    ax.set_ylabel('Forward ΔG‡ / kcal mol⁻¹',fontsize=9)
                    ax.text(.03,.97,f'{panel_id}   n={stat["n"]}\n{label}',transform=ax.transAxes,va='top',fontsize=8)
                    for i,s in enumerate(SYSTEMS):
                        if s in systems:
                            bx.bar(i,rows.loc[s,'work_kcal_mol'],color=SYSTEM_COLORS[s],width=.65,
                                   hatch=None if bool(rows.loc[s,'in_candidate_pool']) else '//')
                        else:bx.text(i,.04,'N/A',transform=bx.get_xaxis_transform(),ha='center',fontsize=8)
                    bx.axhline(0,color='#555555',lw=.7);bx.set_xlim(-.65,3.65)
                    bx.set_xticks(range(4),['MoriPF','R14','R15','R16'],fontsize=8)
                    ratio=stat['std_over_signed_mean']
                    bx.set_title(f'STD/mean = {ratio:.3f}' if np.isfinite(ratio) else 'STD/mean = N/A',fontsize=12)
                    bx.set_ylabel('Atomic work / kcal mol⁻¹' if kind=='atom' else 'Radial work / kcal mol⁻¹',fontsize=9)
                    bx.margins(y=.35)
                    bx.text(.03,.97,f'{panel_id}   n={stat["n"]}\n{label}',transform=bx.transAxes,va='top',fontsize=8,
                            bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=1))
                    mapping='; '.join(f'{s.replace("RAFA","R")}:{rows.loc[s,"atom_or_pair"]}' for s in systems)
                    # Two lines keep full original atom IDs readable in a small panel.
                    pieces=mapping.split('; ')
                    mapping='; '.join(pieces[:2])+'\n'+'; '.join(pieces[2:])
                    for a in [ax,bx]:
                        a.text(.0,-.31,mapping,transform=a.transAxes,fontsize=7,va='top')
                        a.grid(axis='y',alpha=.18);a.tick_params(labelsize=8)
                        a.ticklabel_format(axis='y' if a is bx else 'x',style='sci',scilimits=(-3,4),useOffset=False)
                        if nominated:
                            for spine in a.spines.values():spine.set_edgecolor('#C44444');spine.set_linewidth(1.6)
                for index in range(len(chunk),nrows*4):
                    axs.flat[index].set_visible(False);axb.flat[index].set_visible(False)
                for fig,prefix,desc,pdf in [(fs,25,'scatter',pdfs),(fb,26,'bars',pdfb)]:
                    fig.suptitle(f'{slug} | {desc} | page {page} | C–N 2.50 → 2.10 Å',fontsize=14,y=.995)
                    fig.legend(handles=legends,loc='upper center',bbox_to_anchor=(.5,.975),ncol=4,frameon=False)
                    fig.tight_layout(rect=(0,.015,1,.94),h_pad=4.0,w_pad=2)
                    path=OUT/f'figures/{prefix}_{slug}_{desc}_{page:02d}.png'
                    fig.savefig(path,dpi=160);pdf.savefig(fig);plt.close(fig)
                    figures.append(dict(group=slug,kind=desc,page=page,path=path.name,panel_ids=[p[0] for p in chunk]))
                gallery+=f'#### ページ{page}\n\n'
                gallery+=f'![{japanese}：仕事とΔG‡・ページ{page}](../figures/25_{slug}_scatter_{page:02d}.png)\n\n'
                gallery+=f'![{japanese}：4系の仕事・ページ{page}](../figures/26_{slug}_bars_{page:02d}.png)\n\n'
    pd.DataFrame(fitrows).to_csv(OUT/'csv/work_union_subplot_statistics.csv',index=False)
    pd.DataFrame(plotpoints).to_csv(OUT/'csv/work_union_subplot_points.csv',index=False)
    (OUT/'reports/work_union_subplot_manifest.json').write_text(json.dumps(figures,indent=2)+'\n')
    current_names={f['path'] for f in figures}
    for old in previous:
        # Remove only pages explicitly recorded as generated by this script.
        if old['path'] not in current_names:
            stale=OUT/'figures'/old['path']
            assert stale.parent==OUT/'figures' and stale.name.startswith(('25_','26_'))
            stale.unlink(missing_ok=True)
    gallery+='[全散布図のPDF](../figures/25_work_union_scatter_all.pdf) · [全棒グラフのPDF](../figures/26_work_union_bars_all.pdf)\n'
    print(f'Generated {len(figures)} subplot figures; {len(fitrows)} table rows per chart type.')
    return gallery


if __name__=='__main__':main()
