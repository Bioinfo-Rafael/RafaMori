"""RRHO entropy accounting and a 3 x 3 energy-profile visualization.

Only E is available along the IRC. G and TΔS use a constant TS correction.
The outer connectors are C1 display curves, never computed reaction paths.
"""
from pathlib import Path
import hashlib
import json
import re

import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

from energetics import KCAL, entropy_terms, reverse_barrier_terms
from plotting4 import SYSTEM_COLORS
from reporting import markdown
from figure_explanations import write_annotated_report

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / 'outputs'
QUANTITIES = ('E', 'G', 'T_delta_S')


def connector_width(span, gap, slope):
    """Display rule using E only; reuse the same x positions in all three rows."""
    return float(np.clip(2 * abs(gap) / max(abs(slope), .25), .25 * span, .75 * span))


def connector(x0, y0, slope0, x1, y1, slope1):
    if x1 <= x0:
        raise ValueError('Connector must have a positive horizontal width')
    return CubicHermiteSpline([x0, x1], [y0, y1], [slope0, slope1])


def thermochemistry_tables():
    audit = json.loads((OUT / 'reports/stationary_point_audit.json').read_text())
    lys = next(r for r in audit if r['system'] == 'Lys')
    energy = pd.read_csv(OUT / 'csv/energy_comparison.csv').set_index('system')
    states, checks = [], []
    # Check the direct entropy column independently of the H-G calculation.
    for r in audit:
        path = ROOT / r['path']
        if hashlib.sha256(path.read_bytes()).hexdigest() != r['sha256']:
            raise ValueError(f'Stationary log changed: {path}')
        text = path.read_text(errors='replace')
        matches = re.findall(r'^\s*Total\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$', text, re.M)
        if not matches:
            raise ValueError(f'No thermochemistry Total entropy row: {path}')
        entropy = float(matches[-1][2])
        from_hg = (r['H_hartree'] - r['G_hartree']) * KCAL
        residual = from_hg - r['temperature_K'] * entropy / 1000
        if abs(residual) > .001:
            raise ValueError(f'H-G and printed T*S disagree: {path}')
        checks.append(dict(system=r['system'], role=r['role'], path=str(path.relative_to(ROOT)),
                           freq_thermochemistry_blocks=text.count('Sum of electronic and thermal Free Energies='),
                           entropy_printed_cal_mol_K=entropy, TS_from_H_minus_G_kcal_mol=from_hg,
                           H_minus_G_minus_printed_TS_kcal_mol=residual))
    for path in sorted((ROOT / 'data/completed_20261017/IRC').glob('*.LOG')):
        count = path.read_text(errors='replace').count('Sum of electronic and thermal Free Energies=')
        if count:
            raise ValueError('IRC contains thermochemistry; inspect its point coverage before assuming it absent')
        checks.append(dict(system=path.stem, role='irc', path=str(path.relative_to(ROOT)), freq_thermochemistry_blocks=count))
    for system in energy.index:
        records = {r['role']: r for r in audit if r['system'] == system}
        terms = entropy_terms(records['reactant'], records['ts'], records['product'], lys)
        terms.update(reverse_barrier_terms(records['ts'], records['product']))
        for name, value in terms.items():
            energy.loc[system, name] = value
        ref = {q: records['reactant'][k] + lys[k] for q, k in [('E','energy_hartree'),('H','H_hartree'),('G','G_hartree')]}
        for role in ('separated_reactants', 'ts', 'product'):
            v = ref if role == 'separated_reactants' else {q: records[role][k] for q,k in [('E','energy_hartree'),('H','H_hartree'),('G','G_hartree')]}
            d = {q: (v[q] - ref[q]) * KCAL for q in v}
            states.append(dict(system=system, role=role, temperature_K=lys['temperature_K'], pressure_atm=lys['pressure_atm'],
                               **{f'{q}_hartree': v[q] for q in v}, TS_kcal_mol=(v['H']-v['G'])*KCAL,
                               **{f'delta_{q}_kcal_mol': d[q] for q in d}, T_delta_S_kcal_mol=d['H']-d['G'],
                               G_minus_E_correction_kcal_mol=d['G']-d['E'], RRHO_G_provisional=energy.loc[system,'RRHO_G_provisional']))
    energy.to_csv(OUT / 'csv/energy_comparison.csv')
    pd.DataFrame(states).to_csv(OUT / 'csv/thermochemistry_states.csv', index=False)
    pd.DataFrame(checks).to_csv(OUT / 'csv/thermochemistry_log_audit.csv', index=False)
    columns = ['absolute_electronic_barrier_kcal_mol','enthalpy_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol',
               'T_delta_S_activation_kcal_mol','minus_T_delta_S_activation_kcal_mol','delta_S_activation_cal_mol_K',
               'G_minus_E_activation_kcal_mol','H_minus_E_activation_kcal_mol','RRHO_G_provisional']
    energy[columns].to_csv(OUT / 'csv/activation_thermochemistry.csv')
    return energy


def build_profiles(energy):
    models, curve_rows, join_rows = {}, [], []
    for system, row in energy.iterrows():
        p = pd.read_csv(OUT / f'csv/reaction_force_{system}.csv').sort_values('R')
        x = p.R.to_numpy(); e = p.E_relative_separated_reactants_kcal_mol.to_numpy()
        slope = -p.F_smooth_kcal_mol_per_R.to_numpy()
        # Preserve every calculated IRC energy, including the exact TS value.
        espline = CubicHermiteSpline(x, e, slope)
        offset = row.G_minus_E_activation_kcal_mol
        tds = row.T_delta_S_activation_kcal_mol
        curves = {'E': espline, 'G': CubicHermiteSpline(x,e+offset,slope),
                  'T_delta_S': CubicHermiteSpline(x,np.full_like(x,tds),np.zeros_like(x))}
        products = {'E': row.electronic_reaction_energy_kcal_mol, 'G': row.RRHO_G_reaction_energy_kcal_mol,
                    'T_delta_S': row.T_delta_S_reaction_kcal_mol}
        span = x[-1]-x[0]
        left = x[0]-connector_width(span, e[0], slope[0])
        right = x[-1]+connector_width(span, products['E']-e[-1], slope[-1])
        barwidth = .12*span
        models[system] = dict(x=x, e=e, curves=curves, products=products, left=left, right=right,
                              barwidth=barwidth, connectors={})
        for q, curve in curves.items():
            kind = 'computed_IRC_E_interpolation' if q == 'E' else 'constant_TS_thermal_correction_approximation'
            dense = np.unique(np.r_[np.linspace(x[0],x[-1],501),x])
            for r,v in zip(dense,curve(dense)):
                curve_rows.append(dict(system=system, quantity=q, segment='IRC', R_display=r, value_kcal_mol=v,
                                       provenance=kind, display_coordinate_only=False))
            for side, a,b,ya,yb,sa,sb in [
                ('reactant',left,x[0],0.,float(curve(x[0])),0.,float(curve.derivative()(x[0]))),
                ('product',x[-1],right,float(curve(x[-1])),products[q],float(curve.derivative()(x[-1])),0.)]:
                c = connector(a,ya,sa,b,yb,sb)
                models[system]['connectors'][q,side] = c
                xx = np.linspace(a,b,161)
                roots = c.derivative().roots(extrapolate=False)
                turns = roots[(roots>a+1e-8)&(roots<b-1e-8)]
                yy = c(np.r_[a,b,turns])
                join_rows.append(dict(system=system, quantity=q, side=side, x_start=a, x_end=b,
                                      y_start_kcal_mol=ya, y_end_kcal_mol=yb, slope_start=sa, slope_end=sb,
                                      n_interior_turns=len(turns), connector_min_kcal_mol=yy.min(), connector_max_kcal_mol=yy.max(),
                                      overshoot_below_endpoints_kcal_mol=max(0,min(ya,yb)-yy.min()),
                                      overshoot_above_endpoints_kcal_mol=max(0,yy.max()-max(ya,yb)),
                                      provenance='visual_connector_NOT_computed_IRC'))
                for r,v in zip(xx,c(xx)):
                    curve_rows.append(dict(system=system, quantity=q, segment=side+'_connector', R_display=r, value_kcal_mol=v,
                                           provenance='visual_connector_NOT_computed_IRC', display_coordinate_only=True))
    pd.DataFrame(curve_rows).to_csv(OUT / 'csv/thermochemistry_display_profiles.csv', index=False)
    pd.DataFrame(join_rows).to_csv(OUT / 'csv/thermochemistry_connector_audit.csv', index=False)
    return models


def draw_matrix(models):
    fig, axes = plt.subplots(3,3,figsize=(18,12.5), sharey='row', sharex='col')
    rowlabels = ['Electronic energy ΔE\n(kcal/mol)', 'Gibbs free energy ΔG\nTS correction fixed (kcal/mol)',
                 'Entropy term TΔS = ΔH − ΔG\nTS value fixed (kcal/mol)']
    titles = ['IRC interval only', 'IRC + separate reactant/product levels', 'IRC + smooth illustrative connections']
    for i,q in enumerate(QUANTITIES):
        for j in range(3):
            ax = axes[i,j]
            for system, m in models.items():
                color=SYSTEM_COLORS[system];x=m['x'];curve=m['curves'][q];xx=np.linspace(x[0],x[-1],601)
                style='-' if q=='E' else '--'
                ax.plot(xx,curve(xx),color=color,ls=style,lw=1.7)
                if q=='E':ax.plot(x,m['e'],'.',color=color,ms=2.6,alpha=.5)
                ax.plot(0,curve(0),'o',ms=4,color=color,zorder=5)
                if j:
                    # Constant energy, arbitrary display x; E/G/entropy use identical x positions.
                    ax.plot([m['left']-m['barwidth'],m['left']],[0,0],color=color,lw=2.3)
                    ax.plot([m['right'],m['right']+m['barwidth']],[m['products'][q]]*2,color=color,lw=2.3)
                    ax.plot(m['left'],0,'D',color=color,ms=3.5)
                    ax.plot(m['right'],m['products'][q],'D',color=color,ms=3.5)
                if j==2:
                    for side in ('reactant','product'):
                        c=m['connectors'][q,side];xc=np.linspace(c.x[0],c.x[-1],161)
                        ax.plot(xc,c(xc),color=color,ls=':',lw=1.9,alpha=.85)
            ax.axvline(0,color='#555555',lw=.7,ls='--',alpha=.45)
            ax.axhline(0,color='#888888',lw=.6,alpha=.3)
            ax.grid(alpha=.13);ax.margins(x=.035,y=.12)
            ax.text(.015,.96,chr(65+i*3+j),transform=ax.transAxes,va='top',fontweight='bold',fontsize=12)
            if i==0:ax.set_title(titles[j],fontsize=12,pad=12)
            if j==0:ax.set_ylabel(rowlabels[i],fontsize=11)
            if i==2:ax.set_xlabel('IRC reaction coordinate R' if j==0 else 'R inside IRC; arbitrary display coordinate outside',fontsize=10)
    for j in (1,2):
        axes[0,j].text(.04,.87,'Reactant + Lys',transform=axes[0,j].transAxes,fontsize=9,color='#555555')
        axes[0,j].text(.78,.87,'Product',transform=axes[0,j].transAxes,fontsize=9,color='#555555')
    handles=[Line2D([0],[0],color=c,lw=2,label=s) for s,c in SYSTEM_COLORS.items()]
    handles += [Line2D([0],[0],color='#555555',ls='-',label='E along IRC'),
                Line2D([0],[0],color='#555555',ls='--',label='G / TΔS: constant-TS approximation'),
                Line2D([0],[0],color='#555555',ls=':',label='Illustrative connection (not calculated)'),
                Line2D([0],[0],color='#555555',marker='D',label='Stationary-state level')]
    fig.suptitle('Four-system energy and entropy profiles: calculated values, fixed-TS approximations, and display connections',fontsize=14,y=.992)
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.965),ncol=4,frameon=False,fontsize=10)
    fig.text(.5,.015,'Each system uses its own separated reactant + Lys reference; no MoriPF subtraction.  T = 298.15 K, RRHO at 1 atm.\nG − E includes thermal enthalpy corrections and is not TΔS.  RAFA15 thermochemistry remains provisional (small imaginary reactant mode).',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.065,1,.89),w_pad=1.1,h_pad=1.5)
    for ext in ('png','svg'):
        fig.savefig(OUT/f'figures/23_energy_entropy_profiles_3x3.{ext}',dpi=190)
    plt.close(fig)


def write_report(energy):
    cols=['absolute_electronic_barrier_kcal_mol','enthalpy_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol',
          'T_delta_S_activation_kcal_mol','minus_T_delta_S_activation_kcal_mol','G_minus_E_activation_kcal_mol','H_minus_E_activation_kcal_mol']
    shown=energy[cols].reset_index()
    shown.columns=['分子系','ΔE‡','ΔH‡','ΔG‡','TΔS‡','−TΔS‡','ΔG‡−ΔE‡','ΔH‡−ΔE‡']
    text='''# 活性化障壁の熱力学成分と、反応座標に沿った3×3の図

## 活性化障壁とエントロピーの寄与

すべて各系の分離反応物（基質＋Lys）から遷移状態までの差で、単位はkcal/mol。MoriPFとの差ではない。298.15 K・1 atmのGaussian RRHO値を用いた。

'''+markdown(shown)+'''

電子エネルギーEとエンタルピーHは異なるため、**G−EをそのままTΔSとは呼べない**。次の関係で区別する。

```text
TΔS‡ = ΔH‡ − ΔG‡
ΔG‡ = ΔE‡ + (ΔH‡ − ΔE‡) − TΔS‡
ΔG‡ − ΔE‡ = (ΔH‡ − ΔE‡) − TΔS‡
```

ΔH‡−ΔE‡にはゼロ点エネルギーや熱的な補正が含まれる。表では、依頼されたEとGの差も独立の列として残した。TΔS‡が負なら、エントロピーの寄与−TΔS‡は正となり、自由エネルギー障壁を押し上げる。

RAFA15の−TΔS‡は約11.82 kcal/mol、他の3系は約14.38〜14.59 kcal/mol。今回の配座とRRHO近似では、RAFA15の障壁が低いことには電子エネルギーだけでなく、エントロピーによる押し上げが小さいことも寄与している。ただし、この差をH···O接触だけに割り当てることはできない。RAFA15反応物には微小虚振動が残り、熱化学量は暫定値である。

## 3×3の図

![E・G・TΔSと両端準位の比較](../figures/23_energy_entropy_profiles_3x3.png)

[拡大可能なSVG](../figures/23_energy_entropy_profiles_3x3.svg)

### IRC上のGとTΔSに置いた仮定

4つのIRCログに、freq由来の熱化学ブロックはなかった。そこで、指定に従い、TSの熱補正を全IRC点に固定した参考曲線を描いた。

```text
G_approx(R) = E_IRC(R) + [G_TS − E_TS]
H_approx(R) = E_IRC(R) + [H_TS − E_TS]
ΔG_approx(R) = G_approx(R) − (G_reactant + G_Lys)
TΔS_approx(R) = (H_TS − G_TS)
                 − [(H_reactant − G_reactant) + (H_Lys − G_Lys)]
```

したがって、Gの参考曲線はEと同じ形で上下にずれ、TΔSの参考曲線は水平になる。これは「実際のエントロピーが反応中ずっと一定」と分かったという意味ではなく、今回選んだ近似の結果である。TSの点と、反応物・生成物の横棒は、それぞれのfreq計算の値を使う。TS以外のG・TΔSのIRC曲線は直接計算した値ではない。

### 右列を滑らかにつなぐ方法

計算済みIRC点の横位置とEは変えない。Eの表示曲線は既存のIRC点を通り、既存解析の滑らかな微分値を傾きとして使う3次Hermite補間で描く。左列・中央列・右列のIRC区間はすべて同じ曲線である。

反応物・生成物の横棒はIRC区間の外側に置く。接続幅Lは、電子エネルギーの準位差ΔE_gapとIRC末端の傾きmから次の表示規則で決める。

```text
L = clip(2 |ΔE_gap| / max(|m|, 0.25), 0.25 × IRC幅, 0.75 × IRC幅)
```

エネルギーはkcal/mol、傾きはkcal/mol/R。この規則は見やすい横幅を確保するためのもので、分子間距離や追加のIRC計算から決めた値ではない。各系で決めた同じ横位置をE・G・TΔSの3行に使用する。

接続線は3次Hermite曲線で、IRC側では値と傾きをIRC末端に合わせ、横棒側では傾きを0にする。これにより接続点で値と1階微分が連続になり、垂直な接続や傾きの急な折れを避ける。2階微分の連続性は課していない。

**準位差と末端の傾きの向きが合わなければ、滑らかな接続には途中の谷・山が必要になる。** たとえば左の分離反応物の準位がIRC左端より高く、IRC左端ではTSへ向かってエネルギーが上がっている場合、一度下がってから上がる形になる。この図の点線に現れる谷・山は、未知の中間体や障壁を発見したものではない。RAFA15は2026-10-22受領の新生成物ログに更新した。生成物optはIRC末端より低エネルギーとなり、構造差も小さくなったが、接続線の形は機構の根拠には使わない。[逆向き障壁と生成物の再評価](directional_barriers_discussion_ja.md)を参照。

## 保存した値と検証

- [活性化障壁と熱力学成分](../csv/activation_thermochemistry.csv)：ΔE‡、ΔH‡、ΔG‡、TΔS‡、−TΔS‡、ΔS‡、G−E、H−E。
- [反応物・TS・生成物の値](../csv/thermochemistry_states.csv)：生のE/H/G（Hartree）とTS、各系の分離反応物からの差。
- [ログの確認記録](../csv/thermochemistry_log_audit.csv)：H−Gとログに印字されたT×Sの一致、およびIRCの熱化学ブロック数。
- [図の曲線データ](../csv/thermochemistry_display_profiles.csv)：IRC、固定TS近似、表示用接続線を列で区別。
- [接続線の確認記録](../csv/thermochemistry_connector_audit.csv)：横位置・傾き・途中の極値・端点の範囲を超えた量。

H−Gはログに印字されたT×Sと0.001 kcal/mol以内で一致することを確認する。小さな差は印字桁数による。計算式の背景は[Gaussian公式の熱化学解説](https://gaussian.com/thermo/)を参照。

再生成：`python analysis4_force/thermochemistry_profiles.py`。新しい量子化学計算は実行しない。
'''
    write_annotated_report(OUT/'reports/thermochemistry_profiles_ja.md',text)
    # Include the requested component table in the main summary too.
    path=OUT/'reports/summary_ja.md';summary=path.read_text()
    block='''<!-- thermochemistry:start -->
## 活性化障壁のエントロピー項

各系の分離反応物からTSまでの差（kcal/mol）。TΔS‡はfreqのH−Gから計算し、EとGの差は別の列に示す。

'''+markdown(shown)+'''

G−Eには熱的なエンタルピー補正も含まれるため、TΔSとは異なる。RAFA15はエントロピーによる障壁の押し上げ（−TΔS‡）が他の3系より小さい。微小虚振動による暫定値という扱いは維持する。

[3×3のE・G・TΔS図と、滑らかな接続の説明](thermochemistry_profiles_ja.md)
<!-- thermochemistry:end -->
'''
    summary=re.sub(r'<!-- thermochemistry:start -->.*?<!-- thermochemistry:end -->\n*','',summary,flags=re.S)
    summary=summary.replace('## 仮説に対する解釈',block+'\n## 仮説に対する解釈')
    path.write_text(summary)


def main():
    energy=thermochemistry_tables()
    models=build_profiles(energy)
    draw_matrix(models)
    write_report(energy)
    from directional_barriers import main as write_directional_barriers
    write_directional_barriers()
    print(energy[['T_delta_S_activation_kcal_mol','G_minus_E_activation_kcal_mol']].round(4).to_string())


if __name__=='__main__':
    main()
    from figure_explanations import write_guide
    write_guide()
    import hashlib
    manifest_path=OUT/'reports/output_manifest.json'
    manifest=[dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
              for p in sorted(OUT.rglob('*')) if p.is_file() and p!=manifest_path]
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
