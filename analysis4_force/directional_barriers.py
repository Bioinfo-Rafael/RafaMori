"""Reproducible forward/reverse barriers and their Japanese interpretation."""
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import json

import numpy as np
import pandas as pd

from energetics import KCAL, reverse_barrier_terms
from reporting import markdown

OUT = Path(__file__).resolve().parent / 'outputs'


def main():
    audit = json.loads((OUT / 'reports/stationary_point_audit.json').read_text())
    lys = next(r for r in audit if r['system'] == 'Lys')
    energy = pd.read_csv(OUT / 'csv/energy_comparison.csv').set_index('system')
    endpoints = pd.read_csv(OUT / 'csv/endpoint_correspondence.csv').set_index('system')
    rows = []
    for system in energy.index:
        records = {r['role']: r for r in audit if r['system'] == system}
        reactant, ts, product = (records[k] for k in ('reactant', 'ts', 'product'))
        row = dict(system=system, temperature_K=ts['temperature_K'], pressure_atm=ts['pressure_atm'],
                   product_log=product['path'], product_sha256=product['sha256'],
                   reactant_log=reactant['path'], reactant_sha256=reactant['sha256'],
                   TS_log=ts['path'], TS_sha256=ts['sha256'],
                   forward_RRHO_G_provisional=bool(reactant['imaginary_frequencies']),
                   **reverse_barrier_terms(ts, product))
        for label, key in [('E', 'energy_hartree'), ('H', 'H_hartree'), ('G', 'G_hartree')]:
            reference = reactant[key] + lys[key]
            row[f'forward_delta_{label}_kcal_mol'] = (ts[key] - reference) * KCAL
            row[f'reverse_delta_{label}_kcal_mol'] = (ts[key] - product[key]) * KCAL
            row[f'reaction_delta_{label}_kcal_mol'] = (product[key] - reference) * KCAL
            np.testing.assert_allclose(
                row[f'forward_delta_{label}_kcal_mol'] - row[f'reverse_delta_{label}_kcal_mol'],
                row[f'reaction_delta_{label}_kcal_mol'], atol=1e-8, rtol=0)
        for direction in ('forward', 'reverse', 'reaction'):
            row[f'{direction}_T_delta_S_kcal_mol'] = (
                row[f'{direction}_delta_H_kcal_mol'] - row[f'{direction}_delta_G_kcal_mol'])
        np.testing.assert_allclose(row['forward_T_delta_S_kcal_mol'] - row['reverse_T_delta_S_kcal_mol'],
                                   row['reaction_T_delta_S_kcal_mol'], atol=1e-8, rtol=0)
        rows.append(row)
    table = pd.DataFrame(rows).set_index('system')
    table.to_csv(OUT / 'csv/directional_barriers.csv')
    (OUT / 'reports/directional_barriers.json').write_text(json.dumps(rows, indent=2) + '\n')

    def shown(direction):
        selected = table.loc[['RAFA15', 'MoriPF', 'RAFA16', 'RAFA14'], [
            f'{direction}_delta_G_kcal_mol', f'{direction}_T_delta_S_kcal_mol',
            f'{direction}_delta_E_kcal_mol', f'{direction}_delta_H_kcal_mol']].reset_index()
        suffix = '' if direction == 'reaction' else '‡'
        selected.columns = ['分子系', f'ΔG{suffix}', f'TΔS{suffix}', f'ΔE{suffix}', f'ΔH{suffix}']
        return markdown(selected)

    r = table.loc['RAFA15']
    endpoint = endpoints.loc['RAFA15']
    product = next(x for x in audit if x['system'] == 'RAFA15' and x['role'] == 'product')
    profile = pd.read_csv(OUT / 'csv/reaction_force_RAFA15.csv').sort_values('R')
    endpoint_energy = profile.E_relative_separated_reactants_kcal_mol.iloc[-1]
    updated = datetime.now(ZoneInfo('Asia/Tokyo')).date().isoformat()
    text = f'''# 正方向・逆方向の活性化障壁とRAFA15生成物の解釈

更新日：{updated}（JST）。TS・reactant・product・Lysのfreqログから再計算した。RAFA15の現在の生成物は `{Path(product['path']).name}`。追加のDFT計算や図の接続曲線からの値の読み取りは行っていない。

## 計算の定義

全数値の単位は **kcal/mol**。熱化学量は298.15 K・1 atmのGaussian RRHO値を使う。正方向の出発状態は分離した基質とLysの和、逆方向の出発状態は対応するproductである。

```text
正方向：ΔX‡forward = X(TS) − [X(reactant) + X(Lys)]
逆方向：ΔX‡reverse = X(TS) − X(product)
反応量：ΔXreaction = X(product) − [X(reactant) + X(Lys)]
         ΔXreaction = ΔX‡forward − ΔX‡reverse
X = E, H, G
TΔS‡ = ΔH‡ − ΔG‡
ΔG‡ = ΔH‡ − TΔS‡
```

TΔS‡はエントロピー項であり、それ単独を活性化障壁とは呼ばない。自由エネルギー障壁への寄与は **−TΔS‡**。EとHは異なるため、G−EをTΔSとは扱わない。MoriPFをゼロにする差し引きは行っていない。

## 正方向：reactant＋Lys → TS

以下のΔE・ΔH・ΔGは活性化障壁（‡）を表す。

{shown('forward')}

RAFA15反応物を2026-10-25に `RAFA15__REACTANT_OPTFREQ2 (1).LOG` へ更新した。保存先は従来どおり `Reactant/RAFA15_UP.LOG`。新反応物は虚振動なし（最低12.9416 cm⁻¹）で、旧ログの虚振動による暫定フラグは解除した。負のΔE‡はTSが分離反応物より電子エネルギーで低いことを示し、無障壁反応を意味しない。

自由エネルギー障壁は **{' < '.join(table.sort_values('forward_delta_G_kcal_mol').index)}** の順となる。旧反応物基準のRAFA15はΔE‡ = −4.587、ΔG‡ = 8.943だったが、新基準ではそれぞれ{r.forward_delta_E_kcal_mol:.3f}、{r.forward_delta_G_kcal_mol:.3f}へ増加した。旧解析の「RAFA15が最も低障壁」という結論は撤回する。MoriPF・RAFA16・RAFA14の相互の順序は変わらない。

RAFA15の−TΔS‡は約{-r.forward_T_delta_S_kcal_mol:.3f}で、他の3系より小さい。ただし電子障壁は4系で最大となり、現在のΔG‡も4系で最大である。F···HやH···O接触の因果効果はこの差だけから特定できず、反応物の配座・振動にも依存する。

## 逆方向：product → TS

{shown('reverse')}

RAFA15の逆向き障壁は **ΔE‡ = {r.reverse_delta_E_kcal_mol:.3f}、ΔH‡ = {r.reverse_delta_H_kcal_mol:.3f}、ΔG‡ = {r.reverse_delta_G_kcal_mol:.3f}**。TΔS‡は{r.reverse_T_delta_S_kcal_mol:.3f}で、自由エネルギー障壁へのエントロピー寄与−TΔS‡は{-r.reverse_T_delta_S_kcal_mol:.3f}である。

逆方向の差にはreactantの値が入らないため、今回の反応物置換では逆向き障壁は変わらない。表は指定productとTSの差であり、実際の逆反応障壁と解釈するには同じ反応経路で結ばれていることの確認が必要である。

## RAFA15生成物の置換とIRC末端との対応

2026-10-22に旧 `RAFA15_HB_LYS.LOG` を `RAFA15_LYS_OPTFREQ2.LOG` に置換した。旧ログとその圧縮原本は [archive/20261022/Product](../../../data/completed_20261017/archive/20261022/Product/) へ移動した。

旧ログでは逆向きΔE‡ = −1.750、ΔG‡ = −3.291であった。新ログでは両方が正となり、**「指定生成物がTSより高エネルギーである」という旧ログに基づく解釈は解消した。**

新生成物はopt/freq正常終了、57原子・165振動モード・虚振動なし。M06-2X/6-311++G(d,p)、SMD水、298.15 K・1 atmで既存のTSと条件が一致する。

- 分離反応物基準の電子エネルギー：TSは{r.forward_delta_E_kcal_mol:.3f}、IRC生成物側末端は{endpoint_energy:.3f}、新生成物は{r.reaction_delta_E_kcal_mol:.3f}。
- 新生成物とIRC末端の差：E(product)−E(IRC末端) = {endpoint.product_opt_minus_IRC_endpoint_kcal_mol:.3f}。旧ログの+3.708から低下した。
- 原子番号を固定した重ね合わせRMSD：{endpoint.product_opt_vs_IRC_RMSD_A:.3f} Å（旧ログは2.111 Å）。
- C20–N43距離：新生成物は{endpoint.product_opt_CN_A:.3f} Å、IRC末端は{endpoint.IRC_product_CN_A:.3f} Å。

新生成物はIRC末端より低エネルギーで、構造差も旧ログより小さい。今回のIRCから緩和した生成物としての対応を支持する結果である。ただしRMSDだけで経路接続を証明することはできず、IRC末端から新生成物への連続した最適化履歴は今回の解析では検証していない。

図23の点線は両端準位とIRCを結ぶ表示用の曲線であり、新たに計算された反応経路ではない。点線に現れる谷・山から障壁や中間体を読み取らない。

## 生成物の安定性と生成しやすさ

各系の分離反応物に対する反応量：

{shown('reaction')}

RAFA15の新生成物は分離反応物に対してEが{r.reaction_delta_E_kcal_mol:.3f}、Gが{r.reaction_delta_G_kcal_mol:.3f}となる。新しい分離反応物基準ではE・Gともに生成物側が高くなり、反応ΔGは4系で最大となる。これは生成物そのものを再変更したためではなく、反応物の基準エネルギーが下がったためである。

RAFA15は新基準で正方向のΔG‡が最も高く、一方、生成物から指定TSへ戻るΔG‡は小さいままである。新反応物配座から既存TSへの接続・配座変換経路は未検証である。ただし単一配座・RRHOのエネルギー差から実験の収率や寿命は断定できない。配座集団、濃度、標準状態の変換、低振動数の扱い、後続反応は含めていない。

## 根拠と再生成

- [正方向・逆方向・反応量のCSV](../csv/directional_barriers.csv) と [JSON・入力ハッシュ](directional_barriers.json)
- [freqログの監査値](stationary_point_audit.json) と [状態ごとの熱力学量](../csv/thermochemistry_states.csv)
- [生成物optとIRC末端の比較](../csv/endpoint_correspondence.csv)
- [更新した3×3図](../figures/23_energy_entropy_profiles_3x3.png) と [表示・近似の説明](thermochemistry_profiles_ja.md)
- [旧反応物の退避記録](../../../data/completed_20261017/archive/20261025/manifest.json)：旧 `RAFA15_UP.LOG` のLOGとgzipを保存。
- [現行入力manifest](../../../data/completed_20261017/manifest.json) と [旧入力の退避記録](../../../data/completed_20261017/archive/20261022/manifest.json)

`python analysis4_force/run_analysis.py` で入力から全解析・関連図を再生成できる。既存監査値からの本表・CSV・JSONだけの再生成は `python analysis4_force/directional_barriers.py`。各系で「正方向−逆方向＝反応量」をE・H・G・TΔSについて検証し、小数第3位に丸めて表示した。
'''
    (OUT / 'reports/directional_barriers_discussion_ja.md').write_text(text)


if __name__ == '__main__':
    main()
    import hashlib
    manifest_path = OUT / 'reports/output_manifest.json'
    manifest = [dict(path=str(p.relative_to(OUT)), bytes=p.stat().st_size,
                     sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                for p in sorted(OUT.rglob('*')) if p.is_file() and p.name != '.DS_Store' and p != manifest_path]
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
