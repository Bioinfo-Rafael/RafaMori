"""Fixed-fragment SMD energy differences; Gaussian supplies all energies."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from irc_eda.parsing import read_text, stationary, digest, orientation, HARTREE_KCAL, TABLE

OUT = ROOT / 'irc_eda/outputs/smd_distortion_interaction'
REPORT = ROOT / 'irc_eda/reports/smd_distortion_interaction'
REFS = {'RAFA14':'RAFA14.LOG.gz','RAFA15':'RAFA15_UP.LOG.gz','RAFA16':'RAFA16.LOG.gz','MoriPF':'MORIPF.LOG.gz'}


def save(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2) + '\n', encoding='utf8')


def electronic_check(text):
    from irc_eda.parsing import conditions
    c = conditions(text)
    assert c['method'] == 'M062X' and c['basis'].lower() == '6-311++g(d,p)', c
    assert c['charge'] == 0 and c['multiplicity'] == 1, c
    assert 'smd' in c['solvent'] and 'water' in c['solvent'], c
    assert '(included in total energy above)' in text
    return c


def reference_mapping(reference, point, indices):
    """Match covalent graphs with fixed heavy-atom order; equivalent H arbitrary.

    Used only for atom identity auditing, never to change geometries or energies.
    """
    z=reference['atomic_numbers']; a=np.array(reference['xyz'])
    b=np.array([point['xyz'][i-1] for i in indices]); radii={1:.31,6:.76,7:.71,8:.66,9:.57}
    def graph(x):
        return {(i,j) for i in range(len(z)) for j in range(i) if np.linalg.norm(x[i]-x[j])<1.25*(radii[z[i]]+radii[z[j]])}
    ga,gb=graph(a),graph(b)
    mapping={i:i for i in range(len(z)) if z[i]!=1}
    def hydrogens(g,parent):
        return sorted(i for i in range(len(z)) if z[i]==1 and (max(i,parent),min(i,parent)) in g)
    for parent in list(mapping):
        ha,hb=hydrogens(ga,parent),hydrogens(gb,parent)
        assert len(ha)==len(hb), 'Reference hydrogen connectivity differs'
        mapping.update(zip(ha,hb))
    assert len(mapping)==len(z) and len(set(mapping.values()))==len(z)
    transformed={(max(mapping[i],mapping[j]),min(mapping[i],mapping[j])) for i,j in ga}
    assert transformed==gb, 'Reference is not the same molecular connectivity'
    return {'reference_atom_to_original_atom':[indices[mapping[i]] for i in range(len(z))], 'connectivity_matched':True,'equivalent_hydrogens':'sorted bijection within same parent; not a unique stereochemical assignment','bond_detection':'distance < 1.25 * sum(covalent radii); heavy atom indices held fixed'}


def prepare(system):
    points = [p for p in json.loads((ROOT/'irc_eda/outputs/project/points.json').read_text()) if p['system']==system]
    points.sort(key=lambda p:p['reaction_coordinate'])
    selected = [points[0], min(points,key=lambda p:abs(p['reaction_coordinate'])), points[-1]]
    lys = list(range(45,58)) if system=='MoriPF' else list(range(43,56))
    substrate = [i for i in range(1,58) if i not in lys]
    fragments = {'substrate':substrate, 'lys':lys}
    references = {}
    for frag, path in [('substrate',ROOT/'data/completed_20261017/Reactant'/REFS[system]),('lys',ROOT/'data/completed_20261017/Lys.log.gz')]:
        text = read_text(path)
        electronic_check(text)
        assert 'Normal termination of Gaussian' in text and 'Error termination' not in text
        assert 'Optimization completed' in text
        freqs = [float(v) for row in re.findall(r'Frequencies --([^\n]+)',text) for v in row.split()]
        assert freqs and min(freqs)>0
        d = stationary(text,str(path.relative_to(ROOT)))
        assert d['atomic_numbers']==[selected[0]['atomic_numbers'][i-1] for i in fragments[frag]]
        d.update(sha256=digest(path), minimum_frequency_cm1=min(freqs), normal_termination_count=text.count('Normal termination of Gaussian'))
        d['atom_mapping']=reference_mapping(d,selected[0],fragments[frag])
        references[frag]=d
    jobs=[]; audits=[]
    for role,p in zip(['reactant_side','TS','product_side'],selected):
        source=Path(p['source']); assert digest(source)==p['source_sha256']
        text=read_text(source); electronic_check(text)
        lines=text.splitlines(); z,xyz=orientation(lines,p['geometry_line']-1)
        assert z==p['atomic_numbers'] and np.max(np.abs(np.array(xyz)-p['xyz']))<1e-10
        assert abs(float(re.search(r'=\s*([-\d.]+)',lines[p['energy_line']-1])[1])-p['energy_hartree'])<1e-9
        if role!='TS': assert 'Point Number:' in lines[p['source_line']-1]
        a=np.array(p['xyz']); n=lys[0]-1
        nh=[float(np.linalg.norm(a[n]-a[h-1])) for h in lys if p['atomic_numbers'][h-1]==1 and np.linalg.norm(a[n]-a[h-1])<1.3]
        assert len(nh)==2, 'Proton assignment changed; manual review required'
        audits.append({'role':role,'point_id':p['point_id'],'coordinate':p['reaction_coordinate'], 'source':str(source.relative_to(ROOT)), 'source_sha256':p['source_sha256'],'energy_line':p['energy_line'],'geometry_line':p['geometry_line'],'full_energy_hartree':p['energy_hartree'],'normal_termination_in_source':'Normal termination of Gaussian' in text,'accepted_irc_marker':p['source_line'],'geometry_kind':p['geometry_kind'],'C20_N_distance_A':float(np.linalg.norm(a[19]-a[n])),'C20_O22_distance_A':float(np.linalg.norm(a[19]-a[21])),'NH_distances_A':nh,'merged_log_note':text.splitlines()[:3] if text.startswith(' MERGE NOTE') else []})
        for frag,indices in fragments.items():
            folder=ROOT/'irc_eda/outputs/smd_di_jobs'/system/(role[:1]+'_'+frag);folder.mkdir(parents=True,exist_ok=True)
            atoms=[p['atomic_numbers'][i-1] for i in indices]; coords=[p['xyz'][i-1] for i in indices]
            body='%chk=fragment.chk\n%nprocshared=4\n%mem=4GB\n#p M062X/6-311++G(d,p) 5D 7F SCRF=(SMD,Solvent=Water) NoSymm SCF=(Tight,MaxCycle=200) Integral=UltraFine\n\n'+system+' '+role+' '+frag+' fixed fragment SMD difference analysis\n\n0 1\n'
            body+='\n'.join(f'{TABLE.element[z]} '+ ' '.join(f'{v:.10f}' for v in r) for z,r in zip(atoms,coords))+'\n\n'
            inp=folder/'fragment.gjf'
            if inp.exists(): assert inp.read_text()==body, 'Refusing to overwrite different input'
            else: inp.write_text(body,encoding='ascii')
            jobs.append({'role':role,'fragment':frag,'folder':str(folder.relative_to(ROOT)),'input_sha256':digest(inp),'atomic_numbers':atoms,'xyz':coords,'original_atom_indices':indices,'charge':0,'multiplicity':1})
    manifest={'analysis':'SMD fixed-fragment distortion/interaction differences, uncorrected BSSE','system':system,'references':references,'points':audits,'jobs':jobs,'fragments':fragments,'route':body.splitlines()[3]}
    save(OUT/system/'manifest.json',manifest)
    save(REPORT/system/'source_audit.json',manifest)
    print(system,'prepared 6 fragment SPs; reused 3 full + 2 reference energies',flush=True)


def validate(job):
    folder=ROOT/job['folder']; path=folder/'fragment.out';text=read_text(path)
    assert digest(folder/'fragment.gjf')==job['input_sha256'], 'Input hash mismatch'
    assert 'Normal termination of Gaussian' in text and 'Error termination' not in text
    assert 'Convergence failure' not in text
    electronic_check(text)
    d=stationary(text,str(path.relative_to(ROOT)))
    assert d['atomic_numbers']==job['atomic_numbers']
    delta=float(np.max(np.abs(np.array(d['xyz'])-job['xyz'])))
    assert delta<2e-6, delta
    assert text.count('SCF Done:')==1, 'Expected one fixed-geometry SCF'
    ne=sum(job['atomic_numbers'])
    counts=re.search(r'(\d+) alpha electrons\s+(\d+) beta electrons',text)
    assert counts and int(counts[1])==int(counts[2]) and int(counts[1])+int(counts[2])==ne
    return {'status':'SUCCESS','energy_hartree':d['energy_hartree'],'max_coordinate_difference_A':delta,'output_sha256':digest(path),'input_sha256':job['input_sha256'],'scf_done_count':text.count('SCF Done:'),'normal_termination':True,'conditions':d['conditions'],'scf_stability':'not tested; restricted neutral singlet','output':str(path.relative_to(ROOT))}


def runjob(job):
    folder=ROOT/job['folder']; inp=folder/'fragment.gjf'; assert digest(inp)==job['input_sha256']
    status=folder/'validation.json'
    if status.exists():
        old=json.loads(status.read_text())
        if old.get('status')=='SUCCESS':
            result=validate(job);assert result['output_sha256']==old['output_sha256'];return result
    if (folder/'fragment.out').exists(): raise RuntimeError('Existing unvalidated attempt: '+str(folder))
    env=os.environ.copy();env.update(GAUSS_EXEDIR='C:\\g16w',GAUSS_BSDDIR='C:\\g16w\\basis',GAUSS_LEXEDIR='C:\\g16w\\bin',GAUSS_SCRDIR=str(folder/'scratch'))
    env['PATH']='C:\\g16w;'+env['PATH'];(folder/'scratch').mkdir(exist_ok=True)
    save(status,{'status':'RUNNING','started':time.time(),'input_sha256':job['input_sha256']})
    print('START',job['role'],job['fragment'],flush=True)
    started=time.time()
    with (folder/'launcher.txt').open('w') as log:
        process=subprocess.run(['C:\\g16w\\g16.exe','fragment.gjf'],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT)
    try:
        assert process.returncode==0,process.returncode
        result=validate(job)
    except Exception as e:
        result={'status':'FAILED','reason':str(e),'returncode':process.returncode}
    result['elapsed_seconds']=time.time()-started;save(status,result)
    print('END',job['role'],job['fragment'],result['status'],result.get('energy_hartree'),flush=True)
    if result['status']!='SUCCESS':raise RuntimeError(result)
    return result


def analyze(system):
    m=json.loads((OUT/system/'manifest.json').read_text()); rows=[]
    ref=sum(x['energy_hartree'] for x in m['references'].values())
    for p in m['points']:
        energies={}
        for j in m['jobs']:
            if j['role']==p['role']: energies[j['fragment']]=validate(j)['energy_hartree']
        dist=sum(energies.values())-ref; inter=p['full_energy_hartree']-sum(energies.values()); total=p['full_energy_hartree']-ref
        row={k:p[k] for k in ['role','point_id','coordinate','full_energy_hartree','C20_N_distance_A']}
        row.update(substrate_hartree=energies['substrate'],lys_hartree=energies['lys'],reference_sum_hartree=ref,distortion_kcal_mol=dist*HARTREE_KCAL,interaction_difference_kcal_mol=inter*HARTREE_KCAL,total_relative_kcal_mol=total*HARTREE_KCAL,identity_residual_hartree=total-dist-inter)
        for frag in energies:
            row[frag+'_distortion_kcal_mol']=(energies[frag]-m['references'][frag]['energy_hartree'])*HARTREE_KCAL
        assert abs(row['identity_residual_hartree'])<1e-10
        rows.append(row)
    dest=REPORT/system;dest.mkdir(parents=True,exist_ok=True)
    with (dest/'pilot.csv').open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    save(dest/'pilot.json',rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4.5))
    for col,label in [('distortion_kcal_mol','Distortion (SMD)'),('interaction_difference_kcal_mol','Interaction difference (SMD)'),('total_relative_kcal_mol','Total vs separated references')]:
        ax.plot([r['coordinate'] for r in rows],[r[col] for r in rows],'o--',label=label)
    ax.set(xlabel='IRC coordinate / sqrt(amu) bohr',ylabel='Energy / kcal mol$^{-1}$',title=system+' — 3-point fixed-fragment pilot')
    ax.axhline(0,color='grey',linewidth=.5);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(dest/'pilot.png',dpi=180);plt.close(fig)
    delta_dist=rows[1]['distortion_kcal_mol']-rows[0]['distortion_kcal_mol']
    delta_int=rows[1]['interaction_difference_kcal_mol']-rows[0]['interaction_difference_kcal_mol']
    delta_total=rows[1]['total_relative_kcal_mol']-rows[0]['total_relative_kcal_mol']
    text=f'# {system}: SMD固定断片・変形/相互作用差分pilot\n\n'
    text+='3点・断片6件の正常終了出力を検証した。全系3値と分離参照2値は条件照合して再利用。単位はkcal/mol。\n\n'
    text+='| 点 | 変形 | 相互作用差分 | 分離参照からの全系差 |\n|---|---:|---:|---:|\n'
    for r in rows:
        text+=f"| {r['role']} | {r['distortion_kcal_mol']:.6f} | {r['interaction_difference_kcal_mol']:.6f} | {r['total_relative_kcal_mol']:.6f} |\n"
    text+=f'\n反応物側IRC端点→TSの全系差は {delta_total:.6f} kcal/mol。その内訳として変形項の差 {delta_dist:.6f} と相互作用差分の差 {delta_int:.6f} を得た。これは定義したSMD断片参照に対する数値比較であり、純粋な静電安定化やIQA原子対成分を意味しない。\n'
    text+=f"\n恒等式残差の最大値は {max(abs(r['identity_residual_hartree']) for r in rows):.3e} Eh。算術整合性の確認であり、分解の一意性・物理的妥当性の証明ではない。\n"
    text+='\n固定中性一重項の断片は仮想的な変形反応物である。BSSE未補正、SMD空洞変化を含む。SCF安定性・別スピン状態・全配座探索は未検証。ZPE/熱補正は使用していない。端点は利用可能な受理済みIRC端の座標であり、最適化済み反応物/生成物の意味ではない。\n'
    text+='\nsource_audit.jsonに参照構造の対応表、原ログhash・行番号・正常終了の有無を保存した。元の結合IRCに正常終了がない場合、受理済みステップの収束SCFを再利用したもので、元IRC全体の正常終了を認定していない。pilot.pngの破線は3点間の視覚的補助であり、中間点の計算や補間エネルギーを表さない。\n'
    (dest/'report_ja.md').write_text(text,encoding='utf8')
    save(dest/'fragment_validation.json',[{'job':j,'validation':validate(j)} for j in m['jobs']])
    print(json.dumps(rows,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run','analyze']);parser.add_argument('system',choices=list(REFS));args=parser.parse_args()
    if args.action=='prepare':prepare(args.system)
    elif args.action=='analyze':analyze(args.system)
    else:
        manifest=json.loads((OUT/args.system/'manifest.json').read_text())
        # Gaussian 16W pilot: serialize launches; never continue after failure.
        for job in manifest['jobs']:
            runjob(job)
