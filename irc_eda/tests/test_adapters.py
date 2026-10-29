"""Synthetic fixtures test adapters only; they are NEVER reported as QC results."""
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
from irc_eda.parsing import parse_irc,conditions,signature,rmsd,fchk_info,read_text,HARTREE_KCAL
from irc_eda.upstream import (ROOT,parse_eda,parse_iqa,iqa_input,iqa_method,fragment_text,gaussian_input,
                             eda_template,matrix_block,sp_route)
from irc_eda.pipeline import Pipeline,load_config,dump

REPO=ROOT.parent


def orientation(x):
    return f''' Input orientation:
 ---------------------------------------------------------------------
 Center Atomic Atomic Coordinates (Angstroms)
 Number Number Type X Y Z
 ---------------------------------------------------------------------
 1 1 0 0.0 0.0 0.0
 2 1 0 {x:.6f} 0.0 0.0
 ---------------------------------------------------------------------
'''


def gaussian_fixture():
    return ''' Gaussian 16
 #p irc=(forward,maxpoints=2) hf/sto-3g
 -----------
 Charge = 0 Multiplicity = 1
 Integration on MW PES will use step size of 0.1 sqrt(amu)*bohr.
 Energy From Chk = -1.0
 Point Number: 0 Path Number: 1
 Point Number 1 in FORWARD path direction.
'''+orientation(.75)+''' SCF Done: E(RHF) = -1.010 A.U.
'''+orientation(.76)+''' SCF Done: E(RHF) = -1.020 A.U.
 Delta-x Convergence Met
 Point Number: 1 Path Number: 1
 NET REACTION COORDINATE UP TO THIS POINT = 0.100
 Point Number 2 in FORWARD path direction.
'''+orientation(.77)+''' SCF Done: E(RHF) = -1.030 A.U.
 Delta-x Convergence NOT Met
 Error termination
'''


def test_no_scf_iteration_as_point_or_fabricated_ts(tmp_path):
    file=tmp_path/'synthetic.log';file.write_text(gaussian_fixture())
    with patch('irc_eda.parsing.cclib_read',return_value=(None,'synthetic fixture')):
        p,a=parse_irc(file,'fixture')
    assert len(p)==1 and p[0]['energy_hartree']==-1.020
    assert p[0]['xyz'][1][0]==.76 and p[0]['point_index']==1
    assert a['rejected'][0]['point']==2
    assert p[0]['reaction_coordinate']==.1


def test_reverse_separate_file_uses_announced_direction(tmp_path):
    file=tmp_path/'synthetic.log';file.write_text(gaussian_fixture().replace('FORWARD','REVERSE'))
    with patch('irc_eda.parsing.cclib_read',return_value=(None,'synthetic fixture')):
        p,_=parse_irc(file,'fixture')
    assert p[0]['direction']=='reverse' and p[0]['reaction_coordinate']==-.1
    assert p[0]['path_number']==2 and p[0]['gaussian_path_number']==1


def test_conflicting_duplicate_is_error(tmp_path):
    text=gaussian_fixture().split(' Point Number 2')[0]
    file=tmp_path/'synthetic.log';file.write_text(text+orientation(.8)+''' SCF Done: E(RHF) = -1.05 A.U.
 Point Number: 1 Path Number: 1
 NET REACTION COORDINATE = 0.1
''')
    with patch('irc_eda.parsing.cclib_read',return_value=(None,'synthetic fixture')),pytest.raises(ValueError,match='Conflicting duplicate'):
        parse_irc(file,'fixture')


def gas_point():
    return dict(system='fixture',point_id='fixture_p1_0001',irc_index=0,source='synthetic',source_sha256='synthetic',
                conditions=dict(method='B3LYP',basis='6-31G(d)',charge=0,multiplicity=1,solvent='gas',
                                route='#p irc=(rcfc,maxpoints=5) b3lyp/6-31G(d) geom=allcheck guess=read scf=tight'),
                atomic_numbers=[1,1],xyz=[[0,0,0],[0,0,.74]],reaction_coordinate=0.,coordinate_unit='sqrt(amu)*bohr',
                energy_hartree=-1.,wavefunction_matches=[],readiness={'iqa_status':'READY','iqa_reasons':[]})


def test_sp_removes_optimization_preserves_conditions():
    p=gas_point();p['conditions']['route']='#p opt=(ts,calcfc) freq b3lyp/6-31G(d) geom=allcheck guess=read scf=tight scrf=(smd,solvent=water)'
    p['conditions']['solvent']='scrf=(smd,solvent=water)'
    text=gaussian_input(p,{'threads':2});route=text.splitlines()[3]
    assert all(x not in route for x in ('opt=','freq','geom=','guess='))
    assert 'scf=tight' in route and 'scrf=(smd,solvent=water)' in route
    assert ' sp ' in route and 'B3LYP/6-31G(d)' in route


@pytest.mark.parametrize('method',['M062X','M06-2X','TPSS','SCAN'])
def test_meta_gga_never_silently_replaced(method):
    c=gas_point()['conditions']|{'method':method}
    with pytest.raises(ValueError,match='NOT IMPLEMENTED'):iqa_method(c,{'method':method,'family':'HYB_GGA','EXC_FUNCTIONAL':450,'verified_equivalence_source':'invalid override'})


def test_solvent_is_blocked_for_both_schemes():
    p=gas_point();p['conditions']['solvent']='scrf=(smd,solvent=water)'
    with pytest.raises(ValueError):iqa_input(p['conditions'],{})
    with pytest.raises(ValueError):eda_template(p)


def test_fragments_need_partition_charge_spin():
    p=gas_point()
    fr=[dict(atoms=[1],charge=0,multiplicity=2),dict(atoms=[2],charge=0,multiplicity=-2)]
    assert fragment_text(fr,p)=='2\n0 2\n1\n0 -2\n2\n'
    with pytest.raises(ValueError):fragment_text(None,p)
    with pytest.raises(ValueError):fragment_text([dict(atoms=[1],charge=0,multiplicity=2)]*2,p)
    with pytest.raises(ValueError):fragment_text([dict(atoms=[1],charge=0,multiplicity=1),fr[1]],p)
    with pytest.raises(ValueError):fragment_text([fr[0],dict(atoms=[2],charge=0,multiplicity=2)],p)


def test_official_recorded_sobeda_all_terms_and_closure():
    source=ROOT/'vendor/easyEDA/scr/1.log'
    if not source.exists():pytest.skip('bootstrap upstream required')
    rows,v,stages=parse_eda(source.read_text())
    assert len(rows)==9 and abs(v['closure_kcal_mol'])<1e-8
    d={r['component']:r['energy_kcal_mol'] for r in rows}
    assert d['E_rep']==287. and d['E_xrep']==210.17 and d['E_tot']==-127.49
    assert any(s['stage']=='final state' for s in stages)
    with pytest.raises(ValueError):parse_eda(source.read_text().replace('Orbital (E_orb):','Missing:'))


def matrix_text(title,m):
    n=len(m);text=' '+title+'\n'
    for start in range(0,n,6):
        cols=list(range(start,min(n,start+6)))
        text+=' ----------\n          '+''.join(f'  {i+1:3d} H   ' for i in cols)+'\n ----------\n'
        for i in range(n):text+=f'  {i+1:3d} H  '+''.join(f'{m[i,j]:20.13f}' for j in cols)+'\n'
        text+=' ----------\n'
    return text


def test_apost_matrix_multiblock_and_unique_pair_sum():
    n=7;m=np.eye(n)*-1
    m[0,1]=m[1,0]=.2
    text=matrix_text('FUZZY ATOMS Hartree-Fock ENERGY COMPONENTS',m)+''' Total integrated energy : -6.8000000
 Total energy in Fchk file: -6.8000000
 ...Normal Termination of APOST-3D...
'''
    atoms,pairs,v=parse_iqa(text,n)
    assert v['closure_hartree']==pytest.approx(0)
    assert sum(r['energy_hartree'] for r in atoms if r['component']=='intra_atomic')==-7
    assert len([r for r in pairs if r['component']=='interatomic'])==21
    assert all(r['energy_hartree'] is None for r in pairs if r['component']=='xc')
    with pytest.raises(ValueError):parse_iqa(text.replace('Normal Termination','Abnormal termination'),n)


def test_corrupt_matrix_is_rejected():
    with pytest.raises(ValueError):matrix_block(['title',' 1 H 2 H','1 H -1.0'],0,2)


def test_resume_fingerprint_output_integrity_and_failure(tmp_path):
    c=dict(input_root=str(tmp_path/'input'),output_dir=str(tmp_path/'output'),systems={'fixture':{}},parallel_jobs=1)
    (tmp_path/'input').mkdir();pipe=Pipeline(c);p=gas_point();calls=[]
    def execute(point,work):calls.append(1);(work/'analysis.apost').write_text('synthetic successful execution')
    result={'atomic':[],'pairwise':[],'validation':{'closure_hartree':0.}}
    with patch.object(pipe,'execute_iqa',side_effect=execute),patch.object(pipe,'collect',return_value=result):
        a=pipe.run_point(p,'iqa');b=pipe.run_point(p,'iqa')
        assert a['status']=='SUCCESS' and b['resumed'] and len(calls)==1
        (Path(a['work_dir'])/'analysis.apost').write_text('tampered')
        assert not pipe.run_point(p,'iqa')['resumed'] and len(calls)==2
        p2=p|{'energy_hartree':-1.1}
        assert not pipe.run_point(p2,'iqa')['resumed'] and len(calls)==3
    with patch.object(pipe,'execute_iqa',side_effect=RuntimeError('synthetic failure')):
        failed=pipe.run_point(p|{'energy_hartree':-1.2},'iqa')
    assert failed['status']=='FAILED'
    assert 'synthetic failure' in failed['reasons'][0]


def test_output_directory_guard(tmp_path):
    with pytest.raises(ValueError):Pipeline(dict(input_root=str(tmp_path),output_dir=str(tmp_path)))


@pytest.mark.parametrize('system,count',[('MoriPF',64),('RAFA14',60),('RAFA15',57),('RAFA16',71)])
def test_real_irc_cclib_order_and_existing_profile(system,count):
    c=load_config(ROOT/'config/project.yaml');s=c['systems'][system]
    if not (REPO/s['irc']).exists():pytest.skip('project input unavailable')
    ps,a=parse_irc(REPO/s['irc'],system,s['coordinate_sign'],REPO/s['ts'])
    assert len(ps)==count and sum(p['direction']=='TS' for p in ps)==1
    assert all(p['cclib_crosscheck'] for p in ps if p['direction']!='TS')
    # Independently stored existing Gaussian profile, used only as a regression comparator.
    import pandas as pd
    previous=pd.read_csv(REPO/f'analysis4_force/outputs/csv/reaction_force_{system}.csv')
    assert len(previous)==len(ps)
    np.testing.assert_allclose([p['reaction_coordinate'] for p in ps],previous['signed_irc_coordinate'],atol=1e-9,rtol=0)
    np.testing.assert_allclose([p['energy_hartree'] for p in ps],previous['energy_hartree'],atol=1e-9,rtol=0)
    from irc_eda.parsing import TABLE
    for p in (ps[0],min(ps,key=lambda x:abs(x['reaction_coordinate'])),ps[-1]):
        assert p['atomic_numbers']==ps[0]['atomic_numbers']
        assert len(p['xyz'])==len(p['atomic_numbers'])


def test_real_fchk_not_matched_by_filename(tmp_path):
    file=REPO/'data/completed_20261017/TS/RAFA14_TS2_NBO.fchk'
    if not file.exists():pytest.skip('project fchk unavailable')
    c=load_config(ROOT/'config/project.yaml');c['output_dir']=str(tmp_path/'results')
    p=Pipeline(c);wf=fchk_info(file);p.waveforms=[dict(wf,source=str(file),sha256='test')]
    source=REPO/c['systems']['RAFA14']['ts'];ts=__import__('irc_eda.parsing',fromlist=['stationary']).stationary(read_text(source),source)
    hits=p.match_waveforms(ts)
    assert len(hits)==1 and hits[0]['rmsd_angstrom']<1e-6
    assert not hits[0]['conditions_verified'] # chkbasis is not evidence of basis-name identity
    ts['atomic_numbers']=ts['atomic_numbers'][::-1]
    assert not p.match_waveforms(ts)
