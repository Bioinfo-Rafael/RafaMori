"""Separate official recorded-output checks from genuinely executed samples."""
from pathlib import Path
import shutil
from .upstream import ROOT,parse_eda,parse_iqa
from .parsing import fchk_info
from .pipeline import dump,csv_write


def samples(p):
    folder=p.out/'official_samples';folder.mkdir(exist_ok=True)
    records=[]
    original=ROOT/'vendor/easyEDA/scr/1.log'
    if original.exists():
        rows,v,stages=parse_eda(original.read_text())
        csv_write(folder/'easyEDA_recorded_components.csv',rows)
        dump(folder/'easyEDA_recorded_validation.json',v)
        records.append(dict(sample='easyEDA/scr/1.log',status='RECORDED_OUTPUT_PARSE_PASS',
                            quantum_calculation_executed=False,validation=v))
    sample='C2H6-B3LYP';source=ROOT/'vendor/APOST3D/compiler-testset'
    if not p.tool('apost3d'):
        records.append(dict(sample=sample,status='BLOCKED',quantum_calculation_executed=False,reason='APOST3D executable missing'))
    else:
        work=folder/'apost3d';work.mkdir(exist_ok=True)
        for ext in ('.inp','.fchk'):shutil.copy2(source/(sample+ext),work/(sample+ext))
        try:
            p.command([p.tool('apost3d'),sample],work,'sample.apost',env=__import__('os').environ|{'OMP_NUM_THREADS':str(p.config.get('threads',1)),'KMP_STACKSIZE':'100m'})
            a,b,v=parse_iqa((work/'sample.apost').read_text(),fchk_info(work/(sample+'.fchk'))['atom_count'])
            dump(work/'parsed.json',dict(atomic=a,pairwise=b,validation=v))
            records.append(dict(sample=sample,status='EXECUTED_SAMPLE_PARSED',quantum_calculation_executed=True,validation=v))
        except Exception as e:records.append(dict(sample=sample,status='FAILED',quantum_calculation_executed=True,reason=str(e)))
    deps=[n for n in ('g16','Multiwfn','formchk','unfchk','dos2unix','bc') if not p.tool(n)]
    if p.env['platform']!='Linux' or 'GNU' not in (p.env['tools']['sed']['version'] or ''):deps.append('Linux/WSL with GNU sed')
    if deps:
        records.append(dict(sample='easyEDA coordinates/1.xyz',status='BLOCKED',quantum_calculation_executed=False,reason='Missing: '+', '.join(deps)))
    else:
        work=folder/'easyEDA';work.mkdir(exist_ok=True)
        source=ROOT/'vendor/easyEDA'
        for name in ('fragment.txt','template.gjf'):shutil.copy2(source/name,work/name)
        shutil.copy2(source/'coordinates/1.xyz',work/'1.xyz');shutil.copy2(source/'source/EDA_run.sh',work/'EDA_run.sh')
        try:
            p.command([p.tool('bash'),'EDA_run.sh','1.xyz'],work,'sobEDA.log')
            rows,v,stages=parse_eda((work/'sobEDA.log').read_text());dump(work/'parsed.json',dict(components=rows,validation=v,stages=stages))
            records.append(dict(sample='easyEDA coordinates/1.xyz',status='EXECUTED_SAMPLE_PARSED',quantum_calculation_executed=True,validation=v))
        except Exception as e:records.append(dict(sample='easyEDA coordinates/1.xyz',status='FAILED',quantum_calculation_executed=True,reason=str(e)))
    dump(folder/'sample_status.json',records)
    return records
