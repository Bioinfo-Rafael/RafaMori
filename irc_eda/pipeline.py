"""Auditable orchestration; all quantum numerics remain in external programs."""
from __future__ import annotations
import concurrent.futures
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
import yaml
from . import __version__
from .parsing import (read_text,digest,extension,conditions,stationary,fchk_info,parse_irc,xyz_text,
                      rmsd,signature,FLOAT,number,HARTREE_KCAL)
from .upstream import (ROOT,PINS,environment,iqa_method,iqa_input,gaussian_input,fragment_text,
                       eda_template,parse_iqa,parse_eda,git_revision)


def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    tmp.replace(path)


def csv_write(path,rows,fields=None):
    fields=fields or list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader()
        for r in rows:
            w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})


def load_config(path):
    path=Path(path).resolve();c=yaml.safe_load(path.read_text())
    base=path.parent
    for key in ('input_root','output_dir'):
        c[key]=str((base/c[key]).resolve())
    c['config_file']=str(path)
    if int(c.get('parallel_jobs',1))<1 or int(c.get('threads',1))<1:raise ValueError('Positive concurrency required')
    for name in c.get('systems',{}):
        if not re.fullmatch('[A-Za-z0-9_-]+',name):raise ValueError('Unsafe system identifier')
    for name,value in c.get('executables',{}).items():
        if '/' in value:c['executables'][name]=str((base/value).resolve())
    return c


class Pipeline:
    def __init__(self,config):
        self.config=config;self.root=Path(config['input_root']);self.out=Path(config['output_dir'])
        if self.out==self.root or self.out in self.root.parents:
            raise ValueError('Output must not be the input directory or its ancestor')
        marker=self.out/'.irc_eda_output'
        if self.out.exists() and any(self.out.iterdir()) and not marker.exists():
            raise ValueError('Existing output directory is not owned by this adapter; choose a new directory')
        self.out.mkdir(parents=True,exist_ok=True);marker.touch(exist_ok=True)
        self.env=environment(config)
        self.provenance={'adapter_version':__version__,'adapter_git_commit':self.env['adapter_git_commit'],
                         'adapter_source_sha256':hashlib.sha256(b''.join(p.read_bytes() for p in sorted(ROOT.glob('*.py')))).hexdigest(),
                         'environment':self.env,'config':config,
                         'created_at':datetime.now(ZoneInfo('Asia/Tokyo')).isoformat()}
        dump(self.out/'environment.json',self.env)
        dump(self.out/'provenance.json',self.provenance)
        self.points=[];self.inventory=[];self.audits=[];self.waveforms=[]

    def source(self,path):
        p=(self.root/path).resolve()
        if not p.exists() and Path(str(p)+'.gz').exists():p=Path(str(p)+'.gz')
        return p

    def tool(self,name):return self.env['tools'][name]['path']

    def inspect(self):
        paths=[]
        skip={'vendor','.venv','.git','__pycache__','.pytest_cache'}
        excluded={(self.root/d).resolve() for d in self.config.get('inventory_exclude',[])}
        for directory,dirs,files in os.walk(self.root):
            dirs[:]=[d for d in dirs if d not in skip and (Path(directory)/d).resolve()!=self.out and (Path(directory)/d).resolve() not in excluded]
            for name in files:
                p=Path(directory)/name
                if extension(p) in ('.log','.out','.chk','.fchk','.fch','.gjf','.com','.xyz'):paths.append(p)
        # Scan wavefunctions/geometry once. Duplicate compressed originals remain inventory records.
        for path in sorted(paths):
            ext=extension(path); row=dict(source=str(path),extension=ext,bytes=path.stat().st_size,sha256=digest(path),
                                         calculation_type='UNKNOWN',normal_termination=None,conditions=None,
                                         wavefunction_present=ext in ('.chk','.fchk','.fch'),parse_status='OK')
            row['provenance_file']=str(self.out/'provenance.json')
            try:
                if ext=='.chk':
                    row.update(parse_status='OPAQUE_BINARY',note='formchk required; basename is not correspondence evidence')
                    if self.tool('formchk'):
                        converted=self.out/'conversions'/row['sha256']
                        converted.mkdir(parents=True,exist_ok=True)
                        target=converted/'checkpoint.fchk'
                        if not target.exists():
                            # Binary inputs, including compressed checkpoints, stay untouched.
                            local=converted/'checkpoint.chk'
                            if path.suffix.lower()=='.gz':
                                import gzip
                                with gzip.open(path,'rb') as src,local.open('wb') as dst:shutil.copyfileobj(src,dst)
                            else:shutil.copy2(path,local)
                            self.command([self.tool('formchk'),'checkpoint.chk','checkpoint.fchk'],converted,'formchk.log')
                        wf=fchk_info(target)
                        row.update(parse_status='CONVERTED',conversion_file=str(target),conditions=wf['conditions'],atom_count=wf['atom_count'])
                        self.waveforms.append(dict(wf,source=str(target),sha256=digest(target),original_checkpoint=str(path)))
                elif ext in ('.fchk','.fch'):
                    wf=fchk_info(path); row.update({k:v for k,v in wf.items() if k not in ('xyz','atomic_numbers')})
                    row['calculation_type']=read_text(path).splitlines()[1].split()[0]
                    row['atomic_numbers']=wf['atomic_numbers']
                    self.waveforms.append(dict(wf,source=str(path),sha256=row['sha256']))
                elif ext in ('.log','.out'):
                    text=read_text(path);c=conditions(text);row['conditions']=c
                    row['derived_analysis_log']='DERIVED ANALYSIS LOG' in text[:1000]
                    row['composite_log']='MERGE NOTE' in text[:1000] or len(re.findall('Entering Gaussian System',text))>1
                    end=re.findall(r'Normal termination of Gaussian|Error termination',text)
                    row['normal_termination']=bool(end and end[-1].startswith('Normal')) if end else None
                    row['normal_termination_events']=sum(x.startswith('Normal') for x in end)
                    row['error_termination_events']=sum(x.startswith('Error') for x in end)
                    compact=c['route'].lower()
                    row['calculation_type']='IRC' if c['is_irc'] else 'TS' if re.search(r'opt=.*\bts\b',compact) else 'Optimization' if 'opt' in compact else 'Single Point' if c['method'] else 'UNKNOWN'
                    markers=re.findall(r'Point Number:\s*(\d+)\s+Path Number:\s*(\d+)\s*\n(?:[^\n]*\n){0,3}?\s*NET REACTION COORDINATE[^\n]*',text)
                    row['completed_irc_points']=len(set(markers))
                    row['forward_points']=len({n for n,b in markers if b=='1' and n!='0'})
                    row['reverse_points']=len({n for n,b in markers if b=='2' and n!='0'})
                    row['gaussian_version']=next(iter(re.findall(r'Gaussian 16:\s*([^\n]+)',text)),None)
                    if c['method']:
                        s=stationary(text,path)
                        row.update(atom_count=len(s['atomic_numbers']),atomic_numbers=s['atomic_numbers'],energy_hartree=s['energy_hartree'])
                elif ext in ('.gjf','.com'):
                    text=read_text(path);from .parsing import routes
                    r=routes(text,True);c=conditions(text,r[0] if r else None)
                    cm=re.search(r'^\s*(-?\d+)\s+(\d+)\s*$',text,re.M)
                    if cm:c.update(charge=int(cm[1]),multiplicity=int(cm[2]))
                    row['conditions']=c;row['calculation_type']='IRC' if c['is_irc'] else 'Optimization' if 'opt' in c['route'].lower() else 'Single Point'
                else:
                    row['atom_count']=int(read_text(path).splitlines()[0]);row['note']='XYZ alone contains no electronic conditions'
            except Exception as e:row.update(parse_status='ERROR',error=f'{type(e).__name__}: {e}')
            self.inventory.append(row)
        # Required file inventory is persisted BEFORE path extraction.
        dump(self.out/'inventory.json',self.inventory);csv_write(self.out/'inventory.csv',self.inventory)
        for system,spec in self.config['systems'].items():
            files=spec.get('irc_files',[spec.get('irc')]);allpoints=[];audits=[]
            for index,path in enumerate(files):
                p,a=parse_irc(self.source(path),system,spec.get('coordinate_sign',1),self.source(spec['ts']) if spec.get('ts') and index==0 else None)
                allpoints+=p;audits.append(a)
            unique={}
            for p in allpoints:
                key=(p['path_number'],p['point_index'])
                if key in unique:
                    old=unique[key]
                    if signature(old['conditions'])!=signature(p['conditions']) or abs(old['energy_hartree']-p['energy_hartree'])>1e-7 or rmsd(old['xyz'],p['xyz'])>2e-5 or abs(old['reaction_coordinate']-p['reaction_coordinate'])>1e-6:
                        raise ValueError(f'Inconsistent merged files: {system} {key}')
                else:unique[key]=p
            ps=sorted(unique.values(),key=lambda p:p['reaction_coordinate'])
            if len({signature(p['conditions']) for p in ps})!=1:raise ValueError('Mixed conditions across IRC files')
            for i,p in enumerate(ps):p['irc_index']=i
            ref=spec.get('reference_point',ps[0]['point_id'])
            rp=next((p for p in ps if p['point_id']==ref),None)
            if rp is None:raise ValueError('Reference point not found')
            for p in ps:
                p['provenance_file']=str(self.out/'provenance.json')
                p['reference_point']=ref;p['relative_energy_kcal_mol']=(p['energy_hartree']-rp['energy_hartree'])*HARTREE_KCAL
                p['wavefunction_matches']=self.match_waveforms(p)
                p['readiness']=self.readiness(p)
            self.points+=ps;self.audits.append({'system':system,'files':audits,'reference_point':ref})
        dump(self.out/'points.json',self.points);dump(self.out/'parse_audit.json',self.audits)
        # Per-file correspondence listing complements raw inventory.
        for row in self.inventory:
            row['irc_point_matches']=[p['point_id'] for p in self.points for w in p['wavefunction_matches'] if w['source']==row['source']]
        dump(self.out/'inventory.json',self.inventory);csv_write(self.out/'inventory.csv',self.inventory)
        self.write_profile()
        ready=[dict(system=p['system'],point_id=p['point_id'],irc_index=p['irc_index'],source=p['source'],
                    conditions=p['conditions'],provenance_file=str(self.out/'provenance.json'),**p['readiness']) for p in self.points]
        csv_write(self.out/'readiness.csv',ready);dump(self.out/'readiness.json',ready)
        return self.points

    def match_waveforms(self,p):
        result=[]
        for w in self.waveforms:
            if w['atomic_numbers']!=p['atomic_numbers']:continue
            distance=rmsd(w['xyz'],p['xyz'])
            if distance>self.config.get('geometry_tolerance_angstrom',2e-4):continue
            ec=None if w['energy_hartree'] is None else w['energy_hartree']-p['energy_hartree']
            c=w['conditions'].copy()
            # An explicit basis declaration is an auditable assertion, never guessed from a filename.
            assertion=self.config.get('wavefunction_conditions',{}).get(os.path.relpath(w['source'],self.root),{})
            if assertion:
                if not assertion.get('evidence'):raise ValueError('Wavefunction basis assertion needs evidence')
                c.update({k:v for k,v in assertion.items() if k!='evidence'})
            verified=signature(c)==signature(p['conditions']) and ec is not None and abs(ec)<self.config.get('energy_tolerance_hartree',2e-6)
            result.append(dict(source=w['source'],sha256=w['sha256'],rmsd_angstrom=distance,energy_difference_hartree=ec,
                               conditions_verified=verified,conditions=c,basis_signature=w['basis_signature'],
                               has_mo=w['has_mo'],has_density=w['has_density'],has_energy_terms=w['has_energy_terms'],
                               status='VERIFIED' if verified else 'GEOMETRY_MATCH_CONDITIONS_UNVERIFIED'))
        return result

    def readiness(self,p):
        iqa=[];eda=[];data_status='RECALC_REQUIRED'
        matches=[w for w in p['wavefunction_matches'] if w['conditions_verified'] and w['has_mo'] and w['has_density']]
        if matches:data_status='CONVERSION_REQUIRED' if matches[0]['source'].endswith('.gz') else 'READY'
        elif p['wavefunction_matches']:
            data_status='BLOCKED'
            iqa.append('Geometry-matched wavefunction exists, but electronic conditions are unverified')
        try:iqa_method(p['conditions'],self.config.get('iqa',{}).get('libxc'))
        except ValueError as e:iqa.append(str(e))
        if p['conditions'].get('solvent')!='gas':iqa.append('Solvation/reaction-field energy not validated in ENPART')
        if not self.tool('apost3d'):iqa.append('APOST3D executable missing (ifort build required)')
        if data_status=='RECALC_REQUIRED' and not self.tool('g16'):iqa.append('Gaussian 16 missing for fixed-geometry SP')
        if not matches and not self.tool('formchk'):iqa.append('formchk missing')
        if (not matches or not matches[0]['has_energy_terms']) and not self.tool('get_energy_g16'):
            iqa.append('get_energy_g16 or pre-augmented fchk required for zero-error scheme')
        if not matches:
            try:gaussian_input(p,self.config)
            except ValueError as e:iqa.append(str(e))
        try:fragment_text(self.config['systems'][p['system']].get('fragments'),p)
        except ValueError as e:eda.append(str(e))
        try:eda_template(p)
        except ValueError as e:eda.append(str(e))
        for name in ('g16','formchk','unfchk','Multiwfn','bc','dos2unix','bash'):
            if not self.tool(name):eda.append(name+' missing')
        sed=self.env['tools']['sed']
        if not sed['version'] or 'GNU' not in sed['version']:eda.append('GNU sed required by upstream sed -i syntax')
        if self.env['platform'] not in ('Linux',):eda.append('Unmodified sobEDA shell execution supported on Linux/WSL')
        for name in ('APOST3D','easyEDA'):
            if self.env['upstream'][name]['actual_commit']!=PINS[name][1] or not self.env['upstream'][name]['tracked_source_unmodified']:
                (iqa if name=='APOST3D' else eda).append('Upstream missing or not pinned: '+name)
        return {'iqa_status':'BLOCKED' if iqa else data_status,'iqa_data_status':data_status,'iqa_reasons':iqa,
                'eda_status':'BLOCKED' if eda else 'RECALC_REQUIRED','eda_data_status':'RECALC_REQUIRED','eda_reasons':eda}

    def write_profile(self):
        rows=[]
        for p in self.points:
            r={k:v for k,v in p.items() if k not in ('xyz','atomic_numbers','wavefunction_matches','readiness')}
            r['atom_count']=len(p['xyz']);r['xyz_file']=str(self.out/'xyz'/f'{p["point_id"]}.xyz');rows.append(r)
        csv_write(self.out/'irc_profile.csv',rows);dump(self.out/'irc_profile.json',rows)

    def ensure_points(self):
        if not self.points:self.inspect()

    def prepare(self):
        self.ensure_points();(self.out/'xyz').mkdir(exist_ok=True)
        preparation=[]
        for p in self.points:
            (self.out/'xyz'/f'{p["point_id"]}.xyz').write_text(xyz_text(p))
            folder=self.out/'prepared'/p['point_id'];folder.mkdir(parents=True,exist_ok=True)
            r={'point_id':p['point_id'],'system':p['system'],'gaussian_input':None,'iqa_input':None,'eda_template':None,'errors':[]}
            # Generate missing-wavefunction jobs even when analysis is blocked, but do not run futile SPs.
            if not any(w['conditions_verified'] for w in p['wavefunction_matches']):
                try:
                    file=folder/'single_point.gjf';file.write_text(gaussian_input(p,self.config));r['gaussian_input']=str(file)
                except ValueError as e:r['errors'].append('SP: '+str(e))
            try:
                inp=iqa_input(p['conditions'],self.config.get('iqa',{}));(folder/'analysis.inp').write_text(inp);r['iqa_input']=str(folder/'analysis.inp')
            except ValueError as e:r['errors'].append('IQA: '+str(e))
            try:
                template=eda_template(p);(folder/'template.gjf').write_text(template);r['eda_template']=str(folder/'template.gjf')
            except ValueError as e:r['errors'].append('EDA: '+str(e))
            try:
                fragments=fragment_text(self.config['systems'][p['system']].get('fragments'),p);(folder/'fragment.txt').write_text(fragments)
            except ValueError as e:r['errors'].append('fragment: '+str(e))
            dump(folder/'point.json',p);preparation.append(r)
        # Config templates are edit-once; user edits survive prepare/run-all.
        template=self.out/'fragments.template.yaml'
        if not template.exists():
            data={'instructions':'Copy fragments into systems.<name>.fragments in your config. 1-based original Gaussian order. No partition is inferred.',
                  'systems':{s:{'atom_count':len(next(p for p in self.points if p['system']==s)['xyz']),
                                'fragments':[{'atoms':[],'charge':None,'multiplicity':None},{'atoms':[],'charge':None,'multiplicity':None}]} for s in self.config['systems']}}
            template.write_text(yaml.safe_dump(data,sort_keys=False,allow_unicode=True))
        dump(self.out/'preparation.json',preparation);csv_write(self.out/'preparation.csv',preparation)
        return preparation

    def pilot(self,system):
        ps=[p for p in self.points if p['system']==system]
        selected={ps[0]['point_id'],ps[-1]['point_id'],min(ps,key=lambda p:abs(p['reaction_coordinate']))['point_id']}
        selected.update(p['point_id'] for p in ps if any(w['conditions_verified'] for w in p['wavefunction_matches']))
        return selected

    def run(self,kind):
        self.prepare()
        scope=self.config.get('scope','pilot')
        if scope not in ('pilot','all'):raise ValueError('scope must be pilot or all')
        results=[]
        if scope=='all':
            for s in self.config['systems']:
                for p in self.points:
                    if p['point_id'] in self.pilot(s):
                        m=self.out/'runs'/p['point_id']/kind/'latest.json'
                        if not m.exists() or not self.resume_valid(json.loads(m.read_text()),p,kind):
                            raise ValueError('Full calculation requires successful pilot validation first')
        selected=[p for p in self.points if scope=='all' or p['point_id'] in self.pilot(p['system'])]
        with concurrent.futures.ThreadPoolExecutor(max_workers=int(self.config.get('parallel_jobs',1))) as pool:
            for r in pool.map(lambda p:self.run_point(p,kind),selected):results.append(r)
        dump(self.out/f'run_{kind}.json',results)
        return results

    def fingerprint(self,p,kind):
        calculation_config={k:v for k,v in self.config.items() if k!='scope'}
        material={'point':p,'kind':kind,'config':calculation_config,'environment':self.env,
                  'adapter_source_sha256':self.provenance['adapter_source_sha256']}
        return hashlib.sha256(json.dumps(material,sort_keys=True).encode()).hexdigest()

    def resume_valid(self,state,p,kind):
        return (state.get('status')=='SUCCESS' and state.get('fingerprint')==self.fingerprint(p,kind)
                and bool(state.get('output_hashes'))
                and all(Path(f).exists() and digest(f)==h for f,h in state['output_hashes'].items()))

    def run_point(self,p,kind):
        folder=self.out/'runs'/p['point_id']/kind;folder.mkdir(parents=True,exist_ok=True)
        latest=folder/'latest.json';fingerprint=self.fingerprint(p,kind)
        if latest.exists():
            old=json.loads(latest.read_text())
            if self.resume_valid(old,p,kind):
                return dict(old,resumed=True)
        state={'point_id':p['point_id'],'system':p['system'],'kind':kind,'fingerprint':fingerprint,'status':'BLOCKED',
               'source':p['source'],'source_sha256':p['source_sha256'],'conditions':p['conditions'],'partition':'TFVC' if kind=='iqa' else None,
               'upstream_commit':self.env['upstream']['APOST3D' if kind=='iqa' else 'easyEDA']['actual_commit'],
               'provenance_file':str(self.out/'provenance.json'),'resumed':False,'started_at':datetime.now(ZoneInfo('Asia/Tokyo')).isoformat()}
        if p['readiness'][kind+'_status']=='BLOCKED':
            state['reasons']=p['readiness'][kind+'_reasons'];dump(latest,state);return state
        # Atomic per-point lease protects scratch files from concurrent CLI instances.
        lock=folder/'running.lock'
        try:fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
        except FileExistsError:return dict(state,status='BLOCKED',reasons=['Point already locked; verify prior process before removing running.lock'])
        os.write(fd,str(os.getpid()).encode());os.close(fd)
        attempt=folder/f'attempt_{time.time_ns()}';attempt.mkdir();state['work_dir']=str(attempt)
        try:
            state['status']='RUNNING';dump(latest,state)
            if kind=='iqa':self.execute_iqa(p,attempt)
            else:self.execute_eda(p,attempt)
            payload=self.collect(p,kind,attempt)
            state.update(status='SUCCESS',validation=payload['validation'])
            dump(attempt/'result.json',payload)
            state['output_hashes']={str(f):digest(f) for f in attempt.iterdir() if f.is_file() and f.suffix in ('.json','.log','.apost','.fchk','.inp','.gjf','.txt')}
        except Exception as e:state.update(status='FAILED',reasons=[f'{type(e).__name__}: {e}'])
        finally:
            state['ended_at']=datetime.now(ZoneInfo('Asia/Tokyo')).isoformat();dump(latest,state);dump(attempt/'status.json',state);lock.unlink()
        return state

    def command(self,args,cwd,stdout,stdin=None,env=None):
        with (cwd/stdout).open('w') as out:
            r=subprocess.run(args,cwd=cwd,input=stdin,text=True,stdout=out,stderr=subprocess.STDOUT,
                             env=env,timeout=self.config.get('timeout_seconds',86400))
        if r.returncode:raise RuntimeError(f'Exit {r.returncode}: {args}; see {cwd/stdout}')

    def execute_iqa(self,p,work):
        matches=[w for w in p['wavefunction_matches'] if w['conditions_verified'] and w['has_mo'] and w['has_density']]
        wf=work/'analysis.fchk';log=None
        if matches:
            wf.write_text(read_text(matches[0]['source']))
            # Log association also needs actual geometry/energy/conditions, never just basename.
            for r in self.inventory:
                if r['extension'] not in ('.log','.out') or not r.get('normal_termination'):continue
                try:
                    q=stationary(read_text(r['source']),r['source'])
                    if q['atomic_numbers']==p['atomic_numbers'] and signature(q['conditions'])==signature(p['conditions']) and rmsd(q['xyz'],p['xyz'])<2e-4 and abs(q['energy_hartree']-p['energy_hartree'])<2e-6:
                        candidate=read_text(r['source'])
                        if re.search(r'^ N-N=.*E-N=.*KE=',candidate,re.M) and candidate.count('Normal termination of Gaussian')==1:
                            log=work/'gaussian.log';log.write_text(candidate);break
                except ValueError:pass
        else:
            inp=gaussian_input(p,self.config);(work/'single_point.gjf').write_text(inp)
            self.command([self.tool('g16')],work,'gaussian.log',inp)
            log=work/'gaussian.log'
            if 'Normal termination of Gaussian' not in log.read_text():raise ValueError('Gaussian SP did not terminate normally')
            sp=stationary(log.read_text(),log)
            if signature(sp['conditions'])!=signature(p['conditions']) or sp['atomic_numbers']!=p['atomic_numbers'] or rmsd(sp['xyz'],p['xyz'])>2e-4 or abs(sp['energy_hartree']-p['energy_hartree'])>self.config.get('energy_tolerance_hartree',2e-6):
                raise ValueError('SP condition/geometry/IRC energy validation failed')
            self.command([self.tool('formchk'),'wavefunction.chk','analysis.fchk'],work,'formchk.log')
        info=fchk_info(wf)
        if info['atomic_numbers']!=p['atomic_numbers'] or rmsd(info['xyz'],p['xyz'])>2e-4:raise ValueError('Converted fchk geometry mismatch')
        if not info['has_energy_terms']:
            if not log or not self.tool('get_energy_g16'):raise ValueError('No matched normal single-job #p log / get_energy_g16 for energy augmentation')
            self.command([self.tool('get_energy_g16'),log.name],work,'energy_terms.txt')
            additions=(work/'energy_terms.txt').read_text()
            if not all(x in additions for x in ('Kinetic Energy','Electron-Nuclei Energy','Electron-Electron Energy')):raise ValueError('Incomplete official energy augmentation')
            wf.write_text(wf.read_text()+'\n'+additions)
        (work/'analysis.inp').write_text(iqa_input(p['conditions'],self.config.get('iqa',{})))
        env=os.environ.copy();env.update(OMP_NUM_THREADS=str(self.config.get('threads',1)),KMP_STACKSIZE='100m')
        self.command([self.tool('apost3d'),'analysis'],work,'analysis.apost',env=env)

    def execute_eda(self,p,work):
        (work/'point.xyz').write_text(xyz_text(p))
        (work/'fragment.txt').write_text(fragment_text(self.config['systems'][p['system']].get('fragments'),p))
        (work/'template.gjf').write_text(eda_template(p))
        source=ROOT/'vendor/easyEDA/source/EDA_run.sh'
        shutil.copy2(source,work/'EDA_run.sh')
        # Invoke official per-point kernel; no hard-coded easyEDA GUI/scan path or shared scratch.
        env=os.environ.copy();bins=work/'bin';bins.mkdir()
        for name in ('g16','formchk','unfchk','Multiwfn','bc','dos2unix','sed'):
            (bins/name).symlink_to(self.tool(name))
        env['PATH']=str(bins)+os.pathsep+env.get('PATH','');env['OMP_NUM_THREADS']=str(self.config.get('threads',1))
        self.command([self.tool('bash'),'EDA_run.sh','point.xyz'],work,'sobEDA.log',env=env)
        nfrag=len(self.config['systems'][p['system']]['fragments'])
        for suffix in ['_final.out','_promol.out']+[f'_fragment{i}.out' for i in range(1,nfrag+1)]:
            file=work/('point'+suffix)
            if not file.exists() or 'Normal termination of Gaussian' not in file.read_text():raise ValueError('Missing normal Gaussian stage '+suffix)
            if not all(label in file.read_text() for label in ('SCF Done','ET=','ENTVJ=')):
                raise ValueError('Missing Gaussian L608 energy components in '+suffix)
        # sobEDA output itself may look complete despite an unchecked subordinate command failure.
        final=stationary((work/'point_final.out').read_text(),work/'point_final.out')
        # Multiwfn may concatenate fragments in fragment order. Verify the explicit
        # permutation and retain original indices in every public adapter artifact.
        original=list(range(1,len(p['xyz'])+1))
        grouped=[i for f in self.config['systems'][p['system']]['fragments'] for i in f['atoms']]
        mapping=None
        for order in (original,grouped):
            z=[p['atomic_numbers'][i-1] for i in order];xyz=[p['xyz'][i-1] for i in order]
            if final['atomic_numbers']==z and rmsd(final['xyz'],xyz)<2e-4:mapping=order;break
        if mapping is None or signature(final['conditions'])!=signature(p['conditions']):
            raise ValueError('sobEDA final state conditions/geometry differ from IRC')
        dump(work/'eda_atom_mapping.json',{'final_output_index_to_original_gaussian_index':mapping})
        if abs(final['energy_hartree']-p['energy_hartree'])>self.config.get('energy_tolerance_hartree',2e-6):raise ValueError('sobEDA final total energy differs from IRC')

    def collect(self,p,kind,work):
        if kind=='iqa':
            atomic,pairs,v=parse_iqa((work/'analysis.apost').read_text(),len(p['xyz']))
            v['irc_energy_difference_hartree']=v['target_hartree']-p['energy_hartree']
            if abs(v['closure_hartree'])*HARTREE_KCAL>self.config.get('closure_tolerance_kcal_mol',0.25) or abs(v['matrix_sum_error_hartree'])>1e-5 or abs(v['irc_energy_difference_hartree'])>self.config.get('energy_tolerance_hartree',2e-6):
                dump(work/'failed_validation.json',v);raise ValueError('IQA closure/IRC energy mismatch')
            return {'atomic':atomic,'pairwise':pairs,'validation':v}
        rows,v,stages=parse_eda((work/'sobEDA.log').read_text())
        expected_stages=['promolecular state','final state']+[f'fragment {i}' for i in range(1,len(self.config['systems'][p['system']]['fragments'])+1)]
        for stage in expected_stages:
            present={r['component'] for r in stages if r['stage']==stage}
            if not {'E_tot','E_T','E_els','E_x','E_c','E_disp'}<=present:
                raise ValueError('Incomplete sobEDA stage components: '+stage)
        final=next((r['energy_hartree'] for r in stages if r['stage']=='final state' and r['component']=='E_tot'),None)
        v['final_minus_irc_hartree']=final-p['energy_hartree'] if final is not None else None
        if any(abs(v[k])>0.05 for k in ('closure_kcal_mol','grouped_closure_kcal_mol','xrep_identity_error','c_identity_error')) or final is None or abs(v['final_minus_irc_hartree'])>self.config.get('energy_tolerance_hartree',2e-6):
            dump(work/'failed_validation.json',v);raise ValueError('sobEDA closure/IRC total mismatch')
        return {'components':rows,'stages':stages,'validation':v}

    def analyze(self):
        self.ensure_points();atomic=[];pairs=[];eda=[];validation=[];stages=[]
        for p in self.points:
            common={k:p[k] for k in ('system','point_id','irc_index','reaction_coordinate','source','source_sha256','conditions','reference_point')}
            common['provenance_file']=str(self.out/'provenance.json')
            validation.append(dict(common,scheme='Gaussian IRC',status='ACCEPTED',check='cclib_evaluated_geometry_and_energy',
                                   cclib_crosscheck=p['cclib_crosscheck'],corrector_minus_scf_hartree=p['corrector_minus_scf_hartree']))
            for kind in ('iqa','eda'):
                latest=self.out/'runs'/p['point_id']/kind/'latest.json'
                state=json.loads(latest.read_text()) if latest.exists() else {'status':'NOT_RUN','reasons':p['readiness'][kind+'_reasons']}
                if state['status']=='SUCCESS' and not self.resume_valid(state,p,kind):
                    state=dict(state,status='STALE',reasons=['Input/config/output fingerprints changed'])
                common_kind=dict(common,status=state['status'],upstream_commit=self.env['upstream']['APOST3D' if kind=='iqa' else 'easyEDA']['actual_commit'],partition='TFVC' if kind=='iqa' else 'sobEDA')
                validation.append(dict(common_kind,scheme=kind,check='decomposition_closure',reasons=state.get('reasons',[]),**state.get('validation',{})))
                if state['status']!='SUCCESS':continue
                result=json.loads((Path(state['work_dir'])/'result.json').read_text())
                if kind=='iqa':
                    for key,target in (('atomic',atomic),('pairwise',pairs)):
                        for r in result[key]:target.append(dict(common_kind,**r,energy_kcal_mol=r['energy_hartree']*HARTREE_KCAL if r['energy_hartree'] is not None else None))
                else:
                    eda.extend(dict(common_kind,**r) for r in result['components']);stages.extend(dict(common_kind,**r) for r in result['stages'])
        def deltas(rows):
            reference={(r['system'],r['point_id'],r.get('atom_i'),r.get('atom_j'),r['component']):r for r in rows}
            for r in rows:
                ref=reference.get((r['system'],r['reference_point'],r.get('atom_i'),r.get('atom_j'),r['component']))
                for unit in ('hartree','kcal_mol'):
                    r['delta_'+unit]=r['energy_'+unit]-ref['energy_'+unit] if ref and r['energy_'+unit] is not None and ref['energy_'+unit] is not None else None
        basic=['system','point_id','irc_index','reaction_coordinate','source','source_sha256','conditions','reference_point','status','upstream_commit','partition','provenance_file']
        for name,rows,extra in [('iqa_atomic',atomic,['atom_i','atom_j']),('iqa_pairwise',pairs,['atom_i','atom_j']),('eda_components',eda,[])]:
            deltas(rows);csv_write(self.out/(name+'.csv'),rows,basic+extra+['component','energy_hartree','energy_kcal_mol','delta_hartree','delta_kcal_mol','source_precision_kcal_mol','hartree_conversion']);dump(self.out/(name+'.json'),rows)
        csv_write(self.out/'validation_report.csv',validation);dump(self.out/'validation_report.json',validation)
        csv_write(self.out/'eda_stages.csv',stages,basic+['stage','component','energy_hartree'])
        self.plots(pairs,eda)
        summary={'irc_points':len(self.points),'systems':{s:sum(p['system']==s for p in self.points) for s in self.config['systems']},
                 'iqa_success_points':len({r['point_id'] for r in atomic}),'eda_success_points':len({r['point_id'] for r in eda}),
                 'prepared_sp_inputs':len(list((self.out/'prepared').glob('*/single_point.gjf'))),
                 'executed_gaussian_sp_jobs':len(list((self.out/'runs').glob('*/iqa/attempt_*/gaussian.log'))),
                 'failed_points':[{'point_id':r['point_id'],'scheme':r['scheme'],'reasons':r.get('reasons')} for r in validation if r['status']=='FAILED']}
        dump(self.out/'execution_summary.json',summary)
        return summary

    def plots(self,pairs,eda):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        systems=list(self.config['systems']);fig,axes=plt.subplots(1,len(systems),figsize=(4.3*len(systems),3.7),squeeze=False)
        for ax,s in zip(axes[0],systems):
            ps=[p for p in self.points if p['system']==s]
            ax.plot([p['reaction_coordinate'] for p in ps],[p['relative_energy_kcal_mol'] for p in ps],'.-',ms=3)
            ax.axvline(0,color='gray',lw=.7);ax.set(title=s,xlabel=r'IRC / $\sqrt{amu}$ bohr',ylabel='SCF energy from reference / kcal mol$^{-1}$')
            ax.grid(alpha=.2)
        fig.tight_layout();fig.savefig(self.out/'01_irc_profile.png',dpi=170);plt.close(fig)
        for name,rows,label in [('02_iqa_interactions',pairs,'IQA / TFVC'),('03_eda_components',eda,'sobEDA')]:
            fig,axes=plt.subplots(1,len(systems),figsize=(4.3*len(systems),3.7),squeeze=False)
            for ax,s in zip(axes[0],systems):
                rs=[r for r in rows if r['system']==s and r['delta_kcal_mol'] is not None]
                if not rs:
                    ax.text(.5,.5,'No validated numerical results\nSee readiness.csv',ha='center',va='center',transform=ax.transAxes);ax.set_xticks([]);ax.set_yticks([])
                else:
                    groups={}
                    for r in rs:
                        key=f'{r.get("atom_i","")}-{r.get("atom_j","")} {r["component"]}' if 'atom_i' in r else r['component']
                        groups.setdefault(key,[]).append(r)
                    if name.startswith('02'):
                        requested=self.config.get('plot_pairs',[])
                        if requested:groups={k:v for k,v in groups.items() if [v[0]['atom_i'],v[0]['atom_j']] in requested}
                        groups=dict(sorted(groups.items(),key=lambda kv:np.ptp([r['delta_kcal_mol'] for r in kv[1]]),reverse=True)[:self.config.get('plot_top_pairs',6)])
                    for k,vs in groups.items():
                        vs.sort(key=lambda r:r['reaction_coordinate']);ax.plot([r['reaction_coordinate'] for r in vs],[r['delta_kcal_mol'] for r in vs],'.-',label=k)
                    ax.legend(fontsize=6);ax.set(xlabel=r'IRC / $\sqrt{amu}$ bohr',ylabel='Change from reference / kcal mol$^{-1}$')
                ax.set_title(s+' — '+label)
            fig.tight_layout();fig.savefig(self.out/(name+'.png'),dpi=170);plt.close(fig)
