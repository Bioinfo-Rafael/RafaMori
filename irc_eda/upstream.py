"""Thin interfaces and output readers for the pinned, unmodified upstream codes."""
from __future__ import annotations
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
from pathlib import Path
import numpy as np
from .parsing import FLOAT, HARTREE_KCAL, number, signature

ROOT=Path(__file__).resolve().parent
PINS={'APOST3D':('https://github.com/mgimferrer/APOST3D.git','3ba2d95b58831d7a6adaaedc40cc4ed950eb16c5'),
      'easyEDA':('https://github.com/kangmg/easyEDA.git','08d7a7dfe4775e35c2fc55e321d1262c3afd33e2')}


def git_revision(path):
    r=subprocess.run(['git','-C',str(path),'rev-parse','HEAD'],capture_output=True,text=True)
    return r.stdout.strip() if r.returncode==0 else None


def upstream_clean(path):
    if not (Path(path)/'.git').exists():return False
    return subprocess.run(['git','-C',str(path),'diff','--quiet','HEAD'],capture_output=True).returncode==0


def environment(config):
    names=['g16','formchk','unfchk','Multiwfn','ifort','ifx','icx','gfortran','make','bc','dos2unix','sed','bash','apost3d','get_energy_g16']
    overrides=config.get('executables',{})
    found={}
    for name in names:
        candidate=overrides.get(name) or shutil.which(name)
        if not candidate and name in ('apost3d','get_energy_g16'):
            candidate=ROOT/'vendor/APOST3D'/('utils/get_energy_g16' if name=='get_energy_g16' else name)
        if candidate and not (Path(candidate).is_file() and os.access(candidate,os.X_OK)):
            candidate=shutil.which(str(candidate))
        info={'path':str(Path(candidate).resolve()) if candidate else None,'version':None,'shared_libraries':None}
        if candidate:
            if name in ('ifort','ifx','icx','gfortran','make','bc','sed','bash'):
                try:
                    r=subprocess.run([str(candidate),'--version'],capture_output=True,text=True,timeout=10)
                    info['version']=(r.stdout+r.stderr)[:1500].strip()
                except (OSError,subprocess.TimeoutExpired) as exc: info['version']=str(exc)
            else:
                info['version']=config.get('software_versions',{}).get(name,'unknown; executable present, version not verified')
            if name in ('g16','Multiwfn','apost3d'):
                inspector=shutil.which('otool' if platform.system()=='Darwin' else 'ldd')
                if inspector:
                    args=[inspector]+(['-L'] if platform.system()=='Darwin' else [])+[str(candidate)]
                    r=subprocess.run(args,capture_output=True,text=True)
                    info['shared_libraries']=(r.stdout+r.stderr)[:5000]
            try:
                with open(candidate,'rb') as f: info['binary_sha256']=hashlib.file_digest(f,'sha256').hexdigest()
            except OSError: pass
        found[name]=info
    system=platform.system()
    return {'platform':system,'environment':'WSL' if system=='Linux' and 'microsoft' in platform.release().lower() else system,
            'architecture':platform.machine(),'python':platform.python_version(),'tools':found,
            'upstream':{n:{'url':url,'expected_commit':rev,'actual_commit':git_revision(ROOT/'vendor'/n),
                           'tracked_source_unmodified':upstream_clean(ROOT/'vendor'/n)} for n,(url,rev) in PINS.items()},
            'adapter_git_commit':git_revision(ROOT),
            'search_scope':'PATH, explicit executables configuration, local vendor builds; no remote compute host configured',
            'known_install_locations':{p:Path(p).exists() for p in ['/opt/intel/oneapi','/usr/local/g16','/opt/g16','/Applications/g16']},
            'libxc_archive_present':(ROOT/'vendor/APOST3D/libxc-4.2.3.tar.gz').is_file()}


def iqa_method(c, library=None):
    method=(c.get('method') or '').upper()
    if method in ('M062X','M06-2X','M06','M06L','M06-L','M052X','M05','TPSS','TPSSTPSS','MN15','SCAN'):
        raise ValueError(f'{method}: meta-GGA evaluation is NOT IMPLEMENTED in pinned APOST3D enpart_dft.f:1113-1123')
    if c.get('solvent')!='gas':
        raise ValueError('SMD/PCM energy and reaction-field treatment unvalidated for ENPART; gas-phase substitution prohibited')
    if method in ('HF','B3LYP','BP86'):
        return method
    if library and library.get('method','').upper()==method:
        if library.get('family') not in ('LDA','GGA','HYB_GGA') or not library.get('verified_equivalence_source'):
            raise ValueError('LIBRARY requires verified Gaussian/Libxc equivalence and supported family')
        ids={k:int(v) for k,v in library.items() if k in ('EX_FUNCTIONAL','EC_FUNCTIONAL','EXC_FUNCTIONAL')}
        if not ('EXC_FUNCTIONAL' in ids or {'EX_FUNCTIONAL','EC_FUNCTIONAL'}<=ids.keys()):
            raise ValueError('LIBRARY requires combined XC ID or both exchange/correlation IDs')
        return 'LIBRARY\n'+'\n'.join(f'{k}={v}' for k,v in ids.items())
    raise ValueError(f'No verified APOST3D functional mapping for {method}')


def iqa_input(c,config):
    method=iqa_method(c,config.get('libxc'))
    if config.get('partition','TFVC')!='TFVC':
        raise ValueError('Only verified TFVC interface enabled')
    return '# METHOD\nTFVC\nENPART\nFULLPRECISION\n#\n# ENPART\n'+method+'\nTHREBOD=-1\nTWOELTOLER=0.25\n#\n'


def split_route(route):
    return re.findall(r'[^\s()]+(?:\((?:[^()]|\([^()]*\))*\))?[^\s]*',route)


def sp_route(c):
    if any(c.get(k) is None for k in ('method','basis','charge','multiplicity','solvent')):
        raise ValueError('Incomplete electronic conditions')
    if c['basis'].lower() in ('gen','genecp','chkbasis'):
        raise ValueError('Explicit basis/ECP cards or verified named basis required')
    tokens=split_route(c['route'])
    keep=[]; replaced=False
    for token in tokens:
        low=token.lower()
        if low.startswith('#'):continue
        key=re.split('[=(]',low)[0]
        if key in ('irc','opt','freq','geom','guess','restart','sp','pop','density','prop'): continue
        if low in (c['method'].lower(),c['basis'].lower()): continue
        if '/' in token and not replaced:
            keep.append(f'{c["method"]}/{c["basis"]}');replaced=True
        else: keep.append(token)
    if not replaced:
        if c['method'].lower() not in c['route'].lower() or c['basis'].lower() not in c['route'].lower():
            raise ValueError('Cannot preserve original route safely')
        keep.insert(0,f'{c["method"]}/{c["basis"]}')
    # ECP/ONIOM/external field/general basis need additional data, not a guessed SP.
    if re.search(r'\b(oniom|pseudo|field|external)\b|scrf=[^\s]*read', ' '.join(keep),re.I):
        raise ValueError('Additional route-specific cards unsupported')
    return '#p '+' '.join(keep)+' sp nosymm'


def gaussian_input(p,config):
    from .parsing import TABLE
    route=sp_route(p['conditions'])
    coords='\n'.join(f'{TABLE.element[z]} {x:.10f} {y:.10f} {zz:.10f}' for z,(x,y,zz) in zip(p['atomic_numbers'],p['xyz']))
    return (f'%chk=wavefunction.chk\n%nprocshared={int(config.get("threads",1))}\n%mem={config.get("memory","4GB")}\n'
            f'{route}\n\nFixed IRC geometry {p["point_id"]}\n\n{p["conditions"]["charge"]} {p["conditions"]["multiplicity"]}\n{coords}\n\n')


def fragment_text(fragments,p):
    if not fragments or len(fragments)<2:
        raise ValueError('User-defined fragment atoms/charges/spin multiplicities required')
    indices=[];charges=0; spin_projections=[];lines=[str(len(fragments))]
    for f in fragments:
        atoms=f.get('atoms');charge=f.get('charge');multi=f.get('multiplicity')
        if not atoms or type(charge) is not int or type(multi) is not int or not multi:
            raise ValueError('Fragment atoms, integer charge and nonzero multiplicity required')
        if any(type(i) is not int or i<1 or i>len(p['atomic_numbers']) for i in atoms):
            raise ValueError('Fragment indices must be original Gaussian 1-based indices')
        electrons=sum(p['atomic_numbers'][i-1] for i in atoms)-charge
        if electrons<0 or abs(multi)-1>electrons or (electrons-(abs(multi)-1))%2:
            raise ValueError('Fragment electron count / spin multiplicity incompatible')
        indices+=atoms;charges+=charge;spin_projections.append((abs(multi)-1)*(1 if multi>0 else -1))
        lines.extend([f'{charge} {multi}',','.join(map(str,atoms))])
    if sorted(indices)!=list(range(1,len(p['atomic_numbers'])+1)):
        raise ValueError('Fragments must cover all atoms once, without overlap')
    if charges!=p['conditions']['charge']:
        raise ValueError('Fragment charges do not sum to total charge')
    if abs(sum(spin_projections))!=p['conditions']['multiplicity']-1:
        raise ValueError('Fragment spin projections do not reproduce total multiplicity; negative multiplicity flips spin in sobEDA')
    return '\n'.join(lines)+'\n'


def eda_template(p):
    c=p['conditions']
    if c.get('solvent')!='gas':
        raise ValueError('sobEDA solvent decomposition not validated; SMD/PCM cannot be silently included')
    route=sp_route(c).replace(' sp nosymm','')
    # Source script relies on lower-case nosymm for substitutions and L608 energy output.
    if 'extralinks' in route.lower(): raise ValueError('Existing ExtraLinks needs manual review')
    return f'{route} ExtraLinks=L608 nosymm\n\nFixed IRC geometry sobEDA\n\n{c["charge"]} {c["multiplicity"]}\n[geometry]\n\n-5 5\n\n'


def matrix_block(lines,start,n):
    cells={}; columns=[];started=False
    for line in lines[start+1:]:
        if re.match(r'^\s*(?:\d+\s+[A-Z][a-z]?\s*)+$',line):
            columns=[int(x) for x in re.findall(r'(\d+)\s+[A-Z][a-z]?',line)];continue
        m=re.match(r'^\s*(\d+)\s+[A-Z][a-z]?\s+(.+)$',line)
        if m and columns:
            vals=re.findall(FLOAT,m[2])
            if len(vals)!=len(columns): raise ValueError('Truncated/overflowing APOST matrix')
            row=int(m[1])
            for col,v in zip(columns,vals): cells[(row,col)]=number(v)
            started=True
        elif started and line.strip() and not re.fullmatch(r'[-\s]+',line):break
    if len(cells)!=n*n or any((i,j) not in cells for i in range(1,n+1) for j in range(1,n+1)):
        raise ValueError(f'APOST matrix incomplete: {len(cells)}/{n*n}')
    return cells


def parse_iqa(text,n):
    if 'Normal Termination of APOST-3D' not in text:
        raise ValueError('APOST3D did not terminate normally')
    patterns={
      'total':r'^FUZZY ATOMS (?:KS-DFT|Hartree-Fock) ENERGY COMPONENTS$',
      'electron_nuclear':r'^ELECTRON-NUCLEAR ATTRACTION$',
      'nuclear_nuclear':r'^NUCLEAR-NUCLEAR REPULSION$',
      'coulomb':r'^(?:INTERPOLATED )?COULOMB \(ELECTRON-ELECTRON\) ENERGY TERMS$',
      'xc':r'^(?:(?:INTERPOLATED )?HYBRID KS-DFT XC TERMS|FINAL PURE KS-DFT EXCHANGE-CORRELATION ENERGY COMPONENTS|(?:INTERPOLATED )?HARTREE-FOCK-TYPE EXCHANGE ENERGY TERMS)$'}
    matrices={}; lines=text.splitlines()
    # Last XC block after zero-error interpolation supersedes preceding iterations.
    for i,line in enumerate(lines):
        for name,pattern in patterns.items():
            if re.match(pattern,line.strip()): matrices[name]=matrix_block(lines,i,n)
    if 'total' not in matrices: raise ValueError('No final ENPART energy matrix')
    atom=[];pair=[]
    for (i,j),v in sorted(matrices['total'].items()):
        if j<i:continue  # upstream sums diagonal plus each unordered pair ONCE
        row={'atom_i':i,'atom_j':j,'component':'intra_atomic' if i==j else 'interatomic','energy_hartree':v}
        (atom if i==j else pair).append(row)
        for component in ('xc','classical_electrostatic'):
            energy=None
            if component=='xc' and 'xc' in matrices:energy=matrices['xc'][(i,j)]
            if component=='classical_electrostatic' and {'coulomb','electron_nuclear','nuclear_nuclear'}<=matrices.keys():
                energy=matrices['coulomb'][(i,j)]+matrices['electron_nuclear'][(i,j)]
                if i!=j:energy+=matrices['nuclear_nuclear'][(i,j)]
            (atom if i==j else pair).append(dict(row,component=component,energy_hartree=energy))
    totals=re.findall(r'Total (?:integrated|KS-DFT) energy\s*:\s*('+FLOAT+')',text)
    targets=re.findall(r'Total energy in Fchk file:\s*('+FLOAT+')',text)
    if not totals or not targets: raise ValueError('Missing IQA molecular totals')
    reconstructed=sum(v for (i,j),v in matrices['total'].items() if i<=j)
    validation=dict(reconstructed_hartree=reconstructed,reported_hartree=number(totals[-1]),target_hartree=number(targets[-1]),
                    closure_hartree=reconstructed-number(targets[-1]),matrix_sum_error_hartree=reconstructed-number(totals[-1]))
    return atom,pair,validation


EDA_LABELS={
 'Total interaction energy':'E_tot','Electrostatic (E_els)':'E_els','Exchange (E_x)':'E_x',
 'Pauli repulsion (E_rep)':'E_rep','Exchange-repulsion (E_xrep = E_x + E_rep)':'E_xrep',
 'Orbital (E_orb)':'E_orb','DFT correlation (E_DFTc)':'E_DFTc','Dispersion correction (E_dc)':'E_dc',
 'Coulomb correlation (E_c = E_DFTc + E_dc)':'E_c'}


def parse_eda(text):
    components={}
    for label,component in EDA_LABELS.items():
        values=re.findall(re.escape(label)+r':\s*('+FLOAT+r')\s+kcal/mol',text)
        if not values: raise ValueError(f'sobEDA component missing: {component}')
        components[component]=number(values[-1])
    # E_xrep and E_c are grouped terms. Never add them again to their constituents.
    closure=sum(components[k] for k in ('E_els','E_x','E_rep','E_orb','E_DFTc','E_dc'))-components['E_tot']
    grouped=sum(components[k] for k in ('E_els','E_xrep','E_orb','E_c'))-components['E_tot']
    # Keep all optional upstream sobEDAw printed fields, with distinct names.
    if 'Variant of sobEDA for weak interaction (sobEDAw):' in text:
        tail=text.split('Variant of sobEDA for weak interaction (sobEDAw):')[-1]
        for m in re.finditer(r'^(.+?):\s*('+FLOAT+r')\s*kcal/mol',tail,re.M):
            components['sobEDAw:'+m[1].strip()]=number(m[2])
    rows=[dict(component=k,energy_kcal_mol=v,energy_hartree=v/627.51,source_precision_kcal_mol=0.01,
               hartree_conversion=627.51) for k,v in components.items()]
    stages=[];stage=None
    for line in text.splitlines():
        m=re.match(r'Energy components of (.+):',line)
        if m:stage=m[1]
        m=re.match(r'(E_\w+) = ('+FLOAT+r') Hartree',line)
        if m and stage:stages.append(dict(stage=stage,component=m[1],energy_hartree=number(m[2])))
    return rows,dict(closure_kcal_mol=closure,grouped_closure_kcal_mol=grouped,
                    xrep_identity_error=components['E_xrep']-components['E_x']-components['E_rep'],
                    c_identity_error=components['E_c']-components['E_DFTc']-components['E_dc'],
                    rounding_tolerance_kcal_mol=0.05),stages
