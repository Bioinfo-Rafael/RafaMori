"""Read-only Gaussian parsing. cclib supplies the electronic/geometry cross-check;
completion markers supplement its missing IRC point/path/coordinate identifiers.
Never infer a TS geometry or accept an unfinished corrector point.
"""
from __future__ import annotations
import gzip
import hashlib
import io
import logging
import re
from pathlib import Path
import cclib
import numpy as np
from cclib.parser.utils import PeriodicTable, convertor

FLOAT = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'
SCF = re.compile(r'SCF Done:.*?=\s*(' + FLOAT + ')')
POINT = re.compile(r'Point Number:\s*(\d+)\s+Path Number:\s*(\d+)')
NET = re.compile(r'NET REACTION COORDINATE[^=]*=\s*(' + FLOAT + ')')
HARTREE_KCAL = 627.5094740631
BOHR_A = 0.529177210903
TABLE = PeriodicTable()


def number(s):
    return float(s.replace('D', 'E').replace('d', 'e'))


def read_text(path):
    path = Path(path)
    if path.suffix.lower() == '.gz':
        with gzip.open(path, 'rt', errors='replace') as f:
            return f.read()
    return path.read_text(errors='replace')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def extension(path):
    value=str(path)
    return Path(value[:-3] if value.lower().endswith('.gz') else value).suffix.lower()


def routes(text, input_file=False):
    lines = text.splitlines()
    result = []
    for i, line in enumerate(lines):
        if re.match(r'^\s*#[pPnNtT ]', line):
            if re.match(r'\s*#\s+OF\s',line):
                continue
            parts = [line.lstrip()]
            for s in lines[i + 1:i + 30]:
                if not s.strip() or re.match(r'\s*-{4,}', s):
                    break
                parts.append(s.strip() if input_file else s[1:])
            result.append((' ' if input_file else '').join(parts).strip())
    return result


def conditions(text, route=None):
    rs = routes(text)
    route = route or (rs[0] if rs else '')
    compact = re.sub(r'\s+', '', route).lower()
    mb = re.search(r'([a-zA-Z0-9+_-]+)/([^\s]+)', route)
    charge = re.search(r'Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)', text)
    solvent = re.search(r'scrf=(\([^)]*\)|[^\s]+)', route, re.I)
    basis = re.findall(r'Standard basis:\s*(\S+)', text)
    method = mb[1].upper() if mb else None
    if method is None:
        scf_methods = re.findall(r'SCF Done:\s*E\(([^)]+)\)',text)
        method = scf_methods[-1].upper() if scf_methods else None
    if method and method[0] in 'RU' and method[1:] in ('HF','B3LYP','M062X','PBE1PBE','PBEPBE','BP86'):
        method = method[1:]
    return dict(route=route, method=method, basis=basis[-1] if basis else (mb[2] if mb else None),
                charge=int(charge[1]) if charge else None, multiplicity=int(charge[2]) if charge else None,
                solvent=solvent[0].lower() if solvent else ('gas' if route else None),
                is_irc='irc' in compact, routes=rs)


def signature(c):
    return tuple(re.sub(r'\s+', '', str(c.get(k))).upper() for k in ('method','basis','charge','multiplicity','solvent'))


def orientation(lines, start):
    i, separators = start + 1, 0
    while i < len(lines) and i < start + 12:
        if re.match(r'\s*-{5,}', lines[i]):
            separators += 1
            if separators == 2:
                i += 1
                break
        i += 1
    z, xyz, ids = [], [], []
    while i < len(lines):
        m = re.match(r'\s*(\d+)\s+(\d+)\s+-?\d+\s+('+FLOAT+r')\s+('+FLOAT+r')\s+('+FLOAT+r')\s*$', lines[i])
        if not m:
            break
        ids.append(int(m[1])); z.append(int(m[2])); xyz.append([number(m[k]) for k in (3,4,5)])
        i += 1
    if ids != list(range(1, len(ids)+1)):
        raise ValueError('Gaussian atom ordering is nonsequential')
    return z, xyz


def cclib_read(text):
    """Gaussian parser directly: filename sniffing is unreliable for merged logs."""
    try:
        data = cclib.parser.Gaussian(io.StringIO(text), loglevel=logging.CRITICAL).parse()
        return data, None
    except Exception as exc:
        return None, f'{type(exc).__name__}: {exc}'


def scan(text):
    lines = text.splitlines()
    geoms, energies = [], []
    for i, line in enumerate(lines):
        if re.search(r'(Input|Standard|Z-Matrix) orientation:', line):
            z, xyz = orientation(lines, i)
            if z:
                geoms.append((i, z, xyz))
        m = SCF.search(line)
        if m:
            energies.append((i, number(m[1])))
    return lines, geoms, energies


def stationary(text, source):
    lines, geoms, energies = scan(text)
    if not energies:
        raise ValueError(f'No SCF energy: {source}')
    e = energies[-1]
    gs = [g for g in geoms if g[0] < e[0]]
    if not gs:
        raise ValueError(f'No evaluated geometry: {source}')
    g = gs[-1]
    return dict(source=str(source), energy_hartree=e[1], energy_line=e[0]+1,
                geometry_line=g[0]+1, atomic_numbers=g[1], xyz=g[2], conditions=conditions(text))


def rmsd(a, b):
    a, b = np.asarray(a), np.asarray(b)
    ac, bc = a-a.mean(0), b-b.mean(0)
    u, _, vt = np.linalg.svd(ac.T @ bc)
    d = np.eye(3); d[-1,-1] = np.sign(np.linalg.det(u @ vt))
    return float(np.sqrt(np.mean(np.sum((ac @ u @ d @ vt - bc)**2, axis=1))))


def parse_irc(path, system, sign=1, ts_path=None):
    text = read_text(path)
    if 'DERIVED ANALYSIS LOG' in text[:1000]:
        raise ValueError('Analysis-only oriented log is not an original Gaussian calculation')
    if sign not in (-1,1):
        raise ValueError('coordinate_sign must be +1 or -1')
    c = conditions(text)
    if not c['is_irc']:
        raise ValueError(f'No IRC route in {path}')
    # Multiple jobs may differ in IRC settings, but not electronic conditions.
    for r in c['routes']:
        if 'irc' in r.lower() and signature(conditions(text,r)) != signature(c):
            raise ValueError('Mixed electronic conditions in merged IRC log')
    lines, geoms, energies = scan(text)
    data, error = cclib_read(text)
    audit = dict(source=str(path), cclib_version=cclib.__version__, cclib_error=error,
                 cclib_geometry_count=len(getattr(data,'atomcoords',[])),
                 rejected=[], duplicate_count=0, warnings=[])
    points = {}
    lower = 0
    announced_direction = None
    for i, line in enumerate(lines):
        announcement = re.search(r'Point Number\s+\d+ in (FORWARD|REVERSE)',line)
        if announcement:
            announced_direction = announcement[1]
        m = POINT.search(line)
        if not m:
            continue
        n, original_branch = int(m[1]), int(m[2])
        branch = (1 if announced_direction=='FORWARD' else 2) if announced_direction else original_branch
        if n == 0:
            lower = i
            continue
        following = '\n'.join(lines[i+1:i+6])
        coord = NET.search(following)
        if not coord:
            continue
        es = [e for e in energies if lower < e[0] < i]
        # The last SCF and the orientation preceding THAT SCF are one evaluated pair.
        if not es:
            audit['rejected'].append(dict(point=n,path=branch,reason='missing SCF'))
            lower = i; continue
        e = es[-1]
        gs = [g for g in geoms if lower < g[0] < e[0]]
        if not gs:
            audit['rejected'].append(dict(point=n,path=branch,reason='missing evaluated geometry'))
            lower = i; continue
        g = gs[-1]
        if branch not in (1,2):
            raise ValueError('Unrecognized Gaussian path number')
        checked = False
        if data is not None and hasattr(data,'atomcoords') and hasattr(data,'scfenergies'):
            if list(data.atomnos) != g[1]:
                raise ValueError('cclib / orientation atom order disagreement')
            ce = convertor(np.asarray(data.scfenergies),'eV','hartree')
            candidates = np.flatnonzero(np.abs(ce-e[1]) < 1e-7)
            # Some Gaussian jobs print input plus standard frames. Compare rotation-invariantly.
            checked = any(j < len(data.atomcoords) and rmsd(g[2],data.atomcoords[j]) < 2e-5 for j in candidates)
            if not checked:
                raise ValueError(f'cclib geometry/SCF cross-check failed: {path}:{i+1}')
        p = dict(system=system, point_index=n, path_number=branch,
                 gaussian_path_number=original_branch,
                 direction='forward' if branch==1 else 'reverse',
                 reaction_coordinate=abs(number(coord[1]))*(1 if branch==1 else -1)*sign,
                 coordinate_unit='sqrt(amu)*bohr', coordinate_sign=sign,
                 energy_hartree=e[1], atomic_numbers=g[1], xyz=g[2], source=str(Path(path).resolve()),
                 source_sha256=digest(path), source_line=i+1, energy_line=e[0]+1, geometry_line=g[0]+1,
                 geometry_kind='last_evaluated_SCF_geometry_at_accepted_IRC_step',
                 conditions=c, cclib_crosscheck=checked)
        # Corrected endpoints can be interpolants, rather than separately evaluated geometries.
        corrected = re.findall(r'Corrected End Point Energy\s*=\s*('+FLOAT+')','\n'.join(lines[max(lower,i-60):i]))
        p['corrector_energy_hartree'] = number(corrected[-1]) if corrected else None
        p['corrector_minus_scf_hartree'] = p['corrector_energy_hartree']-e[1] if corrected else None
        key = (branch,n)
        if key in points:
            old = points[key]
            if abs(old['energy_hartree']-e[1])>1e-7 or rmsd(old['xyz'],g[2])>2e-5 or abs(old['reaction_coordinate']-p['reaction_coordinate'])>1e-6:
                raise ValueError(f'Conflicting duplicate point {key}')
            audit['duplicate_count'] += 1
        else:
            points[key] = p
        lower = i
    if not points:
        raise ValueError('No completed IRC points')
    attempted = {(1 if m[2]=='FORWARD' else 2,int(m[1])) for m in re.finditer(r'Point Number\s+(\d+) in (FORWARD|REVERSE)',text)}
    for key in sorted(attempted-set(points)):
        audit['rejected'].append(dict(point=key[1],path=key[0],reason='no completion + NET coordinate'))
    values = list(points.values())
    if any(p['atomic_numbers'] != values[0]['atomic_numbers'] for p in values):
        raise ValueError('Atom identities change along IRC')
    if ts_path:
        ts = stationary(read_text(ts_path), Path(ts_path).resolve())
        checkpoints = [number(x) for x in re.findall(r'Energy From Chk\s*=\s*('+FLOAT+')',text)]
        neighbors = sorted(values,key=lambda p:abs(p['reaction_coordinate']))[:2]
        if (signature(ts['conditions']) != signature(c) or ts['atomic_numbers'] != values[0]['atomic_numbers']
            or not checkpoints or max(abs(e-ts['energy_hartree']) for e in checkpoints)>1e-6
            or max(rmsd(ts['xyz'],p['xyz']) for p in neighbors) > 0.1):
            raise ValueError('Explicit TS failed condition/order/checkpoint energy/neighbor geometry checks')
        values.append(dict(values[0], **ts, point_index=0,path_number=0,direction='TS',reaction_coordinate=0.,
                           source_line=ts['geometry_line'],source_sha256=digest(ts_path),
                           geometry_kind='explicit_TS_log_checked_against_IRC_checkpoint_energy',
                           cclib_crosscheck=None,corrector_energy_hartree=None,corrector_minus_scf_hartree=None))
    else:
        # A real orientation BEFORE the initial point-0 marker is acceptable;
        # the first displaced point is not a TS, even if its energy is nearby.
        first_zero = next((i for i,line in enumerate(lines) if (m:=POINT.search(line)) and int(m[1])==0),None)
        initial_geometries = [g for g in geoms if first_zero is not None and g[0]<first_zero]
        checkpoints = re.findall(r'Energy From Chk\s*=\s*('+FLOAT+')','\n'.join(lines[:first_zero])) if first_zero is not None else []
        if initial_geometries and checkpoints:
            g=initial_geometries[-1]
            if g[1]!=values[0]['atomic_numbers']:raise ValueError('Initial TS atom order mismatch')
            values.append(dict(values[0],point_index=0,path_number=0,gaussian_path_number=0,direction='TS',reaction_coordinate=0.,
                               energy_hartree=number(checkpoints[-1]),xyz=g[2],geometry_line=g[0]+1,source_line=first_zero+1,
                               energy_line=None,geometry_kind='explicit_initial_IRC_orientation_and_checkpoint_energy',
                               cclib_crosscheck=None,corrector_energy_hartree=None,corrector_minus_scf_hartree=None))
        else:
            audit['warnings'].append('No explicit TS geometry; no TS fabricated')
    values.sort(key=lambda p:p['reaction_coordinate'])
    for index,p in enumerate(values):
        p['irc_index'] = index
        p['point_id'] = f'{system}_p{p["path_number"]}_{p["point_index"]:04d}'
    audit['accepted_points'] = len(values)
    audit['path_counts'] = {str(k):sum(p['path_number']==k for p in values) for k in (0,1,2)}
    for branch in (1,2):
        ns = sorted(p['point_index'] for p in values if p['path_number']==branch)
        audit.setdefault('missing_indices',{})[str(branch)] = sorted(set(range(1,max(ns,default=0)+1))-set(ns))
    return values,audit


def fchk_info(path):
    """Read metadata/geometry/basis signatures, skipping huge orbital matrices."""
    text = read_text(path); lines = text.splitlines()
    scalars, arrays = {}, {}
    wanted = {'Atomic numbers','Current cartesian coordinates','Shell types','Number of primitives per shell',
              'Shell to atom map','Primitive exponents','Contraction coefficients','P(S=P) Contraction coefficients'}
    i=2; route=''
    while i<len(lines):
        line=lines[i]; i+=1
        if len(line)<45: continue
        label=line[:43].strip(); rest=line[43:].strip()
        m=re.match(r'([IRCL])\s+(?:N=\s*(\d+)|(.+))',rest)
        if not m: continue
        if m[2]:
            n=int(m[2])
            if label=='Route':
                chunks=[]
                while i<len(lines) and not re.match(r'^.{43}[IRCL]\s',lines[i]):
                    chunks.append(lines[i]);i+=1
                route=''.join(chunks).strip()
            elif label in wanted:
                vals=[]
                while len(vals)<n and i<len(lines):
                    vals.extend(lines[i].split());i+=1
                arrays[label]=[number(v) for v in vals[:n]]
        else:
            try: scalars[label]=number(m[3])
            except ValueError: pass
    c=conditions('',route=route)
    header=lines[1].split()
    if not c['method'] and len(header)>=2:
        c['method']=header[1].removeprefix('R').removeprefix('U')
    c['charge']=int(scalars['Charge']);c['multiplicity']=int(scalars['Multiplicity'])
    coords=np.asarray(arrays['Current cartesian coordinates']).reshape(-1,3)*BOHR_A
    return dict(conditions=c, atomic_numbers=[int(z) for z in arrays['Atomic numbers']],xyz=coords.tolist(),
                atom_count=int(scalars['Number of atoms']),energy_hartree=scalars.get('Total Energy'),
                nbasis=int(scalars.get('Number of basis functions',0)),
                basis_signature=hashlib.sha256(repr({k:v for k,v in arrays.items() if k not in ('Atomic numbers','Current cartesian coordinates')}).encode()).hexdigest(),
                has_mo='Alpha MO coefficients' in text,has_density='Total SCF Density' in text,
                has_energy_terms=all(k in scalars for k in ('Kinetic Energy','Electron-Nuclei Energy','Electron-Electron Energy')))


def xyz_text(p):
    return f'{len(p["atomic_numbers"])}\n{p["point_id"]}; R={p["reaction_coordinate"]} {p["coordinate_unit"]}; E={p["energy_hartree"]} Eh\n'+''.join(
        f'{TABLE.element[z]:2s} {x:.10f} {y:.10f} {zz:.10f}\n' for z,(x,y,zz) in zip(p['atomic_numbers'],p['xyz']))
