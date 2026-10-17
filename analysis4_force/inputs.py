"""Auditable Gaussian input extraction and analysis-only oriented LOG export."""
from __future__ import annotations
import collections, gzip, hashlib, json, re
from pathlib import Path
import numpy as np
from parse_gaussian_irc import FLOAT, _orientation, parse_log
from io_parsers import _route, SCF, _number


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def restore_inputs(root):
    """Expand byte-identical compressed originals and verify every input hash."""
    for record in json.loads((root/'manifest.json').read_text()):
        path = root / Path(record['path']).relative_to('data/completed_20261017')
        if not path.exists():
            with gzip.open(str(path)+'.gz','rb') as src:
                path.write_bytes(src.read())
        if sha256(path) != record['sha256']:
            raise ValueError(f'Input checksum mismatch: {path}')


def parse_stationary(path):
    text=path.read_text(errors='replace'); lines=text.splitlines()
    scfs=[(i,_number(m[1])) for i,line in enumerate(lines) if (m:=SCF.search(line))]
    if not scfs: raise ValueError(f'No SCF energy: {path}')
    energy_line,energy=scfs[-1]
    orientations=[i for i,line in enumerate(lines[:energy_line+1]) if re.search(r'(Input|Standard|Z-Matrix) orientation:',line)]
    geometry_line=orientations[-1]
    _,atoms,z,xyz=_orientation(lines,geometry_line)
    # Frequency output belongs to the final harmonic analysis, if present.
    harmonic=text.rfind('Harmonic frequencies')
    freqtext=text[harmonic:] if harmonic>=0 else ''
    frequencies=[float(x) for m in re.finditer(r'Frequencies --\s+([^\n]+)',freqtext) for x in m[1].split()]
    def last(pattern):
        matches=re.findall(pattern,text);return float(matches[-1]) if matches else None
    route=_route(lines);compact=re.sub(r'\s+','',route).upper()
    method='M062X' if 'M062X' in compact else ('B3LYP' if 'B3LYP' in compact else 'UNKNOWN')
    basis='6-311++G(D,P)' if '6-311++G(D,P)' in compact else ('6-311+G(2D,P)' if '6-311+G(2D,P)' in compact else 'UNKNOWN')
    charge=re.findall(r'Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)',text)[0]
    termination=[m.group(0) for m in re.finditer(r'Normal termination of Gaussian|Error termination',text)]
    tp=re.findall(r'Temperature\s+([\d.]+) Kelvin\.\s+Pressure\s+([\d.]+) Atm',text)
    return dict(path=str(path),sha256=sha256(path),energy_hartree=energy,energy_source_line=energy_line+1,geometry_source_line=geometry_line+1,
        atom_numbers=atoms,atomic_numbers=z,xyz=xyz,atom_count=len(atoms),composition=dict(collections.Counter(map(int,z))),
        route=route,method_basis=f'{method}/{basis}',solvent='SMD(water)' if 'SMD' in compact and 'SOLVENT=WATER' in compact else ('PCM(water)' if 'SOLVENT=WATER' in compact else 'UNKNOWN'),
        charge=int(charge[0]),multiplicity=int(charge[1]),optimization_completed='Optimization completed' in text,
        normal_termination=bool(termination and termination[-1].startswith('Normal')),frequency_count=len(frequencies),
        imaginary_frequencies=[f for f in frequencies if f<0],G_hartree=last(r'Sum of electronic and thermal Free Energies=\s*('+FLOAT+')'),
        H_hartree=last(r'Sum of electronic and thermal Enthalpies=\s*('+FLOAT+')'),
        ZPE_corrected_hartree=last(r'Sum of electronic and zero-point Energies=\s*('+FLOAT+')'),
        temperature_K=float(tp[-1][0]) if tp else None,pressure_atm=float(tp[-1][1]) if tp else None)


def serializable(record):
    return {k:v for k,v in record.items() if k not in ('xyz','atom_numbers','atomic_numbers')}


def aligned(mobile, reference):
    mc=mobile-mobile.mean(0);rc=reference-reference.mean(0)
    u,_,vt=np.linalg.svd(mc.T@rc);d=np.eye(3);d[-1,-1]=np.sign(np.linalg.det(u@vt))
    result=mc@(u@d@vt)+reference.mean(0)
    return result,float(np.sqrt(np.mean(np.sum((result-reference)**2,axis=1))))


def bond(frame,a,b):
    x=frame.set_index('atom_number')[['x','y','z']]
    return float(np.linalg.norm(x.loc[a]-x.loc[b]))


def write_oriented_log(frame,path,source,multiplier):
    """Canonical parsed view, not a newly executed or verbatim Gaussian output."""
    lines=['# DERIVED ANALYSIS LOG — NOT A GAUSSIAN CALCULATION',f'# Original: {source.name}',f'# Original SHA256: {sha256(source)}',
           f'# R_analysis = {multiplier:+d} * R_parser; C-N long to short.',
           '# Converged coordinate records only; TS geometry/energy from matching TS frequency LOG.',
           '# Cartesian forces unchanged under path reorientation; TS force unavailable.',
           '# Original point identifiers are recorded below. No new quantum calculation performed.']
    for _,p in frame.groupby('signed_irc_coordinate',sort=True):
        p=p.sort_values('atom_number');r=float(p.signed_irc_coordinate.iloc[0]);pathnum=0 if r==0 else (1 if r>0 else 2)
        point=int(p.point_number.iloc[0]);original_path=int(p.path_number.iloc[0])
        lines += [f'# Source path {original_path}, point {point}',f'Point Number: {point} Path Number: {pathnum}',
                  ' Standard orientation:',' ---------------------------------------------------------------------',
                  ' Center Atomic Atomic Coordinates (Angstroms)',' Number Number Type X Y Z',' ---------------------------------------------------------------------']
        for a in p.itertuples(): lines.append(f'{a.atom_number:6d} {a.atomic_number:6d} 0 {a.x:.12f} {a.y:.12f} {a.z:.12f}')
        lines += [' ---------------------------------------------------------------------',f' SCF Done: E(RM062X) = {p.energy_hartree.iloc[0]:.12f} A.U.']
        if p.fx.notna().all():
            lines += [' Forces (Hartrees/Bohr)',' ---------------------------------------------------------------------',' Center Atomic X Y Z']
            for a in p.itertuples():lines.append(f'{a.atom_number:6d} {a.atomic_number:6d} {a.fx:.12f} {a.fy:.12f} {a.fz:.12f}')
            lines.append(' ---------------------------------------------------------------------')
        lines.append(f' NET REACTION COORDINATE = {r:.12f}')
    path.write_text('\n'.join(lines)+'\n')


def prepare_irc(source,ts,system,n_atom,output):
    frame,meta=parse_log(source,system,1,ts['atom_count'],627.509474)
    dropped=frame.loc[frame.irc_coordinate_estimated,['path_number','point_number','signed_irc_coordinate']].drop_duplicates().to_dict('records')
    frame=frame[~frame.irc_coordinate_estimated].copy()
    # Validate atom identity before using separately supplied TS coordinates.
    midpoint=frame[frame.signed_irc_coordinate.abs()<1e-12].sort_values('atom_number')
    if not np.array_equal(midpoint.atomic_number.to_numpy(),ts['atomic_numbers']):raise ValueError(f'{system}: TS atom identities differ')
    checkpoint_E=float(midpoint.energy_hartree.iloc[0]);energy_error=(ts['energy_hartree']-checkpoint_E)*627.509474
    if abs(energy_error)>.01:raise ValueError(f'{system}: TS/IRC energy mismatch {energy_error}')
    xyz,rmsd=aligned(ts['xyz'],midpoint[['x','y','z']].to_numpy())
    if rmsd>.10:raise ValueError(f'{system}: TS/IRC structure mismatch {rmsd}')
    frame.loc[midpoint.index,['x','y','z']]=xyz
    frame.loc[midpoint.index,'energy_hartree']=ts['energy_hartree']
    frame['energy_relative_kcal']=(frame.energy_hartree-ts['energy_hartree'])*627.509474
    lo=frame[frame.signed_irc_coordinate==frame.signed_irc_coordinate.min()];hi=frame[frame.signed_irc_coordinate==frame.signed_irc_coordinate.max()]
    dlo,dhi=bond(lo,20,n_atom),bond(hi,20,n_atom)
    if abs(dlo-dhi)<.1:raise ValueError(f'{system}: ambiguous CN direction')
    multiplier=1 if dlo>dhi else -1
    frame['signed_irc_coordinate']*=multiplier
    frame=frame.sort_values(['signed_irc_coordinate','atom_number']).reset_index(drop=True)
    write_oriented_log(frame,output,source,multiplier)
    # The analysis uses the exported LOG, with a round-trip assertion.
    reread,_=parse_log(output,system,1,ts['atom_count'],627.509474)
    for cols in (['signed_irc_coordinate','energy_hartree','x','y','z'],['fx','fy','fz']):
        np.testing.assert_allclose(reread[cols],frame[cols],atol=1e-10,rtol=0,equal_nan=True)
    cn=[bond(p,20,n_atom) for _,p in reread.groupby('signed_irc_coordinate',sort=True)]
    audit=dict(system=system,source=str(source),oriented_log=str(output),coordinate_multiplier=multiplier,
        original_negative_CN_A=dlo,original_positive_CN_A=dhi,reactant_side_CN_A=cn[0],product_side_CN_A=cn[-1],
        CN_monotonically_decreasing=bool(np.all(np.diff(cn)<=1e-5)),max_CN_increase_A=float(max(0,np.diff(cn).max())),
        point_count=len(cn),TS_midpoint_alignment_RMSD_A=rmsd,TS_energy_difference_kcal_mol=energy_error,
        excluded_unconverged_points=dropped,original_parser_warnings=meta['warnings'])
    return reread,audit
