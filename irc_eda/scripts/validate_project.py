"""Real-data round-trip validation, without invoking any quantum chemistry engine."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from irc_eda.pipeline import load_config,dump,csv_write
from irc_eda.parsing import digest,TABLE
from irc_eda.upstream import ROOT
c=load_config(ROOT/'config/project.yaml');out=Path(c['output_dir'])
points=json.loads((out/'points.json').read_text());rows=[]
for system in c['systems']:
    ps=[p for p in points if p['system']==system]
    previous=ROOT.parent/'analysis4_force/outputs/csv'/f'reaction_force_{system}.csv'
    old=pd.read_csv(previous)
    assert len(old)==len(ps)
    maxr=float(np.max(np.abs(old.signed_irc_coordinate-np.array([p['reaction_coordinate'] for p in ps]))))
    maxe=float(np.max(np.abs(old.energy_hartree-np.array([p['energy_hartree'] for p in ps]))))
    assert maxr<1e-9 and maxe<1e-9
    xyzerror=0.
    for p in ps:
        file=out/'xyz'/f'{p["point_id"]}.xyz';lines=file.read_text().splitlines()
        assert int(lines[0])==len(p['xyz'])
        atoms=[x.split() for x in lines[2:]]
        assert [a[0] for a in atoms]==[TABLE.element[z] for z in p['atomic_numbers']]
        xyz=np.array([[float(x) for x in a[1:]] for a in atoms])
        err=float(np.max(np.abs(xyz-np.array(p['xyz']))));xyzerror=max(xyzerror,err)
        assert err<1e-9 and digest(p['source'])==p['source_sha256']
    rows.append(dict(system=system,point_count=len(ps),existing_profile=str(previous),
                     max_coordinate_difference=maxr,max_energy_difference_hartree=maxe,
                     xyz_max_difference_angstrom=xyzerror,source_hashes_unchanged=True,
                     cclib_checked_non_TS_points=sum(p['cclib_crosscheck'] is True for p in ps),
                     max_corrector_minus_scf_hartree=max((abs(p['corrector_minus_scf_hartree']) for p in ps
                                                        if p['corrector_minus_scf_hartree'] is not None),default=None)))
csv_write(ROOT/'reports/project_validation.csv',rows);dump(ROOT/'reports/project_validation.json',rows)
print(json.dumps(rows,indent=2))
