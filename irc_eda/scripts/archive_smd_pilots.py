"""Archive or verify successful pilot inputs/outputs; never launch Gaussian."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from irc_eda.scripts.smd_distortion_pilot import OUT,REPORT,validate,save,analyze

SYSTEMS=['RAFA14','MoriPF','RAFA15','RAFA16']
ARCHIVE=ROOT/'irc_eda/data/smd_pilots'

def sha(data):return hashlib.sha256(data).hexdigest()

def archive():
    entries=[]
    for system in SYSTEMS:
        manifest=json.loads((OUT/system/'manifest.json').read_text())
        for job in manifest['jobs']:
            validation=validate(job)
            folder=ROOT/job['folder']
            recorded=json.loads((folder/'validation.json').read_text())
            assert recorded['status']=='SUCCESS' and recorded['output_sha256']==validation['output_sha256']
            target=ARCHIVE/system/(job['role']+'_'+job['fragment'])
            target.mkdir(parents=True,exist_ok=True)
            inp=(folder/'fragment.gjf').read_bytes();output=(folder/'fragment.out').read_bytes()
            for name,data in [('fragment.gjf',inp),('fragment.out.gz',gzip.compress(output,mtime=0))]:
                path=target/name
                if path.exists():assert path.read_bytes()==data, 'Archive differs: '+str(path)
                else:path.write_bytes(data)
            entries.append({'system':system,'role':job['role'],'fragment':job['fragment'],
                'input':(target/'fragment.gjf').relative_to(ARCHIVE).as_posix(),
                'output_gzip':(target/'fragment.out.gz').relative_to(ARCHIVE).as_posix(),
                'input_sha256':sha(inp),'uncompressed_output_sha256':sha(output),
                'gzip_sha256':sha((target/'fragment.out.gz').read_bytes()),'uncompressed_output_bytes':len(output)})
        save(ARCHIVE/system/'manifest.json',manifest)
    save(ARCHIVE/'index.json',{'format_version':1,'successful_jobs':len(entries),'entries':entries})

def verify(reanalyze=False):
    index=json.loads((ARCHIVE/'index.json').read_text())
    assert index['successful_jobs']==24 and len(index['entries'])==24
    OUT.mkdir(parents=True,exist_ok=True)
    checks=[]
    # Isolated reconstruction: does not use original job folders or checkpoints.
    with tempfile.TemporaryDirectory(prefix='archive_verify_',dir=OUT) as scratch:
        scratch=Path(scratch)
        for system in SYSTEMS:
            manifest=json.loads((ARCHIVE/system/'manifest.json').read_text())
            entries=[e for e in index['entries'] if e['system']==system]
            assert len(entries)==6
            for job in manifest['jobs']:
                entry=next(e for e in entries if e['role']==job['role'] and e['fragment']==job['fragment'])
                inp=(ARCHIVE/entry['input']).read_bytes();packed=(ARCHIVE/entry['output_gzip']).read_bytes()
                output=gzip.decompress(packed)
                assert sha(inp)==entry['input_sha256']==job['input_sha256']
                assert sha(packed)==entry['gzip_sha256']
                assert sha(output)==entry['uncompressed_output_sha256']
                folder=scratch/system/(job['role']+'_'+job['fragment']);folder.mkdir(parents=True)
                (folder/'fragment.gjf').write_bytes(inp);(folder/'fragment.out').write_bytes(output)
                job['folder']=folder.relative_to(ROOT).as_posix()
                result=validate(job)
                checks.append({'system':system,'role':job['role'],'fragment':job['fragment'],
                    'energy_hartree':result['energy_hartree'],'max_coordinate_difference_A':result['max_coordinate_difference_A'],
                    'input_sha256':entry['input_sha256'],'output_sha256':entry['uncompressed_output_sha256'],'status':'SUCCESS'})
            if reanalyze:
                # Redirect only manifest input; reports use the usual tracked path.
                import irc_eda.scripts.smd_distortion_pilot as pilot
                save(scratch/'manifests'/system/'manifest.json',manifest)
                original_out=pilot.OUT
                try:
                    pilot.OUT=scratch/'manifests';analyze(system)
                finally:pilot.OUT=original_out
                # Persist original provenance, not temporary verification paths.
                original=json.loads((ARCHIVE/system/'manifest.json').read_text())
                save(REPORT/system/'source_audit.json',original)
                records=json.loads((REPORT/system/'fragment_validation.json').read_text())
                for record,job in zip(records,original['jobs']):
                    record['job']=job
                    record['validation']['output']=job['folder'].replace('\\','/')+'/fragment.out'
                save(REPORT/system/'fragment_validation.json',records)
    save(REPORT/'archive_validation.json',{'successful_fragment_outputs':len(checks),'verified_from':'tracked .gjf and .out.gz only','checks':checks})
    print('Verified 24 archived successful Gaussian outputs; Gaussian launches: 0')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['pack','verify','reanalyze']);args=parser.parse_args()
    if args.action=='pack':archive()
    verify(reanalyze=args.action=='reanalyze')
