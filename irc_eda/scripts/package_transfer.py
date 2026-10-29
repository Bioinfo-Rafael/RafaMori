"""Build a portable, hash-checked source/data ZIP; no licensed binaries or venv."""
from pathlib import Path
import hashlib
import json
import zipfile

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
DATA = REPO / 'data/completed_20261017'
REQUIRED = [
    'IRC/MORIPF_IRC1.LOG.gz', 'IRC/RAFA14_IRC_MERGED.LOG.gz',
    'IRC/RAFA15_IRC2.LOG.gz', 'IRC/RAFA16_IRC_COMPLETE_MERGED.LOG.gz',
    'TS/MORIPF_TS.LOG.gz', 'TS/RAFA14_TS2_FREQ.LOG.gz',
    'TS/RAFA15_TS6.LOG.gz', 'TS/RAFA16_TS.LOG.gz',
    'TS/RAFA14_TS2_NBO.LOG.gz', 'TS/RAFA15_TS6_NBO.LOG.gz',
    'TS/RAFA14_TS2_NBO.fchk.gz', 'TS/RAFA15_TS6_NBO.fchk.gz',
]
SKIP = {'.venv', 'vendor', 'outputs', '__pycache__', '.pytest_cache', 'reports'}

def build():
    files = [p for p in ROOT.rglob('*') if p.is_file() and not (set(p.relative_to(ROOT).parts) & SKIP)
             and p.name != 'local.yaml']
    files += [DATA / name for name in REQUIRED]
    for p in files:
        if not p.is_file():
            raise FileNotFoundError(p)
    manifest = {str(p.relative_to(REPO)).replace('\\','/'): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(files)}
    destination = ROOT / 'outputs/transfer'
    destination.mkdir(parents=True, exist_ok=True)
    target = destination / 'RafaMori-irc-eda-transfer.zip'
    with zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED) as archive:
        for name in manifest:
            archive.write(REPO / name, name)
        archive.writestr('transfer_manifest.json', json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'zip':str(target), 'file_count':len(files), 'bytes':target.stat().st_size,
                      'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}, indent=2))
    return target

if __name__ == '__main__':
    build()
