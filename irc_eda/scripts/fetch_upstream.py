"""Clone exact upstream revisions, preserving local checkouts and their licenses."""
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from irc_eda.upstream import PINS,ROOT,git_revision
for name,(url,rev) in PINS.items():
    target=ROOT/'vendor'/name
    if not target.exists():
        subprocess.run(['git','clone',url,str(target)],check=True)
        subprocess.run(['git','-C',str(target),'checkout','--detach',rev],check=True)
    if git_revision(target)!=rev:
        raise SystemExit(f'{target}: existing revision differs; checkout is preserved. Expected {rev}')
    if subprocess.check_output(['git','-C',str(target),'status','--porcelain'],text=True).strip():
        raise SystemExit(f'{target}: local modifications detected; preserved, please review')
(ROOT/'upstream-lock.json').write_text(json.dumps({n:dict(url=u,commit=r) for n,(u,r) in PINS.items()},indent=2)+'\n')
