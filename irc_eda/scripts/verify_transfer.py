"""Verify an extracted transfer ZIP using only Python's standard library."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parents[2]
manifest_file = root / 'transfer_manifest.json'
if not manifest_file.exists():
    raise SystemExit('No transfer_manifest.json: this verifier is for the extracted transfer ZIP.')
manifest = json.loads(manifest_file.read_text())
bad = [name for name, expected in manifest.items()
       if not (root/name).is_file() or hashlib.sha256((root/name).read_bytes()).hexdigest()!=expected]
if bad:
    raise SystemExit('Missing or modified files: '+', '.join(bad))
print(f'OK: {len(manifest)} source/data files match their SHA-256 hashes.')
