#!/usr/bin/env bash
# The official scripts are used unchanged; no replacement compiler or DFT engine.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export APOST3D_PATH="$(cd -- "$script_dir/../vendor/APOST3D" && pwd)"
if ! command -v ifort >/dev/null 2>&1; then
  echo 'BLOCKED: ifort missing. Load a licensed/available Intel oneAPI 2024 Linux environment.' >&2
  exit 2
fi
if ! command -v icx >/dev/null 2>&1; then
  echo 'BLOCKED: icx missing; official Libxc build requires it.' >&2
  exit 2
fi
cd "$APOST3D_PATH"
if [[ ! -f libxc-4.2.3.tar.gz && ! -f libxc-4.2.3/lib/libxcf90.a ]]; then
  echo 'BLOCKED: pinned upstream does not include libxc-4.2.3.tar.gz; obtain the exact upstream-compatible archive. Do not use a newer Libxc silently.' >&2
  exit 2
fi
if [[ ! -f libxc-4.2.3/lib/libxcf90.a ]]; then bash -e ./compile_libxc.sh; fi
bash -e ./make_compile.sh
[[ -x apost3d && -x utils/get_energy_g16 ]]
