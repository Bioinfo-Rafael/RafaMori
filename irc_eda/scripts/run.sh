#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd -- "$script_dir/../.." && pwd)"
cd "$repo_dir"
exec "$repo_dir/irc_eda/.venv/bin/python" -m irc_eda "$@"
