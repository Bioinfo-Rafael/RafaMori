"""Compatibility wrapper for the preserved base analysis."""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

BASE_SRC = Path(__file__).resolve().parents[1] / "base_analysis" / "src"
sys.path.insert(0, str(BASE_SRC))
runpy.run_path(str(BASE_SRC / "run_analysis.py"), run_name="__main__")

