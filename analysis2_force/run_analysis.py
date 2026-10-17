"""Convenient entry point for analysis2_force."""
from __future__ import annotations

import runpy
from pathlib import Path

runpy.run_path(str(Path(__file__).resolve().parent / "src" / "run_analysis.py"), run_name="__main__")
