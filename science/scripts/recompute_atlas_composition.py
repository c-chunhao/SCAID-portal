#!/usr/bin/env python3
"""Current entry: coverage-aware complete specimen-grid recomputation."""
from pathlib import Path
import runpy

if __name__ == '__main__':
    target = Path(__file__).resolve().parents[1] / 'analysis_vNext/statistics/recompute_atlas_composition.py'
    runpy.run_path(str(target), run_name='__main__')
