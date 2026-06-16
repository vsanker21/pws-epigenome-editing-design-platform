"""
Run Phase 3 curation pipeline end-to-end.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPTS = [
    "curate_editing_screens.py",
    "curate_expression_methylation.py",
    "curate_hypothalamic_locus.py",
    "build_digital_twin.py",
]

ROOT = Path(__file__).resolve().parents[1]


def main():
    for script in SCRIPTS:
        path = ROOT / "scripts" / script
        print(f"\n{'='*60}\nRunning {script}\n{'='*60}")
        result = subprocess.run([sys.executable, str(path)], cwd=str(ROOT))
        if result.returncode != 0:
            sys.exit(result.returncode)
    print("\nPhase 3 curation complete.")


if __name__ == "__main__":
    main()
