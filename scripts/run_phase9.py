"""Run Phase 9: validation suite + final catalog synthesis."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    "validate_held_out.py",
    "validate_organoid_concordance.py",
    "assess_collateral_imprinting.py",
    "build_final_catalog.py",
]


def main():
    for script in SCRIPTS:
        print(f"\n{'='*60}\nRunning {script}\n{'='*60}")
        r = subprocess.run([sys.executable, str(ROOT / "scripts" / script)], cwd=str(ROOT))
        if r.returncode != 0:
            sys.exit(r.returncode)
    print("\nPhase 9 validation complete.")


if __name__ == "__main__":
    main()
