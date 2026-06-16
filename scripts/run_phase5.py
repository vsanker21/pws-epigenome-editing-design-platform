"""Run Phase 5: organoid parsing + forward model training."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ["parse_organoid_scrna.py", "train_forward_model.py"]


def main():
    for script in SCRIPTS:
        print(f"\n{'='*60}\nRunning {script}\n{'='*60}")
        r = subprocess.run([sys.executable, str(ROOT / "scripts" / script)], cwd=str(ROOT))
        if r.returncode != 0:
            sys.exit(r.returncode)
    print("\nPhase 5 complete.")


if __name__ == "__main__":
    main()
