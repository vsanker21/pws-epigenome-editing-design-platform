"""Phase 12: Advanced regulatory + off-target + methylation DMR integration."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SCRIPTS = [
    "cas_offinder_scan.py",
    "download_encode_reference.py",
    "reprocess_gse152098_raw.py",
    "map_methylation_dmrs.py",
    "sensitivity_analysis.py",
    "build_final_catalog.py",
    "build_novelty_synthesis.py",
]


def run_script(name: str) -> int:
    path = ROOT / "scripts" / name
    print(f"\n{'='*60}\nRunning {name}\n{'='*60}")
    return subprocess.run([sys.executable, str(path)], cwd=str(ROOT)).returncode


def main():
    for script in SCRIPTS:
        if run_script(script) != 0:
            sys.exit(1)
    print("\nPhase 12 complete.")


if __name__ == "__main__":
    main()
