"""Run Phase 11: advanced analyses (CRISPRepi, Pareto, GNN, ODE, Bayesian, novelty)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

PHASE_11 = [
    "curate_crisprepi.py",
    "curate_neuron_expression.py",
    "parse_organoid_celltypes.py",
    "train_graph_propagation_model.py",
    "pareto_optimize_designs.py",
    "bayesian_optimize_designs.py",
    "simulate_circuit_ode.py",
    "map_circuit_impact.py",
    "validate_neuron_concordance.py",
    "build_novelty_synthesis.py",
    "build_final_catalog.py",
]


def main():
    for script in PHASE_11:
        print(f"\n{'='*60}\nRunning {script}\n{'='*60}")
        r = subprocess.run([sys.executable, str(ROOT / "scripts" / script)], cwd=str(ROOT))
        if r.returncode != 0:
            sys.exit(r.returncode)
    print("\nPhase 11 advanced analyses complete.")


if __name__ == "__main__":
    main()
