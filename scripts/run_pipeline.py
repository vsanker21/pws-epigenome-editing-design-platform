"""Master pipeline: run full PWS epigenome-editing design platform end-to-end."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

PHASES = {
    "3": ["run_phase3.py", "integrate_and_merge.py"],
    "5": ["run_phase5.py"],
    "6": ["run_phase6.py"],
    "9": ["run_phase9.py"],
    "10": [
        "score_hypothalamic_accessibility.py",
        "assess_locus_offtarget.py",
        "genome_wide_offtarget.py",
        "generate_experimental_protocol.py",
        "build_final_catalog.py",
    ],
    "11": [
        "curate_crisprepi.py",
        "curate_neuron_expression.py",
        "parse_organoid_celltypes.py",
        "train_graph_propagation_model.py",
        "pareto_optimize_designs.py",
        "bayesian_optimize_designs.py",
        "simulate_circuit_ode.py",
        "map_circuit_impact.py",
        "sensitivity_analysis.py",
        "build_novelty_synthesis.py",
        "build_final_catalog.py",
    ],
    "12": [
        "cas_offinder_scan.py",
        "download_encode_reference.py",
        "reprocess_gse152098_raw.py",
        "map_methylation_dmrs.py",
        "sensitivity_analysis.py",
        "build_final_catalog.py",
        "build_novelty_synthesis.py",
    ],
}

ALL_SCRIPTS = [
    "run_phase3.py",
    "integrate_and_merge.py",
    "run_phase5.py",
    "run_phase6.py",
    "run_phase9.py",
    "score_hypothalamic_accessibility.py",
    "assess_locus_offtarget.py",
    "genome_wide_offtarget.py",
    "generate_experimental_protocol.py",
    "curate_crisprepi.py",
    "curate_neuron_expression.py",
    "parse_organoid_celltypes.py",
    "train_graph_propagation_model.py",
    "pareto_optimize_designs.py",
    "bayesian_optimize_designs.py",
    "simulate_circuit_ode.py",
    "map_circuit_impact.py",
    "sensitivity_analysis.py",
    "cas_offinder_scan.py",
    "download_encode_reference.py",
    "reprocess_gse152098_raw.py",
    "map_methylation_dmrs.py",
    "build_novelty_synthesis.py",
    "build_final_catalog.py",
]


def run_script(name: str) -> int:
    path = ROOT / "scripts" / name
    print(f"\n{'='*60}\nRunning {name}\n{'='*60}")
    return subprocess.run([sys.executable, str(path)], cwd=str(ROOT)).returncode


def main():
    parser = argparse.ArgumentParser(description="PWS epigenome-editing design pipeline")
    parser.add_argument(
        "--phase",
        choices=["3", "5", "6", "9", "10", "11", "12", "all"],
        default="11",
        help="Which phase(s) to run (default: 11 = advanced analyses)",
    )
    args = parser.parse_args()

    if args.phase == "all":
        scripts = ALL_SCRIPTS
    else:
        scripts = PHASES[args.phase]

    for script in scripts:
        if run_script(script) != 0:
            sys.exit(1)

    print(f"\nPipeline phase '{args.phase}' complete.")
    if args.phase in ("10", "11", "12", "all"):
        print("\nKey outputs:")
        print("  data/models/final_catalog/pws_therapeutic_design_catalog.csv")
        print("  data/models/final_catalog/project_synthesis_report.json")
        print("  data/models/experimental_protocol/experimental_validation_protocol.json")
        print("  data/models/validation/held_out_benchmark.json")


if __name__ == "__main__":
    main()
