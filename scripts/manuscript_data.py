"""Load and expose all project data for manuscript generation."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(r"g:\Prader Willi Syndrome")

GITHUB_REPO = "https://github.com/vsanker21/pws-epigenome-editing-design-platform"
ZENODO_DOI = "10.5281/zenodo.20723805"  # Replace after Zenodo archives GitHub release v1.0.0
ZENODO_URL = f"https://doi.org/{ZENODO_DOI}"
CODE_AVAILABILITY = (
    f"{GITHUB_REPO} (archived on Zenodo: {ZENODO_URL})"
)


def load_json(rel: str) -> dict:
    with open(ROOT / rel, encoding="utf-8") as f:
        return json.load(f)


def load_csv(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel)


# Core reports
held_out = load_json("data/models/validation/held_out_benchmark.json")
organoid = load_json("data/models/validation/organoid_concordance.json")
forward = load_json("data/models/forward_model_report.json")
sensitivity = load_json("data/models/sensitivity/sensitivity_report.json")
optimization = load_json("data/models/optimization/optimization_report.json")
synthesis = load_json("data/models/final_catalog/project_synthesis_report.json")
collateral = load_json("data/models/validation/collateral_imprinting_risk.json")
gse152098 = load_json("data/curated/gse152098_reprocessing_report.json")
dmr_map = load_json("data/curated/methylation_dmr_mapping_report.json")
encode = load_json("data/encode_reference/summary.json")
cas_off = load_json("data/models/validation/cas_offinder_offtarget.json")
integration = load_json("data/integrated/integration_summary.json")
locus_summary = load_json("data/locus_annotation/locus_summary.json")
neuron_conc = load_json("data/models/validation/neuron_ipsc_concordance.json")
locus_ot = load_json("data/models/validation/locus_offtarget_assessment.json")
pareto = load_json("data/models/pareto/pareto_report.json")
bayesian = load_json("data/models/bayesian_optimization/bayesian_dose_optimization.json")
graph = load_json("data/models/graph_model/graph_propagation_report.json")
circuit_ode = load_json("data/models/circuit_ode/circuit_ode_report.json")
crisprepi = load_json("data/curated/crisprepi_transfer_priors.json")
accessibility = load_json("data/models/cell_type_scoring/accessibility_report.json")
novelty = load_json("data/models/novelty_synthesis/novelty_synthesis_report.json")
protocol = load_json("data/models/experimental_protocol/experimental_validation_protocol.json")
genome_ot = load_json("data/models/validation/genome_wide_offtarget.json")

# DataFrames
catalog = load_csv("data/models/final_catalog/pws_therapeutic_design_catalog.csv")
tet1_ranked = load_csv("data/models/therapeutic_candidates_ranked.csv")
pareto_df = load_csv("data/models/pareto/pareto_therapeutic_designs.csv")
circuit_sim = load_csv("data/models/circuit_ode/circuit_ode_simulations.csv")
organoid_fc = load_csv("data/curated/organoid_pws_foldchange.csv")
sensitivity_df = load_csv("data/models/sensitivity/catalog_window_probability.csv")
gp_samples = load_csv("data/models/bayesian_optimization/gp_training_samples.csv")
gse28525 = load_csv("data/curated/gse28525_pws_probe_methylation.csv")
forward_pred = load_csv("data/models/forward_model_predictions.csv")

TOP = synthesis["top_recommended_design"]
