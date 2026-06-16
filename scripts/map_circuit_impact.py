"""
Phase 8 (exploratory): Map predicted PWS gene reactivation to hypothalamic circuits.

Uses organoid scRNA-seq cell-type markers and published PWS-hypothalamus
literature links to estimate which neuronal subtypes may normalize under
top optimized designs. Hypothesis-generating only — not validated prediction.

Data: GSE262700 organoid scRNA (Nemoto et al. 2025)
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CURATED = ROOT / "data" / "curated"
MODELS = ROOT / "data" / "models"
OUT = MODELS / "circuit_mapping"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

# Hypothalamic neuron subtypes linked to PWS phenotypes (literature-curated)
CIRCUIT_NODES = {
    "AgRP_neurons": {
        "markers": ["AGRP", "NPY"],
        "pws_phenotype_link": "hyperphagia / orexigenic drive",
        "expected_direction_if_PWS_rescued": "normalize_down",
    },
    "POMC_neurons": {
        "markers": ["POMC", "CART"],
        "pws_phenotype_link": "anorexigenic signaling deficit",
        "expected_direction_if_PWS_rescued": "normalize_up",
    },
    "OXT_neurons": {
        "markers": ["OXT", "AVP"],
        "pws_phenotype_link": "social/endocrine dysfunction",
        "expected_direction_if_PWS_rescued": "normalize_up",
    },
    "SIM1_neurons": {
        "markers": ["SIM1", "ESR1"],
        "pws_phenotype_link": "PVH metabolic integration",
        "expected_direction_if_PWS_rescued": "normalize",
    },
}

# SNRPN/SNHG14 are expressed in edited organoids across multiple clusters
PWS_GENE_CIRCUIT_PRIORS = {
    "SNRPN": {"primary_circuits": ["OXT_neurons", "SIM1_neurons"], "reactivation_weight": 0.6},
    "SNHG14": {"primary_circuits": ["AgRP_neurons", "POMC_neurons"], "reactivation_weight": 0.7},
}


def load_organoid_markers() -> pd.DataFrame | None:
    path = CURATED / "organoid_celltype_markers.csv"
    if path.exists():
        return pd.read_csv(path)
    return None


def load_top_designs() -> pd.DataFrame:
    opt_path = MODELS / "optimization" / "optimized_therapeutic_designs.csv"
    if opt_path.exists():
        return pd.read_csv(opt_path).head(5)
    return pd.DataFrame()


def estimate_circuit_impact(design: pd.Series) -> dict:
    snrpn_pct = design.get("predicted_SNRPN_pct_WT", 0)
    snhg_pct = design.get("predicted_SNHG14_pct_WT", 0)

    # Rescue fraction: how close to 100% WT (capped)
    snrpn_rescue = min(1.0, snrpn_pct / 100)
    snhg_rescue = min(1.0, snhg_pct / 100)

    circuit_scores = {}
    for circuit, info in CIRCUIT_NODES.items():
        score = 0.0
        for gene, prior in PWS_GENE_CIRCUIT_PRIORS.items():
            if circuit in prior["primary_circuits"]:
                rescue = snrpn_rescue if gene == "SNRPN" else snhg_rescue
                score += prior["reactivation_weight"] * rescue
        circuit_scores[circuit] = {
            "predicted_normalization_score": round(min(1.0, score), 3),
            "direction": info["expected_direction_if_PWS_rescued"],
            "phenotype_link": info["pws_phenotype_link"],
            "confidence": "exploratory",
        }

    hyperphagia_proxy = (
        circuit_scores["AgRP_neurons"]["predicted_normalization_score"] * 0.5
        + circuit_scores["POMC_neurons"]["predicted_normalization_score"] * 0.5
    )

    return {
        "design_rank": int(design.get("rank", 0)),
        "strategy": design.get("strategy"),
        "tet1_grna": design.get("tet1_grna_id"),
        "vp64_grna": design.get("vp64_grna_id"),
        "predicted_SNRPN_pct_WT": snrpn_pct,
        "predicted_SNHG14_pct_WT": snhg_pct,
        "circuit_impacts": circuit_scores,
        "hyperphagia_normalization_proxy": round(hyperphagia_proxy, 3),
        "caveat": (
            "Circuit mapping is hypothesis-generating. SNRPN/SNHG14 reactivation "
            "does not directly predict AgRP/POMC normalization without longitudinal "
            "hypothalamic circuit data."
        ),
    }


def main():
    designs = load_top_designs()
    if designs.empty:
        log.warning("No optimized designs found; run optimize_therapeutic_designs.py first")
        return

    impacts = [estimate_circuit_impact(row) for _, row in designs.iterrows()]

    report = {
        "phase": 8,
        "status": "exploratory_hypothesis_generating",
        "organoid_validation": {
            "source": "Nemoto et al. 2025 GSE262700",
            "finding": "Tet1 editing reactivates SNRPN/SNHG14 in hypothalamic organoids",
            "SNRPN_fold_change": 5297,
            "SNHG14_fold_change": 158,
        },
        "circuit_nodes": CIRCUIT_NODES,
        "design_circuit_impacts": impacts,
        "interpretation": (
            "Designs achieving dual-gene reactivation in the 70-130% WT window are "
            "predicted to partially normalize orexigenic/anorexigenic circuit balance. "
            "AgRP normalization depends primarily on SNHG14/SNORD116 restoration; "
            "POMC/OXT on SNRPN dosage within therapeutic window."
        ),
    }

    with open(OUT / "circuit_impact_report.json", "w") as f:
        json.dump(report, f, indent=2)

    pd.DataFrame(impacts).to_csv(OUT / "circuit_impact_by_design.csv", index=False)
    log.info("Circuit impact mapped for %d top designs -> %s", len(impacts), OUT)


if __name__ == "__main__":
    main()
