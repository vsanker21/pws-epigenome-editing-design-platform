"""
Phase 10c: Generate structured experimental validation protocol.

Produces a wet-lab-ready protocol outline for testing top in-silico designs
in PWS UPD15 hypothalamic neurons/organoids, aligned with published methods
(Nemoto 2025 organoids; Rohm 2025 iPSC editing).
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "data" / "models"
OUT = MODELS / "experimental_protocol"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)


def load_top_designs(n: int = 3) -> list[dict]:
    catalog = pd.read_csv(MODELS / "final_catalog" / "pws_therapeutic_design_catalog.csv")
    return catalog.head(n).to_dict(orient="records")


def build_protocol() -> dict:
    designs = load_top_designs(3)

    return {
        "protocol_version": "1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "title": "Experimental validation of in-silico PWS epigenome-editing designs (maternal UPD15)",
        "objective": (
            "Test whether computationally prioritized hybrid Tet1+VP64 designs "
            "achieve dual-gene reactivation (SNRPN + SNHG14) within 70-130% of "
            "paternal-equivalent expression in hypothalamic neurons, with acceptable "
            "collateral imprinting profile."
        ),
        "model_system": {
            "primary": "PWS UPD15 iPSC-derived hypothalamic neurons or organoids",
            "rationale": "Clinically relevant cell type; Nemoto 2025 demonstrated Tet1-class reactivation in this system",
            "reference": "Nemoto et al. Nat Commun 2025 (GSE262700)",
            "alternative": "PWS UPD15 iPSC neural progenitors (Rohm 2025 protocol)",
        },
        "test_arms": [
            {
                "arm_id": "A",
                "label": "Top hybrid (bisulfite-validated Tet1 + minimal VP64)",
                "design": designs[0] if designs else None,
                "editors": ["dCas9-SunTag-TET1 or dCas9-Tet1", "dCas9-VP64 (minimal dose ratio ~6%)"],
                "priority": "primary",
            },
            {
                "arm_id": "B",
                "label": "Tet1 monotherapy comparator",
                "design": next((d for d in designs if d.get("strategy") == "Tet1_monotherapy"), designs[0] if designs else None),
                "editors": ["dCas9-SunTag-TET1 or dCas9-Tet1"],
                "priority": "comparator",
            },
            {
                "arm_id": "C",
                "label": "Non-targeting / unedited PWS UPD15 control",
                "design": None,
                "editors": ["non-targeting gRNA or no editor"],
                "priority": "control",
            },
            {
                "arm_id": "D",
                "label": "Healthy (WT) iPSC comparator",
                "design": None,
                "editors": ["none"],
                "priority": "reference",
            },
        ],
        "delivery": {
            "options": ["AAV (hypothalamic tropism)", "LNP-mRNA (transient)", "RNP electroporation (iPSC)"],
            "recommended_first": "RNP electroporation in iPSC (Rohm 2025) or mRNA/LNP for dose titration",
            "dose_titration": "Titrate VP64 to avoid SNRPN >130% WT; start at in-silico w_vp64 ~0.06 relative to Tet1",
        },
        "primary_endpoints": [
            {
                "endpoint": "PWS-ICR CpG methylation",
                "method": "Targeted bisulfite amplicon sequencing (GSE285300 loci)",
                "success_criterion": "Edited methylation approaches WT (~40%) from PWS baseline (~95%)",
                "guides": ["PWS_G.4363 amplicon region"],
            },
            {
                "endpoint": "SNRPN expression",
                "method": "qPCR and/or bulk RNA-seq",
                "success_criterion": "70-130% of WT iPSC SNRPN (paternal-equivalent window)",
            },
            {
                "endpoint": "SNHG14 / SNORD116 expression",
                "method": "qPCR and/or RNA-seq",
                "success_criterion": "70-130% of WT; SNHG14 reactivation confirms ICR demethylation",
            },
        ],
        "secondary_endpoints": [
            {"endpoint": "UBE3A expression", "rationale": "maternally expressed — collateral imprinting check"},
            {"endpoint": "ATP10A expression", "rationale": "maternally expressed — collateral check"},
            {"endpoint": "GABRB3/GABRA5 expression", "rationale": "GABAA cluster preservation"},
            {"endpoint": "Genome-wide methylation (optional)", "method": "Illumina EPIC array or WGBS"},
        ],
        "sample_size_guidance": {
            "minimum_per_arm": "n=3 biological replicates",
            "rationale": "Consistent with Rohm 2025 and Nemoto 2025 replicate structure",
            "power_note": "Formal power calculation requires pilot variance estimate from Arm A",
        },
        "analysis_plan": [
            "Compare measured SNRPN/SNHG14 %WT to forward model predictions",
            "Compute Pearson/Spearman correlation with in-silico ranked designs (model refinement)",
            "Report collateral gene expression as fold-change vs unedited PWS control",
            "Update phenomenological hybrid dose model with empirical dose-response",
        ],
        "go_no_go_criteria": {
            "go": "Both SNRPN and SNHG14 in 70-130% WT; no >2-fold collateral disruption at UBE3A/GABRB3",
            "no_go": "SNRPN >200% WT (overshoot) or SNHG14 <50% WT without SNRPN rescue",
            "iterate": "SNRPN 50-70% WT with SNHG14 in window — reduce VP64 dose and retest",
        },
        "computational_feedback_loop": (
            "Measured outcomes feed back into train_forward_model.py editor priors "
            "and optimize_therapeutic_designs.py dose parameters for v2 model."
        ),
    }


def main():
    protocol = build_protocol()
    with open(OUT / "experimental_validation_protocol.json", "w") as f:
        json.dump(protocol, f, indent=2, default=str)

    log.info("Protocol generated: %d test arms, %d primary endpoints",
             len(protocol["test_arms"]), len(protocol["primary_endpoints"]))
    log.info("Primary arm A: %s + %s",
             protocol["test_arms"][0]["design"].get("tet1_grna_id") if protocol["test_arms"][0]["design"] else "?",
             protocol["test_arms"][0]["design"].get("vp64_grna_id") if protocol["test_arms"][0]["design"] else "?")
    log.info("Output -> %s", OUT)


if __name__ == "__main__":
    main()
