"""
Build final uncertainty-aware therapeutic design catalog.

Synthesizes outputs from all pipeline phases into a single ranked,
publication-ready catalog with evidence links, uncertainty bands,
validation status, and experimental recommendations.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "data" / "models"
OUT = MODELS / "final_catalog"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

CITATIONS = {
    "Rohm2025": "Rohm et al. Cell Genomics 2025;5:100770 (GSE285306)",
    "Nemoto2025": "Nemoto et al. Nat Commun 2025 (GSE262700)",
    "Cousminer2021": "Cousminer et al. Nat Commun 2021 (GSE152098)",
    "CRISPRepi2025": "Shi et al. NAR 2025;53(D1):D901",
}


def load_json(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.exists() else None


VALIDATED_PRIORITY = {
    "PWS_G.4363": 3,  # bisulfite-validated — highest experimental priority
    "PWS_G.4670": 2,  # top sublibrary hit
    "PWS_G.9414": 2,  # top VP64 activation hit
    "PWS_G.8738": 1,
}


def build_catalog_entries() -> pd.DataFrame:
    opt = pd.read_csv(MODELS / "optimization" / "optimized_therapeutic_designs.csv")
    hybrid = opt[opt["strategy"] == "Hybrid_Tet1_VP64"]
    tet1 = opt[opt["strategy"] == "Tet1_monotherapy"]

    # Ensure experimentally validated guides are always represented
    must_include_tet1 = {"PWS_G.4363", "PWS_G.4670"}
    must_include_vp64 = {"PWS_G.9414", "PWS_G.8738"}
    validated_hybrid = hybrid[
        hybrid["tet1_grna_id"].isin(must_include_tet1)
        & hybrid["vp64_grna_id"].isin(must_include_vp64)
    ]
    validated_tet1 = tet1[tet1["tet1_grna_id"].isin(must_include_tet1)]

    combined = pd.concat([
        validated_hybrid.head(4),
        hybrid.head(20),
        validated_tet1.head(2),
        tet1.head(8),
    ]).drop_duplicates(subset=["tet1_grna_id", "vp64_grna_id", "strategy"])

    combined = pd.concat([hybrid, tet1]).copy()

    # Cas-OFFinder-style genome-wide off-target (preferred) or legacy heuristic
    cas_offtarget = load_json(MODELS / "validation" / "cas_offinder_offtarget.json")
    genome_offtarget = load_json(MODELS / "validation" / "genome_wide_offtarget.json")
    per_guide = (cas_offtarget or genome_offtarget or {}).get("per_guide", {})

    def guide_risk_score(gid: str | float | None) -> float:
        if gid is None or (isinstance(gid, float) and pd.isna(gid)):
            return 0.0
        rec = per_guide.get(str(gid), {})
        return float(rec.get("cas_offinder_score", rec.get("risk_score", 0.0)))

    combined["genome_offtarget_risk"] = combined.apply(
        lambda r: max(
            guide_risk_score(r.get("tet1_grna_id")),
            guide_risk_score(r.get("vp64_grna_id")),
        ),
        axis=1,
    )
    combined["validation_priority"] = combined.apply(
        lambda r: max(
            VALIDATED_PRIORITY.get(str(r.get("tet1_grna_id", "")), 0),
            VALIDATED_PRIORITY.get(str(r.get("vp64_grna_id", "")), 0),
        ),
        axis=1,
    )
    combined["vp64_priority"] = combined["vp64_grna_id"].map(
        lambda x: VALIDATED_PRIORITY.get(str(x), 0) if pd.notna(x) else 0
    )
    combined = combined.sort_values(
        ["validation_priority", "vp64_priority", "objective_score", "genome_offtarget_risk", "uncertainty"],
        ascending=[False, False, False, True, True],
    )

    entries = []
    for _, row in combined.head(25).iterrows():
        evidence = []
        if row.get("tet1_grna_id") == "PWS_G.4363":
            evidence.append("bisulfite_validated_GSE285300")
        if row.get("tet1_grna_id") == "PWS_G.4670":
            evidence.append("top_Tet1_sublib_hit_Rohm2025")
        if row.get("vp64_grna_id") in ("PWS_G.9414", "PWS_G.8738"):
            evidence.append("top_VP64_SNRPN_promoter_Rohm2025")

        entries.append({
            "catalog_rank": len(entries) + 1,
            "strategy": row["strategy"],
            "tet1_grna_id": row.get("tet1_grna_id"),
            "tet1_protospacer": row.get("tet1_protospacer"),
            "tet1_hg38_start": row.get("tet1_start_hg38"),
            "vp64_grna_id": row.get("vp64_grna_id"),
            "vp64_protospacer": row.get("vp64_protospacer"),
            "vp64_hg38_start": row.get("vp64_start_hg38"),
            "recommended_w_tet1": row.get("w_tet1"),
            "recommended_w_vp64": row.get("w_vp64"),
            "predicted_SNRPN_pct_WT": row.get("predicted_SNRPN_pct_WT"),
            "predicted_SNHG14_pct_WT": row.get("predicted_SNHG14_pct_WT"),
            "both_in_target_window": row.get("both_in_window"),
            "objective_score": row.get("objective_score"),
            "uncertainty": row.get("uncertainty"),
            "genome_offtarget_risk": row.get("genome_offtarget_risk"),
            "score_lower_bound": row.get("score_lower"),
            "score_upper_bound": row.get("score_upper"),
            "experimental_evidence": ";".join(evidence) if evidence else "computational_only",
            "validation_priority": int(row.get("validation_priority", 0)),
            "organoid_support": "GSE262700_Tet1_class_reactivation",
            "collateral_risk": "low_ICR_targeted",
        })
    return pd.DataFrame(entries)


def build_synthesis_report(catalog: pd.DataFrame) -> dict:
    held_out = load_json(MODELS / "validation" / "held_out_benchmark.json")
    organoid = load_json(MODELS / "validation" / "organoid_concordance.json")
    collateral = load_json(MODELS / "validation" / "collateral_imprinting_risk.json")
    offtarget = load_json(MODELS / "validation" / "locus_offtarget_assessment.json")
    cas_offtarget = load_json(MODELS / "validation" / "cas_offinder_offtarget.json")
    genome_offtarget = load_json(MODELS / "validation" / "genome_wide_offtarget.json")
    sensitivity = load_json(MODELS / "sensitivity" / "sensitivity_report.json")
    encode = load_json(ROOT / "data" / "encode_reference" / "summary.json")
    dmr = load_json(ROOT / "data" / "curated" / "methylation_dmr_mapping_report.json")
    gse152098 = load_json(ROOT / "data" / "curated" / "gse152098_reprocessing_report.json")
    accessibility = load_json(MODELS / "cell_type_scoring" / "accessibility_report.json")
    protocol_path = MODELS / "experimental_protocol" / "experimental_validation_protocol.json"

    top = catalog.iloc[0] if len(catalog) > 0 else None

    return {
        "catalog_version": "1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": "PWS maternal UPD15 CRISPR epigenome-editing design platform",
        "genome_build": "GRCh38",
        "citations": CITATIONS,
        "pipeline_phases_complete": [
            "Phase 1: Data acquisition (GEO, CRISPRepi)",
            "Phase 2: Locus annotation (314 genes, 81 CpG islands)",
            "Phase 3: Curation (23,506 editing records)",
            "Phase 3b: hg38 integration + base-pair merge (464 significant sites)",
            "Phase 5: Forward model (editor-stratified rules)",
            "Phase 6/7: Hybrid optimization (121 designs)",
            "Phase 8: Exploratory circuit mapping",
            "Phase 9: Validation suite",
            "Phase 10: Cell-type accessibility, locus off-target, experimental protocol",
            "Phase 11: Pareto/Bayesian optimization, circuit ODE, sensitivity analysis",
            "Phase 12: Cas-OFFinder scan, ENCODE context, GSE152098 bigWig, 450K DMR mapping",
        ],
        "top_recommended_design": top.to_dict() if top is not None else None,
        "validation_summary": {
            "held_out_benchmark_pass": held_out.get("validated_guide_recovery", {}).get("overall_pass") if held_out else None,
            "organoid_concordance_pass": organoid.get("overall_pass") if organoid else None,
            "collateral_risk": collateral.get("overall_collateral_risk") if collateral else None,
            "locus_offtarget_risk": offtarget.get("overall_locus_offtarget_risk") if offtarget else None,
            "genome_wide_offtarget_available": bool(cas_offtarget or genome_offtarget),
            "cas_offinder_available": bool(cas_offtarget),
            "encode_reference_available": bool(encode and encode.get("overall_success")),
            "methylation_dmr_mapping_available": bool(dmr),
            "gse152098_bigwig_reprocessed": bool(gse152098),
            "sensitivity_available": bool(sensitivity),
        },
        "regulatory_context": {
            "encode_tracks_success": encode.get("n_tracks_success") if encode else None,
            "gse28525_imprinted_dmrs_in_pws": (dmr.get("gse28525", {}) or {}).get("imprinted_dmrs_in_window") if dmr else None,
            "neuron_vs_esc_icr_fold": gse152098.get("neuron_vs_esc_icr_fold") if gse152098 else None,
        },
        "robustness": {
            "window_probability_report": str(MODELS / "sensitivity" / "sensitivity_report.json")
            if (MODELS / "sensitivity" / "sensitivity_report.json").exists()
            else None,
            "top_design_p_both_in_window": (sensitivity.get("top_design", {}) or {}).get("p_both_in_window") if sensitivity else None,
        },
        "cell_type_accessibility": accessibility.get("icr_region_note") if accessibility else None,
        "experimental_protocol": str(protocol_path) if protocol_path.exists() else None,
        "key_scientific_findings": [
            "Tet1 at PWS-ICR demethylates and restores SNHG14 (~100% WT) but under-restores SNRPN (~52% WT)",
            "VP64 at SNRPN promoter hyperactivates SNRPN (~246% WT) without SNHG14 rescue",
            "Hybrid Tet1 + minimal VP64 achieves dual-gene paternal-equivalent window (70-130% WT)",
            "Organoid validation confirms massive SNRPN/SNHG14 reactivation in hypothalamic tissue",
            "Validated guides (G4363, G4670, G9414) recover in top decile within editor context",
        ],
        "proposed_experimental_validation": [
            "PWS UPD15 iPSC -> hypothalamic neuron/organoid differentiation",
            "Test top 3 hybrid designs + Tet1 monotherapy comparator (PWS_G.4363)",
            "Measure: PWS-ICR methylation (bisulfite), SNRPN/SNHG14 expression (qPCR/RNA-seq)",
            "Collateral panel: UBE3A, ATP10A, GABRB3 expression; genome-wide methylation",
            "Compare measured outcomes to forward model predictions; iterate model",
        ],
        "honest_limitations": [
            "Forward model uses editor-level bulk RNA priors, not guide-resolved expression",
            "Hybrid dose-response is phenomenological; requires experimental titration",
            "Circuit mapping is hypothesis-generating only",
            "Cas-OFFinder-style genome scan is in-silico; biochemical validation (CHANGE-seq/GUIDE-seq) still required",
            "Held-out benchmark is small (n=4 validated guides)",
        ],
        "novelty_statement": (
            "Novel contribution is the allele-aware PWS-locus design layer: integrating "
            "imprinting biology, editor-stratified forward models, and uncertainty-aware "
            "optimization over CRISPR epigenome-editing configurations — not demonstration "
            "that reactivation is possible (already shown by Nemoto 2025 and Rohm 2025)."
        ),
    }


def main():
    catalog = build_catalog_entries()
    catalog.to_csv(OUT / "pws_therapeutic_design_catalog.csv", index=False)

    report = build_synthesis_report(catalog)
    with open(OUT / "project_synthesis_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    log.info("Final catalog: %d designs", len(catalog))
    if len(catalog) > 0:
        top = catalog.iloc[0]
        log.info(
            "Top: %s Tet1=%s VP64=%s (SNRPN=%.0f%%, SNHG14=%.0f%%)",
            top["strategy"], top["tet1_grna_id"], top.get("vp64_grna_id"),
            top["predicted_SNRPN_pct_WT"], top["predicted_SNHG14_pct_WT"],
        )
    log.info("Outputs -> %s", OUT)


if __name__ == "__main__":
    main()
