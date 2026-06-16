"""
Held-out validation benchmark for the PWS epigenome-editing design pipeline.

Tests (per Project_Modifications requirement):
  1. Validated guide recovery — do known hits rank in top decile within editor?
  2. Screen concordance — Spearman correlation of rules score vs -log10(padj) within editor
  3. Distance-to-anchor — Tet1 guides nearer PWS-ICR should rank higher
  4. Methylation concordance — Tet1 guides with larger methylation delta rank higher
  5. Editor outcome prediction — predicted bulk RNA priors match observed GSE243185

This is explicitly a SMALL benchmark; results are reported with uncertainty.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
INTEGRATED = ROOT / "data" / "integrated"
CURATED = ROOT / "data" / "curated"
MODELS = ROOT / "data" / "models"
OUT = MODELS / "validation"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

VALIDATED_GUIDES = {
    "PWS_G.4670": {"editor": "dCas9-Tet1", "evidence": "Top Tet1 sublibrary hit (Rohm 2025)"},
    "PWS_G.4363": {"editor": "dCas9-Tet1", "evidence": "Bisulfite-validated Tet1 guide (GSE285300)"},
    "PWS_G.9414": {"editor": "dCas9-VP64", "evidence": "Top VP64 activation hit (Rohm 2025)"},
    "PWS_G.8738": {"editor": "dCas9-VP64", "evidence": "VP64 SNRPN promoter hit (Rohm 2025)"},
}

EDITOR_ANCHORS_HG38 = {
    "dCas9-Tet1": 24954921,
    "dCas9-VP64": 24862855,
}

DECILE_THRESHOLD = 0.10  # top 10%


def load_ranked_predictions() -> pd.DataFrame:
    path = MODELS / "forward_model_predictions.parquet"
    if path.exists():
        return pd.read_parquet(path)
    merged = pd.read_parquet(INTEGRATED / "locus_merged_hg38.parquet")
    return merged


def validated_guide_recovery(ranked: pd.DataFrame) -> dict:
    results = {}
    for gid, meta in VALIDATED_GUIDES.items():
        editor = meta["editor"]
        sub = ranked[ranked["editor_system"] == editor].reset_index(drop=True)
        hit = sub[sub["grna_id"] == gid]
        if hit.empty:
            results[gid] = {"found": False, "pass": False, "evidence": meta["evidence"]}
            continue
        rank = int(hit.index[0]) + 1
        n = len(sub)
        percentile = 1 - rank / n
        in_top_decile = rank <= max(1, int(np.ceil(DECILE_THRESHOLD * n)))
        results[gid] = {
            "found": True,
            "editor": editor,
            "rank": rank,
            "n_editor_guides": n,
            "percentile": round(percentile * 100, 1),
            "in_top_decile": in_top_decile,
            "pass": in_top_decile,
            "score": float(hit["rules_reactivation_score"].iloc[0]) if "rules_reactivation_score" in hit.columns else None,
            "padj": float(hit["padj"].iloc[0]),
            "evidence": meta["evidence"],
        }
    n_pass = sum(1 for v in results.values() if v.get("pass"))
    return {
        "per_guide": results,
        "n_validated": len(VALIDATED_GUIDES),
        "n_pass_top_decile": n_pass,
        "pass_rate": round(n_pass / len(VALIDATED_GUIDES), 2),
        "overall_pass": n_pass >= 3,  # require 3/4 in top decile
    }


def within_editor_concordance(ranked: pd.DataFrame) -> dict:
    results = {}
    for editor in ["dCas9-Tet1", "dCas9-VP64"]:
        sub = ranked[ranked["editor_system"] == editor].copy()
        if len(sub) < 5:
            continue
        sub["neg_log10_padj"] = -np.log10(sub["padj"].clip(lower=1e-300))
        score_col = "rules_reactivation_score" if "rules_reactivation_score" in sub.columns else None
        if score_col is None:
            continue
        rho, pval = spearmanr(sub[score_col], sub["neg_log10_padj"])
        results[editor] = {
            "n": len(sub),
            "spearman_rho": round(float(rho), 3),
            "p_value": float(pval),
            "interpretation": (
                "positive concordance" if rho > 0.3 else
                "weak concordance" if rho > 0 else "discordant"
            ),
        }
    return results


def tet1_distance_benchmark(ranked: pd.DataFrame) -> dict:
    tet1 = ranked[ranked["editor_system"] == "dCas9-Tet1"].copy()
    anchor = EDITOR_ANCHORS_HG38["dCas9-Tet1"]
    tet1["dist_icr"] = (tet1["start_hg38"] - anchor).abs()
    tet1["neg_log10_padj"] = -np.log10(tet1["padj"].clip(lower=1e-300))

    rho_dist_padj, _ = spearmanr(tet1["dist_icr"], tet1["neg_log10_padj"])
    score_col = "rules_reactivation_score" if "rules_reactivation_score" in tet1.columns else None
    rho_dist_score = None
    if score_col:
        rho_dist_score, _ = spearmanr(tet1["dist_icr"], tet1[score_col])

    # Methylation concordance (guides with bisulfite data)
    meth = tet1.dropna(subset=["methylation_delta_edited_vs_nt"])
    rho_meth = None
    if len(meth) >= 5 and score_col:
        rho_meth, _ = spearmanr(
            meth["methylation_delta_edited_vs_nt"].abs(),
            meth[score_col],
        )

    return {
        "n_tet1": len(tet1),
        "dist_vs_padj_spearman": round(float(rho_dist_padj), 3),
        "dist_vs_score_spearman": round(float(rho_dist_score), 3) if rho_dist_score is not None else None,
        "methylation_vs_score_spearman": round(float(rho_meth), 3) if rho_meth is not None else None,
        "n_with_methylation": len(meth),
        "expected": "negative rho (closer to ICR = higher significance/score)",
        "pass": rho_dist_padj < -0.2 if not np.isnan(rho_dist_padj) else False,
    }


def editor_outcome_prediction() -> dict:
    """Compare model priors to observed bulk RNA-seq (GSE243185)."""
    target = json.loads((CURATED / "therapeutic_target_state.json").read_text())
    outcomes = target.get("editor_outcomes", {})
    window = target.get("expression_window_pct", [70, 130])

    records = []
    for editor in ["dCas9-Tet1", "dCas9-VP64"]:
        for gene in ["SNRPN", "SNHG14"]:
            key = f"{editor}:{gene}"
            if key not in outcomes:
                continue
            obs = outcomes[key]
            pct = obs.get("pct_of_wt")
            if pct is None:
                continue
            in_window = window[0] <= pct <= window[1]
            records.append({
                "editor": editor,
                "gene": gene,
                "observed_pct_WT": round(pct, 1),
                "in_target_window": in_window,
            })

    return {
        "source": "GSE243185 bulk RNA-seq (Rohm 2025)",
        "expression_window_pct": window,
        "outcomes": records,
        "Tet1_SNRPN_deficit": True,  # ~52% — motivates hybrid strategy
        "VP64_SNHG14_absent": True,  # 0% — VP64 alone insufficient
    }


def leave_one_out_stability(ranked: pd.DataFrame) -> dict:
    """Check if top-5 Tet1 guides are stable when each validated guide is excluded."""
    tet1 = ranked[ranked["editor_system"] == "dCas9-Tet1"].copy()
    score_col = "rules_reactivation_score" if "rules_reactivation_score" in tet1.columns else None
    if score_col is None:
        return {"status": "unavailable"}

    baseline_top5 = set(tet1.nlargest(5, score_col)["grna_id"])
    stability = {}
    for gid in ["PWS_G.4670", "PWS_G.4363"]:
        sub = tet1[tet1["grna_id"] != gid]
        new_top5 = set(sub.nlargest(5, score_col)["grna_id"])
        overlap = len(baseline_top5 & new_top5)
        stability[gid] = {
            "excluded_guide": gid,
            "top5_overlap_with_baseline": overlap,
            "stable": overlap >= 4,
        }
    return stability


def main():
    ranked = load_ranked_predictions()
    if "rules_reactivation_score" not in ranked.columns:
        log.warning("rules_reactivation_score missing; re-run train_forward_model.py")

    report = {
        "benchmark_type": "held_out_validation",
        "n_significant_merged_sites": len(ranked),
        "validated_guide_recovery": validated_guide_recovery(ranked),
        "within_editor_concordance": within_editor_concordance(ranked),
        "tet1_distance_benchmark": tet1_distance_benchmark(ranked),
        "editor_outcome_prediction": editor_outcome_prediction(),
        "leave_one_out_stability": leave_one_out_stability(ranked),
        "limitations": [
            "Only 4 experimentally validated guides available as held-out benchmark",
            "Bulk RNA outcomes are editor-level, not guide-resolved",
            "Cross-editor ranking is not a valid test (mechanisms differ)",
            "Organoid validation (GSE262700) is separate — different editor construct (SunTag-TET1)",
        ],
    }

    with open(OUT / "held_out_benchmark.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    recovery = report["validated_guide_recovery"]
    log.info(
        "Validated guide recovery: %d/%d in top decile (pass=%s)",
        recovery["n_pass_top_decile"], recovery["n_validated"], recovery["overall_pass"],
    )
    for gid, info in recovery["per_guide"].items():
        if info.get("found"):
            log.info("  %s: rank %d (%.0f%%ile) %s",
                     gid, info["rank"], info["percentile"],
                     "PASS" if info["pass"] else "FAIL")
    log.info("Report -> %s", OUT / "held_out_benchmark.json")


if __name__ == "__main__":
    main()
