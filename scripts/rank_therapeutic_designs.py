"""
Phase 4/7 bridge: Rank therapeutic editor/gRNA designs using curated data.

Composite score (maximize):
  - Predicted SNRPN/SNHG14 reactivation toward 70-130% WT window
  - Element-editor compatibility (Tet1 -> demethylation sites; VP64 -> activation sites)
  - Subregion relevance (pws_critical > outside)
  - Statistical confidence (-log10 padj)

Outputs uncertainty-aware ranked designs, NOT point-optimal predictions.

Output: data/curated/ranked_therapeutic_designs.csv
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "locus.yaml"
CURATED = ROOT / "data" / "curated"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

# Empirical editor-element compatibility from Rohm 2025
# Tet1 demethylates at PWS-ICR (~25.2M hg19); VP64 activates at SNRPN promoter (~25.1M hg19)
EDITOR_ELEMENT_PRIORS = {
    "dCas9-Tet1": {"preferred_start_range": (25190000, 25220000), "mechanism": "demethylation"},
    "dCas9-VP64": {"preferred_start_range": (25070000, 25115000), "mechanism": "transcriptional_activation"},
    "dCas9-KRAB": {"preferred_start_range": (25190000, 25230000), "mechanism": "repression_mapping"},
}


def load_target_state() -> dict:
    path = CURATED / "therapeutic_target_state.json"
    if path.exists():
        return json.loads(path.read_text())
    return {}


def element_compatibility_score(editor: str, start: int) -> float:
    prior = EDITOR_ELEMENT_PRIORS.get(editor)
    if not prior:
        return 0.5
    lo, hi = prior["preferred_start_range"]
    if lo <= start <= hi:
        return 1.0
    dist = min(abs(start - lo), abs(start - hi))
    return max(0.0, 1.0 - dist / 500_000)


def reactivation_potential(editor: str, target_state: dict) -> float:
    """Empirical reactivation strength from bulk RNA-seq (normalized 0-1)."""
    outcomes = target_state.get("editor_outcomes", {})
    snrpn_key = f"{editor}:SNRPN"
    snhg_key = f"{editor}:SNHG14"
    scores = []
    genes = target_state.get("genes", {})
    window = target_state.get("expression_window_pct", [70, 130])

    for key, gene in [(snrpn_key, "SNRPN"), (snhg_key, "SNHG14")]:
        if key not in outcomes or gene not in genes:
            continue
        pct = outcomes[key].get("pct_of_wt")
        if pct is None:
            continue
        # Score peaks at 100% WT, penalizes under and over
        if window[0] <= pct <= window[1]:
            scores.append(1.0)
        elif pct < window[0]:
            scores.append(pct / window[0])
        else:
            scores.append(max(0.0, 1.0 - (pct - window[1]) / window[1]))
    return float(np.mean(scores)) if scores else 0.3


def rank_designs(screens: pd.DataFrame, target_state: dict) -> pd.DataFrame:
    sig = screens[screens["significant"]].copy()
    sig = sig[sig["editor_system"].isin(["dCas9-Tet1", "dCas9-VP64"])]

    # Deduplicate: best padj per (editor, start)
    sig = sig.sort_values("padj").groupby(["editor_system", "start"], as_index=False).first()

    records = []
    for editor in sig["editor_system"].unique():
        emp_reactivation = reactivation_potential(editor, target_state)
        sub = sig[sig["editor_system"] == editor]
        for _, r in sub.iterrows():
            elem_score = element_compatibility_score(editor, int(r["start"]))
            subregion_bonus = 1.0 if r["subregion"] == "pws_critical" else 0.5
            confidence = min(float(r["neg_log10_padj"]), 50) / 50
            composite = (
                0.30 * emp_reactivation
                + 0.25 * elem_score
                + 0.25 * confidence
                + 0.20 * subregion_bonus
            )
            records.append({
                "grna_id": r["grna_id"],
                "protospacer": r["protospacer"],
                "editor_system": editor,
                "chrom": r["chrom"],
                "start": int(r["start"]),
                "end": int(r["end"]),
                "subregion": r["subregion"],
                "genome_build": r["genome_build"],
                "padj": r["padj"],
                "mechanism": EDITOR_ELEMENT_PRIORS.get(editor, {}).get("mechanism", "unknown"),
                "element_compatibility": round(elem_score, 3),
                "empirical_reactivation_score": round(emp_reactivation, 3),
                "confidence_score": round(confidence, 3),
                "composite_score": round(composite, 4),
                "uncertainty_note": "Scores are hypothesis-generating; based on bulk RNA-seq and screen statistics",
            })

    ranked = pd.DataFrame(records).sort_values("composite_score", ascending=False)
    return ranked


def main():
    screens = pd.read_parquet(CURATED / "editing_screens.parquet")
    target_state = load_target_state()
    ranked = rank_designs(screens, target_state)

    out = CURATED / "ranked_therapeutic_designs.csv"
    ranked.to_csv(out, index=False)

    top = ranked.head(15)
    log.info("Top 15 therapeutic designs:")
    for _, r in top.iterrows():
        log.info(
            "  %s %s @ %d (score=%.3f, elem=%.2f, padj=%.2e)",
            r["editor_system"], r["grna_id"], r["start"],
            r["composite_score"], r["element_compatibility"], r["padj"],
        )
    log.info("Wrote %d ranked designs -> %s", len(ranked), out)


if __name__ == "__main__":
    main()
