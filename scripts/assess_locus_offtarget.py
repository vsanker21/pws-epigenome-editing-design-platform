"""
Phase 10b: Locus-level off-target and protospacer uniqueness assessment.

This is NOT genome-wide CRISPR off-target prediction (would require Cas-OFFinder,
BLAST, or CRISPRepi off-target modules). Instead we assess:

  1. Protospacer uniqueness within the PWS critical region (chr15:24.8-32.7M)
  2. Duplicate protospacer sites in the editing screen library
  3. Near-duplicate seeds (12 bp PAM-proximal) within the locus

Scientific rationale: collateral imprinting risk at neighboring sites within
15q11-q13 is the primary local off-target concern for PWS therapy.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INTEGRATED = ROOT / "data" / "integrated"
MODELS = ROOT / "data" / "models"
OUT = MODELS / "validation"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

SEED_LEN = 12


def hamming(s1: str, s2: str) -> int:
    return sum(a != b for a, b in zip(s1, s2))


def assess_protospacer(seq: str, all_screens: pd.DataFrame, guide_id: str) -> dict:
    seq = seq.upper()
    seed = seq[-SEED_LEN:]

    exact_matches = all_screens[all_screens["protospacer"].str.upper() == seq]
    other_sites = exact_matches[exact_matches["grna_id"] != guide_id]

    seed_matches = []
    for _, r in all_screens.iterrows():
        other = str(r["protospacer"]).upper()
        if len(other) >= SEED_LEN and other != seq:
            if hamming(seed, other[-SEED_LEN:]) <= 1:
                seed_matches.append({
                    "grna_id": r["grna_id"],
                    "start_hg38": int(r["start_hg38"]),
                    "editor": r["editor_system"],
                    "hamming": hamming(seed, other[-SEED_LEN:]),
                })

    n_exact = len(other_sites)
    n_seed = len(seed_matches)
    risk = "low"
    if n_exact > 0:
        risk = "moderate"
    if n_exact > 2 or n_seed > 5:
        risk = "elevated"

    return {
        "protospacer": seq,
        "exact_duplicates_in_locus": n_exact,
        "near_seed_matches_1mm": n_seed,
        "offtarget_locus_risk": risk,
        "duplicate_sites": other_sites[["grna_id", "start_hg38", "editor_system"]].head(5).to_dict(orient="records"),
        "seed_near_matches": seed_matches[:5],
    }


def main():
    merged = pd.read_parquet(INTEGRATED / "locus_merged_hg38.parquet")
    catalog = pd.read_csv(MODELS / "final_catalog" / "pws_therapeutic_design_catalog.csv")

    top = catalog[catalog["catalog_rank"] == 1].iloc[0]
    guides_to_assess = []
    if pd.notna(top.get("tet1_grna_id")):
        guides_to_assess.append(("tet1", str(top["tet1_grna_id"]), str(top["tet1_protospacer"])))
    if pd.notna(top.get("vp64_grna_id")):
        guides_to_assess.append(("vp64", str(top["vp64_grna_id"]), str(top["vp64_protospacer"])))

    assessments = {}
    for role, gid, proto in guides_to_assess:
        assessments[f"{role}_{gid}"] = assess_protospacer(proto, merged, gid)

    report = {
        "assessment_type": "locus_level_protospacer_uniqueness",
        "scope": "PWS critical region significant screen sites only (not genome-wide)",
        "top_recommended_design": {
            "tet1_grna_id": top.get("tet1_grna_id"),
            "vp64_grna_id": top.get("vp64_grna_id"),
        },
        "guide_assessments": assessments,
        "limitations": [
            "Does not predict genome-wide Cas9 binding or cutting",
            "Limited to sites tiled in Rohm 2025 screen library",
            "Recommend Cas-OFFinder or CHANGE-seq for preclinical off-target validation",
        ],
        "overall_locus_offtarget_risk": max(
            (a["offtarget_locus_risk"] for a in assessments.values()),
            key=lambda x: ["low", "moderate", "elevated"].index(x),
        ) if assessments else "unknown",
    }

    with open(OUT / "locus_offtarget_assessment.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    log.info("Top design locus off-target risk: %s", report["overall_locus_offtarget_risk"])
    for key, a in assessments.items():
        log.info("  %s: %s (%d exact dup, %d seed near-matches)",
                 key, a["offtarget_locus_risk"], a["exact_duplicates_in_locus"], a["near_seed_matches_1mm"])
    log.info("Report -> %s", OUT / "locus_offtarget_assessment.json")


if __name__ == "__main__":
    main()
