"""
Phase 12a: Cas-OFFinder-style genome-wide off-target scan with coordinates.

Implements the canonical SpCas9 NGG search used by Cas-OFFinder:
  - Enumerate all PAM sites on both strands
  - Score 20nt protospacer mismatches (0-4)
  - Report mismatch counts + top genomic hit coordinates per guide

This is a deterministic in-silico screen using hg38 reference FASTA.
It complements (and supersedes for ranking) the faster seed-only heuristic.
"""

from __future__ import annotations

import argparse
import json
import logging
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from offtarget_core import (
    DEFAULT_CHROMS,
    Guide,
    build_seed_to_guides,
    cas_offinder_score,
    download_chr_if_missing,
    load_chrom_seq,
    risk_label,
    scan_chromosome,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
GENOME = DATA / "genome" / "hg38"
MODELS = DATA / "models"
VALIDATION = MODELS / "validation"
VALIDATION.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)


def load_guides() -> list[Guide]:
    catalog = pd.read_csv(MODELS / "final_catalog" / "pws_therapeutic_design_catalog.csv")
    guides: dict[str, str] = {}
    for _, r in catalog.iterrows():
        for role in ("tet1", "vp64"):
            gid = r.get(f"{role}_grna_id")
            seq = r.get(f"{role}_protospacer")
            if pd.notna(gid) and pd.notna(seq) and len(str(seq)) == 20:
                guides[str(gid)] = str(seq).upper()
    if not guides:
        raise RuntimeError("No guides found in catalog.")
    return [Guide(k, v) for k, v in sorted(guides.items())]


def main():
    parser = argparse.ArgumentParser(description="Cas-OFFinder-style genome scan")
    parser.add_argument("--chromosomes", default="all", help="Comma list or 'all'")
    args = parser.parse_args()

    if str(args.chromosomes).lower() == "all":
        chromosomes = list(DEFAULT_CHROMS)
    else:
        chromosomes = [c.strip() for c in str(args.chromosomes).split(",") if c.strip()]

    guides = load_guides()
    seed_to_guides = build_seed_to_guides(guides, max_seed_mm=1)
    log.info("Scanning %d guides across %d chromosomes", len(guides), len(chromosomes))

    per_guide_counts: dict[str, dict[str, int]] = defaultdict(
        lambda: {f"mm{i}": 0 for i in range(5)}
    )
    per_guide_hits: dict[str, list] = defaultdict(list)
    chrom_stats = {}

    for chrom in chromosomes:
        try:
            fa = download_chr_if_missing(GENOME, chrom)
            seq = load_chrom_seq(fa)
        except Exception as e:
            log.warning("Skip %s: %s", chrom, e)
            chrom_stats[chrom] = {"scanned": False, "error": str(e)}
            continue

        log.info("Scanning %s (%d bp)", chrom, len(seq))
        counts, hits = scan_chromosome(
            chrom, seq, seed_to_guides, max_mismatches=4, collect_hits=True, max_hits_per_guide=40
        )
        chrom_stats[chrom] = {"scanned": True, "length": len(seq)}

        for gid, c in counts.items():
            for k, v in c.items():
                per_guide_counts[gid][k] += int(v)
        for gid, hlist in hits.items():
            per_guide_hits[gid].extend([asdict(h) for h in hlist])

    per_guide_summary = {}
    for gid, c in per_guide_counts.items():
        mm0 = int(c.get("mm0", 0))
        score = cas_offinder_score(c)
        hit_list = sorted(per_guide_hits.get(gid, []), key=lambda h: (h["mismatches"], h["chrom"], h["start"]))[:40]
        per_guide_summary[gid] = {
            "counts": c,
            "cas_offinder_score": round(score, 3),
            "risk_label": risk_label(score, mm0),
            "top_hits": hit_list,
        }

    report = {
        "assessment_type": "cas_offinder_style_genome_scan",
        "reference": {"genome": "hg38", "source": "UCSC hgdownload", "chromosomes": chromosomes, "chromosome_stats": chrom_stats},
        "parameters": {"pam": "SpCas9 NGG", "guide_len": 20, "max_mismatches": 4, "seed_prefilter_mm": 1},
        "per_guide": per_guide_summary,
        "limitations": [
            "Mismatch-only model; no RNA-DNA bulges",
            "Does not model chromatin accessibility",
            "mm0 hits include intended on-target sites",
        ],
    }

    out = VALIDATION / "cas_offinder_offtarget.json"
    out.write_text(json.dumps(report, indent=2))
    log.info("Cas-OFFinder report -> %s", out)

    # Also refresh genome_wide_offtarget.json for backward compatibility
    legacy = {
        "assessment_type": "genome_wide_offtarget_seed_scan",
        "reference": report["reference"],
        "parameters": {"pam": "SpCas9 NGG", "guide_len": 20, "seed_len": 12, "max_seed_mismatches": 1, "max_total_mismatches": 4},
        "per_guide": {
            gid: {
                "counts": rec["counts"],
                "risk_score": rec["cas_offinder_score"],
                "risk_label": rec["risk_label"],
            }
            for gid, rec in per_guide_summary.items()
        },
        "limitations": report["limitations"] + ["Derived from Cas-OFFinder-style scan"],
    }
    (VALIDATION / "genome_wide_offtarget.json").write_text(json.dumps(legacy, indent=2))


if __name__ == "__main__":
    main()
