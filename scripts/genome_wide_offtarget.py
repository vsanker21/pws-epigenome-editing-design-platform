"""
Phase 10c: Genome-wide CRISPR off-target *potential* scan (seed+PAM heuristic).

Goal
----
Provide a reproducible, offline-capable approximation of genome-wide off-target
liability for the catalog's protospacers by scanning reference genome FASTA.

This is intentionally conservative and does NOT replace experimental assays
(CHANGE-seq, GUIDE-seq) or high-fidelity aligners; it is a computational
screening layer to help rank designs by safety.

Method (SpCas9 NGG)
-------------------
We scan both strands for canonical SpCas9 PAMs:
  - Forward strand:  N{20} + NGG
  - Reverse strand:  CCN + N{20}  (equivalently NGG on reverse)

To keep runtime tractable, we:
  1) Load the final catalog and extract unique protospacers.
  2) For each guide, compute the PAM-proximal seed (last 12 nt).
  3) Generate all seed variants within <=1 mismatch (Hamming distance).
  4) Scan chromosomes and only evaluate sites whose seed is in the variant map.
  5) For each candidate site, compute full 20nt mismatches vs each matching guide
     and count hits at mismatch <= MAX_MISMATCHES.

Inputs
------
 - data/models/final_catalog/pws_therapeutic_design_catalog.csv
 - (downloaded if missing) data/genome/hg38/chr*.fa.gz  (UCSC hgdownload)

Output
------
 - data/models/validation/genome_wide_offtarget.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
GENOME = DATA / "genome" / "hg38"
MODELS = DATA / "models"
VALIDATION = MODELS / "validation"
VALIDATION.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

# SpCas9 defaults
PAM_FWD_1 = "G"
PAM_FWD_2 = "G"  # ... NGG
PAM_REV_0 = "C"
PAM_REV_1 = "C"  # ... CCN (reverse-strand PAM equivalent)

SEED_LEN = 12
GUIDE_LEN = 20
MAX_SEED_MM = 1
MAX_MISMATCHES = 3

UCSC_CHR_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz"

DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
PIPELINE_DEFAULT_CHROMS = ["chr15", "chr14", "chr16"]


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return 10**9
    return sum(x != y for x, y in zip(a, b))


def download_chr_if_missing(chrom: str) -> Path:
    GENOME.mkdir(parents=True, exist_ok=True)
    out = GENOME / f"{chrom}.fa.gz"
    if out.exists():
        return out

    url = UCSC_CHR_URL.format(chrom=chrom)
    log.info("Downloading %s -> %s", url, out)
    with requests.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(out, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    return out


def iter_fasta_sequence(path_gz: Path) -> Iterable[str]:
    """Yield concatenated FASTA sequence lines (uppercase), header stripped."""
    with gzip.open(path_gz, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                continue
            yield line.strip().upper()


def load_chrom_seq(path_gz: Path) -> str:
    # Join in-memory; ok for a chromosome at a time on typical RAM.
    return "".join(iter_fasta_sequence(path_gz))


def seed_variants(seed: str, max_mm: int = 1) -> set[str]:
    seed = seed.upper()
    alphabet = ["A", "C", "G", "T"]
    out = {seed}
    if max_mm <= 0:
        return out
    # Only implement <=1 mismatch (kept explicit for performance/reproducibility).
    for i, ch in enumerate(seed):
        for a in alphabet:
            if a != ch:
                out.add(seed[:i] + a + seed[i + 1 :])
    return out


@dataclass(frozen=True)
class Guide:
    guide_id: str
    seq: str  # 20nt protospacer

    @property
    def seed(self) -> str:
        return self.seq[-SEED_LEN:]


def load_guides_from_catalog() -> list[Guide]:
    catalog_path = MODELS / "final_catalog" / "pws_therapeutic_design_catalog.csv"
    df = pd.read_csv(catalog_path)

    guides: dict[str, str] = {}
    for _, r in df.iterrows():
        for role in ("tet1", "vp64"):
            gid = r.get(f"{role}_grna_id")
            seq = r.get(f"{role}_protospacer")
            if pd.notna(gid) and pd.notna(seq):
                seq = str(seq).upper()
                if len(seq) == GUIDE_LEN:
                    guides[str(gid)] = seq
    out = [Guide(guide_id=k, seq=v) for k, v in sorted(guides.items())]
    if not out:
        raise RuntimeError("No guides found in final catalog.")
    return out


def build_seed_to_guides(guides: list[Guide]) -> dict[str, list[Guide]]:
    seed_to_guides: dict[str, list[Guide]] = defaultdict(list)
    for g in guides:
        for s in seed_variants(g.seed, MAX_SEED_MM):
            seed_to_guides[s].append(g)
    return seed_to_guides


def scan_chromosome_for_guides(chrom: str, seq: str, seed_to_guides: dict[str, list[Guide]]) -> dict[str, dict[str, int]]:
    """
    Returns counts per guide:
      {guide_id: {"mm0": n, "mm1": n, "mm2": n, "mm3": n}}
    """
    counts: dict[str, dict[str, int]] = defaultdict(lambda: {"mm0": 0, "mm1": 0, "mm2": 0, "mm3": 0})
    n = len(seq)

    # Forward-strand PAM: N{20} + NGG
    for i in range(0, n - (GUIDE_LEN + 3) + 1):
        pam2 = seq[i + GUIDE_LEN + 1 : i + GUIDE_LEN + 3]
        if pam2 != (PAM_FWD_1 + PAM_FWD_2):
            continue
        proto = seq[i : i + GUIDE_LEN]
        if "N" in proto:
            continue
        seed = proto[-SEED_LEN:]
        if seed not in seed_to_guides:
            continue
        for g in seed_to_guides[seed]:
            mm = hamming(proto, g.seq)
            if mm <= MAX_MISMATCHES:
                counts[g.guide_id][f"mm{mm}"] += 1

    # Reverse-strand PAM (as seen on forward sequence): CCN + N{20}
    for i in range(0, n - (3 + GUIDE_LEN) + 1):
        if seq[i] != PAM_REV_0 or seq[i + 1] != PAM_REV_1:
            continue
        proto_fwd = seq[i + 3 : i + 3 + GUIDE_LEN]
        if "N" in proto_fwd:
            continue
        proto = revcomp(proto_fwd)
        seed = proto[-SEED_LEN:]
        if seed not in seed_to_guides:
            continue
        for g in seed_to_guides[seed]:
            mm = hamming(proto, g.seq)
            if mm <= MAX_MISMATCHES:
                counts[g.guide_id][f"mm{mm}"] += 1

    return counts


def summarize_risk(per_guide_counts: dict[str, dict[str, int]]) -> dict[str, dict]:
    """
    Convert mismatch counts into a coarse risk label + numeric risk score.
    Risk score weights are heuristic:
      mm0: 1.0, mm1: 0.3, mm2: 0.1, mm3: 0.03
    """
    weights = {"mm0": 1.0, "mm1": 0.3, "mm2": 0.1, "mm3": 0.03}
    out = {}
    for gid, c in per_guide_counts.items():
        score = sum(weights[k] * float(c.get(k, 0)) for k in weights)
        # Excluding on-target is not possible here without guide coordinates; we label as "potential".
        risk = "low"
        if c.get("mm0", 0) >= 5 or score >= 10:
            risk = "elevated"
        elif c.get("mm0", 0) >= 2 or score >= 3:
            risk = "moderate"
        out[gid] = {
            "counts": c,
            "risk_score": score,
            "risk_label": risk,
        }
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Genome-wide off-target potential scan (seed+PAM heuristic)."
    )
    parser.add_argument(
        "--chromosomes",
        default=",".join(PIPELINE_DEFAULT_CHROMS),
        help=(
            "Comma-separated chromosomes to scan (e.g. chr15,chr14) or 'all'. "
            "Default is a PWS-relevant subset for pipeline speed."
        ),
    )
    args = parser.parse_args()

    guides = load_guides_from_catalog()
    seed_to_guides = build_seed_to_guides(guides)

    if str(args.chromosomes).strip().lower() == "all":
        chromosomes = list(DEFAULT_CHROMS)
    else:
        chromosomes = [c.strip() for c in str(args.chromosomes).split(",") if c.strip()]
        if not chromosomes:
            chromosomes = list(PIPELINE_DEFAULT_CHROMS)

    log.info("Guides: %d unique protospacers", len(guides))
    log.info("Chromosomes to scan: %s", ",".join(chromosomes))

    per_guide_total: dict[str, dict[str, int]] = defaultdict(lambda: {"mm0": 0, "mm1": 0, "mm2": 0, "mm3": 0})
    chrom_stats = {}

    for chrom in chromosomes:
        try:
            fa = download_chr_if_missing(chrom)
            seq = load_chrom_seq(fa)
        except Exception as e:
            log.warning("Skipping %s (genome download/read failed): %s", chrom, e)
            chrom_stats[chrom] = {"scanned": False, "error": str(e)}
            continue

        log.info("Scanning %s (len=%d)", chrom, len(seq))
        chrom_counts = scan_chromosome_for_guides(chrom, seq, seed_to_guides)
        chrom_stats[chrom] = {"scanned": True, "length": len(seq)}

        for gid, c in chrom_counts.items():
            for k in ("mm0", "mm1", "mm2", "mm3"):
                per_guide_total[gid][k] += int(c.get(k, 0))

    per_guide_summary = summarize_risk(per_guide_total)
    report = {
        "assessment_type": "genome_wide_offtarget_seed_scan",
        "reference": {
            "genome": "hg38",
            "source": "UCSC hgdownload (chromosome FASTA)",
            "local_cache": str(GENOME),
            "chromosomes_attempted": chromosomes,
            "chromosome_stats": chrom_stats,
        },
        "parameters": {
            "pam": "SpCas9 NGG",
            "guide_len": GUIDE_LEN,
            "seed_len": SEED_LEN,
            "max_seed_mismatches": MAX_SEED_MM,
            "max_total_mismatches": MAX_MISMATCHES,
        },
        "per_guide": per_guide_summary,
        "limitations": [
            "Counts potential binding sites; does not model chromatin accessibility or editor-specific binding kinetics",
            "Does not subtract known on-target locus coordinates (reports include on-target-like matches)",
            "Does not handle bulges/indels; mismatches only",
            "Default scan is limited to a subset of chromosomes unless additional chr*.fa.gz are cached",
        ],
    }

    out_path = VALIDATION / "genome_wide_offtarget.json"
    out_path.write_text(json.dumps(report, indent=2))
    log.info("Report -> %s", out_path)


if __name__ == "__main__":
    main()

