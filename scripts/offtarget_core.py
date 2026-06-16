"""Shared SpCas9 off-target scanning utilities (Cas-OFFinder-compatible heuristic)."""

from __future__ import annotations

import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import requests

GUIDE_LEN = 20
SEED_LEN = 12
MAX_MISMATCHES = 4  # Cas-OFFinder default upper bound for screening

UCSC_CHR_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz"
DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return 10**9
    return sum(x != y for x, y in zip(a, b))


def seed_variants(seed: str, max_mm: int = 1) -> set[str]:
    seed = seed.upper()
    alphabet = ["A", "C", "G", "T"]
    out = {seed}
    if max_mm <= 0:
        return out
    for i, ch in enumerate(seed):
        for a in alphabet:
            if a != ch:
                out.add(seed[:i] + a + seed[i + 1 :])
    return out


@dataclass(frozen=True)
class Guide:
    guide_id: str
    seq: str

    @property
    def seed(self) -> str:
        return self.seq[-SEED_LEN:]


@dataclass
class OffTargetHit:
    chrom: str
    start: int
    end: int
    strand: str
    mismatches: int
    protospacer: str
    pam: str


def download_chr_if_missing(genome_dir: Path, chrom: str) -> Path:
    genome_dir.mkdir(parents=True, exist_ok=True)
    out = genome_dir / f"{chrom}.fa.gz"
    if out.exists():
        return out
    url = UCSC_CHR_URL.format(chrom=chrom)
    with requests.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(out, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    return out


def load_chrom_seq(path_gz: Path) -> str:
    parts: list[str] = []
    with gzip.open(path_gz, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                continue
            parts.append(line.strip().upper())
    return "".join(parts)


def build_seed_to_guides(guides: list[Guide], max_seed_mm: int = 1) -> dict[str, list[Guide]]:
    seed_to_guides: dict[str, list[Guide]] = defaultdict(list)
    for g in guides:
        for s in seed_variants(g.seed, max_seed_mm):
            seed_to_guides[s].append(g)
    return seed_to_guides


def scan_chromosome(
    chrom: str,
    seq: str,
    seed_to_guides: dict[str, list[Guide]],
    *,
    max_mismatches: int = MAX_MISMATCHES,
    collect_hits: bool = False,
    max_hits_per_guide: int = 50,
) -> tuple[dict[str, dict[str, int]], dict[str, list[OffTargetHit]]]:
    counts: dict[str, dict[str, int]] = defaultdict(
        lambda: {f"mm{i}": 0 for i in range(max_mismatches + 1)}
    )
    hits: dict[str, list[OffTargetHit]] = defaultdict(list)
    n = len(seq)

    def maybe_add_hit(gid: str, hit: OffTargetHit) -> None:
        if not collect_hits:
            return
        bucket = hits[gid]
        bucket.append(hit)
        bucket.sort(key=lambda h: (h.mismatches, h.chrom, h.start))
        if len(bucket) > max_hits_per_guide:
            del bucket[max_hits_per_guide:]

    # Forward NGG
    for i in range(0, n - (GUIDE_LEN + 3) + 1):
        if seq[i + GUIDE_LEN + 1 : i + GUIDE_LEN + 3] != "GG":
            continue
        proto = seq[i : i + GUIDE_LEN]
        if "N" in proto:
            continue
        seed = proto[-SEED_LEN:]
        if seed not in seed_to_guides:
            continue
        for g in seed_to_guides[seed]:
            mm = hamming(proto, g.seq)
            if mm <= max_mismatches:
                counts[g.guide_id][f"mm{mm}"] += 1
                maybe_add_hit(
                    g.guide_id,
                    OffTargetHit(chrom, i, i + GUIDE_LEN, "+", mm, proto, seq[i + GUIDE_LEN : i + GUIDE_LEN + 3]),
                )

    # Reverse CCN
    for i in range(0, n - (3 + GUIDE_LEN) + 1):
        if seq[i] != "C" or seq[i + 1] != "C":
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
            if mm <= max_mismatches:
                counts[g.guide_id][f"mm{mm}"] += 1
                maybe_add_hit(
                    g.guide_id,
                    OffTargetHit(chrom, i + 3, i + 3 + GUIDE_LEN, "-", mm, proto, seq[i : i + 3]),
                )

    return counts, hits


def cas_offinder_score(counts: dict) -> float:
    """Cas-OFFinder-inspired mismatch weighting (higher = more liability)."""
    weights = {0: 1.0, 1: 0.5, 2: 0.25, 3: 0.125, 4: 0.0625}
    total = 0.0
    for k, v in counts.items():
        if isinstance(k, int):
            mm = k
        else:
            mm = int(str(k).replace("mm", ""))
        total += weights.get(mm, 0.0) * float(v)
    return total


def risk_label(score: float, mm0: int) -> str:
    if mm0 >= 5 or score >= 12:
        return "elevated"
    if mm0 >= 2 or score >= 4:
        return "moderate"
    return "low"
