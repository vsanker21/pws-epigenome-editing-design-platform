"""
Coordinate harmonization utilities: hg19 → GRCh38 liftOver.

Uses pyliftover with UCSC chain files (downloaded on first use).
Scientific rationale: Gersbach 2025 data is hg19; locus annotation and
hypothalamic atlas are hg38. Base-pair merging requires a single reference.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from typing import NamedTuple

import pandas as pd

log = logging.getLogger(__name__)

REFERENCE_BUILD = "GRCh38"
SOURCE_BUILD_GERSBACH = "hg19"


class LiftResult(NamedTuple):
    chrom: str
    start: int
    end: int
    strand: str
    mapped: bool


@lru_cache(maxsize=1)
def get_liftover(source: str = "hg19", target: str = "hg38"):
    from pyliftover import LiftOver
    return LiftOver(source, target)


def lift_position(chrom: str, pos: int, lo=None) -> tuple[str, int, str] | None:
    """Lift single 0-based or 1-based position; returns (chrom, pos, strand) or None."""
    if lo is None:
        lo = get_liftover()
    c = chrom if chrom.startswith("chr") else f"chr{chrom}"
    results = lo.convert_coordinate(c, int(pos))
    if not results:
        return None
    r = results[0]
    return r[0], int(r[1]), r[2]


def lift_interval(
    chrom: str,
    start: int,
    end: int,
    lo=None,
) -> LiftResult:
    """Lift interval by lifting start and end independently; require same chrom."""
    if lo is None:
        lo = get_liftover()
    c = chrom if chrom.startswith("chr") else f"chr{chrom}"
    rs = lo.convert_coordinate(c, int(start))
    re = lo.convert_coordinate(c, int(end))
    if not rs or not re:
        return LiftResult(c, start, end, "+", False)
    if rs[0][0] != re[0][0]:
        # Use start lift only + preserve width as fallback
        new_start = rs[0][1]
        width = max(1, int(end) - int(start))
        return LiftResult(rs[0][0], new_start, new_start + width, rs[0][2], True)
    new_start, new_end = sorted([rs[0][1], re[0][1]])
    if new_start == new_end:
        new_end = new_start + max(1, int(end) - int(start))
    return LiftResult(rs[0][0], new_start, new_end, rs[0][2], True)


def lift_dataframe(
    df: pd.DataFrame,
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str | None = "end",
    source_build: str = "hg19",
) -> pd.DataFrame:
    """Add hg38 coordinate columns to a dataframe."""
    if source_build == REFERENCE_BUILD.lower() or source_build == "hg38":
        out = df.copy()
        out["chrom_hg38"] = out[chrom_col]
        out["start_hg38"] = out[start_col]
        out["end_hg38"] = out[end_col] if end_col else out[start_col]
        out["liftover_mapped"] = True
        out["genome_build_hg38"] = REFERENCE_BUILD
        return out

    lo = get_liftover("hg19", "hg38")
    chroms, starts, ends, strands, mapped = [], [], [], [], []

    for _, row in df.iterrows():
        chrom = str(row[chrom_col])
        start = int(row[start_col])
        end = int(row[end_col]) if end_col and end_col in row and pd.notna(row[end_col]) else start + 20
        result = lift_interval(chrom, start, end, lo)
        chroms.append(result.chrom)
        starts.append(result.start)
        ends.append(result.end)
        strands.append(result.strand)
        mapped.append(result.mapped)

    out = df.copy()
    out["chrom_hg38"] = chroms
    out["start_hg38"] = starts
    out["end_hg38"] = ends
    out["strand_hg38"] = strands
    out["liftover_mapped"] = mapped
    out["genome_build_hg38"] = REFERENCE_BUILD
    return out


def interval_overlap(s1: int, e1: int, s2: int, e2: int) -> bool:
    return s1 <= e2 and s2 <= e1


def overlap_length(s1: int, e1: int, s2: int, e2: int) -> int:
    if not interval_overlap(s1, e1, s2, e2):
        return 0
    return min(e1, e2) - max(s1, s2) + 1
