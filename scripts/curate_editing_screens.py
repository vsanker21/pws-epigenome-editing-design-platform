"""
Phase 3: Curate PWS epigenome-editing screen data into unified training tables.

Scientific rationale (Rohm et al. Cell Genomics 2025):
- Tet1 demethylation and VP64 activation reactivate maternal SNRPN at DISTINCT elements
- KRAB screen maps repressive elements on paternal allele
- Sublibrary screens tile the PWS locus; full-genome screens provide genome-wide context

Output: data/curated/editing_screens.parquet
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "locus.yaml"
DATA = ROOT / "data"
OUT = DATA / "curated"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

SCREEN_FILES = [
    {
        "path": DATA / "gse285289" / "GSE285289_Tet1_sublib_results.csv.gz",
        "editor": "dCas9-Tet1",
        "screen_type": "sublib_tiling",
        "study": "Rohm2025",
        "genome_build": "hg19",
        "coord_col": "chr15_coord",
    },
    {
        "path": DATA / "gse285295" / "GSE285295_VP64_sublib_results.csv.gz",
        "editor": "dCas9-VP64",
        "screen_type": "sublib_tiling",
        "study": "Rohm2025",
        "genome_build": "hg19",
        "coord_col": "chr15_coord",
    },
    {
        "path": DATA / "gse285293" / "GSE285293_VP64_full_results_with_coords.csv.gz",
        "editor": "dCas9-VP64",
        "screen_type": "full_locus",
        "study": "Rohm2025",
        "genome_build": "hg19",
        "coord_col": "chr15.coordinate",
    },
    {
        "path": DATA / "gse285285" / "GSE285285_KRAB_full_results_with_counts.csv.gz",
        "editor": "dCas9-KRAB",
        "screen_type": "full_locus",
        "study": "Rohm2025",
        "genome_build": "hg19",
        "coord_col": "chr15.coordinate",
    },
]


def load_subregions() -> list[dict]:
    with open(CONFIG, encoding="utf-8") as f:
        locus = yaml.safe_load(f)
    regions = []
    for name, reg in locus["regions"].items():
        regions.append({"name": name, "start": reg["start"], "end": reg["end"]})
    return sorted(regions, key=lambda r: r["start"])


def annotate_subregion(pos: int, regions: list[dict]) -> str:
    for reg in regions:
        if reg["start"] <= pos <= reg["end"]:
            return reg["name"]
    return "outside_pws_critical"


def load_sublib(path: Path, editor: str, screen_type: str, study: str,
                genome_build: str, coord_col: str, regions: list[dict]) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f)

    df = df[df[coord_col].notna()].copy()
    df[coord_col] = df[coord_col].astype(int)
    low_cols = [c for c in df.columns if "_L_" in c]
    high_cols = [c for c in df.columns if "_H_" in c]

    out = pd.DataFrame({
        "grna_id": df["gRNA"],
        "protospacer": df["Protospacer"],
        "editor_system": editor,
        "screen_type": screen_type,
        "study": study,
        "genome_build": genome_build,
        "chrom": "chr15",
        "start": df[coord_col],
        "end": df[coord_col] + 20,
        "strand": ".",
        "log2fc": None,
        "padj": df["padj"],
        "neg_log10_padj": -np.log10(df["padj"].clip(lower=1e-300)),
        "mean_low_bin_counts": df[low_cols].mean(axis=1) if low_cols else None,
        "mean_high_bin_counts": df[high_cols].mean(axis=1) if high_cols else None,
        "screen_score": df[high_cols].mean(axis=1) - df[low_cols].mean(axis=1) if high_cols else None,
        "source_file": path.name,
    })
    out["subregion"] = out["start"].apply(lambda x: annotate_subregion(int(x), regions))
    out["significant"] = out["padj"] < 0.05
    return out


def load_deseq_screen(path: Path, editor: str, screen_type: str, study: str,
                      genome_build: str, coord_col: str, regions: list[dict]) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f)

    out = pd.DataFrame({
        "grna_id": df["gRNA"],
        "protospacer": df["Protospacer"],
        "editor_system": editor,
        "screen_type": screen_type,
        "study": study,
        "genome_build": genome_build,
        "chrom": "chr15",
        "start": pd.to_numeric(df[coord_col], errors="coerce"),
        "end": pd.to_numeric(df.get("end_coord", df[coord_col]), errors="coerce"),
        "strand": df.get("Strand", "."),
        "log2fc": df["log2FoldChange"],
        "padj": df["padj"],
        "neg_log10_padj": -np.log10(df["padj"].clip(lower=1e-300)),
        "base_mean": df.get("baseMean"),
        "screen_score": df["log2FoldChange"],
        "source_file": path.name,
    })
    valid = out["start"].notna()
    out = out[valid].copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].fillna(out["start"] + 20).astype(int)
    out["subregion"] = out["start"].apply(lambda x: annotate_subregion(int(x), regions))
    out["significant"] = out["padj"] < 0.05
    return out


def main():
    regions = load_subregions()
    frames = []

    for spec in SCREEN_FILES:
        if not spec["path"].exists():
            log.warning("Missing %s", spec["path"])
            continue
        log.info("Loading %s", spec["path"].name)
        if spec["screen_type"] == "sublib_tiling":
            df = load_sublib(**spec, regions=regions)
        else:
            df = load_deseq_screen(**spec, regions=regions)
        log.info("  %d guides, %d significant (padj<0.05)", len(df), df["significant"].sum())
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.sort_values(["editor_system", "padj"])

    out_parquet = OUT / "editing_screens.parquet"
    out_csv = OUT / "editing_screens.csv"
    combined.to_parquet(out_parquet, index=False)
    combined.to_csv(out_csv, index=False)

    # Summary by editor and subregion
    summary = (
        combined.groupby(["editor_system", "subregion", "significant"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    summary.to_csv(OUT / "editing_screens_summary.csv", index=False)

    # Top hits per editor (scientifically actionable designs)
    tops = (
        combined[combined["significant"]]
        .sort_values("padj")
        .groupby("editor_system", group_keys=False)
        .head(20)
    )
    tops.to_csv(OUT / "top_editing_hits.csv", index=False)

    log.info("Wrote %d editing records -> %s", len(combined), out_parquet)


if __name__ == "__main__":
    main()
