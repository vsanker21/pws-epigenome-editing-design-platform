"""Parse Nemoto 2025 GSE262700 organoid scRNA-seq for PWS gene validation."""

from __future__ import annotations

import gzip
import json
import logging
from pathlib import Path

import pandas as pd
import scipy.io

ROOT = Path(__file__).resolve().parents[1]
EXTRACTED = ROOT / "data" / "gse262700" / "extracted"
OUT = ROOT / "data" / "curated"

PWS_GENES = ["SNRPN", "SNHG14", "MAGEL2", "NDN", "MKRN3", "NPAP1", "UBE3A", "ATP10A"]

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)


def load_mtx_sample(prefix: str) -> tuple[pd.DataFrame, list[str], object]:
    with gzip.open(EXTRACTED / f"{prefix}_features.tsv.gz", "rt") as f:
        features = [line.strip().split("\t") for line in f]
    feat_df = pd.DataFrame(features, columns=["ensembl", "symbol", "type"])
    with gzip.open(EXTRACTED / f"{prefix}_barcodes.tsv.gz", "rt") as f:
        barcodes = [line.strip() for line in f]
    mat = scipy.io.mmread(EXTRACTED / f"{prefix}_matrix.mtx.gz").tocsr()
    return feat_df, barcodes, mat


def summarize_pws_genes(feat_df: pd.DataFrame, mat, sample_name: str) -> list[dict]:
    records = []
    for sym in PWS_GENES:
        rows = feat_df[feat_df["symbol"] == sym]
        if rows.empty:
            continue
        idx = rows.index[0]
        counts = mat[idx, :].toarray().flatten()
        records.append({
            "sample": sample_name,
            "gene_symbol": sym,
            "ensembl": rows.iloc[0]["ensembl"],
            "n_cells": len(counts),
            "mean_counts": float(counts.mean()),
            "median_counts": float(pd.Series(counts).median()),
            "pct_expressing": float((counts > 0).mean() * 100),
            "total_counts": int(counts.sum()),
        })
    return records


def main():
    samples = {
        "GSM8174019_guide_plus": "edited_organoid",
        "GSM8174020_guide_minus": "control_organoid",
    }
    all_records = []
    for prefix, label in samples.items():
        feat_df, barcodes, mat = load_mtx_sample(prefix)
        log.info("%s: %d cells, %d genes", label, len(barcodes), len(feat_df))
        all_records.extend(summarize_pws_genes(feat_df, mat, label))

    df = pd.DataFrame(all_records)
    df.to_parquet(OUT / "organoid_pws_expression.parquet", index=False)
    df.to_csv(OUT / "organoid_pws_expression.csv", index=False)

    # Fold-change edited vs control
    pivot = df.pivot(index="gene_symbol", columns="sample", values="mean_counts")
    if "edited_organoid" in pivot.columns and "control_organoid" in pivot.columns:
        pivot["fold_change_edited_vs_control"] = pivot["edited_organoid"] / pivot["control_organoid"].replace(0, 1)
        pivot.to_csv(OUT / "organoid_pws_foldchange.csv")

    summary = {
        "source": "Nemoto et al. 2025 GSE262700",
        "edited_sample": "guide_plus (epigenome-edited hypothalamic organoids)",
        "control_sample": "guide_minus (unedited PWS organoids)",
        "genes": df.to_dict(orient="records"),
    }
    with open(OUT / "organoid_validation_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    log.info("Organoid PWS expression summary:\n%s", pivot if "fold_change_edited_vs_control" in pivot.columns else df)


if __name__ == "__main__":
    main()
