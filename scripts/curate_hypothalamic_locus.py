"""
Extract 15q11-q13 cis-regulatory features from hypothalamic neuron atlas (GSE152098).

Rationale: Digital twin requires cell-type-matched baseline chromatin architecture
(ATAC accessibility, promoter-anchored loops, expression) at the PWS locus in
hypothalamic neurons — the clinically relevant cell type (Cousminer 2021).

Output: data/curated/hypothalamic_atac_chr15.parquet
        data/curated/hypothalamic_loops_chr15.parquet
        data/curated/hypothalamic_expression_chr15.parquet
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "locus.yaml"
DATA = ROOT / "data"
OUT = DATA / "curated"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)


def load_locus_window() -> tuple[str, int, int]:
    with open(CONFIG, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    r = cfg["regions"]["pws_critical"]
    return cfg["chromosome"], r["start"], r["end"]


def in_window(chrom: str, start: int, end: int, target_chrom: str, win_start: int, win_end: int) -> bool:
    c = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return c == target_chrom and start <= win_end and end >= win_start


def curate_atac(chrom: str, start: int, end: int) -> pd.DataFrame:
    path = DATA / "gse152090" / "GSE152090_HypothalamusESC_atac_counts_fpkm.txt.gz"
    log.info("Filtering ATAC peaks from %s", path.name)
    chunks = []
    with gzip.open(path, "rt") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            peak_chr, peak_start, peak_end = parts[2], int(parts[3]), int(parts[4])
            if in_window(peak_chr, peak_start, peak_end, chrom, start, end):
                chunks.append(parts)

    cols = ["id_num", "id", "chrom", "start", "end", "peak_len", "sample", "count", "lib_size", "fpkm", "log2fpkm"]
    df = pd.DataFrame(chunks, columns=cols)
    for c in ["start", "end", "peak_len", "count", "fpkm", "log2fpkm"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["genome_build"] = "hg38"
    log.info("  %d ATAC peaks in PWS window", len(df))
    return df


def curate_loops(chrom: str, start: int, end: int) -> pd.DataFrame:
    files = {
        "ESC": DATA / "gse152092" / "GSE152092_HumanESCells_4frag_chicagoRes.ibed.txt.gz",
        "HypothalamicProg": DATA / "gse152092" / "GSE152092_HypothalamicProg_4frag_chicagoRes.ibed.txt.gz",
        "HypothalamicNeurons": DATA / "gse152092" / "GSE152092_HypothalamicNeurons_4frag_chicagoRes.ibed.txt.gz",
    }
    frames = []
    for cell_type, path in files.items():
        if not path.exists():
            continue
        log.info("Filtering Capture-C loops: %s (%s)", path.name, cell_type)
        rows = []
        with gzip.open(path, "rt") as f:
            f.readline()  # header
            for line in f:
                p = line.strip().split("\t")
                if len(p) < 9:
                    continue
                b_chr, b_start, b_end = p[0], int(p[1]), int(p[2])
                o_chr, o_start, o_end = p[4], int(p[5]), int(p[6])
                if (in_window(b_chr, b_start, b_end, chrom, start, end) or
                        in_window(o_chr, o_start, o_end, chrom, start, end)):
                    rows.append({
                        "cell_type": cell_type,
                        "bait_chrom": b_chr, "bait_start": b_start, "bait_end": b_end, "bait_name": p[3],
                        "other_chrom": o_chr, "other_start": o_start, "other_end": o_end, "other_name": p[7],
                        "n_reads": int(p[8]) if p[8].isdigit() else p[8],
                        "score": float(p[9]) if len(p) > 9 else None,
                        "genome_build": "hg38",
                    })
        frames.append(pd.DataFrame(rows))
        log.info("  %d interactions involving PWS locus", len(rows))

    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def curate_expression(chrom: str, start: int, end: int) -> pd.DataFrame:
    path = DATA / "gse152097" / "GSE152097_HypothalamusESC_rnaseq_tpm_perReplicate.csv.gz"
    log.info("Loading RNA TPM: %s", path.name)
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f)

    # Filter by PWS gene symbols
    pws_symbols = ["SNRPN", "SNHG14", "UBE3A", "MAGEL2", "NDN", "MKRN3", "NPAP1", "ATP10A",
                   "GABRB3", "GABRA5", "GABRG3", "PWRN1", "IPW"]
    gene_col = "gene_name" if "gene_name" in df.columns else df.columns[0]
    sub = df[df[gene_col].astype(str).isin(pws_symbols)].copy()

    sub["genome_build"] = "hg38"
    log.info("  %d expression records for PWS genes", len(sub))
    return sub


def main():
    chrom, start, end = load_locus_window()

    atac = curate_atac(chrom, start, end)
    atac.to_parquet(OUT / "hypothalamic_atac_chr15.parquet", index=False)

    loops = curate_loops(chrom, start, end)
    loops.to_parquet(OUT / "hypothalamic_loops_chr15.parquet", index=False)

    expr = curate_expression(chrom, start, end)
    expr.to_parquet(OUT / "hypothalamic_expression_chr15.parquet", index=False)

    log.info("Hypothalamic locus curation complete -> %s", OUT)


if __name__ == "__main__":
    main()
