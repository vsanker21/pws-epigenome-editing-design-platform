"""
Curate neuron-specific RNA outcomes from GSE285305 (Rohm 2025 neuronal differentiation).

Provides cell-type-matched validation beyond iPSC bulk RNA (GSE243185).
Maps chr15 features to annotated PWS genes via coordinate overlap.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "data" / "curated"
CONFIG = ROOT / "config" / "locus.yaml"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

PWS_GENES = ["SNRPN", "SNHG14", "MAGEL2", "NDN", "MKRN3", "NPAP1", "UBE3A", "ATP10A"]


def load_gene_bed() -> pd.DataFrame:
    rows = []
    with open(ROOT / "data" / "locus_annotation" / "genes.bed") as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 4:
                continue
            rows.append({"chrom": p[0], "start": int(p[1]), "end": int(p[2]), "gene": p[3]})
    return pd.DataFrame(rows).drop_duplicates(subset=["gene"])


PWS_GENE_REGIONS_HG19 = {
    "SNRPN": (25100000, 25110000),
    "SNHG14": (24820000, 25450000),
    "UBE3A": (25600000, 25800000),
    "ATP10A": (25900000, 26200000),
    "MAGEL2": (23700000, 23900000),
    "NDN": (23950000, 24100000),
    "MKRN3": (23600000, 23750000),
    "NPAP1": (24800000, 24900000),
}


def overlap_gene_hg19(chrom: str, start: int, end: int) -> str | None:
    if chrom != "chr15":
        return None
    mid = (start + end) // 2
    for gene, (lo, hi) in PWS_GENE_REGIONS_HG19.items():
        if lo <= mid <= hi:
            return gene
    return None


def parse_featurecounts(path: Path, editor: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    records = []
    for _, row in df.iterrows():
        if row["Chr"] != "chr15":
            continue
        gene = overlap_gene_hg19(row["Chr"], int(row["Start"]), int(row["End"]))
        if gene not in PWS_GENES:
            continue
        for col in df.columns:
            if col in ("Geneid", "Chr", "Start", "End", "Strand", "Length"):
                continue
            records.append({
                "gene_symbol": gene,
                "editor_system": editor,
                "sample": col,
                "raw_count": float(row[col]),
                "cell_type": "neuron",
                "source_file": path.name,
            })
    return pd.DataFrame(records)


def summarize_outcomes(expr: pd.DataFrame) -> dict:
    outcomes = {}
    for editor in expr["editor_system"].unique():
        for gene in PWS_GENES:
            sub = expr[(expr["editor_system"] == editor) & (expr["gene_symbol"] == gene)]
            if sub.empty:
                continue
            edited = sub[sub["sample"].str.contains("mat", case=False)]["raw_count"].mean()
            wt = sub[sub["sample"].str.contains("WT", case=False)]["raw_count"].mean()
            pws_nt = sub[sub["sample"].str.contains("PWS") & sub["sample"].str.contains("NT", case=False)]["raw_count"].mean()
            outcomes[f"{editor}:{gene}"] = {
                "neuron_edited_mean": float(edited) if pd.notna(edited) else None,
                "neuron_wt_mean": float(wt) if pd.notna(wt) else None,
                "neuron_pws_nt_mean": float(pws_nt) if pd.notna(pws_nt) else None,
                "pct_of_wt_neuron": float(100 * edited / wt) if wt and wt > 0 and pd.notna(edited) else None,
            }
    return outcomes


def main():
    tet1_path = DATA / "gse285305" / "GSE285305_Tet1-featurecounts.csv.gz"
    vp64_path = DATA / "gse285305" / "GSE285305_VP64-featurecounts.csv.gz"

    frames = []
    if tet1_path.exists():
        frames.append(parse_featurecounts(tet1_path, "dCas9-Tet1"))
    if vp64_path.exists():
        frames.append(parse_featurecounts(vp64_path, "dCas9-VP64"))

    if not frames:
        log.warning("GSE285305 not found")
        return

    expr = pd.concat(frames, ignore_index=True)
    expr.to_parquet(OUT / "neuron_expression_outcomes.parquet", index=False)
    expr.to_csv(OUT / "neuron_expression_outcomes.csv", index=False)

    outcomes = summarize_outcomes(expr)
    report = {
        "source": "GSE285305 neuron featureCounts (Rohm 2025)",
        "cell_type": "differentiated neuron (not iPSC)",
        "outcomes": outcomes,
        "comparison_to_ipsc": "Enables neuron-specific held-out validation vs GSE243185 iPSC bulk RNA",
    }
    with open(OUT / "neuron_target_state.json", "w") as f:
        json.dump(report, f, indent=2)

    snrpn = outcomes.get("dCas9-Tet1:SNRPN", {})
    log.info("Neuron Tet1 SNRPN: %.1f%% WT", snrpn.get("pct_of_wt_neuron", 0) or 0)
    log.info("Output -> %s", OUT / "neuron_target_state.json")


if __name__ == "__main__":
    main()
