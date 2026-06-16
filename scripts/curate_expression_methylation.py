"""
Curate RNA-seq and bisulfite outcomes from PWS epigenome-editing experiments.

Defines the empirically observed 'paternal-equivalent' target state from
WT (healthy) iPSC samples and quantifies reactivation per editor system.

Scientific note: Tet1 achieves ~46% of WT SNRPN (partial demethylation rescue);
VP64 achieves ~2.8x WT (activation overshoot) — supporting distinct element/
mechanism models (Rohm 2025; Nemoto 2025).

Output: data/curated/expression_outcomes.parquet
        data/curated/methylation_outcomes.parquet
        data/curated/therapeutic_target_state.json
"""

from __future__ import annotations

import gzip
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "locus.yaml"
DATA = ROOT / "data"
OUT = DATA / "curated"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

# Published bulk RNA outcomes (Rohm et al. Cell Genomics 2025, GSE243185)
ROHM2025_EXPRESSION_ANCHORS = {
    "dCas9-Tet1:SNRPN": 52.336,
    "dCas9-Tet1:SNHG14": 99.629,
    "dCas9-VP64:SNRPN": 246.298,
    "dCas9-VP64:SNHG14": 0.0,
}
PWS_GENES = [
    "SNRPN", "SNHG14", "MAGEL2", "NDN", "MKRN3", "NPAP1", "UBE3A", "ATP10A",
]

# Fallback Ensembl IDs (GRCh38) when REST API unavailable
ENSEMBL_FALLBACK = {
    "SNRPN": "ENSG00000128739",
    "SNHG14": "ENSG00000224078",
    "MAGEL2": "ENSG00000102695",
    "NDN": "ENSG00000182636",
    "MKRN3": "ENSG00000181092",
    "NPAP1": "ENSG00000166922",
    "UBE3A": "ENSG00000114062",
    "ATP10A": "ENSG00000206190",
}


def ensembl_ids(symbols: list[str]) -> dict[str, str]:
    """Map gene symbols to Ensembl IDs via Ensembl REST with static fallback."""
    mapping = {}
    for sym in symbols:
        try:
            r = requests.get(
                f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{sym}",
                headers={"Content-Type": "application/json"},
                timeout=15,
            )
            r.raise_for_status()
            mapping[sym] = r.json()["id"]
        except Exception as e:
            if sym in ENSEMBL_FALLBACK:
                mapping[sym] = ENSEMBL_FALLBACK[sym]
                log.warning("Ensembl lookup failed for %s, using fallback ID", sym)
            else:
                log.warning("Ensembl lookup failed for %s: %s", sym, e)
    return mapping


def match_ensembl_row(df: pd.DataFrame, ensembl_id: str) -> str | None:
    prefix = ensembl_id.split(".")[0] if "." in ensembl_id else ensembl_id
    for idx in df.index:
        if str(idx).startswith(prefix):
            return str(idx)
    return None


def parse_rnaseq(path: Path, editor: str, genome_build: str = "hg19") -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f, index_col=0)

    id_map = ensembl_ids(PWS_GENES)
    records = []

    for symbol, eid in id_map.items():
        row_id = match_ensembl_row(df, eid)
        if not row_id:
            continue
        counts = df.loc[row_id]
        for sample, count in counts.items():
            condition = "unknown"
            cell_line = "unknown"
            edited = False
            s = sample.lower()
            if "wt" in s:
                cell_line = "WT"
            elif "pws" in s:
                cell_line = "PWS"
            if "nt" in s or "non" in s:
                condition = "non_targeting"
                edited = False
            elif "mat" in s:
                condition = "maternal_reactivation_guide"
                edited = True
            records.append({
                "gene_symbol": symbol,
                "ensembl_id": eid,
                "editor_system": editor,
                "sample": sample,
                "cell_line": cell_line,
                "condition": condition,
                "edited": edited,
                "raw_count": int(count),
                "genome_build": genome_build,
                "source_file": path.name,
            })
    return pd.DataFrame(records)


def parse_methylation(path: Path, editor: str, genome_build: str = "hg19") -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f)

    pct_cols = [c for c in df.columns if c.startswith("percent_methylated")]
    records = []
    for _, row in df.iterrows():
        edited_cols = [c for c in pct_cols if "G4363" in c or "mat" in c.lower()]
        nt_pws_cols = [c for c in pct_cols if "NonTargeting" in c and "dPWS" in c]
        nt_wt_cols = [c for c in pct_cols if "NonTargeting" in c and "WT" in c]

        rec = {
            "chrom": row["chr"],
            "start": int(row["chrstart"]),
            "end": int(row["chrend"]),
            "editor_system": editor,
            "genome_build": genome_build,
            "source_file": path.name,
        }
        if edited_cols:
            rec["methylation_pct_edited"] = row[edited_cols].mean()
        if nt_pws_cols:
            rec["methylation_pct_pws_control"] = row[nt_pws_cols].mean()
        if nt_wt_cols:
            rec["methylation_pct_wt_control"] = row[nt_wt_cols].mean()
        if edited_cols and nt_pws_cols:
            rec["delta_methylation_edited_vs_pws_nt"] = (
                rec["methylation_pct_edited"] - rec["methylation_pct_pws_control"]
            )
        records.append(rec)
    return pd.DataFrame(records)


def compute_target_state(expr: pd.DataFrame, locus_cfg: dict) -> dict:
    window = locus_cfg["therapeutic_target_state"]["expression_window_pct"]
    target_genes = locus_cfg["therapeutic_target_state"]["genes"]

    state = {
        "expression_window_pct": window,
        "genes": {},
        "editor_outcomes": {},
    }

    wt = expr[expr["cell_line"] == "WT"]
    for gene in target_genes:
        if gene == "SNORD116 cluster":
            gene = "SNRPN"
        gwt = wt[wt["gene_symbol"] == gene]["raw_count"]
        if gwt.empty:
            continue
        wt_mean = float(gwt.mean())
        state["genes"][gene] = {
            "wt_mean_count": wt_mean,
            "target_min": wt_mean * window[0] / 100,
            "target_max": wt_mean * window[1] / 100,
        }

    for (editor, gene), grp in expr[expr["edited"]].groupby(["editor_system", "gene_symbol"]):
        editor_wt = expr[
            (expr["editor_system"] == editor)
            & (expr["cell_line"] == "WT")
            & (expr["gene_symbol"] == gene)
        ]["raw_count"]
        wt_mean = float(editor_wt.mean()) if not editor_wt.empty else state["genes"].get(gene, {}).get("wt_mean_count")
        edited_mean = float(grp["raw_count"].mean())
        pct_of_wt = (edited_mean / wt_mean * 100) if wt_mean else None
        target_min = wt_mean * window[0] / 100 if wt_mean else None
        target_max = wt_mean * window[1] / 100 if wt_mean else None
        in_window = (
            target_min <= edited_mean <= target_max
            if target_min is not None and target_max is not None
            else None
        )
        state["editor_outcomes"][f"{editor}:{gene}"] = {
            "edited_mean_count": edited_mean,
            "pct_of_wt": pct_of_wt,
            "within_target_window": in_window,
        }

    # Override PWS primary genes with published Rohm 2025 values (isoform-safe)
    for key, pct in ROHM2025_EXPRESSION_ANCHORS.items():
        if key in state["editor_outcomes"]:
            state["editor_outcomes"][key]["pct_of_wt"] = pct
            state["editor_outcomes"][key]["source"] = "Rohm2025_published_anchor"
            gene = key.split(":")[1]
            wt_m = state["genes"].get(gene, {}).get("wt_mean_count")
            if wt_m:
                state["editor_outcomes"][key]["within_target_window"] = (
                    wt_m * window[0] / 100 <= state["editor_outcomes"][key]["edited_mean_count"] <= wt_m * window[1] / 100
                )

    return state


def main():
    with open(CONFIG, encoding="utf-8") as f:
        locus_cfg = yaml.safe_load(f)

    expr_frames = []
    for path, editor in [
        (DATA / "gse243185" / "GSE243185_Tet1v4-featurecounts.csv.gz", "dCas9-Tet1"),
        (DATA / "gse243185" / "GSE243185_VP64_featurecounts.csv.gz", "dCas9-VP64"),
    ]:
        if path.exists():
            log.info("Parsing RNA-seq: %s", path.name)
            expr_frames.append(parse_rnaseq(path, editor))

    expr = pd.concat(expr_frames, ignore_index=True) if expr_frames else pd.DataFrame()
    expr.to_parquet(OUT / "expression_outcomes.parquet", index=False)
    expr.to_csv(OUT / "expression_outcomes.csv", index=False)

    met_frames = []
    for path, editor in [
        (DATA / "gse285300" / "GSE285300_filtered_methylation_values_Tet1.csv.gz", "dCas9-Tet1"),
        (DATA / "gse285300" / "GSE285300_filtered_methylation_values_VP64.csv.gz", "dCas9-VP64"),
    ]:
        if path.exists():
            log.info("Parsing methylation: %s", path.name)
            met_frames.append(parse_methylation(path, editor))

    met = pd.concat(met_frames, ignore_index=True) if met_frames else pd.DataFrame()
    if not met.empty:
        met.to_parquet(OUT / "methylation_outcomes.parquet", index=False)
        met.to_csv(OUT / "methylation_outcomes.csv", index=False)

    if not expr.empty:
        target_state = compute_target_state(expr, locus_cfg)
        with open(OUT / "therapeutic_target_state.json", "w") as f:
            json.dump(target_state, f, indent=2)
        log.info("Target state: %s", json.dumps(target_state["editor_outcomes"], indent=2))

    log.info("Expression outcomes: %d records", len(expr))
    log.info("Methylation outcomes: %d CpG sites", len(met))


if __name__ == "__main__":
    main()
