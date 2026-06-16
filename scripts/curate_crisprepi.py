"""
Curate CRISPRepi atlas for transfer-learning priors on epigenome-editing outcomes.

Integrates 2,520 human CRISPRepi records to provide editor-class effect priors
beyond the PWS-specific Rohm 2025 dataset (n=57 Tet1 guides).

Scientific rationale: CRISPRepi captures demethylation vs activation outcomes
across loci; PWS-specific model uses these as informative priors on mechanism
class, not as direct PWS outcome predictors.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CRISPREPI = ROOT / "data" / "crisprepi"
OUT = ROOT / "data" / "curated"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

EDITOR_MAP = {
    "dCas9-Tet1": ["dCas9_Tet1", "dCas9-Tet1", "SunTag14aa_TET1", "SunTag22aa_TET1", "SunTag-TET1"],
    "dCas9-VP64": ["dCas9_VP64", "dCas9-VP64", "VP64"],
    "dCas9-KRAB": ["dCas9_KRAB", "dCas9-KRAB", "KRAB"],
    "dCas9-p300": ["dCas9_p300", "p300", "SunTag_p300"],
}


def load_crisprepi() -> pd.DataFrame:
    path = CRISPREPI / "Homo_sapiens.tsv.gz"
    return pd.read_csv(path, sep="\t")


def normalize_editor(dtype: str) -> str:
    if pd.isna(dtype):
        return "other"
    s = str(dtype)
    for canonical, patterns in EDITOR_MAP.items():
        if any(p.lower() in s.lower() for p in patterns):
            return canonical
    return "other"


def summarize_editor_priors(df: pd.DataFrame) -> dict:
    df = df.copy()
    df["editor_canonical"] = df["dCas9_type"].apply(normalize_editor)

    priors = {}
    for editor in ["dCas9-Tet1", "dCas9-VP64", "dCas9-KRAB", "dCas9-p300"]:
        sub = df[df["editor_canonical"] == editor]
        if sub.empty:
            continue
        epi = sub["EpiEffect"].fillna("unknown").value_counts(normalize=True)
        priors[editor] = {
            "n_crisprepi_records": len(sub),
            "primary_epigenetic_effects": epi.head(5).to_dict(),
            "demethylation_fraction": float((sub["EpiEffect"] == "demethylation").mean()),
            "activation_fraction": float(
                sub["EpiEffect"].astype(str).str.contains("transcription|activation|H3K27ac", case=False, na=False).mean()
            ),
            "cell_types": sub["Cell_type"].dropna().value_counts().head(5).to_dict(),
            "tissues": sub["Tissue_type"].dropna().value_counts().head(5).to_dict(),
        }
    return priors


def pws_locus_coordinate_hits(df: pd.DataFrame) -> pd.DataFrame:
    """Records with sgRNA coordinates on chr15."""
    if "sgRNA_coordinate" not in df.columns:
        return pd.DataFrame()
    coords = df["sgRNA_coordinate"].fillna("").astype(str)
    mask = coords.str.contains(r"chr15[:/]", case=False, regex=True)
    return df.loc[mask]


def main():
    df = load_crisprepi()
    priors = summarize_editor_priors(df)
    chr15_hits = pws_locus_coordinate_hits(df)

    # Per-file editor summaries
    file_summaries = {}
    for f in CRISPREPI.glob("*.tsv.gz"):
        if f.name == "Homo_sapiens.tsv.gz":
            continue
        sub = pd.read_csv(f, sep="\t")
        file_summaries[f.stem] = {
            "n_records": len(sub),
            "unique_genes": int(sub["Target_gene"].nunique()) if "Target_gene" in sub.columns else 0,
        }

    report = {
        "source": "CRISPRepi (Shi et al. NAR 2025;53(D1):D901)",
        "n_human_records": len(df),
        "editor_class_priors": priors,
        "chr15_coordinate_hits": len(chr15_hits),
        "chr15_records": chr15_hits[["Target_gene", "dCas9_type", "EpiEffect", "sgRNA_coordinate"]].head(10).to_dict(orient="records") if len(chr15_hits) else [],
        "file_summaries": file_summaries,
        "transfer_learning_use": (
            "Tet1 demethylation prior: {:.0%} of CRISPRepi Tet1 records report demethylation. "
            "VP64 activation prior: {:.0%} report transcriptional effects. "
            "These inform mechanism-class confidence in forward model, not locus-specific outcomes."
        ).format(
            priors.get("dCas9-Tet1", {}).get("demethylation_fraction", 0),
            priors.get("dCas9-VP64", {}).get("activation_fraction", 0),
        ),
    }

    with open(OUT / "crisprepi_transfer_priors.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    df.to_parquet(OUT / "crisprepi_human_records.parquet", index=False)
    log.info("CRISPRepi: %d human records, %d chr15 hits", len(df), len(chr15_hits))
    log.info("Tet1 demethylation prior: %.0f%%", 100 * priors.get("dCas9-Tet1", {}).get("demethylation_fraction", 0))
    log.info("Output -> %s", OUT / "crisprepi_transfer_priors.json")


if __name__ == "__main__":
    main()
