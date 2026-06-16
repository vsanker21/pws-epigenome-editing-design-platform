"""
Phase 12c: Probe-level DMR mapping for GSE28525 and GSE298378 at PWS locus.

Uses the Illumina HumanMethylation450 manifest (hg19 coordinates) from UCSC ENCODE
supplemental, lifts probes to hg38, and maps imprinted DMRs / editing-associated
methylation changes within chr15:24.8-32.7M.

GSE28525 (Kelsey 2011): genome-wide UPD imprinting screen — identifies canonical
imprinted DMRs via pUPD vs mUPD vs biparental controls.

GSE298378: processed 450K beta values for editing-related samples (Rohm cohort).
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
DATA = ROOT / "data"
MANIFEST_DIR = DATA / "manifests"
CURATED = DATA / "curated"
OUT = DATA / "curated"
MODELS = DATA / "models"

MANIFEST_URL = (
    "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/"
    "wgEncodeHaibMethyl450/supplemental/wgEncodeHaibMethyl450CpgIslandDetails.txt"
)
MANIFEST_PATH = MANIFEST_DIR / "HumanMethylation450_manifest.txt"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)


def load_pws_window() -> tuple[str, int, int]:
    with open(ROOT / "config" / "locus.yaml", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    r = cfg["regions"]["pws_critical"]
    return cfg["chromosome"], int(r["start"]), int(r["end"])


def download_manifest() -> Path:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    if MANIFEST_PATH.exists() and MANIFEST_PATH.stat().st_size > 1_000_000:
        return MANIFEST_PATH
    log.info("Downloading 450K manifest (~188MB)...")
    with requests.get(MANIFEST_URL, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(MANIFEST_PATH, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    return MANIFEST_PATH


def parse_manifest_chr15(path: Path) -> pd.DataFrame:
    """Stream-parse manifest; keep chr15 probes only. Skips Illumina BPM preamble."""
    rows = []
    header: list[str] | None = None
    idx: dict[str, int] = {}
    need = ["IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island"]

    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if header is None:
                if line.startswith("IlmnID,"):
                    header = line.strip().split(",")
                    idx = {h: i for i, h in enumerate(header)}
                continue
            parts = line.strip().split(",")
            if len(parts) < len(header or []):
                continue
            chrom = parts[idx["CHR"]].replace('"', "")
            if chrom not in ("15", "chr15"):
                continue
            rows.append({k: parts[idx[k]].replace('"', "") for k in need if k in idx})
    df = pd.DataFrame(rows)
    df["chrom_hg19"] = "chr15"
    df["pos_hg19"] = pd.to_numeric(df["MAPINFO"], errors="coerce")
    df = df.dropna(subset=["pos_hg19"])
    df["probe_id"] = df["IlmnID"]
    return df


def lift_manifest_to_hg38(manifest15: pd.DataFrame) -> pd.DataFrame:
    from liftover_utils import lift_position

    chroms, starts, mapped = [], [], []
    for _, r in manifest15.iterrows():
        res = lift_position("chr15", int(r["pos_hg19"]))
        if res:
            chroms.append(res[0])
            starts.append(res[1])
            mapped.append(True)
        else:
            chroms.append("chr15")
            starts.append(int(r["pos_hg19"]))
            mapped.append(False)
    out = manifest15.copy()
    out["chrom_hg38"] = chroms
    out["pos_hg38"] = starts
    out["liftover_mapped"] = mapped
    return out


def beta_from_signals(df: pd.DataFrame, sample: str) -> pd.Series:
    u_col = f"{sample} Unmethylated Signal"
    m_col = f"{sample} Methylated Signal"
    if u_col not in df.columns or m_col not in df.columns:
        return pd.Series(np.nan, index=df.index)
    u = pd.to_numeric(df[u_col], errors="coerce")
    m = pd.to_numeric(df[m_col], errors="coerce")
    return m / (u + m + 100.0)


def analyze_gse28525(manifest: pd.DataFrame, chrom: str, start: int, end: int) -> pd.DataFrame:
    path = DATA / "gse28525" / "GSE28525_signal_intensities.txt.gz"
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f, sep="\t")

    df = df.merge(manifest, left_on="ID_REF", right_on="probe_id", how="inner")
    df = df[(df["chrom_hg38"] == chrom) & (df["pos_hg38"] >= start) & (df["pos_hg38"] <= end)]

    pupd_cols = ["pUPD1", "pUPD2", "pUPD3"]
    biparental_cols = ["Lymphocyte1", "Lymphocyte2", "Lymphocyte3", "Muscle", "Brain", "Buccal DNA"]

    for s in pupd_cols + ["mUPD"] + biparental_cols:
        df[f"beta_{s}"] = beta_from_signals(df, s)

    df["beta_pUPD_mean"] = df[[f"beta_{s}" for s in pupd_cols]].mean(axis=1)
    df["beta_biparental_mean"] = df[[f"beta_{c}" for c in biparental_cols if f"beta_{c}" in df.columns]].mean(axis=1)

    # Imprinted DMR criterion (Kelsey 2011): maternal and paternal UPD deviate from biparental in opposite directions
    df["delta_mUPD"] = df["beta_mUPD"] - df["beta_biparental_mean"]
    df["delta_pUPD"] = df["beta_pUPD_mean"] - df["beta_biparental_mean"]
    df["imprinted_dmr"] = (
        (df["delta_mUPD"].abs() >= 0.25)
        & (df["delta_pUPD"].abs() >= 0.25)
        & (np.sign(df["delta_mUPD"]) != np.sign(df["delta_pUPD"]))
    )
    df["maternally_methylated"] = df["imprinted_dmr"] & (df["delta_mUPD"] > 0)
    df["paternally_methylated"] = df["imprinted_dmr"] & (df["delta_pUPD"] > 0)
    return df


def analyze_gse298378(manifest: pd.DataFrame, chrom: str, start: int, end: int) -> pd.DataFrame:
    path = DATA / "gse298378" / "GSE298378_processed_data.txt.gz"
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f, sep="\t", low_memory=False)

    probe_col = "TargetID" if "TargetID" in df.columns else "NAME"
    beta_cols = [c for c in df.columns if c.endswith(".AVG_Beta")]
    slim = df[[probe_col] + beta_cols].copy()
    slim = slim.rename(columns={probe_col: "probe_id"})
    slim = slim.merge(manifest, on="probe_id", how="inner")
    slim = slim[(slim["chrom_hg38"] == chrom) & (slim["pos_hg38"] >= start) & (slim["pos_hg38"] <= end)]

    for c in beta_cols:
        slim[c] = pd.to_numeric(slim[c], errors="coerce")
    slim["beta_mean"] = slim[beta_cols].mean(axis=1)
    slim["beta_std"] = slim[beta_cols].std(axis=1)
    return slim


def main():
    chrom, start, end = load_pws_window()
    manifest_path = download_manifest()
    manifest15 = parse_manifest_chr15(manifest_path)
    manifest = lift_manifest_to_hg38(manifest15)
    manifest.to_parquet(MANIFEST_DIR / "methylation450_chr15_hg38.parquet", index=False)
    log.info("Manifest chr15 probes: %d (lifted %d)", len(manifest), int(manifest["liftover_mapped"].sum()))

    g28525 = analyze_gse28525(manifest, chrom, start, end)
    g298378 = analyze_gse298378(manifest, chrom, start, end)

    g28525.to_csv(CURATED / "gse28525_pws_probe_methylation.csv", index=False)
    g298378.to_csv(CURATED / "gse298378_pws_probe_methylation.csv", index=False)

    # Summarize imprinted DMRs near PWS-ICR
    icr_lo, icr_hi = 24940000, 25010000
    icr_dmrs = g28525[
        g28525["imprinted_dmr"]
        & (g28525["pos_hg38"] >= icr_lo)
        & (g28525["pos_hg38"] <= icr_hi)
    ]

    report = {
        "manifest_source": MANIFEST_URL,
        "manifest_probes_chr15": len(manifest),
        "pws_window": f"{chrom}:{start}-{end}",
        "gse28525": {
            "probes_in_window": len(g28525),
            "imprinted_dmrs_in_window": int(g28525["imprinted_dmr"].sum()),
            "icr_imprinted_dmrs": int(len(icr_dmrs)),
            "pws_icr_maternally_methylated_probes": int(g28525.loc[
                (g28525["pos_hg38"] >= icr_lo) & (g28525["pos_hg38"] <= icr_hi), "maternally_methylated"
            ].sum()),
        },
        "gse298378": {
            "probes_in_window": len(g298378),
            "mean_beta_range": [
                float(g298378["beta_mean"].min()) if len(g298378) else None,
                float(g298378["beta_mean"].max()) if len(g298378) else None,
            ],
            "n_samples": len([c for c in g298378.columns if c.endswith(".AVG_Beta")]),
        },
        "scientific_rationale": (
            "Probe-level mapping anchors methylation array data to hg38 PWS coordinates, "
            "enabling collateral imprinting risk assessment at known imprinted DMRs."
        ),
        "outputs": {
            "manifest": str(MANIFEST_DIR / "methylation450_chr15_hg38.parquet"),
            "gse28525": str(CURATED / "gse28525_pws_probe_methylation.csv"),
            "gse298378": str(CURATED / "gse298378_pws_probe_methylation.csv"),
        },
    }
    (CURATED / "methylation_dmr_mapping_report.json").write_text(json.dumps(report, indent=2))
    log.info(
        "DMR mapping: GSE28525 %d imprinted DMRs in PWS window; GSE298378 %d probes",
        report["gse28525"]["imprinted_dmrs_in_window"],
        report["gse298378"]["probes_in_window"],
    )


if __name__ == "__main__":
    main()
