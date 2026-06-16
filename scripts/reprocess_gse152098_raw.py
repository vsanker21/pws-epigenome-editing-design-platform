"""
Phase 12b: Reprocess GSE152098 RAW bigWig ATAC signal at the PWS locus.

GSE152098 RAW contains hypothalamic neuron / progenitor / ESC ATAC bigWig files
(Cousminer 2021). We extract them and quantify accessibility at:
  - PWS subregions (ICR, SNRPN-SNHG14, GABAA cluster)
  - Top therapeutic guide coordinates

Uses UCSC bigWigSummary (Linux binary via WSL) for scientifically standard
signal quantification when pyBigWig is unavailable on Windows/Python 3.13.
"""

from __future__ import annotations

import json
import logging
import shlex
import subprocess
import tarfile
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW_TAR = DATA / "gse152098" / "GSE152098_RAW.tar"
BW_DIR = DATA / "gse152098" / "bigwig"
OUT = DATA / "curated"
MODELS = ROOT / "data" / "models"
BW_SUMMARY = ROOT / "tools" / "ucsc" / "bigWigSummary"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)


def load_locus() -> dict:
    with open(ROOT / "config" / "locus.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)


def extract_bigwigs() -> list[Path]:
    BW_DIR.mkdir(parents=True, exist_ok=True)
    existing = sorted(BW_DIR.glob("*.bigwig"))
    if existing:
        return existing
    if not RAW_TAR.exists():
        raise FileNotFoundError(RAW_TAR)
    log.info("Extracting bigWig files from %s", RAW_TAR.name)
    with tarfile.open(RAW_TAR, "r") as tar:
        for member in tar.getmembers():
            if member.name.endswith(".bigwig"):
                tar.extract(member, path=BW_DIR)
    return sorted(BW_DIR.glob("**/*.bigwig"))


def to_wsl_path(path: Path) -> str:
    s = str(path.resolve()).replace("\\", "/")
    if len(s) >= 2 and s[1] == ":":
        drive = s[0].lower()
        return f"/mnt/{drive}{s[2:]}"
    return s


def bigwig_mean(path: Path, chrom: str, start: int, end: int) -> float | None:
    """Return mean signal for interval via WSL bigWigSummary."""
    if not BW_SUMMARY.exists():
        raise FileNotFoundError(f"Missing {BW_SUMMARY}")
    wsl_bw = to_wsl_path(path)
    wsl_tool = to_wsl_path(BW_SUMMARY)
    script = f"{shlex.quote(wsl_tool)} {shlex.quote(wsl_bw)} {chrom} {start} {end} 1"
    cmd = ["wsl", "bash", "-lc", script]
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True, timeout=120)
        parts = out.strip().split()
        if len(parts) == 1:
            return float(parts[0])
        if len(parts) >= 3:
            return float(parts[2])
    except subprocess.CalledProcessError as e:
        log.warning("bigWigSummary failed for %s: %s", path.name, e.output[:200])
    return None


def cell_type_from_name(name: str) -> str:
    if "HypothalamicNeurons" in name:
        return "HypothalamicNeurons"
    if "HypothalamicProg" in name:
        return "HypothalamicProg"
    if "HumanESCells" in name:
        return "HumanESCells"
    return "unknown"


def main():
    cfg = load_locus()
    chrom = cfg["chromosome"]
    regions = cfg["regions"]

    bigwigs = extract_bigwigs()
    log.info("Found %d bigWig files", len(bigwigs))

    region_rows = []
    for rname, rdef in regions.items():
        start, end = int(rdef["start"]), int(rdef["end"])
        for bw in bigwigs:
            mean_sig = bigwig_mean(bw, chrom, start, end)
            region_rows.append({
                "region": rname,
                "chrom": chrom,
                "start": start,
                "end": end,
                "cell_type": cell_type_from_name(bw.name),
                "replicate": bw.name,
                "mean_signal": mean_sig,
                "genome_build": "hg38",
            })

    region_df = pd.DataFrame(region_rows)
    region_df.to_parquet(OUT / "gse152098_bigwig_regions.parquet", index=False)

    # Guide-level: top 3 catalog designs, one replicate per cell type for speed
    catalog = pd.read_csv(MODELS / "final_catalog" / "pws_therapeutic_design_catalog.csv")
    rep_by_type: dict[str, Path] = {}
    for bw in bigwigs:
        ct = cell_type_from_name(bw.name)
        rep_by_type.setdefault(ct, bw)

    guide_rows = []
    window = 250
    for _, row in catalog.head(3).iterrows():
        for role in ("tet1", "vp64"):
            gid = row.get(f"{role}_grna_id")
            pos = row.get(f"{role}_hg38_start")
            if pd.isna(gid) or pd.isna(pos):
                continue
            pos = int(pos)
            for ct, bw in rep_by_type.items():
                mean_sig = bigwig_mean(bw, chrom, max(0, pos - window), pos + window)
                guide_rows.append({
                    "catalog_rank": int(row["catalog_rank"]),
                    "grna_id": gid,
                    "role": role,
                    "chrom": chrom,
                    "start": pos - window,
                    "end": pos + window,
                    "cell_type": cell_type_from_name(bw.name),
                    "replicate": bw.name,
                    "mean_signal": mean_sig,
                })

    guide_df = pd.DataFrame(guide_rows)
    guide_df.to_parquet(OUT / "gse152098_bigwig_guides.parquet", index=False)

    # Aggregate neuron vs ESC at ICR for top design
    icr = regions["pws_icr"]
    icr_df = region_df[region_df["region"] == "pws_icr"].copy()
    summary = (
        icr_df.groupby("cell_type", as_index=False)["mean_signal"]
        .mean()
        .rename(columns={"mean_signal": "icr_mean_signal"})
    )
    neuron = summary.loc[summary["cell_type"] == "HypothalamicNeurons", "icr_mean_signal"]
    esc = summary.loc[summary["cell_type"] == "HumanESCells", "icr_mean_signal"]
    fold = None
    if len(neuron) and len(esc) and pd.notna(esc.iloc[0]) and esc.iloc[0] != 0:
        fold = float(neuron.iloc[0] / esc.iloc[0])

    report = {
        "dataset": "GSE152098_RAW",
        "method": "UCSC bigWigSummary via WSL",
        "n_bigwig_files": len(bigwigs),
        "regions_quantified": list(regions.keys()),
        "icr_celltype_means": summary.to_dict(orient="records"),
        "neuron_vs_esc_icr_fold": fold,
        "outputs": {
            "regions": str(OUT / "gse152098_bigwig_regions.parquet"),
            "guides": str(OUT / "gse152098_bigwig_guides.parquet"),
        },
        "rationale": "Direct bigWig reprocessing provides replicate-resolved ATAC signal at PWS-ICR in hypothalamic neurons vs ESC.",
    }
    (OUT / "gse152098_reprocessing_report.json").write_text(
        json.dumps(report, indent=2, default=lambda x: None if (isinstance(x, float) and x != x) else x)
    )
    log.info("GSE152098 reprocessing complete -> %s", OUT)


if __name__ == "__main__":
    main()
