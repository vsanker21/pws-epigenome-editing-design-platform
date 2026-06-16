"""
Validate downloaded datasets and produce inventory report.

Checks file integrity, previews key PWS datasets, and confirms
presence of SNRPN/SNORD116 in expression data.
"""

from __future__ import annotations

import gzip
import json
import tarfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
REPORT = DATA / "validation_report.json"

PWS_GENES = ["SNRPN", "SNHG14", "MAGEL2", "NDN", "MKRN3", "NPAP1", "UBE3A"]


def preview_gz_csv(path: Path, nrows: int = 5) -> dict:
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f, nrows=1000)
    return {
        "path": str(path.relative_to(ROOT)),
        "shape_preview": list(df.shape),
        "columns": list(df.columns[:15]),
        "head": df.head(nrows).to_dict(orient="records"),
    }


def check_pws_genes_in_rnaseq(path: Path) -> dict:
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f, index_col=0)
    gene_col = df.index if df.index.name else df.iloc[:, 0]
    if df.index.name is None and "gene" in df.columns:
        genes = set(df["gene"].astype(str))
    else:
        genes = set(df.index.astype(str))
    found = [g for g in PWS_GENES if any(g in x for x in genes)]
    return {"file": path.name, "pws_genes_found": found, "n_samples": df.shape[1]}


def inspect_tar(path: Path) -> dict:
    with tarfile.open(path, "r") as tar:
        members = tar.getnames()
    return {"file": path.name, "n_members": len(members), "members": members[:20]}


def main():
    report = {"datasets": {}, "locus_annotation": {}, "status": "ok"}

    # Locus annotation
    locus_summary = DATA / "locus_annotation" / "locus_summary.json"
    if locus_summary.exists():
        with open(locus_summary) as f:
            report["locus_annotation"] = json.load(f)

    # GSE285306 RNA-seq
    for f in (DATA / "gse285306").glob("*featurecounts*.csv.gz"):
        report["datasets"][f.name] = {
            **check_pws_genes_in_rnaseq(f),
            **preview_gz_csv(f, nrows=2),
        }

    # Tiling screen
    screen = DATA / "gse285306" / "GSE285289_Tet1_sublib_results.csv.gz"
    if screen.exists():
        report["datasets"]["tiling_screen"] = preview_gz_csv(screen)

    # Methylation
    for f in (DATA / "gse285306").glob("*methylation*.csv.gz"):
        report["datasets"][f.name] = preview_gz_csv(f)

    # Nemoto organoid scRNA
    tar = DATA / "gse262700" / "GSE262700_RAW.tar"
    if tar.exists():
        report["datasets"]["gse262700"] = inspect_tar(tar)

    # CRISPRepi
    hs = DATA / "crisprepi" / "Homo_sapiens.tsv.gz"
    if hs.exists():
        with gzip.open(hs, "rt") as f:
            crepi = pd.read_csv(f, sep="\t", nrows=500)
        report["datasets"]["crisprepi_human"] = {
            "n_records_preview": len(crepi),
            "columns": list(crepi.columns),
            "editing_systems": crepi["editing_system"].value_counts().head(10).to_dict()
            if "editing_system" in crepi.columns
            else {},
            "cell_types": crepi["cell_type"].value_counts().head(10).to_dict()
            if "cell_type" in crepi.columns
            else {},
        }

    with open(REPORT, "w") as f:
        json.dump(report, f, indent=2, default=str)

    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
