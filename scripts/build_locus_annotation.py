"""
Build annotated regulatory map of 15q11-q13 PWS locus (GRCh38).

Fetches from UCSC Table Browser API and Ensembl REST:
- Gene annotations
- CpG islands
- Known imprinting control elements
- CTCF binding sites (ENCODE)
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]
LOCUS_CONFIG = ROOT / "config" / "locus.yaml"
OUT_DIR = ROOT / "data" / "locus_annotation"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

UCSC_API = "https://api.genome.ucsc.edu"


def load_locus() -> dict:
    with open(LOCUS_CONFIG, encoding="utf-8") as f:
        return yaml.safe_load(f)


def ucsc_track(chrom: str, start: int, end: int, track: str, params: dict | None = None) -> list[dict]:
    """Query UCSC Genome Browser API for a track region."""
    url = f"{UCSC_API}/getData/track"
    q = {
        "genome": "hg38",
        "track": track,
        "chrom": chrom,
        "start": start,
        "end": end,
        **(params or {}),
    }
    log.info("UCSC query: track=%s region=%s:%d-%d", track, chrom, start, end)
    r = requests.get(url, params=q, timeout=60)
    r.raise_for_status()
    data = r.json()
    return data.get(track, data.get("cpgIslandExt", []))


def ucsc_refgene(chrom: str, start: int, end: int) -> list[dict]:
    """Fetch gene annotations from UCSC refGene track."""
    url = f"{UCSC_API}/getData/track"
    q = {"genome": "hg38", "track": "refGene", "chrom": chrom, "start": start, "end": end}
    log.info("UCSC refGene query: %s:%d-%d", chrom, start, end)
    r = requests.get(url, params=q, timeout=120)
    r.raise_for_status()
    return r.json().get("refGene", [])


def ensembl_genes(chrom: str, start: int, end: int, max_span: int = 500_000) -> list[dict]:
    """Fetch genes via Ensembl REST; chunks large regions to avoid API limits."""
    ens_chrom = chrom.replace("chr", "")
    genes = []
    for chunk_start in range(start, end, max_span):
        chunk_end = min(chunk_start + max_span, end)
        region = f"{ens_chrom}:{chunk_start + 1}-{chunk_end}:1"
        url = f"https://rest.ensembl.org/overlap/region/human/{region}"
        params = {"feature": "gene", "content-type": "application/json"}
        log.info("Ensembl gene query: %s", region)
        r = requests.get(url, params=params, headers={"Content-Type": "application/json"}, timeout=60)
        r.raise_for_status()
        genes.extend(r.json())
    return genes


def fetch_jaspar_motifs() -> dict:
    """ZNF274 and CTCF are key repressors at PWS-ICR."""
    motifs = {}
    for motif_id, name in [("MA1145.1", "ZNF274"), ("MA0139.1", "CTCF"), ("MA0506.1", "G9a/EHMT2_target")]:
        try:
            url = f"https://jaspar.genereg.net/api/v1/matrix/{motif_id}/"
            r = requests.get(url, timeout=30)
            if r.ok:
                motifs[name] = r.json()
        except Exception as e:
            log.warning("JASPAR fetch %s failed: %s", name, e)
    return motifs


def write_bed(records: list[dict], path: Path, fields: list[str]):
    """Write simple BED file from record list."""
    with open(path, "w", encoding="utf-8") as f:
        for rec in records:
            if "txStart" in rec:
                chrom = rec["chrom"]
                start, end = rec["txStart"], rec["txEnd"]
                name = rec.get("name2", rec.get("name", "."))
            elif "chromStart" in rec:
                chrom, start, end = rec["chrom"], rec["chromStart"], rec["chromEnd"]
                name = rec.get("name", ".")
            else:
                chrom = rec.get("seq_region_name", "chr15")
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"
                start = rec.get("start", 0) - 1
                end = rec.get("end", 0)
                name = rec.get("external_name", rec.get("id", "."))
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    locus = load_locus()
    chrom = locus["chromosome"]
    region = locus["regions"]["pws_critical"]
    start, end = region["start"], region["end"]

    annotation = {
        "genome_build": locus["genome_build"],
        "region": f"{chrom}:{start}-{end}",
        "fetched_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
        "subregions": locus["regions"],
        "key_genes": locus["key_genes"],
    }

    # Genes (UCSC refGene primary; Ensembl supplementary)
    genes = ucsc_refgene(chrom, start, end)
    annotation["gene_count"] = len(genes)
    with open(OUT_DIR / "genes_ucsc.json", "w", encoding="utf-8") as f:
        json.dump(genes, f, indent=2)
    write_bed(genes, OUT_DIR / "genes.bed", ["chrom", "chromStart", "chromEnd", "name"])
    log.info("Wrote %d UCSC refGene entries", len(genes))

    try:
        ens_genes = ensembl_genes(chrom, start, end)
        with open(OUT_DIR / "genes_ensembl.json", "w", encoding="utf-8") as f:
            json.dump(ens_genes, f, indent=2)
        log.info("Wrote %d Ensembl genes", len(ens_genes))
    except Exception as e:
        log.warning("Ensembl gene fetch failed: %s", e)

    # CpG islands
    try:
        cpgs = ucsc_track(chrom, start, end, "cpgIslandExt")
        annotation["cpg_island_count"] = len(cpgs)
        with open(OUT_DIR / "cpg_islands.json", "w", encoding="utf-8") as f:
            json.dump(cpgs, f, indent=2)
        write_bed(cpgs, OUT_DIR / "cpg_islands.bed", ["chrom", "chromStart", "chromEnd", "name"])
        log.info("Wrote %d CpG islands", len(cpgs))
    except Exception as e:
        log.warning("CpG island fetch failed: %s", e)

    # CTCF ENCODE peaks (wgEncodeRegTfbsClusteredV3)
    try:
        ctcf = ucsc_track(chrom, start, end, "wgEncodeRegTfbsClusteredV3", {"factor": "CTCF"})
        annotation["ctcf_peak_count"] = len(ctcf)
        with open(OUT_DIR / "ctcf_peaks.json", "w", encoding="utf-8") as f:
            json.dump(ctcf, f, indent=2)
        log.info("Wrote %d CTCF peaks", len(ctcf))
    except Exception as e:
        log.warning("CTCF fetch failed: %s", e)

    # Subregion BED for PWS-ICR, SNORD116
    with open(OUT_DIR / "pws_subregions.bed", "w", encoding="utf-8") as f:
        for name, reg in locus["regions"].items():
            f.write(f"{chrom}\t{reg['start']}\t{reg['end']}\t{name}\n")

    with open(OUT_DIR / "locus_summary.json", "w", encoding="utf-8") as f:
        json.dump(annotation, f, indent=2)

    log.info("Locus annotation complete -> %s", OUT_DIR)


if __name__ == "__main__":
    main()
