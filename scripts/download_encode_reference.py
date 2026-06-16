"""
Download ENCODE / Roadmap regulatory context for chr15 PWS locus via UCSC REST API.

Fixes prior failures by using correct hg38 track internal names discovered via
/list/tracks (legacy wgEncode* container names are not queryable directly).

Scientific rationale: neuronal DNase/H3K27ac/H3K4me3/CTCF at PWS locus informs
whether epigenome-editing targets lie in active regulatory chromatin in relevant
cell types (brain, embryo/iPSC proxy).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]
LOCUS = ROOT / "config" / "locus.yaml"
OUT = ROOT / "data" / "encode_reference"
UCSC = "https://api.genome.ucsc.edu"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

# Leaf tracks with data at hg38 (validated)
TRACKS = [
    ("encRegTfbsClustered", "ENCODE TFBS clusters (CTCF etc.)"),
    ("wgEncodeReg4DnaseAllBrain", "ENCODE DNase-seq — brain"),
    ("wgEncodeReg4DnaseAllEmbryo", "ENCODE DNase-seq — embryo (iPSC proxy)"),
    ("wgEncodeReg4MarkH3k27acAllBrain", "H3K27ac — brain"),
    ("wgEncodeReg4MarkH3k27acAllEmbryo", "H3K27ac — embryo"),
    ("wgEncodeReg4MarkH3k4me3AllBrain", "H3K4me3 — brain"),
    ("wgEncodeReg4MarkH3k4me3AllEmbryo", "H3K4me3 — embryo"),
    ("wgEncodeReg4MarkCtcfAllBrain", "CTCF ChIP-seq — brain"),
]


def load_locus() -> tuple[str, int, int]:
    with open(LOCUS, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    r = cfg["regions"]["pws_critical"]
    return cfg["chromosome"], int(r["start"]), int(r["end"])


def fetch_track(track: str, chrom: str, start: int, end: int, max_items: int = 5000) -> dict:
    params = {
        "genome": "hg38",
        "track": track,
        "chrom": chrom,
        "start": start,
        "end": end,
        "maxItemsOutput": max_items,
    }
    log.info("Fetching %s (%s:%d-%d)", track, chrom, start, end)
    r = requests.get(f"{UCSC}/getData/track", params=params, timeout=180)
    if r.status_code not in (200, 206):
        return {"error": r.text, "status_code": r.status_code}
    data = r.json()
    # Track key is usually the track name
    items = data.get(track, [])
    if not isinstance(items, list):
        items = []
    return {"n_items": len(items), "status_code": r.status_code, "items": items}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    chrom, start, end = load_locus()
    summary = {
        "region": f"{chrom}:{start}-{end}",
        "genome_build": "GRCh38",
        "api": UCSC,
        "tracks": {},
    }

    for track, desc in TRACKS:
        try:
            result = fetch_track(track, chrom, start, end)
            n = result.get("n_items", 0)
            summary["tracks"][track] = {
                "description": desc,
                "n_items": n,
                "status_code": result.get("status_code"),
                "error": result.get("error"),
            }
            out_path = OUT / f"{track}.json"
            payload = result.get("items", [])
            out_path.write_text(json.dumps(payload, indent=2))
            log.info("  %s: %d items -> %s", track, n, out_path.name)
        except Exception as e:
            log.warning("Track %s failed: %s", track, e)
            summary["tracks"][track] = {"description": desc, "error": str(e)}

    n_ok = sum(1 for t in summary["tracks"].values() if t.get("n_items", 0) > 0)
    summary["n_tracks_success"] = n_ok
    summary["overall_success"] = n_ok >= 4
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    log.info("ENCODE reference: %d/%d tracks with data -> %s", n_ok, len(TRACKS), OUT)


if __name__ == "__main__":
    main()
