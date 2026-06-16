"""
Discover and download all GEO supplementary files for listed accessions.

Queries NCBI GEO metadata for each accession and downloads every
supplementary_file entry. Skips duplicates already in manifest.
"""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path

import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "datasets.yaml"
DATA_DIR = ROOT / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Reuse download infrastructure
sys.path.insert(0, str(Path(__file__).resolve().parent))
from download_data import download_file, normalize_url, update_manifest

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(DATA_DIR / "download.log", encoding="utf-8"),
    ],
)
log = logging.getLogger(__name__)

GEO_META_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


def fetch_supplementary_urls(accession: str) -> list[str]:
    """Parse GEO SOFT format for supplementary_file lines."""
    r = requests.get(
        GEO_META_URL,
        params={"acc": accession, "targ": "self", "form": "text", "view": "quick"},
        timeout=60,
    )
    r.raise_for_status()
    urls = []
    for line in r.text.splitlines():
        if line.startswith("!Series_supplementary_file"):
            match = re.search(r"=\s*(ftp://\S+|https://\S+)", line)
            if match:
                urls.append(match.group(1).strip())
    return urls


def accession_output_dir(accession: str) -> Path:
    return DATA_DIR / accession.lower()


def load_geo_accessions(include_raw: bool = False) -> list[str]:
    with open(CONFIG, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    accessions = set()
    for ds_id, ds in cfg.get("datasets", {}).items():
        if ds.get("geo_accession"):
            accessions.add(ds["geo_accession"])
        for sub in ds.get("subseries", []):
            accessions.add(sub)
        for extra in ds.get("extra_geo", []):
            accessions.add(extra)

    return sorted(accessions)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Download all GEO supplementary files")
    parser.add_argument("--include-raw", action="store_true", help="Include *_RAW.tar archives")
    parser.add_argument("--accession", action="append", help="Limit to specific GEO accession(s)")
    args = parser.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    accessions = args.accession or load_geo_accessions()
    results = []
    seen_urls: set[str] = set()

    for acc in accessions:
        log.info("Discovering supplementary files for %s", acc)
        try:
            urls = fetch_supplementary_urls(acc)
        except Exception as e:
            log.error("Failed to query %s: %s", acc, e)
            continue

        if not urls:
            log.warning("No supplementary files found for %s", acc)
            continue

        log.info("  Found %d file(s)", len(urls))
        out_dir = accession_output_dir(acc)

        for url in urls:
            norm = normalize_url(url)
            if norm in seen_urls:
                log.info("  Skip duplicate URL: %s", Path(url).name)
                continue
            if not args.include_raw and "_RAW.tar" in url:
                log.info("  Skip RAW archive (use --include-raw): %s", Path(url).name)
                continue
            seen_urls.add(norm)
            dest = out_dir / Path(urlparse_name(url))
            results.append(download_file(url, dest))

    update_manifest(results)
    ok = sum(1 for r in results if r["status"] in ("downloaded", "skipped_existing", "complete"))
    fail = sum(1 for r in results if r["status"] == "failed")
    log.info("GEO discovery complete: %d ok, %d failed", ok, fail)


def urlparse_name(url: str) -> str:
    from urllib.parse import urlparse
    return Path(urlparse(url).path).name


if __name__ == "__main__":
    main()
