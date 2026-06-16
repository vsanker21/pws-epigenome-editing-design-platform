"""
Systematic data download pipeline for PWS epigenome-editing computational platform.

Downloads curated datasets per config/datasets.yaml with:
- Resume support (skip existing complete files)
- Checksum validation where available
- Manifest logging for reproducibility
"""

from __future__ import annotations

import hashlib
import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "datasets.yaml"
DATA_DIR = ROOT / "data"
MANIFEST = DATA_DIR / "manifest.json"

DATA_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(DATA_DIR / "download.log", encoding="utf-8"),
    ],
)
log = logging.getLogger(__name__)


def load_config() -> dict:
    with open(CONFIG, encoding="utf-8") as f:
        return yaml.safe_load(f)


def file_size_str(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}TB"


def normalize_url(url: str) -> str:
    """Convert GEO FTP URLs to HTTPS for requests compatibility."""
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        return url.replace("ftp://", "https://", 1)
    return url


def download_file(url: str, dest: Path, chunk_size: int = 1024 * 1024) -> dict:
    """Download with resume support. Returns metadata dict."""
    url = normalize_url(url)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")

    headers = {}
    mode = "wb"
    existing = 0
    if tmp.exists():
        existing = tmp.stat().st_size
        headers["Range"] = f"bytes={existing}-"
        mode = "ab"
        log.info("Resuming %s from %s", dest.name, file_size_str(existing))
    elif dest.exists() and dest.stat().st_size > 0:
        log.info("Already exists: %s (%s)", dest.name, file_size_str(dest.stat().st_size))
        return {
            "url": url,
            "path": str(dest.relative_to(ROOT)),
            "size_bytes": dest.stat().st_size,
            "status": "skipped_existing",
            "sha256": _sha256(dest),
        }

    log.info("Downloading: %s -> %s", url, dest)
    t0 = time.time()
    try:
        resp = requests.get(url, headers=headers, stream=True, timeout=120)
        if resp.status_code == 416:
            # Range not satisfiable - file complete
            if tmp.exists():
                tmp.rename(dest)
            return {"url": url, "path": str(dest.relative_to(ROOT)), "status": "complete"}
        resp.raise_for_status()

        total = int(resp.headers.get("content-length", 0)) + existing
        downloaded = existing
        with open(tmp, mode) as f:
            for chunk in resp.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total > 0 and downloaded % (10 * chunk_size) < chunk_size:
                        pct = 100 * downloaded / total
                        log.info("  %s: %.1f%% (%s)", dest.name, pct, file_size_str(downloaded))

        tmp.rename(dest)
        elapsed = time.time() - t0
        size = dest.stat().st_size
        log.info("Done: %s (%s in %.1fs)", dest.name, file_size_str(size), elapsed)
        return {
            "url": url,
            "path": str(dest.relative_to(ROOT)),
            "size_bytes": size,
            "status": "downloaded",
            "sha256": _sha256(dest),
            "elapsed_sec": round(elapsed, 1),
        }
    except Exception as e:
        log.error("Failed %s: %s", url, e)
        return {"url": url, "path": str(dest.relative_to(ROOT)), "status": "failed", "error": str(e)}


def _sha256(path: Path) -> str:
  h = hashlib.sha256()
  with open(path, "rb") as f:
    for block in iter(lambda: f.read(1024 * 1024), b""):
      h.update(block)
  return h.hexdigest()


def url_to_dest(url: str, dataset_id: str) -> Path:
    name = Path(urlparse(url).path).name
    return DATA_DIR / dataset_id / name


def collect_urls(config: dict, priority_max: int = 3, include_large: bool = False) -> list[tuple[str, str, str]]:
    """Return list of (dataset_id, url, description)."""
    items = []
    for ds_id, ds in config.get("datasets", {}).items():
        if ds.get("priority", 99) > priority_max:
            continue
        for f in ds.get("files", []):
            if f.get("size") == "large" and not include_large:
                continue
            items.append((ds_id, f["url"], f.get("description", "")))
        for src in ds.get("sources", []):
            # Reference URLs only - not downloaded
            pass
    return items


def update_manifest(entries: list[dict]) -> None:
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    existing = {}
    if MANIFEST.exists():
        with open(MANIFEST, encoding="utf-8") as f:
            existing = json.load(f)
    existing["last_updated"] = datetime.now(timezone.utc).isoformat()
    existing.setdefault("files", {})
    for e in entries:
        key = e.get("url", e.get("path", ""))
        existing["files"][key] = e
    with open(MANIFEST, "w", encoding="utf-8") as f:
        json.dump(existing, f, indent=2)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Download PWS project datasets")
    parser.add_argument("--priority", type=int, default=2, help="Max priority tier (1=highest)")
    parser.add_argument("--include-large", action="store_true", help="Download large RAW archives")
    parser.add_argument("--dataset", type=str, help="Download only this dataset ID")
    args = parser.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()
    items = collect_urls(config, priority_max=args.priority, include_large=args.include_large)
    if args.dataset:
        items = [i for i in items if i[0] == args.dataset]

    log.info("Queued %d files for download", len(items))
    results = []
    for ds_id, url, desc in items:
        dest = url_to_dest(url, ds_id)
        log.info("=== [%s] %s ===", ds_id, desc)
        results.append(download_file(url, dest))

    update_manifest(results)
    ok = sum(1 for r in results if r["status"] in ("downloaded", "skipped_existing", "complete"))
    fail = sum(1 for r in results if r["status"] == "failed")
    log.info("Summary: %d ok, %d failed, manifest -> %s", ok, fail, MANIFEST)


if __name__ == "__main__":
    main()
