"""
Generate complete data inventory report after all downloads.

Summarizes every file in data/ with size, dataset group, and download status.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
INVENTORY = DATA / "inventory.json"


def human_size(n: int) -> str:
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} TB"


def main():
    files = []
    total = 0
    by_dataset: dict[str, dict] = {}

    for path in sorted(DATA.rglob("*")):
        if not path.is_file():
            continue
        if path.name in ("download.log",):
            continue
        size = path.stat().st_size
        total += size
        rel = str(path.relative_to(DATA))
        parts = rel.split("\\")
        ds = parts[0] if parts else "root"

        entry = {
            "path": rel.replace("\\", "/"),
            "size_bytes": size,
            "size_human": human_size(size),
        }
        files.append(entry)

        if ds not in by_dataset:
            by_dataset[ds] = {"file_count": 0, "total_bytes": 0, "files": []}
        by_dataset[ds]["file_count"] += 1
        by_dataset[ds]["total_bytes"] += size
        by_dataset[ds]["files"].append(path.name)

    inventory = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "total_files": len(files),
        "total_bytes": total,
        "total_human": human_size(total),
        "datasets": {
            k: {**v, "total_human": human_size(v["total_bytes"])}
            for k, v in sorted(by_dataset.items())
        },
        "files": files,
    }

    with open(INVENTORY, "w", encoding="utf-8") as f:
        json.dump(inventory, f, indent=2)

    print(f"Inventory: {len(files)} files, {human_size(total)} total")
    print(f"Written -> {INVENTORY}")
    for ds, info in inventory["datasets"].items():
        print(f"  {ds}: {info['file_count']} files, {info['total_human']}")


if __name__ == "__main__":
    main()
