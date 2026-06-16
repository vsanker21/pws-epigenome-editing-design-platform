"""Poll Zenodo for the archived GitHub release DOI and update manuscript constants."""
from __future__ import annotations

import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_DATA = ROOT / "scripts" / "manuscript_data.py"
ZENODO_QUERY = "pws-epigenome-editing-design-platform"


def fetch_zenodo_doi() -> str | None:
    url = (
        "https://zenodo.org/api/records/?"
        + urllib.parse.urlencode(
            {"q": ZENODO_QUERY, "all_versions": "true", "size": 25, "sort": "mostrecent"}
        )
    )
    with urllib.request.urlopen(url, timeout=30) as resp:
        payload = json.load(resp)

    for hit in payload.get("hits", {}).get("hits", []):
        metadata = hit.get("metadata", {})
        related = metadata.get("related_identifiers", [])
        for item in related:
            identifier = item.get("identifier", "")
            if "pws-epigenome-editing-design-platform" in identifier:
                doi = metadata.get("doi") or hit.get("doi")
                if doi:
                    return doi
        title = metadata.get("title", "").lower()
        if "pws" in title and "epigenome" in title:
            doi = metadata.get("doi") or hit.get("doi")
            if doi:
                return doi
    return None


def update_manuscript_data(doi: str) -> None:
    text = MANUSCRIPT_DATA.read_text(encoding="utf-8")
    new_text, n = re.subn(
        r'ZENODO_DOI = "10\.5281/zenodo\.[^"]*"',
        f'ZENODO_DOI = "{doi}"',
        text,
        count=1,
    )
    if n != 1:
        raise RuntimeError("Could not update ZENODO_DOI in manuscript_data.py")
    MANUSCRIPT_DATA.write_text(new_text, encoding="utf-8")

    readme = ROOT / "README.md"
    readme_text = readme.read_text(encoding="utf-8")
    readme_text = readme_text.replace(
        "https://doi.org/10.5281/zenodo.XXXXXXX *(update after Zenodo archives release v1.0.0)*",
        f"https://doi.org/{doi}",
    )
    readme.write_text(readme_text, encoding="utf-8")


def main() -> int:
    doi = fetch_zenodo_doi()
    if not doi:
        print(
            "No Zenodo record found yet. Enable GitHub integration at "
            "https://zenodo.org/account/settings/github/ and ensure release v1.0.0 is archived."
        )
        return 1

    print(f"Found Zenodo DOI: {doi}")
    update_manuscript_data(doi)
    print(f"Updated {MANUSCRIPT_DATA.name} and README.md")

    import subprocess

    subprocess.run([sys.executable, str(ROOT / "scripts" / "generate_nar_manuscript.py")], check=True)
    print("Regenerated manuscript with Zenodo DOI.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
