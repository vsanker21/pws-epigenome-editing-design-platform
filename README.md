# PWS Epigenome-Editing Design Platform

Computational design-prioritization platform for CRISPR epigenome-editing therapy in **maternal UPD15 Prader–Willi syndrome** (PWS).

Manuscript target: *Nucleic Acids Research* — Methods  
Title: *An allele-aware computational design platform for CRISPR epigenome-editing therapy in maternal UPD15 Prader–Willi syndrome*

## Overview

This repository integrates public multi-omic datasets (Rohm 2025, Nemoto 2025, Cousminer 2021, ENCODE, CRISPRepi) into an allele-aware GRCh38 PWS-locus digital twin, applies editor-stratified forward models, optimizes hybrid Tet1+VP64 therapeutic designs, and exports a ranked uncertainty-aware catalog with retrospective validation.

**Top recommended design:** Hybrid PWS_G.4670 (Tet1) + PWS_G.9414 (VP64) at w_Tet1≈0.987, w_VP64≈0.062 — predicted SNRPN 70.0% WT, SNHG14 98.3% WT.

## Repository contents

| Path | Description |
|------|-------------|
| `scripts/` | Full analysis pipeline (phases 3–12) |
| `scripts/run_pipeline.py` | Master orchestrator |
| `data/models/` | Processed outputs, validation reports, therapeutic catalog |
| `data/integrated/` | hg38-merged locus data |
| `data/curated/` | Curated summaries and probe mappings |
| `data/encode_reference/` | ENCODE regulatory context (UCSC API) |
| `manuscript/` | NAR Methods manuscript and figures |

Raw GEO files (~10+ GB) are **not** included; download with `scripts/download_data.py` (see below).

## Requirements

- Python 3.10+
- See `requirements.txt`

```bash
pip install -r requirements.txt
```

Optional: WSL + UCSC `bigWigSummary` for GSE152098 reprocessing (Phase 12).

## Quick start (reproduce processed outputs)

```bash
# Run full pipeline (requires GEO data downloaded first)
python scripts/run_pipeline.py --phase all

# Or run individual phases:
python scripts/run_pipeline.py --phase 3    # curation + hg38 integration
python scripts/run_pipeline.py --phase 5    # forward model
python scripts/run_pipeline.py --phase 6    # hybrid optimization
python scripts/run_pipeline.py --phase 9    # validation suite
python scripts/run_pipeline.py --phase 10   # accessibility + protocol
python scripts/run_pipeline.py --phase 11   # Pareto, Bayesian GP, sensitivity
python scripts/run_phase12.py               # Cas-OFFinder, ENCODE, bigWig, DMR
```

## Key outputs

- `data/models/final_catalog/pws_therapeutic_design_catalog.csv` — 25 ranked designs
- `data/models/final_catalog/project_synthesis_report.json` — pipeline audit
- `data/models/validation/held_out_benchmark.json` — guide recovery (3/4 pass)
- `data/models/experimental_protocol/experimental_validation_protocol.json` — wet-lab protocol

## Data sources (public)

| Accession | Reference | Use |
|-----------|-----------|-----|
| GSE285306 | Rohm et al. 2025 Cell Genomics | gRNA screen |
| GSE243185 | Rohm 2025 | Bulk RNA |
| GSE285300 | Rohm 2025 | Bisulfite |
| GSE262700 | Nemoto 2025 Nat Commun | Organoid scRNA |
| GSE285305 | Rohm 2025 | Neuron RNA |
| GSE152098 | Cousminer 2021 | Hypothalamic ATAC |
| GSE28525, GSE298378 | Imprinting methylation | 450K DMR mapping |
| CRISPRepi | Shi et al. 2025 NAR | Editor priors |

## Citation

If you use this code, please cite the accompanying NAR Methods manuscript and this repository:

- **GitHub:** https://github.com/vsanker21/pws-epigenome-editing-design-platform
- **Zenodo:** https://doi.org/10.5281/zenodo.XXXXXXX *(update after Zenodo archives release v1.0.0)*

## License

MIT License — see LICENSE.

## Contact

[vsanker21](https://github.com/vsanker21)
