# Release v1.0.0 — NAR Methods manuscript deposit

First public release accompanying the *Nucleic Acids Research* Methods submission:

**An allele-aware computational design platform for CRISPR epigenome-editing therapy in maternal UPD15 Prader-Willi syndrome**

## Contents

- Full analysis pipeline (`scripts/run_pipeline.py`, phases 3–12)
- Processed models, validation outputs, and therapeutic design catalog
- Manuscript figure assets (`manuscript/figures/`)
- Configuration and requirements for reproduction

## Reproduce

```bash
pip install -r requirements.txt
python scripts/download_data.py   # fetch public GEO sources (~10+ GB)
python scripts/run_pipeline.py --phase all
```

## Key outputs

- `data/models/final_catalog/pws_therapeutic_design_catalog.csv`
- `data/models/validation/held_out_benchmark.json`
- `data/models/final_catalog/project_synthesis_report.json`

## Citation

Cite the GitHub repository and Zenodo DOI (assigned upon Zenodo archive of this release).
