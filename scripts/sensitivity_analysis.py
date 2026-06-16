"""
Phase 11b: Robustness / sensitivity analysis for top therapeutic designs.

Purpose
-------
Turn heuristic "uncertainty" into an actionable probability that a design lands
in the desired dual-gene paternal-equivalent window (70–130% WT), under a simple
stochastic perturbation model.

This supports publication-quality claims like:
  - "Top design has X% probability of meeting both gene targets under model uncertainty"
  - "Ranking is stable / unstable under uncertainty"

Model
-----
For each design, we perturb predicted gene outcomes with Gaussian noise:
  SNRPN'  ~ Normal(mu = predicted_SNRPN_pct_WT,  sigma = K * uncertainty)
  SNHG14' ~ Normal(mu = predicted_SNHG14_pct_WT, sigma = K * uncertainty)

Where uncertainty comes from Phase 6/7 (guide confidence heuristics).

This is a pragmatic sensitivity layer; it does not claim mechanistic fidelity.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "data" / "models"
OUT = MODELS / "sensitivity"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

WINDOW = (70.0, 130.0)


def in_window(x: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return (x >= lo) & (x <= hi)


def main():
    catalog_path = MODELS / "final_catalog" / "pws_therapeutic_design_catalog.csv"
    df = pd.read_csv(catalog_path).copy()

    # Analyze top N only (publication-oriented).
    top_n = min(25, len(df))
    df = df.sort_values("catalog_rank").head(top_n).reset_index(drop=True)

    # Monte Carlo settings
    n_samples = 50_000
    # Scale: convert heuristic uncertainty (~0.12–0.18) into %-WT sigma.
    # Chosen so typical designs have ~3–5% WT 1-sigma, matching "dose uncertainty" scale.
    K = 25.0

    lo, hi = WINDOW
    results = []
    rng = np.random.default_rng(7)

    for _, r in df.iterrows():
        mu_snrpn = float(r["predicted_SNRPN_pct_WT"])
        mu_snhg14 = float(r["predicted_SNHG14_pct_WT"])
        unc = float(r.get("uncertainty", 0.15))
        sigma = max(1e-6, K * unc)

        snrpn = rng.normal(mu_snrpn, sigma, size=n_samples)
        snhg14 = rng.normal(mu_snhg14, sigma, size=n_samples)

        p_snrpn = float(np.mean(in_window(snrpn, lo, hi)))
        p_snhg14 = float(np.mean(in_window(snhg14, lo, hi)))
        p_both = float(np.mean(in_window(snrpn, lo, hi) & in_window(snhg14, lo, hi)))

        results.append(
            {
                "catalog_rank": int(r["catalog_rank"]),
                "strategy": str(r["strategy"]),
                "tet1_grna_id": r.get("tet1_grna_id"),
                "vp64_grna_id": r.get("vp64_grna_id"),
                "uncertainty": unc,
                "sigma_pct": sigma,
                "p_SNRPN_in_window": p_snrpn,
                "p_SNHG14_in_window": p_snhg14,
                "p_both_in_window": p_both,
            }
        )

    out_df = pd.DataFrame(results).sort_values("catalog_rank")
    out_df.to_csv(OUT / "catalog_window_probability.csv", index=False)

    report = {
        "analysis": "window_probability_monte_carlo",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "window": {"SNRPN_pct_WT": [lo, hi], "SNHG14_pct_WT": [lo, hi]},
        "n_designs": int(top_n),
        "n_samples": int(n_samples),
        "sigma_model": {"sigma_pct": f"{K} * uncertainty", "K": K},
        "top_design": out_df.iloc[0].to_dict() if len(out_df) else None,
        "top5": out_df.head(5).to_dict(orient="records"),
        "limitations": [
            "Assumes Gaussian perturbations and shared sigma across genes per design",
            "Uncertainty is heuristic (guide confidence), not a calibrated posterior",
            "Does not include delivery variability or cell-type heterogeneity explicitly",
        ],
    }
    (OUT / "sensitivity_report.json").write_text(json.dumps(report, indent=2, default=str))
    log.info("Sensitivity outputs -> %s", OUT)
    if len(out_df):
        log.info(
            "Top design p(both in window)=%.3f (sigma=%.2f%%WT)",
            float(out_df.iloc[0]["p_both_in_window"]),
            float(out_df.iloc[0]["sigma_pct"]),
        )


if __name__ == "__main__":
    main()

