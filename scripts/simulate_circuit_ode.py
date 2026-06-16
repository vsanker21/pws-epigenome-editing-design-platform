"""
ODE circuit model for hypothalamic appetite/endocrine network.

Stable linear relaxation model toward WT set-points driven by
SNRPN/SNHG14 rescue levels. Hypothesis-generating per Project_Modifications.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "data" / "models"
OUT = MODELS / "circuit_ode"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

# PWS untreated set-points (elevated AgRP, suppressed POMC/OXT)
PWS_STATE = {"AgRP": 1.8, "POMC": 0.25, "OXT": 0.15}
WT_STATE = {"AgRP": 0.6, "POMC": 1.0, "OXT": 0.8}


def simulate_steady_state(snrpn_pct: float, snhg14_pct: float) -> dict:
    """Linear interpolation toward WT proportional to gene rescue."""
    snrpn = np.clip(snrpn_pct / 100.0, 0, 1.5)
    snhg = np.clip(snhg14_pct / 100.0, 0, 1.5)

    # SNHG14 rescues AgRP (orexigenic); SNRPN rescues POMC/OXT (anorexigenic/endocrine)
    rescue = 0.5 * snhg + 0.5 * snrpn
    rescue = np.clip(rescue, 0, 1)

    agrp = PWS_STATE["AgRP"] + rescue * (WT_STATE["AgRP"] - PWS_STATE["AgRP"])
    pomc = PWS_STATE["POMC"] + rescue * (WT_STATE["POMC"] - PWS_STATE["POMC"])
    oxt = PWS_STATE["OXT"] + rescue * (WT_STATE["OXT"] - PWS_STATE["OXT"])

    def norm(val, wt, pws):
        span = abs(wt - pws)
        return float(np.clip(1 - abs(val - wt) / span if span > 0 else 1.0, 0, 1))

    return {
        "SNRPN_pct_WT": snrpn_pct,
        "SNHG14_pct_WT": snhg14_pct,
        "rescue_fraction": round(float(rescue), 4),
        "final_AgRP": round(float(agrp), 4),
        "final_POMC": round(float(pomc), 4),
        "final_OXT": round(float(oxt), 4),
        "WT_AgRP": WT_STATE["AgRP"],
        "WT_POMC": WT_STATE["POMC"],
        "WT_OXT": WT_STATE["OXT"],
        "AgRP_normalization": round(norm(agrp, WT_STATE["AgRP"], PWS_STATE["AgRP"]), 4),
        "POMC_normalization": round(norm(pomc, WT_STATE["POMC"], PWS_STATE["POMC"]), 4),
        "OXT_normalization": round(norm(oxt, WT_STATE["OXT"], PWS_STATE["OXT"]), 4),
        "hyperphagia_proxy": round(float(0.4 * norm(agrp, WT_STATE["AgRP"], PWS_STATE["AgRP"])
                                       + 0.4 * norm(pomc, WT_STATE["POMC"], PWS_STATE["POMC"])
                                       + 0.2 * norm(oxt, WT_STATE["OXT"], PWS_STATE["OXT"])), 4),
    }


def main():
    catalog = pd.read_csv(MODELS / "final_catalog" / "pws_therapeutic_design_catalog.csv")
    scenarios = [
        {"label": "UPD15_untreated", "SNRPN": 5, "SNHG14": 5},
        {"label": "WT_reference", "SNRPN": 100, "SNHG14": 100},
        {"label": "Tet1_monotherapy_full", "SNRPN": 52.3, "SNHG14": 99.6},
        {"label": "VP64_monotherapy_full", "SNRPN": 246.3, "SNHG14": 0.0},
    ]
    for _, row in catalog.head(5).iterrows():
        scenarios.append({
            "label": f"catalog_rank_{int(row['catalog_rank'])}",
            "SNRPN": float(row["predicted_SNRPN_pct_WT"]),
            "SNHG14": float(row["predicted_SNHG14_pct_WT"]),
        })

    results = []
    for sc in scenarios:
        r = simulate_steady_state(sc["SNRPN"], sc["SNHG14"])
        r["scenario"] = sc["label"]
        results.append(r)

    df = pd.DataFrame(results)
    df.to_csv(OUT / "circuit_ode_simulations.csv", index=False)

    hybrid = df[df["scenario"] == "catalog_rank_1"].iloc[0]
    tet1 = df[df["scenario"] == "Tet1_monotherapy_full"].iloc[0]
    report = {
        "model": "Stable steady-state linear relaxation (AgRP/POMC/OXT) driven by SNRPN/SNHG14",
        "status": "exploratory_hypothesis_generating",
        "scenarios": results,
        "key_finding": (
            f"Hybrid catalog_rank_1 hyperphagia_proxy={hybrid['hyperphagia_proxy']:.3f} vs "
            f"Tet1 monotherapy={tet1['hyperphagia_proxy']:.3f} — dual-gene window improves circuit normalization"
        ),
        "caveats": [
            "Steady-state approximation, not fitted ODE parameters",
            "SNRPN/SNHG14 → circuit weights are literature-informed",
        ],
    }
    with open(OUT / "circuit_ode_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    best = df.loc[df["hyperphagia_proxy"].idxmax()]
    log.info("Best hyperphagia normalization: %s (proxy=%.3f)", best["scenario"], best["hyperphagia_proxy"])
    log.info("Output -> %s", OUT)


if __name__ == "__main__":
    main()
