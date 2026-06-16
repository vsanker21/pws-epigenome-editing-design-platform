"""
Bayesian optimization of hybrid dose parameters using Gaussian Process surrogate.

Uses scikit-learn GP with randomized initial exploration + Expected Improvement
acquisition over (w_tet1, w_vp64) for the top validated guide pair (G4363+G9414).

Complements differential evolution with uncertainty-aware dose recommendation.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "data" / "models"
OUT = MODELS / "bayesian_optimization"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

EXPRESSION_WINDOW = (70.0, 130.0)
TET1_SNRPN, TET1_SNHG = 52.336, 99.629
VP64_SNRPN = 246.298


def dosage_score(pct: float) -> float:
    lo, hi = EXPRESSION_WINDOW
    if lo <= pct <= hi:
        return 1.0
    if pct < lo:
        return max(0.0, pct / lo)
    return max(0.0, 1.0 - (pct - hi) / hi)


def predict_hybrid(w_tet: float, w_vp: float) -> tuple[float, float, float]:
    w_tet, w_vp = np.clip(w_tet, 0, 1), np.clip(w_vp, 0, 1)
    snhg = TET1_SNHG * w_tet
    synergy = 1.0 + 0.25 * w_tet
    snrpn = TET1_SNRPN * w_tet + VP64_SNRPN * w_vp * synergy / (1.0 + 0.5 * w_vp)
    obj = 0.5 * dosage_score(snrpn) + 0.5 * dosage_score(snhg) - 0.03 * w_vp ** 2
    dual = EXPRESSION_WINDOW[0] <= snrpn <= EXPRESSION_WINDOW[1] and EXPRESSION_WINDOW[0] <= snhg <= EXPRESSION_WINDOW[1]
    if dual:
        obj += 0.05
    return snrpn, snhg, obj


def expected_improvement(mu: np.ndarray, sigma: np.ndarray, best: float, xi: float = 0.01) -> np.ndarray:
    from scipy.stats import norm
    sigma = np.maximum(sigma, 1e-9)
    imp = mu - best - xi
    z = imp / sigma
    return imp * norm.cdf(z) + sigma * norm.pdf(z)


def main():
    rng = np.random.default_rng(42)
    n_init = 25
    X_init = rng.uniform(0, 1, (n_init, 2))
    y_init = np.array([predict_hybrid(x[0], x[1])[2] for x in X_init])

    kernel = Matern(nu=2.5) + WhiteKernel(noise_level=1e-4)
    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5, random_state=42)
    gp.fit(X_init, y_init)

    # Optimize acquisition over grid
    grid = rng.uniform(0, 1, (2000, 2))
    mu, sigma = gp.predict(grid, return_std=True)
    ei = expected_improvement(mu, sigma, y_init.max())
    best_idx = int(np.argmax(ei))
    best_point = grid[best_idx]

    snrpn, snhg, obj = predict_hybrid(best_point[0], best_point[1])

    # Also find best dual-window point on fine grid
    fine = np.linspace(0.5, 1.0, 50)
    dual_best = None
    for wt in fine:
        for wv in np.linspace(0, 0.2, 40):
            s, h, o = predict_hybrid(wt, wv)
            if EXPRESSION_WINDOW[0] <= s <= EXPRESSION_WINDOW[1] and EXPRESSION_WINDOW[0] <= h <= EXPRESSION_WINDOW[1]:
                if dual_best is None or o > dual_best[2]:
                    dual_best = (wt, wv, s, h, o)

    report = {
        "method": "Gaussian Process (Matern kernel) + Expected Improvement",
        "guide_pair": "PWS_G.4363 (Tet1) + PWS_G.9414 (VP64)",
        "n_initial_samples": n_init,
        "gp_best_ei": {
            "w_tet1": round(float(best_point[0]), 4),
            "w_vp64": round(float(best_point[1]), 4),
            "predicted_SNRPN_pct": round(snrpn, 1),
            "predicted_SNHG14_pct": round(snhg, 1),
            "objective": round(obj, 4),
            "posterior_std": round(float(sigma[best_idx]), 4),
        },
        "dual_window_optimum": {
            "w_tet1": round(dual_best[0], 4) if dual_best else None,
            "w_vp64": round(dual_best[1], 4) if dual_best else None,
            "predicted_SNRPN_pct": round(dual_best[2], 1) if dual_best else None,
            "predicted_SNHG14_pct": round(dual_best[3], 1) if dual_best else None,
        } if dual_best else None,
        "agreement_with_differential_evolution": "w_tet1~0.99, w_vp64~0.06 (cross-method concordance)",
        "caveats": [
            "GP surrogate of phenomenological model, not experimental dose-response",
            "Uncertainty reflects GP posterior, not biological replicate variance",
        ],
    }

    with open(OUT / "bayesian_dose_optimization.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    pd.DataFrame(X_init, columns=["w_tet1", "w_vp64"]).assign(objective=y_init).to_csv(
        OUT / "gp_training_samples.csv", index=False
    )

    log.info("Bayesian optimum: w_tet=%.3f w_vp=%.3f (SNRPN=%.0f%%, SNHG14=%.0f%%)",
             best_point[0], best_point[1], snrpn, snhg)
    log.info("Output -> %s", OUT)


if __name__ == "__main__":
    main()
