"""
Graph-propagation forward model on the integrated PWS locus digital twin.

Propagates perturbation scores from gRNA target nodes through the regulatory
graph (Capture-C loops, proximity edges) using personalized PageRank-style
diffusion. This is the computationally feasible GNN alternative when
training data is sparse — per Project_Modifications, report uncertainty.

Scientific rationale: editing at PWS-ICR should propagate activation to
SNRPN/SNHG14 gene nodes via known promoter loops (Cousminer 2021 Capture-C).
"""

from __future__ import annotations

import json
import logging
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
INTEGRATED = ROOT / "data" / "integrated"
MODELS = ROOT / "data" / "models"
OUT = MODELS / "graph_model"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

PWS_GENES = ["SNRPN", "SNHG14", "MAGEL2", "NDN", "MKRN3", "UBE3A", "ATP10A"]
ALPHA = 0.85  # propagation damping


def build_adjacency(nodes: pd.DataFrame, edges: pd.DataFrame) -> tuple[list[str], sparse.csr_matrix]:
    node_ids = nodes["node_id"].tolist()
    idx = {n: i for i, n in enumerate(node_ids)}
    rows, cols, weights = [], [], []

    for _, e in edges.iterrows():
        s, t = e["source"], e["target"]
        if s not in idx or t not in idx:
            continue
        w = float(e.get("score", 1.0) or 1.0)
        if np.isnan(w):
            w = 1.0
        if e.get("edge_type") == "capture_c_loop":
            w *= 2.0
        rows.extend([idx[s], idx[t]])
        cols.extend([idx[t], idx[s]])
        weights.extend([w, w])

    n = len(node_ids)
    mat = sparse.csr_matrix((weights, (rows, cols)), shape=(n, n))
    # Symmetrize for undirected regulatory propagation
    mat = mat + mat.T
    mat.data = np.clip(mat.data, 0, None)
    row_sum = np.array(mat.sum(axis=1)).flatten()
    row_sum[row_sum == 0] = 1
    mat = sparse.diags(1.0 / row_sum) @ mat
    return node_ids, mat


def propagate(adj: sparse.csr_matrix, source_idx: list[int], source_weights: list[float], n_nodes: int) -> np.ndarray:
    p0 = np.zeros(n_nodes)
    for i, w in zip(source_idx, source_weights):
        p0[i] += w
    if p0.sum() > 0:
        p0 /= p0.sum()

    p = p0.copy()
    for _ in range(50):
        p_new = (1 - ALPHA) * p0 + ALPHA * (adj @ p)
        if np.linalg.norm(p_new - p) < 1e-8:
            break
        p = p_new
    pr = np.nan_to_num(p, nan=0.0, posinf=0.0, neginf=0.0)
    return np.clip(pr, 0, None)


def score_guide_propagation(
    grna_row: pd.Series,
    nodes: pd.DataFrame,
    adj: sparse.csr_matrix,
    node_ids: list[str],
    idx: dict[str, int],
) -> dict:
    gid = grna_row["grna_id"]
    editor = grna_row["editor_system"]
    grna_node = f"grna:{editor}:{gid}"
    if grna_node not in idx:
        matches = [n for n in node_ids if gid in n]
        grna_node = matches[0] if matches else None
    if grna_node is None or grna_node not in idx:
        return {"status": "no_graph_node"}

    sign = 1.0 if "Tet1" in editor or "VP64" in editor else -0.5
    screen_w = min(-np.log10(max(grna_row["padj"], 1e-300)), 50) / 50

    pr = propagate(adj, [idx[grna_node]], [sign * screen_w], len(node_ids))

    gene_scores = {}
    for gene in PWS_GENES:
        matches = [
            i for i, n in enumerate(node_ids)
            if f"gene:{gene}" == n or n.endswith(f":{gene}")
        ]
        if not matches and "name" in nodes.columns:
            matches = [
                i for i, row in nodes.iterrows()
                if row.get("name") == gene and i < len(node_ids)
            ]
        if matches:
            gene_scores[gene] = float(max(pr[i] for i in matches if i < len(pr)))
        elif gene == "SNRPN":
            # SNRPN often annotated within SNHG14 transcriptional unit
            gene_scores[gene] = gene_scores.get("SNHG14", 0) * 0.52

    return {
        "grna_id": gid,
        "editor_system": editor,
        "propagation_scores": gene_scores,
        "SNRPN_propagated": gene_scores.get("SNRPN", 0),
        "SNHG14_propagated": gene_scores.get("SNHG14", 0),
        "combined_propagation": 0.5 * gene_scores.get("SNRPN", 0) + 0.5 * gene_scores.get("SNHG14", 0),
    }


def main():
    nodes = pd.read_parquet(INTEGRATED / "integrated_nodes_hg38.parquet")
    edges = pd.read_parquet(INTEGRATED / "integrated_edges_hg38.parquet")
    merged = pd.read_parquet(INTEGRATED / "locus_merged_hg38.parquet")

    node_ids, adj = build_adjacency(nodes, edges)
    idx = {n: i for i, n in enumerate(node_ids)}

    tet1 = merged[merged["editor_system"] == "dCas9-Tet1"]
    vp64 = merged[merged["editor_system"] == "dCas9-VP64"]
    test_guides = pd.concat([
        tet1[tet1["grna_id"].isin(["PWS_G.4363", "PWS_G.4670"])],
        vp64[vp64["grna_id"].isin(["PWS_G.9414", "PWS_G.8738"])],
    ])

    results = [score_guide_propagation(r, nodes, adj, node_ids, idx) for _, r in test_guides.iterrows()]
    results_df = pd.DataFrame([r for r in results if r.get("status") != "no_graph_node"])

    # Rank all Tet1 guides by propagation
    all_tet1 = []
    for _, r in tet1.iterrows():
        s = score_guide_propagation(r, nodes, adj, node_ids, idx)
        if s.get("status") != "no_graph_node":
            all_tet1.append(s)
    tet1_ranked = pd.DataFrame(all_tet1).sort_values("combined_propagation", ascending=False)

    results_df.to_csv(OUT / "validated_guide_propagation.csv", index=False)
    tet1_ranked.to_csv(OUT / "tet1_propagation_ranked.csv", index=False)

    report = {
        "method": "Personalized PageRank propagation on integrated locus graph",
        "n_nodes": len(node_ids),
        "n_edges": len(edges),
        "damping_alpha": ALPHA,
        "validated_guide_propagation": results_df.to_dict(orient="records"),
        "top5_tet1_by_propagation": tet1_ranked.head(5)[["grna_id", "SNRPN_propagated", "SNHG14_propagated", "combined_propagation"]].to_dict(orient="records"),
        "caveats": [
            "Graph is static; no learned edge weights",
            "Propagation magnitude is relative, not absolute %WT",
            "Complements rules-based forward model, does not replace it",
        ],
    }
    with open(OUT / "graph_propagation_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    with open(OUT / "graph_propagation_model.pkl", "wb") as f:
        pickle.dump({"node_ids": node_ids, "alpha": ALPHA}, f)

    log.info("Graph: %d nodes, %d edges", len(node_ids), len(edges))
    if not results_df.empty:
        log.info("G4363 propagation SNRPN=%.4f SNHG14=%.4f",
                 results_df[results_df.grna_id == "PWS_G.4363"]["SNRPN_propagated"].iloc[0] if "PWS_G.4363" in results_df.grna_id.values else 0,
                 results_df[results_df.grna_id == "PWS_G.4363"]["SNHG14_propagated"].iloc[0] if "PWS_G.4363" in results_df.grna_id.values else 0)
    log.info("Output -> %s", OUT)


if __name__ == "__main__":
    main()
