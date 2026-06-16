"""
Assemble the PWS locus digital twin: integrated regulatory graph + feature matrix.

Nodes: genes, CpG islands, significant gRNA target sites, ATAC peaks
Edges: Capture-C promoter loops (hypothalamic neurons)
Features per node: expression (WT/edited), methylation, accessibility, screen scores

Scientific framing: static annotated map + empirical perturbation outcomes.
NOT a dynamical simulator — uncertainty-aware design substrate (per Project_Modifications).

Output: data/curated/digital_twin_nodes.parquet
        data/curated/digital_twin_edges.parquet
        data/curated/digital_twin_summary.json
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "locus.yaml"
DATA = ROOT / "data"
CURATED = DATA / "curated"
LOCUS = DATA / "locus_annotation"
OUT = CURATED

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)


def load_locus_regions() -> dict:
    with open(CONFIG, encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_gene_nodes(locus_cfg: dict) -> pd.DataFrame:
    bed = LOCUS / "genes.bed"
    if not bed.exists():
        return pd.DataFrame()
    rows = []
    with open(bed) as f:
        for line in f:
            chrom, start, end, name = line.strip().split("\t")
            pe = name in locus_cfg["key_genes"]["paternally_expressed"]
            me = name in locus_cfg["key_genes"]["maternally_expressed"]
            rows.append({
                "node_id": f"gene:{name}",
                "node_type": "gene",
                "name": name,
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "imprinting": "paternal" if pe else ("maternal" if me else "other"),
                "genome_build": "hg38",
            })
    return pd.DataFrame(rows).drop_duplicates(subset=["node_id"])


def build_cpg_nodes() -> pd.DataFrame:
    path = LOCUS / "cpg_islands.bed"
    if not path.exists():
        return pd.DataFrame()
    rows = []
    with open(path) as f:
        for i, line in enumerate(f):
            chrom, start, end, name = line.strip().split("\t")
            rows.append({
                "node_id": f"cpg:{chrom}:{start}-{end}",
                "node_type": "cpg_island",
                "name": name or f"CpG_{i}",
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "genome_build": "hg38",
            })
    return pd.DataFrame(rows)


def build_grna_nodes(screens: pd.DataFrame, significant_only: bool = True) -> pd.DataFrame:
    df = screens[screens["significant"]].copy() if significant_only else screens.copy()
    if df.empty:
        return pd.DataFrame()
    # Deduplicate by position — keep best padj per (editor, start)
    df = df.sort_values("padj").groupby(["editor_system", "start"], as_index=False).first()
    rows = []
    for _, r in df.iterrows():
        rows.append({
            "node_id": f"grna:{r['editor_system']}:{r['grna_id']}",
            "node_type": "grna_target",
            "name": r["grna_id"],
            "chrom": r["chrom"],
            "start": int(r["start"]),
            "end": int(r["end"]),
            "editor_system": r["editor_system"],
            "subregion": r["subregion"],
            "padj": r["padj"],
            "log2fc": r.get("log2fc"),
            "screen_score": r.get("screen_score"),
            "protospacer": r.get("protospacer"),
            "genome_build": r["genome_build"],
        })
    return pd.DataFrame(rows)


def build_atac_nodes(atac: pd.DataFrame) -> pd.DataFrame:
    if atac.empty:
        return pd.DataFrame()
    peaks = atac.groupby(["chrom", "start", "end", "id"]).agg(
        mean_fpkm=("fpkm", "mean"),
        max_fpkm=("fpkm", "max"),
        samples=("sample", "nunique"),
    ).reset_index()
    rows = []
    for _, r in peaks.iterrows():
        rows.append({
            "node_id": f"atac:{r['chrom']}:{r['start']}-{r['end']}",
            "node_type": "atac_peak",
            "name": r["id"],
            "chrom": r["chrom"],
            "start": int(r["start"]),
            "end": int(r["end"]),
            "mean_fpkm": r["mean_fpkm"],
            "max_fpkm": r["max_fpkm"],
            "genome_build": "hg38",
        })
    return pd.DataFrame(rows)


def build_edges(loops: pd.DataFrame) -> pd.DataFrame:
    if loops.empty:
        return pd.DataFrame()
    hn = loops[loops["cell_type"] == "HypothalamicNeurons"].copy()
    edges = []
    for _, r in hn.iterrows():
        bait_id = f"loop_anchor:{r['bait_chrom']}:{r['bait_start']}-{r['bait_end']}"
        other_id = f"loop_anchor:{r['other_chrom']}:{r['other_start']}-{r['other_end']}"
        edges.append({
            "source": bait_id,
            "target": other_id,
            "edge_type": "capture_c_loop",
            "cell_type": r["cell_type"],
            "n_reads": r["n_reads"],
            "score": r["score"],
            "bait_name": r["bait_name"],
            "other_name": r["other_name"],
            "genome_build": "hg38",
        })
    return pd.DataFrame(edges)


def attach_expression_features(nodes: pd.DataFrame, expr: pd.DataFrame) -> pd.DataFrame:
    if expr.empty or nodes.empty:
        return nodes
    wt = expr[expr["cell_line"] == "WT"].groupby("gene_symbol")["raw_count"].mean()
    edited = expr[expr["edited"]].groupby(["editor_system", "gene_symbol"])["raw_count"].mean()

    nodes = nodes.copy()
    nodes["wt_expression"] = nodes["name"].map(wt)
    for editor in edited.index.get_level_values(0).unique():
        col = f"edited_expr_{editor.replace('dCas9-', '')}"
        nodes[col] = nodes["name"].map(edited[editor])
    return nodes


def main():
    locus_cfg = load_locus_regions()

    screens = pd.read_parquet(CURATED / "editing_screens.parquet")
    expr = pd.read_parquet(CURATED / "expression_outcomes.parquet") if (CURATED / "expression_outcomes.parquet").exists() else pd.DataFrame()
    atac = pd.read_parquet(CURATED / "hypothalamic_atac_chr15.parquet") if (CURATED / "hypothalamic_atac_chr15.parquet").exists() else pd.DataFrame()
    loops = pd.read_parquet(CURATED / "hypothalamic_loops_chr15.parquet") if (CURATED / "hypothalamic_loops_chr15.parquet").exists() else pd.DataFrame()

    gene_nodes = build_gene_nodes(locus_cfg)
    cpg_nodes = build_cpg_nodes()
    grna_nodes = build_grna_nodes(screens)
    atac_nodes = build_atac_nodes(atac)

    nodes = pd.concat([gene_nodes, cpg_nodes, grna_nodes, atac_nodes], ignore_index=True)
    gene_subset = nodes[nodes["node_type"] == "gene"]
    nodes = attach_expression_features(nodes, expr)

    edges = build_edges(loops)

    nodes.to_parquet(OUT / "digital_twin_nodes.parquet", index=False)
    edges.to_parquet(OUT / "digital_twin_edges.parquet", index=False)

    # Load target state if available
    target_path = CURATED / "therapeutic_target_state.json"
    target_state = json.loads(target_path.read_text()) if target_path.exists() else {}

    summary = {
        "description": "PWS 15q11-q13 digital twin (static annotated map + empirical editing outcomes)",
        "genome_builds": {"locus_annotation": "hg38", "editing_screens": "hg19", "note": "Cross-build integration requires liftOver for base-pair alignment"},
        "node_counts": nodes["node_type"].value_counts().to_dict(),
        "edge_count": len(edges),
        "significant_grna_targets": int(grna_nodes.shape[0]),
        "pws_key_genes": int(gene_subset[gene_subset["imprinting"].isin(["paternal", "maternal"])].shape[0]),
        "therapeutic_target_state": target_state.get("editor_outcomes", {}),
        "key_finding": (
            "Tet1 demethylation and VP64 activation hit distinct elements with different "
            "dosage outcomes — Tet1 partial rescue (~46% WT SNRPN), VP64 overshoot (~281% WT)"
        ),
    }
    with open(OUT / "digital_twin_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    log.info("Digital twin: %d nodes, %d edges", len(nodes), len(edges))
    log.info("Node types: %s", summary["node_counts"])


if __name__ == "__main__":
    main()
