"""Detailed figure and table legends with quantitative inferences for NAR Methods manuscript."""
from __future__ import annotations

import manuscript_data as D

T = D.TOP
ICR_PCT = 100 * D.integration["subregion_hg38_distribution"]["pws_icr"] / 464
G4670 = D.cas_off["per_guide"]["PWS_G.4670"]
G9414 = D.cas_off["per_guide"]["PWS_G.9414"]


def _inf(text: str) -> str:
    return f"\n\nInference: {text}"


def _stats(text: str) -> str:
    return f"\n\nStatistical note: {text}"


def _xref(text: str) -> str:
    return f"\n\nCross-references: {text}"


def main_figure_legends() -> dict[str, str]:
    return {
        "fig1": (
            "Figure 1. End-to-end computational workflow for PWS epigenome-editing therapeutic design prioritization.\n\n"
            "Seven sequential modules are shown. (1) Data curation ingested 23,506 CRISPR epigenome-editing records from "
            "Rohm et al. (2025) GSE285306 and CRISPRepi, yielding 464 padj-significant sites after base-pair merge on GRCh38. "
            "(2) The allele-aware digital twin integrated 897 hypothalamic ATAC peaks, 153 genes, 81 CpG islands, 80 assayed "
            "methylation CpGs, and 464 gRNA targets into 1,675 nodes connected by 2,315 edges. (3) Editor-stratified forward "
            "modeling ranked 57 significant dCas9-Tet1 and 17 dCas9-VP64 guides using rules_reactivation_score (cross-editor ML "
            "rejected; CV R²<0). (4) Hybrid optimization evaluated 121 Tet1×VP64 configurations; 96 (79.3%) achieved dual-gene "
            "expression within the 70–130% WT window. (5) Retrospective validation integrated held-out guides (n=4), organoid "
            "scRNA-seq (GSE262700), neuron RNA (GSE285305), 450K methylation (GSE28525/GSE298378), and Cas-OFFinder-compatible "
            "genome scans. (6) The final catalog ranks 25 uncertainty-aware designs. (7) A prospective four-arm experimental "
            "protocol specifies wet-lab falsification of the top hybrid hypothesis."
            + _inf(
                "The workflow transforms a genome-scale screen into a short, mechanism-ranked list suitable for experimental "
                "follow-up. The bottleneck for therapeutic adequacy is not data volume (23,506 records) but editor-class "
                "complementarity: no monotherapy arm in the Rohm dataset satisfies dual-gene window constraints, motivating "
                "the hybrid optimization module. Every module emits JSON-audited outputs enabling full reproducibility."
            )
            + _xref("Materials and Methods §2.1–2.12; Supplementary Methods S1; Supplementary Table S11.")
        ),
        "fig2": (
            "Figure 2. Editor-class bulk RNA outcomes reveal complementary, non-overlapping therapeutic gaps (GSE243185).\n\n"
            "(A) dCas9-Tet1 editing (ICR demethylation class): SNRPN expression reached 52.3% of WT iPSC mean "
            "(WT mean count 231.8; edited 121.3), falling 17.7 percentage points below the 70% therapeutic lower bound. "
            "SNHG14 reached 99.6% of WT (edited 12,621 vs WT 23,127.5), within the 70–130% green shaded window. "
            "(B) dCas9-VP64 editing (SNRPN promoter activation class): SNRPN reached 246.3% of WT (edited 571.0 counts), "
            "exceeding the 130% upper bound by 116.3 points. SNHG14 remained at 0% of WT (edited 5,487.7 counts vs "
            "near-complete silencing baseline), indicating VP64 cannot substitute for ICR demethylation."
            + _inf(
                "These four numbers—52.3, 99.6, 246.3, and 0.0—are the quantitative foundation of the hybrid hypothesis. "
                "Tet1 solves the SNHG14/SNORD116 deficit that defines PWS imprinting pathology but leaves SNRPN sub-therapeutic; "
                "VP64 closes the SNRPN gap only at the cost of overshoot and complete SNHG14 failure. A linear combination "
                "model predicting SNRPN = (1−w_VP64)×52.3 + w_VP64×246.3 identifies w_VP64≈0.062 as the dose placing SNRPN "
                "at the 70% boundary while preserving Tet1-dominated SNHG14 rescue (~98.3%)."
            )
            + _stats("Expression normalized to WT UPD15 iPSC bulk RNA (Rohm 2025); therapeutic window predefined as 70–130% WT.")
            + _xref("Results §2; Figure 5; Supplementary Table S2; Methods §2.6.")
        ),
        "fig3": (
            f"Figure 3. Spatial architecture of 464 significant gRNA hits within the PWS digital twin (GRCh38).\n\n"
            f"Horizontal bars show site counts per subregion: pws_icr (n=373, {ICR_PCT:.1f}%), snrpn_snh14 (n=75, 16.2%), "
            f"pws_critical other (n=9), snord116_cluster (n=4), outside_pws_critical (n=3). The ICR concentration indicates "
            f"that Tet1-mediated demethylation is the dominant screen signal, consistent with Rohm et al.'s biological conclusion "
            f"that ICR hypomethylation drives SNHG14 rescue. VP64 hits concentrate in snrpn_snh14 (SNRPN promoter/proximal region), "
            f"~92 kb telomeric to the top Tet1 hit PWS_G.4670 (chr15:24,954,921 vs VP64 PWS_G.9414 at 24,862,855)."
            + _inf(
                f"ICR-centric architecture ({ICR_PCT:.1f}% of hits) justifies ICR-targeted Tet1 as the primary therapeutic axis "
                f"and VP64 as a secondary, spatially distinct supplementation. Only 3 guides mapped outside the critical region, "
                f"supporting locus specificity of the screen. The spatial separation between Tet1 ICR guides and VP64 SNRPN "
                f"promoter guides (>90 kb) reduces risk of competing cis-regulatory interference in hybrid delivery."
            )
            + _xref("Supplementary Figure S1; Supplementary Table S7; Methods §2.4.")
        ),
        "fig4": (
            "Figure 4. Tet1 guide ranking as a function of distance to assayed CpG and screen significance.\n\n"
            "Each point is one of 57 significant dCas9-Tet1 guides. X-axis: distance (bp) to nearest bisulfite-assayed CpG "
            "(GSE285300). Y-axis: rules_reactivation_score. Color intensity: −log10(padj). Red marker: PWS_G.4670 "
            "(score=1.02, padj=9.2×10⁻¹⁵, distance=275 bp, edited methylation 37.2%, Δ−54.9%). Orange marker: PWS_G.4363 "
            "(score=1.003, padj=2.95×10⁻⁶, rank 22/57, percentile 61.4%)—bisulfite-validated but mid-ranked by composite score."
            + _inf(
                "ICR-proximal guides cluster at high scores (Spearman ρ=−0.86 for distance vs score, n=57), confirming that "
                "physical proximity to the methylation anchor is the dominant Tet1 design principle. The G4363 outlier demonstrates "
                "that methylation validation and screen ranking capture partially orthogonal properties: G4363 is experimentally "
                "confirmed for ICR demethylation yet ranks below 19 other guides by padj-weighted composite score. This discordance "
                "mandates retaining G4363 as a bisulfite comparator in wet-lab validation (Supplementary Table S1) rather than "
                "treating computational rank as sole ground truth."
            )
            + _stats(
                "Methylation correlation with score: Spearman ρ=0.664 (n=48 guides with methylation data). "
                "Tet1 within-editor Ridge CV R²=−0.17 (supplementary only)."
            )
            + _xref("Supplementary Table S5; Results §2; Methods §2.6.")
        ),
        "fig5": (
            "Figure 5. Phenomenological hybrid dose-response for top guide pair PWS_G.4670 (Tet1) + PWS_G.9414 (VP64).\n\n"
            "Blue curve: predicted SNRPN % WT; red curve: SNHG14 % WT; green band: 70–130% therapeutic window. "
            "Gold point: differential-evolution optimum at w_VP64=0.062, w_Tet1=0.987, predicting SNRPN=70.0% and SNHG14=98.3%. "
            "Blue point: Tet1 monotherapy (w_VP64=0): SNRPN=52.3%, SNHG14=99.6%. Dashed vertical line: optimal VP64 fraction. "
            "SNRPN crosses the 70% lower bound at w_VP64≈0.062; exceeds 130% at w_VP64≈0.38."
            + _inf(
                "The optimizer lands on the SNRPN window boundary rather than the centroid (100% WT), reflecting a conservative "
                "strategy: minimize VP64 exposure while barely satisfying SNRPN adequacy. SNHG14 remains robust (>98%) across "
                "all w_VP64<0.15 because Tet1-class ICR demethylation dominates the SNHG14 axis. Monte Carlo perturbation "
                f"(σ=3.68% for rank-1 design) yields p(both in window)=0.498—meaning the deterministic prediction sits at "
                "a knife-edge under uncertainty and experimental VP64 titration is essential. Bayesian GP on G4363+G9414 "
                f"suggests w_VP64≈0.16 (Supplementary Figure S4), confirming guide-pair dependence."
            )
            + _xref("Table 1; Supplementary Figures S2, S4; Methods §2.7.")
        ),
        "fig6": (
            "Figure 6. Held-out retrospective validation of four experimentally characterized guides (n=4).\n\n"
            "Y-axis: within-editor percentile rank of rules_reactivation_score. Dashed line: 90th percentile (top-decile threshold). "
            "PWS_G.4670 (Tet1): rank 1/57, 98.2nd percentile, pass. PWS_G.4363 (Tet1): rank 22/57, 61.4th percentile, fail. "
            "PWS_G.9414 (VP64): rank 1/17, 94.1st percentile, pass. PWS_G.8738 (VP64): rank 2/17, 88.2nd percentile, pass. "
            "Overall pass rate: 3/4 (75%)."
            + _inf(
                "Within-editor stratification is mandatory: cross-editor ranking is biologically invalid because Tet1 and VP64 "
                "operate on different cis-elements with non-comparable outcome metrics. Tet1 concordance (ρ=0.401, p=0.002, n=57) "
                "and VP64 concordance (ρ=0.630, p=0.007, n=17) support ranking utility within class. The G4363 failure is scientifically "
                "informative: it proves the forward model prioritizes screen significance and ICR proximity over bisulfite-validated "
                "methylation efficacy, defining a specific model limitation to address in v2 with guide-resolved expression data."
            )
            + _xref("Supplementary Table S1; Results §4; Methods §2.9.1.")
        ),
        "fig7": (
            "Figure 7. Convergent validation across three independent experimental contexts.\n\n"
            "(A) Nemoto 2025 organoid scRNA-seq (GSE262700): SNRPN fold-change 5,297× (control 0.00038→edited 2.01 mean counts); "
            "SNHG14 158× (0.034→5.37). SNRPN expressed in 58.5% of edited cells; SNHG14 in 73.6%. "
            "(B) Rohm 2025 bulk RNA: Tet1 SNRPN 52.3% vs neuron GSE285305 47.2% WT (Δ−5.1); Tet1 SNHG14 99.6% vs neuron 56.8% "
            "(Δ−42.8). (C) GSE152098 bigWig ICR signal: ESC 0.057, hypothalamic progenitor 0.332, hypothalamic neuron 0.603 "
            f"({D.gse152098['neuron_vs_esc_icr_fold']:.1f}× neuron/ESC)."
            + _inf(
                "Three independent datasets converge on biological plausibility despite scale mismatch. Organoids confirm massive "
                "reactivation from silenced baseline in hypothalamic tissue—the clinically relevant cell type. Bulk/neuron comparison "
                "shows Tet1-direction concordance for SNRPN but attenuated SNHG14 in neurons, suggesting cell-type-specific calibration "
                "may be needed. The 10.5× ICR accessibility enrichment in hypothalamic neurons supports targeting this locus in "
                "hypothalamic differentiation paradigms, even though sparse ATAC peak calling at the ICR (1 peak in Cousminer atlas) "
                "underestimates accessibility (Nemoto 2025 demonstrated effective Tet1 demethylation regardless)."
            )
            + _xref("Supplementary Figures S6, S9; Supplementary Tables S8–S9; Methods §2.3, §2.9.")
        ),
        "fig8": (
            "Figure 8. Integrated safety assessment: genome-wide off-target and imprinting collateral context.\n\n"
            "(A) Cas-OFFinder-compatible scores (hg38, SpCas9 NGG, mm0–mm4): PWS_G.4670 score=1.438 (mm0=1 on-target, mm3=1, mm4=5); "
            "PWS_G.9414 score=2.625 (mm0=1, mm2=1, mm3=2, mm4=18). Both labeled low risk. "
            f"(B) 450K methylation infrastructure: {D.dmr_map['gse28525']['imprinted_dmrs_in_window']} imprinted DMRs in PWS window "
            f"(GSE28525), {D.dmr_map['gse28525']['icr_imprinted_dmrs']} ICR-associated maternally methylated probes, "
            f"{D.dmr_map['gse298378']['probes_in_window']} probes mapped in GSE298378 (8 samples)."
            + _inf(
                "Safety is multifaceted: in-silico off-target scores are low for both recommended guides, but biochemical validation "
                "(CHANGE-seq/GUIDE-seq) remains mandatory before translation. Locus-level analysis found zero exact protospacer "
                "duplicates and zero 1-mm seed matches among 464 screened sites (Supplementary Table S10). All 25 catalog hybrids "
                "target the ICR >3.5 Mb from the GABAA cluster (chr15:28.5–31.0 Mb), with zero GABAA proximity flags. Collateral "
                "bulk RNA shows minimal UBE3A/ATP10A perturbation under Tet1 (51.8% and 65.2% WT). VP64 dose minimization (6.2%) "
                "is the primary mitigation against SNRPN overshoot and distant methylation risk."
            )
            + _xref("Supplementary Figures S7–S8; Supplementary Tables S2, S10; Methods §2.10.")
        ),
    }


def supplementary_figure_legends() -> dict[str, str]:
    tet1_ode = D.circuit_sim[D.circuit_sim["scenario"] == "Tet1_monotherapy_full"].iloc[0]
    hybrid_ode = D.circuit_sim[D.circuit_sim["scenario"] == "catalog_rank_1"].iloc[0]
    return {
        "sfig1": (
            "Supplementary Figure S1. Composition of the PWS-locus digital twin graph (GRCh38).\n\n"
            f"Pie chart: ATAC peaks 897 (53.6%), gRNA targets 464 (27.7%), genes 153 (9.1%), CpG islands 81 (4.8%), "
            f"methylation CpGs 80 (4.8%). Total nodes: 1,675; edges: 2,315. Integration complete across builds "
            f"(23,506 screen records liftOver-mapped at 100% success rate)."
            + _inf(
                "The twin is ATAC-heavy because hypothalamic neuron accessibility was tiled densely across chr15:24.8–32.7 Mb, "
                "while significant gRNA targets represent the experimentally validated editing space. This composition ensures "
                "guide ranking can be contextualized against chromatin state without conflating screen coverage with regulatory annotation density."
            )
            + _xref("Figure 3; Methods §2.4; Supplementary Table S7.")
        ),
        "sfig2": (
            "Supplementary Figure S2. Monte Carlo window probability for top 5 catalog designs.\n\n"
            f"Bar chart: p(both SNRPN and SNHG14 in 70–130% window) for catalog ranks 1–5. Rank 1: p=0.498 "
            f"(p_SNRPN=0.498, p_SNHG14=1.000). Ranks 2–5: p≈0.500–0.502. Red dashed line: 50% threshold. "
            f"n=50,000 Gaussian samples per design; σ_pct=25×uncertainty (rank 1: σ=3.68%)."
            + _inf(
                "SNHG14 predictions are virtually certain to remain in-window under perturbation because Tet1-dominated dosing "
                "maintains ~99% WT even with ±3.7% noise. Joint probability is entirely governed by SNRPN sitting at the 70% "
                "lower boundary—small upward perturbations fail SNRPN while small downward perturbations pass. This explains why "
                "experimental dose titration should bracket w_VP64 around 0.06 rather than treating 70.0% as a confident point estimate."
            )
            + _xref("Figure 5; Table 1; Methods §2.11.")
        ),
        "sfig3": (
            f"Supplementary Figure S3. Four-objective Pareto frontier (n={D.pareto['n_pareto_optimal']} non-dominated / 121 evaluated).\n\n"
            "Gray points: all 121 hybrid designs. Red points: Pareto-optimal designs (SNRPN window score, SNHG14 window score, "
            "safety proxy, confidence). All Pareto optima include PWS_G.4670 as Tet1 guide with predicted 70.0%/98.3% WT and "
            "objective scores of 1.0 for both gene window terms."
            + _inf(
                f"Only {100*D.pareto['n_pareto_optimal']/121:.1f}% of designs are Pareto-optimal, indicating most guide pairs are "
                "dominated on at least one objective. The tight cluster at maximum window scores confirms that once the correct "
                "Tet1 anchor (G4670) is fixed, VP64 guide choice primarily affects safety and confidence—not the dual-window prediction, "
                "which is dose-dominated."
            )
            + _xref("Results §3; Methods §2.8.")
        ),
        "sfig4": (
            "Supplementary Figure S4. Bayesian Gaussian Process dose optimization vs differential-evolution optimum.\n\n"
            f"Scatter: 25 GP training samples (G4363+G9414 pair); color=objective. Red star: GP best EI "
            f"(w_Tet1={D.bayesian['gp_best_ei']['w_tet1']:.3f}, w_VP64={D.bayesian['gp_best_ei']['w_vp64']:.3f}, "
            f"SNRPN={D.bayesian['gp_best_ei']['predicted_SNRPN_pct']:.1f}%, posterior_std=0.098). "
            "Gold point: DE optimum for G4670+G9414 (w_VP64=0.062, SNRPN=70.0%)."
            + _inf(
                "Cross-method discordance in w_VP64 (0.062 vs 0.160) reflects guide-pair dependence of the phenomenological surrogate, "
                "not optimizer failure. G4363 has lower screen score but validated methylation; the GP explores higher VP64 fractions "
                "to compensate. We report DE on G4670+G9414 as primary because both guides are top within-editor screen hits."
            )
            + _xref("Supplementary Table S4; Methods §2.8.")
        ),
        "sfig5": (
            "Supplementary Figure S5. ODE hypothalamic circuit model: hyperphagia proxy by editing scenario (exploratory).\n\n"
            f"Scenarios: UPD15 untreated (proxy=0.05, AgRP=1.74 vs WT 0.6); Tet1 monotherapy (proxy=0.760, AgRP=0.889); "
            f"hybrid catalog rank 1 (proxy=0.842, AgRP=0.790); WT reference (proxy=1.0). VP64 monotherapy proxy=0.750."
            + _inf(
                f"Hybrid design improves circuit normalization (+8.2 points vs Tet1 mono: 0.842 vs 0.760) by moving SNRPN from "
                f"52.3% toward 70% while maintaining SNHG14 at 98.3%. This model is hypothesis-generating only—not used for "
                f"catalog ranking—and requires physiological validation against PWS feeding phenotypes."
            )
            + _xref("Supplementary Methods S10; Discussion.")
        ),
        "sfig6": (
            "Supplementary Figure S6. Complete organoid reactivation panel (Nemoto 2025, GSE262700).\n\n"
            "Log-scale fold-changes (edited vs control organoid): SNRPN 5,297×, SNHG14 158×, MKRN3 42.8×, ATP10A 2.9×, "
            "UBE3A 1.05×, MAGEL2/NDN/NPAP1 near zero from silenced baseline."
            + _inf(
                "SNRPN and SNHG14 show orders-of-magnitude reactivation, confirming biological feasibility of Tet1-class editing "
                "in hypothalamic organoids. Collateral paternal genes (MKRN3) show unintended reactivation (42.8×), highlighting "
                "the need for collateral monitoring in experimental validation. UBE3A (maternally expressed) is unchanged (1.05×), "
                "supporting low collateral risk at the expression level."
            )
            + _xref("Figure 7A; Supplementary Table S8; Methods §2.9.2.")
        ),
        "sfig7": (
            "Supplementary Figure S7. Cas-OFFinder-compatible mismatch spectrum for recommended guides (full hg38).\n\n"
            f"Tet1 PWS_G.4670: mm0=1, mm1=0, mm2=0, mm3=1 (chr13:108,927,895), mm4=5. Score=1.438. "
            f"VP64 PWS_G.9414: mm0=1, mm1=0, mm2=1 (chr7), mm3=2, mm4=18. Score=2.625. No mm0 off-targets on chr15 for Tet1 beyond on-target."
            + _inf(
                "Mismatch spectra confirm low aggregate off-target burden. The VP64 guide carries more mm4 sites (18 vs 5), "
                "consistent with its higher aggregate score, but both remain in the 'low risk' label category. In-silico scans "
                "do not model bulges, epigenetic modulators, or delivery chemistry; CHANGE-seq/GUIDE-seq is required preclinically."
            )
            + _xref("Figure 8A; Methods §2.10.")
        ),
        "sfig8": (
            "Supplementary Figure S8. ICR imprinted DMR methylation beta values (GSE28525 450K reference).\n\n"
            "Mean beta across ICR maternally methylated imprinted probes: paternal UPD reference, maternal UPD, brain. "
            "Maternal UPD shows elevated beta at ICR probes consistent with biallelic methylation in PWS; paternal UPD shows hypomethylation."
            + _inf(
                "Probe-level mapping anchors collateral methylation monitoring to known imprinted DMRs. Ten imprinted DMRs in the "
                f"PWS window—including {D.dmr_map['gse28525']['icr_imprinted_dmrs']} at the ICR—provide a reference frame for "
                "detecting unintended imprinting disruption in experimental validation (target: edited ICR methylation ~40% from ~95% baseline)."
            )
            + _xref("Figure 8B; Methods §2.5; Supplementary Table S3 (methylation).")
        ),
        "sfig9": (
            "Supplementary Figure S9. Cell-type-specific expression: neuron (GSE285305) minus iPSC bulk (GSE243185).\n\n"
            "Tet1 SNRPN: Δ=−5.2% (52.3→47.2, concordant direction). Tet1 SNHG14: Δ=−42.9% (99.6→56.8, concordant but attenuated). "
            "VP64 SNRPN: Δ=−220.6% (246.3→25.7, discordant magnitude)."
            + _inf(
                "Neuron-specific RNA confirms partial SNRPN rescue under Tet1 but reveals SNHG14 attenuation in neuronal context "
                "(56.8% vs 99.6% in iPSC bulk). This cell-type gap suggests the phenomenological model calibrated on iPSC bulk may "
                "overestimate SNHG14 in neurons and supports hypothalamic neuron/organoid as the primary validation system rather than iPSC bulk."
            )
            + _xref("Figure 7B; Supplementary Table S9; Methods §2.9.3.")
        ),
        "sfig10": (
            "Supplementary Figure S10. Catalog-wide relationship between design uncertainty and genome off-target risk.\n\n"
            f"Scatter of 24 hybrid designs: x=uncertainty (0.147–0.171), y=genome off-target risk (1.44–3.94). "
            f"Color=catalog rank. Rank 1: uncertainty=0.147, off-target=2.625."
            + _inf(
                "Top-ranked designs cluster at low uncertainty (<0.16) with moderate off-target scores driven primarily by VP64 guide choice. "
                "Designs pairing G4670 with G8738 achieve lower off-target (2.438) at similar efficacy, offering an alternative if "
                "G9414 off-target sites warrant concern upon biochemical validation."
            )
            + _xref("Supplementary Table S6; Results §3.")
        ),
    }


def main_table_notes() -> dict[str, str]:
    return {
        "table1": (
            "Table 1. Complete specification of the top-ranked hybrid therapeutic design (catalog rank 1).\n\n"
            f"Tet1 guide PWS_G.4670 targets chr15:24,954,921–24,954,941 (protospacer GCAGGCTGGCGCGCATGCTC, PAM AGG) within the PWS-ICR. "
            f"VP64 guide PWS_G.9414 targets chr15:24,862,855–24,862,875 (protospacer CAGTTCAGGGCATGAGATAA) at the SNRPN promoter region. "
            f"Optimized dosing: w_Tet1=0.987, w_VP64=0.062 (effective VP64 is 6.2% of combined editor load). "
            f"Predicted outcomes: SNRPN=70.0% WT (at lower therapeutic bound), SNHG14=98.3% WT. Objective score=1.0 (both genes in window). "
            f"Uncertainty=0.147; Monte Carlo p(both in window)=0.498. Genome off-target risk score=2.625 (Cas-OFFinder-compatible). "
            f"Experimental evidence: top Tet1 sublibrary hit + top VP64 activation hit (Rohm 2025)."
            + _inf(
                "This design is the primary falsifiable hypothesis for wet-lab validation (protocol Arm A). The 70.0% SNRPN prediction "
                "is intentionally at the window boundary—experimental success requires SNRPN≥70% AND SNHG14≥70% without SNRPN>130%. "
                "If Arm A achieves dual-window expression, the phenomenological model gains credibility; if SNRPN remains <70% at w_VP64=0.06, "
                "the iterate criterion in the protocol mandates VP64 dose escalation."
            )
            + _xref("Figure 5; Supplementary Table S6; Supplementary Table S12.")
        ),
        "table2": (
            "Table 2. Top 10 designs from the uncertainty-aware therapeutic catalog (n=25 total).\n\n"
            "Ranks 1–19: Hybrid_Tet1_VP64 with w_Tet1≈0.987, w_VP64≈0.062, predicted SNRPN=70.0%, SNHG14=98.3%. "
            "Rank 20: Tet1 monotherapy (PWS_G.4670 only): SNRPN=52.3%, SNHG14=99.6%, objective=0.874, fails dual-window. "
            "Off-target scores range 1.44–3.94; top designs with G9414 score 2.625; G4670+G8738 scores 2.438."
            + _inf(
                "The catalog's redundancy (ranks 1–19 share identical predictions) reflects dose-dominated optimization: once optimal "
                "w_VP64 is found, most high-scoring Tet1 ICR guides are interchangeable. Discrimination among hybrids depends on "
                "off-target score and uncertainty. Tet1 monotherapy at rank 20 provides the critical comparator showing SNRPN deficit "
                "without VP64 supplementation."
            )
            + _xref("Supplementary Table S6; Figure 5.")
        ),
    }


def supplementary_table_notes() -> dict[str, str]:
    return {
        "S1": (
            "Supplementary Table S1. Held-out guide validation metrics (n=4 guides, 2 editor classes).\n\n"
            "Pass criterion: within-editor percentile ≥90. Overall pass rate 75% (3/4). G4363 fail is the key discordance."
            + _inf("See Figure 6 for visual representation. G4363 should be prioritized for bisulfite validation despite rank 22.")
            + _xref("Results §4; Methods §2.9.1.")
        ),
        "S2": (
            "Supplementary Table S2. Collateral imprinting gene expression under editor-class editing (% WT, GSE243185).\n\n"
            "Tet1 shows low perturbation at UBE3A (51.8%) and ATP10A (65.2%). NDN/MKRN3 remain at 0% (paternal genes on maternal chr). "
            "VP64 shows moderate NPAP1 reduction (35.3%) and SNRPN overshoot (246.3%)."
            + _inf("Collateral panel defines experimental secondary endpoints in protocol Arms A–B.")
            + _xref("Figure 2; Results §6.")
        ),
        "S6": (
            f"Supplementary Table S6. Complete therapeutic design catalog (n={len(D.catalog)} ranked configurations).\n\n"
            "24 Hybrid_Tet1_VP64 + 1 Tet1_monotherapy. All hybrids share w_Tet1≈0.987, w_VP64≈0.062 when dual-window achieved."
            + _inf("Full machine-readable catalog: data/models/final_catalog/pws_therapeutic_design_catalog.csv")
            + _xref("Table 2; Supplementary Figure S10.")
        ),
    }
