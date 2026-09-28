# Ferroptosis defense collapse in intrauterine adhesions

Analysis code and key results for:

> **Epithelial ferroptosis defense collapse marks a potentially reversible susceptible state in intrauterine adhesions** (manuscript in preparation, Frontiers in Immunology)

All data are publicly available from GEO. No raw sequencing data are redistributed here — this repository provides the analysis pipeline, intermediate result tables, and publication figures.

## Data sources

| GEO accession | Type | Role in study |
|---|---|---|
| GSE224093 | bulk RNA-seq (severe IUA vs controls) | Layer 1 primary analysis: ferroptosis-defense imbalance |
| GSE215968 | single-cell atlas (~106k cells, 17 cell types) | Layer 2: epithelial dual-hit, myeloid iron sink, HMGB1 axis |
| GSE311899 | paired pre/post CD133+ BMDSC therapy biopsies (same phase 1/2 trial as GSE215968) | Layer 3: therapy-response context (mRNA level, honest-negative for epithelial defense recovery) |
| GSE160633 | within-patient pooled RNA-seq | Layer 4: direction-level concordance check |

## Repository structure

```
code/       analysis scripts (Python 3.13 / R 4.6.1), numbered in run order
results/    key output tables (CSV)
figures/    publication figures (300 dpi PNG)
```

### Code index

| Script | Purpose |
|---|---|
| `01_bulk_deg_ferroptosis.py` | GSE224093 differential expression + ferroptosis pathway scoring |
| `03_sc_ferroptosis.py` | GSE215968 single-cell module scoring (activity / defense / iron / inflammation) |
| `17_p2_hmgb1_sc.py` | HMGB1-axis (20-gene) cell-type scoring |
| `18_p2_sample_robust.py` | sample-level robustness (Welch/MWU/Hedges g) |
| `21/22_p2_gse311899_*.R` | GSE311899 preliminary processing |
| `23_p2_pubfigs.py` / `24_p2_fig1_design.py` | publication figures (Arial 7 pt, 180/85 mm, 300 dpi) |
| `28_p2_gse311899_samples.py` | GSE311899 download / streaming pseudobulk / paired statistics (Wilcoxon + paired t) |
| `30_p2_fig7a.py` | Figure 7a: paired pre/post module line plots |
| `34_p2_figs2.py` | Figure S2: per-gene direction heatmap (37-gene panel) |
| `36_en_refs_renumber.py` | citation renumbering for the manuscript |

## Key findings (summary)

1. **Bulk (GSE224093):** ferroptosis activity elevated, defenses not compensating (0.532 vs 0.563, p = 0.0087, padj = 0.0175); only 4 DEGs at padj < 0.05 & \|log2FC\| > 0.5 — the signal is pathway-level, not single-gene.
2. **Single-cell (GSE215968):** epithelium is the only major cell type with a dual hit (activity +0.145 / defense −0.244, both padj ≈ 0); macrophages are the dominant myeloid iron sink (+0.163, padj = 2.5e-5); ferroptosis–inflammation coupling strongest in macrophages (ρ = 0.30) and dendritic cells (ρ = 0.31).
3. **HMGB1 axis:** epithelial HMGB1 down (−0.155) with myeloid activation (DC +0.291, Mac +0.114) — consistent with release-and-uptake redistribution.
4. **Therapy cohort (GSE311899, 236,867 cells, 8 epithelial + 9 stromal pairs):** epithelial defenses do **not** recover at mRNA level (Δ −0.016, p = 0.46); stromal pro-ferroptotic (+0.487, p = 0.027) and defensive (+0.279, p = 0.039) modules rise in parallel, interpreted as regenerative transcriptional rebound (hypothesis-generating only). Reported transparently as a negative result for epithelial defense recovery.

## Environment

- Python 3.13 (pandas, numpy, scipy, statsmodels, matplotlib, scanpy)
- R 4.6.1 (Matrix, data.table) — scripts 21/22
- Some scripts contain absolute local paths from the original Windows run environment; adjust `BASE`/`RES` constants at the top of each script for your setup.

## Data availability

All datasets are public on GEO (accessions above). Code will be archived with a Zenodo DOI upon publication.

## License

MIT
