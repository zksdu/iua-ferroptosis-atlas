# -*- coding: utf-8 -*-
"""
Layer-1 bioinformatics support for NSFC IUA/exosome project (GSE224093).
Severe IUA (n=7) vs normal endometrium (n=7), bulk RNA-seq FPKM.

Outputs:
  results/GSE224093_DEG.csv            full DEG table
  results/GSE224093_pathway_scores.csv sample-level pathway/ferroptosis scores
  results/GSE224093_pathway_stats.csv  group stats + spearman correlations
  figures/fig1_volcano.png
  figures/fig2_pathway_scores_box.png
  figures/fig3_score_correlation.png
  figures/fig4_ferroptosis_heatmap.png
"""
import numpy as np
import pandas as pd
import scipy.stats as st
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = r"C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo"
DATA = ROOT + "/data"
RES = ROOT + "/results"
FIG = ROOT + "/figures"

# ---------------- gene sets (curated public-knowledge symbols) ----------------
FERROPTOSIS = [
    "ACSL4","AIFM2","ALOX12","ALOX15","ALOX5","ATG5","ATG7","BECN1","CDKN1A",
    "CISD1","CPS1","DHODH","DHFR","DPP4","FANCD2","FTH1","FTL","G6PD","GCH1",
    "GLS","GLS2","GPX4","GSS","GCLC","GCLM","HMOX1","HSPB1","IREB2","KEAP1",
    "LPCAT3","METTL14","NCOA4","NFE2L2","NOS2","NOX1","NOX4","NQO1","PHKG2",
    "PTGS2","RRM2","SAT1","SLC3A2","SLC7A11","SLC40A1","SLC11A2","SQSTM1",
    "STEAP3","TFRC","TK1","TPP2","VDAC2","VDAC3","ATF4","HAMP","CP","ABCB6",
    "ACSL3","MANF","SLC39A8","SLC39A14","CAV1","ETS1","FANCD2","ISCU","Pdia6",
]
HMGB1_AXIS = [
    "HMGB1","TLR2","TLR4","TLR9","AGER","NLRP3","AIM2","PYCARD","CASP1","CASP4",
    "CASP5","IL1B","IL18","MYD88","TICAM1","RELA","NFKB1","NFKBIA","CXCL8","CCL2",
]
NOTCH = ["NOTCH1","NOTCH2","NOTCH3","NOTCH4","DLL1","DLL3","DLL4","JAG1","JAG2",
         "RBPJ","MAML1","MAML2","HES1","HES5","HEY1","HEYL"]
WNT = ["WNT1","WNT2","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT7A","WNT7B","WNT10A",
       "WNT11","WNT16","FZD1","FZD2","FZD3","FZD4","FZD5","FZD6","FZD7","FZD8","FZD9",
       "FZD10","LRP5","LRP6","CTNNB1","APC","AXIN1","AXIN2","GSK3B","TCF7L1","TCF7L2",
       "LEF1","MYC","CCND1","DKK1","SFRP1","WIF1","ROR1","ROR2"]
GENE_SETS = {"FERROPTOSIS": FERROPTOSIS, "HMGB1_AXIS": HMGB1_AXIS, "NOTCH": NOTCH, "WNT": WNT}

# ---------------- 1. load matrix ----------------
df = pd.read_csv(f"{DATA}/GSE224093_gene_FPKM.txt.gz", sep="\t", index_col=0)
df = df[~df.index.duplicated(keep="first")]
samples = list(df.columns)
groups = ["Control" if s.startswith("Con") else "IUA" for s in samples]
print("matrix:", df.shape, "| groups:", pd.Series(groups).value_counts().to_dict())

log2 = np.log2(df + 1)
iua_idx = [i for i, g in enumerate(groups) if g == "IUA"]
con_idx = [i for i, g in enumerate(groups) if g == "Control"]

# ---------------- 2. DEG (Welch t-test + BH) ----------------
X_iua = log2.iloc[:, iua_idx].values
X_con = log2.iloc[:, con_idx].values
t_stat, p_val = st.ttest_ind(X_iua, X_con, axis=1, equal_var=False)
ok = ~np.isnan(p_val)
padj = np.full_like(p_val, np.nan)
padj[ok] = multipletests(p_val[ok], method="fdr_bh")[1]  # BH on valid tests only
log2fc = X_iua.mean(axis=1) - X_con.mean(axis=1)
deg = pd.DataFrame({"GeneName": log2.index, "log2FC": log2fc, "t": t_stat,
                    "pval": p_val, "padj": padj})
deg.to_csv(f"{RES}/GSE224093_DEG.csv", index=False)
sig_up = ((deg.padj < 0.05) & (deg.log2FC > 0.5)).sum()
sig_dn = ((deg.padj < 0.05) & (deg.log2FC < -0.5)).sum()
print(f"DEG (padj<0.05, |log2FC|>0.5): up={sig_up}, down={sig_dn}")

# ---------------- 3. sample-level gene-set scores (rank-based singscore) ----------------
# Transparent mean-rank scoring: per sample, rank genes (high expr = high rank),
# score = (mean_rank_within_set - mean_rank_all) / max_rank -> [-1, 1]
scores = {}
for name, genes in GENE_SETS.items():
    genes_in = [g for g in genes if g in log2.index]
    missing = len(genes) - len(genes_in)
    ranked = log2.rank(axis=0)          # per-sample ranks, high expr = high rank
    n_all = ranked.shape[0]
    mr_all = ranked.mean(axis=0)
    mr_set = ranked.loc[genes_in].mean(axis=0)
    sc = (mr_set - mr_all) / (n_all - mr_all)
    scores[name] = sc
    print(f"score {name}: {len(genes_in)}/{len(genes)} genes matched")
score_df = pd.DataFrame(scores)
score_df.insert(0, "group", groups)
score_df.insert(0, "sample", samples)
score_df.to_csv(f"{RES}/GSE224093_pathway_scores.csv", index=False)

# group stats per pathway
rows = []
for name in GENE_SETS:
    a = score_df.loc[score_df.group == "IUA", name]
    b = score_df.loc[score_df.group == "Control", name]
    tt = st.ttest_ind(a, b, equal_var=False)
    rows.append({"pathway": name, "mean_IUA": a.mean(), "mean_Control": b.mean(),
                 "diff": a.mean() - b.mean(), "t": tt.statistic, "pval": tt.pvalue})
pw_stats = pd.DataFrame(rows)
pw_stats["padj"] = multipletests(pw_stats.pval, method="fdr_bh")[1]

# spearman correlation of scores across samples
corr = score_df[list(GENE_SETS)].corr(method="spearman")
corr_p = pd.DataFrame(index=corr.index, columns=corr.columns, dtype=float)
for i in corr.index:
    for j in corr.columns:
        if i != j:
            rho, p = st.spearmanr(score_df[i], score_df[j])
            corr_p.loc[i, j] = p
        else:
            corr_p.loc[i, j] = 0.0
pw_stats.to_csv(f"{RES}/GSE224093_pathway_stats.csv", index=False)
corr.to_csv(f"{RES}/GSE224093_score_correlation.csv")
corr_p.to_csv(f"{RES}/GSE224093_score_correlation_pval.csv")
print(pw_stats.round(4).to_string())
print("\nSpearman corr:\n", corr.round(3).to_string())

# ---------------- 4. figures ----------------
sns.set_theme(style="whitegrid", context="paper")
pal = {"Control": "#5DCAA5", "IUA": "#D85A30"}

# fig1 volcano
plt.figure(figsize=(6, 5))
d = deg.copy()
d["-log10p"] = -np.log10(d.pval.clip(lower=1e-300))
d["sig"] = np.where((d.padj < 0.05) & (d.log2FC > 0.5), "up",
            np.where((d.padj < 0.05) & (d.log2FC < -0.5), "down", "ns"))
for k, c in [("ns", "#B4B2A9"), ("down", "#1D9E75"), ("up", "#D85A30")]:
    sub = d[d.sig == k]
    plt.scatter(sub.log2FC, sub["-log10p"], s=4, c=c, label=k, alpha=0.6, linewidths=0)
for g, dx, dy in [("HMGB1", 1.5, 32), ("NCOA4", 1.5, 26), ("SLC7A11", 1.5, 20),
                  ("PTGS2", 1.5, 14), ("GPX4", 1.5, 8)]:
    if g in d.GeneName.values:
        r = d[d.GeneName == g].iloc[0]
        plt.annotate(g, (r.log2FC, r["-log10p"]), fontsize=8,
                     xytext=(r.log2FC + dx * 0.1, r["-log10p"] + dy * 0.3))
plt.legend(frameon=False, fontsize=8)
plt.xlabel("log2FC (IUA vs Control)"); plt.ylabel("-log10(p)")
plt.title("GSE224093 differential expression", fontsize=11)
plt.tight_layout(); plt.savefig(f"{FIG}/fig1_volcano.png", dpi=300); plt.close()

# fig2 pathway score boxes
plt.figure(figsize=(7, 4.2))
melt = score_df.melt(id_vars=["sample", "group"], value_vars=list(GENE_SETS),
                     var_name="pathway", value_name="score")
ax = sns.boxplot(data=melt, x="pathway", y="score", hue="group", palette=pal,
                 width=0.6, fliersize=0)
sns.stripplot(data=melt, x="pathway", y="score", hue="group", palette=pal,
              dodge=True, size=5, alpha=0.9, ax=ax, linewidth=0)
for i, name in enumerate(GENE_SETS):
    p = pw_stats.loc[pw_stats.pathway == name, "pval"].iloc[0]
    y = melt.score.max() + 0.02
    ax.text(i, y, f"p={p:.2g}", ha="center", fontsize=8)
ax.set_xlabel(""); ax.set_ylabel("rank-based score")
ax.set_title("Endometrial pathway activity: IUA vs control (GSE224093)", fontsize=11)
plt.xticks(rotation=15)
plt.legend(frameon=False, fontsize=8)
plt.tight_layout(); plt.savefig(f"{FIG}/fig2_pathway_scores_box.png", dpi=300); plt.close()

# fig3 correlation heatmap
plt.figure(figsize=(5.2, 4.2))
mask_p = corr_p < 0.05
annot = corr.round(2).astype(str) + np.where(mask_p, " *", "")
sns.heatmap(corr, annot=annot, fmt="", cmap="RdBu_r", vmin=-1, vmax=1,
            square=True, linewidths=0.5, cbar_kws={"shrink": 0.8})
plt.title("Pathway score correlation (Spearman, n=14)\n* p<0.05", fontsize=10)
plt.tight_layout(); plt.savefig(f"{FIG}/fig3_score_correlation.png", dpi=300); plt.close()

# fig4 ferroptosis gene heatmap
ferr_in = [g for g in FERROPTOSIS if g in log2.index]
z = log2.loc[ferr_in].T
z = (z - z.mean()) / z.std().replace(0, 1)
z = z.clip(-2.5, 2.5)
order = np.argsort(groups)
z = z.iloc[order]
col_colors = pd.Series([pal[g] for g in np.array(groups)[order]], index=z.index)
g = sns.clustermap(z, col_colors=col_colors, cmap="vlag", center=0,
                   figsize=(8, 9), dendrogram_ratio=0.12, cbar_pos=(0.02, 0.8, 0.03, 0.15),
                   xticklabels=True, yticklabels=False)
g.ax_heatmap.set_xlabel("")
g.figure.suptitle("Ferroptosis regulator genes (z-scored log2FPKM)", y=1.0, fontsize=11)
g.savefig(f"{FIG}/fig4_ferroptosis_heatmap.png", dpi=300)

print("\nDONE. figures + tables written.")
