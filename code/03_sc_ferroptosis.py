# -*- coding: utf-8 -*-
"""Layer-2 (memory-light): IUA single-cell ferroptosis sensitivity atlas.
GSE215968: AS (40,336) vs WOI Control (66,064), 106,400 cells x 33,538 genes.
Adapted from nsfc-direction-B/scripts/01_pe_ferroptosis_map.py (proven block-scan).

Key outputs:
  results/GSE215968_celltype_ferroptosis.csv   per-celltype x module stats (AS vs CTL)
  results/GSE215968_ferro_inflam_coupling.csv  ferro-inflammation coupling per celltype
  results/GSE215968_per_cell_scores.csv.gz     per-cell scores
  figures/fig5_celltype_box.png                dual-target baseline (Epi vs Stroma focus)
  figures/fig6_module_heatmap.png              median module score heatmap
  figures/fig7_coupling_heatmap.png            ferro-inflam coupling heatmap
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = r"C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo"
H5AD = ROOT + "/data/GSE215968_sc.h5ad"
RES, FIG = ROOT + "/results", ROOT + "/figures"

PRO_FERRO = ["ACSL4","LPCAT3","PTGS2","ALOX15","ALOX5","ALOXE3","TFRC","SAT1",
             "NCOA4","SLC39A14","NOX1","NOX4","PLA2G6"]
DEFENSE   = ["GPX4","SLC7A11","SLC3A2","FTH1","FTL","DHODH","GCLC","GCLM",
             "NQO1","SOD2","PRDX1","TXNRD1"]
IRON      = ["TFRC","FTH1","FTL","HMOX1","NCOA4","SLC40A1","SLC39A14"]
INFLAM    = ["IL1B","IL6","CCL2","CCL3","CCL4","CXCL8","TNF","NLRP3","SPP1","CD83"]
ALL_Genes = list(dict.fromkeys(PRO_FERRO + DEFENSE + IRON + INFLAM))

# ---- metadata via h5py (cheap) ----
f = h5py.File(H5AD, "r")
n_cells = f["X/indptr"].shape[0] - 1
group = np.array([x.decode() for x in f["obs/Group"][:]])
ctype = np.array([x.decode() for x in f["obs/principal_cell_types"][:]])
sampl = np.array([x.decode() for x in f["obs/orig.ident"][:]])
genes = np.array([x.decode() if isinstance(x, bytes) else x
                  for x in f["var/features"][:]])
print("cells", n_cells, "| AS:", (group == "AS").sum(),
      "| WOI CTL:", (group == "WOI Control").sum())

gset = set(genes.tolist())
wanted = [g for g in ALL_Genes if g in gset]
col_idx = np.array([int(np.where(genes == g)[0][0]) for g in wanted])
print("genes found:", len(wanted), "/", len(ALL_Genes),
      "| missing:", [g for g in ALL_Genes if g not in gset])

# ---- block-scan CSR extraction ----
X = f["X"]
indptr = X["indptr"][:]
M = np.zeros((n_cells, len(wanted)), dtype=np.float32)
BLOCK = 5000
posmap = {int(c): i for i, c in enumerate(col_idx)}
for s in range(0, n_cells, BLOCK):
    e = min(s + BLOCK, n_cells)
    i0, i1 = int(indptr[s]), int(indptr[e])
    idx = X["indices"][i0:i1]
    val = X["data"][i0:i1]
    lp = indptr[s:e + 1].astype(np.int64) - i0
    rows = np.repeat(np.arange(s, e, dtype=np.int64), np.diff(lp))
    mask = np.isin(idx, col_idx)
    if mask.any():
        h_idx, h_val, h_row = idx[mask], val[mask], rows[mask]
        pos = np.fromiter((posmap[int(c)] for c in h_idx),
                          dtype=np.int64, count=len(h_idx))
        M[h_row, pos] = h_val
f.close()
print("extracted:", M.shape, "| zero frac:", round(float((M == 0).mean()), 4))

# ---- z-score per gene, module scores ----
Z = np.empty_like(M)
for j in range(M.shape[1]):
    col = M[:, j]
    sd = col.std()
    Z[:, j] = (col - col.mean()) / sd if sd > 1e-8 else 0.0

def mod_score(names):
    ii = [wanted.index(g) for g in names if g in wanted]
    return Z[:, ii].mean(axis=1)

obs = pd.DataFrame(dict(group=group, celltype=ctype, sample=sampl))
obs["ferro_activity"] = mod_score(PRO_FERRO)
obs["ferro_defense"]  = mod_score(DEFENSE)
obs["iron_load"]      = mod_score(IRON)
obs["inflam_score"]   = mod_score(INFLAM)

# ---- per-celltype AS vs CTL ----
rows = []
for ct, g in obs.groupby("celltype"):
    a_, c_ = g[g.group == "AS"], g[g.group == "WOI Control"]
    if len(a_) < 15 or len(c_) < 15:
        continue
    for m in ["ferro_activity", "ferro_defense", "iron_load", "inflam_score"]:
        u, p = stats.mannwhitneyu(a_[m], c_[m], alternative="two-sided")
        rows.append(dict(celltype=ct, module=m, n_AS=len(a_), n_CTL=len(c_),
                         median_AS=round(a_[m].median(), 4),
                         median_CTL=round(c_[m].median(), 4),
                         delta=round(a_[m].median() - c_[m].median(), 4),
                         pval=p))
res = pd.DataFrame(rows)
res["padj"] = multipletests(res.pval, method="fdr_bh")[1]
res.to_csv(RES + "/GSE215968_celltype_ferroptosis.csv", index=False)
print(res.to_string(index=False))

# ---- ferro-inflammation coupling per celltype ----
crows = []
for ct, g in obs.groupby("celltype"):
    if len(g) < 30:
        continue
    r, p = stats.spearmanr(g.ferro_activity, g.inflam_score)
    crows.append(dict(celltype=ct, n=len(g), rho=round(r, 3), pval=p))
cpl = pd.DataFrame(crows)
cpl["padj"] = multipletests(cpl.pval, method="fdr_bh")[1]
cpl.to_csv(RES + "/GSE215968_ferro_inflam_coupling.csv", index=False)
print(cpl.to_string(index=False))

# ---- per-sample robustness (sample-level medians) ----
srows = []
for (ct, sp, gr), g in obs.groupby(["celltype", "sample", "group"]):
    if len(g) < 20:
        continue
    srows.append(dict(celltype=ct, sample=sp, group=gr,
                      ferro_activity=g.ferro_activity.median(),
                      ferro_defense=g.ferro_defense.median(),
                      iron_load=g.iron_load.median()))
smed = pd.DataFrame(srows)
smed.to_csv(RES + "/GSE215968_sample_medians.csv", index=False)
print("sample-level rows:", len(smed))

# ---- figures ----
plt.rcParams.update({"figure.dpi": 150, "font.size": 9})
# order: key contrast first
order = ["Epithelium", "AS-Epithelium", "Stromal", "AS-uAS-uSMC",
         "Macrophages", "Endothelium", "Perivascular", "Ciliated Epithelium",
         "NK", "CD8+ T cells", "B cells", "DC", "Mast cells",
         "Erythrocytes"]
order = [c for c in order if c in set(obs.celltype)]
pal = {"AS": "#c0392b", "WOI Control": "#2f7d4f"}

# fig5: boxplots activity+defense, main cell types
sub = obs[obs.celltype.isin(order)]
fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.6))
for ax, m, ttl in zip(axes, ["ferro_activity", "ferro_defense"],
                      ["Pro-ferroptosis activity", "Ferroptosis defense"]):
    sns.boxplot(data=sub, x="celltype", y=m, hue="group", order=order,
                hue_order=["WOI Control", "AS"], palette=pal,
                showfliers=False, ax=ax, boxprops=dict(alpha=.7), width=.75)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
    ax.set_title(ttl); ax.set_xlabel("")
fig.tight_layout()
fig.savefig(FIG + "/fig5_celltype_box.png", bbox_inches="tight")
plt.close(fig)

# fig6: median module heatmap per celltype (delta AS-CTL)
piv = res.pivot_table(index="celltype", columns="module", values="delta")
piv = piv.reindex(order)[["ferro_activity", "ferro_defense", "iron_load", "inflam_score"]]
fig, ax = plt.subplots(figsize=(5.6, 5.2))
sns.heatmap(piv, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
            cbar_kws={"label": "delta median (AS - CTL)"}, ax=ax)
ax.set_title("Ferroptosis module shift by cell type\n(GSE215968, AS vs WOI control)")
fig.tight_layout()
fig.savefig(FIG + "/fig6_module_heatmap.png", bbox_inches="tight")
plt.close(fig)

# fig7: coupling heatmap
hm = cpl.set_index("celltype").loc[[c for c in order if c in set(cpl.celltype)], ["rho"]]
fig, ax = plt.subplots(figsize=(4.6, 4.6))
sns.heatmap(hm, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
            cbar_kws={"label": "Spearman rho"}, ax=ax)
ax.set_title("Ferroptosis activity x inflammation coupling\n(per cell type, all cells)")
fig.tight_layout()
fig.savefig(FIG + "/fig7_coupling_heatmap.png", bbox_inches="tight")
plt.close(fig)

obs.to_csv(RES + "/GSE215968_per_cell_scores.csv.gz",
           index=False, compression="gzip")
print("DONE")
