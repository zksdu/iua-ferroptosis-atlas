# -*- coding: utf-8 -*-
"""40: bulk deconvolution-style composition control + GSE311899 stromal GSEA.
A. MCP-counter-style single-sample marker enrichment for major compartments in GSE224093
   bulk (log2FPKM, n=14): do IUA bulk ferroptosis-score decreases persist after
   adjusting for estimated epithelial fraction? (linear regression score ~ group + epi)
B. GSEA prerank (gseapy, MSigDB_Hallmark_2020 + KEGG_2021) on GSE311899 paired
   per-gene mean delta (Post-Pre) in Stroma (the double-rise compartment).
Outputs:
  results/P2_bulk_deconv.csv, results/P2_bulk_deconv_adjusted.csv
  results/P2_GSE311899_GSEA_stroma.csv (+ epithelium), figures/figS4_deconv_gsea.png
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = r"C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo"
RES, FIG = ROOT + "/results", ROOT + "/figures"

MARKERS = {
 "Epithelium": ["EPCAM","KRT8","KRT18","KRT19","KRT5","PAEP","KRT15","CLU","MUC1"],
 "Stroma":     ["DCN","LUM","COL1A1","COL3A1","COL6A1","C7","PDGFRA","AEBP1"],
 "Endothelium":["PECAM1","VWF","CLDN5","CDH5","EGFL7","ADGRF5"],
 "Myeloid":    ["LYZ","CD68","CD14","FCGR3A","CST3","CD163","MSR1","MRC1"],
 "Lymphoid":   ["CD3D","CD3E","NKG7","GNLY","MS4A1","CD79A","CD2","CD7"],
}

# ---- A. bulk composition ----
bulk = pd.read_csv(ROOT + "/data/GSE224093_gene_FPKM.txt.gz", sep="\t", index_col=0)
bulk = bulk.groupby(level=0).mean()
l2 = np.log2(bulk + 1)
def zrow(g):
    v = l2.loc[g]
    return (v - v.mean()) / v.std()
comp = {}
for ct, genes in MARKERS.items():
    gs = [g for g in genes if g in l2.index]
    comp[ct] = pd.Series(np.mean([zrow(g) for g in gs], axis=0), index=l2.columns)
comp = pd.DataFrame(comp)
scores = pd.read_csv(RES + "/GSE224093_pathway_scores.csv").set_index("sample")
meta = comp.join(scores, how="inner")
meta["group"] = np.where(meta.index.str.startswith("Con"), "Control", "IUA")
print(meta[["Epithelium","Stroma","Myeloid","FERROPTOSIS","group"]].to_string())

comp.to_csv(RES + "/P2_bulk_deconv.csv")
# group test on composition itself
rows = []
for ct in MARKERS:
    a = meta[meta.group == "IUA"][ct]; c = meta[meta.group == "Control"][ct]
    u, p = stats.mannwhitneyu(a, c, alternative="two-sided")
    rows.append(dict(compartment=ct, median_IUA=a.median(), median_CTL=c.median(),
                     delta=a.median()-c.median(), pval=p))
ctab = pd.DataFrame(rows)
ctab.to_csv(RES + "/P2_bulk_deconv_groups.csv", index=False)
print(ctab.to_string(index=False))

# adjusted regression: score ~ group + epithelial score (and each pathway)
def adj_test(scorecol, covars):
    y = meta[scorecol].values.astype(float)
    Xd = []
    Xd.append((meta.group == "IUA").astype(float).values)  # group effect
    for cv in covars:
        Xd.append(meta[cv].values.astype(float))
    Xd.append(np.ones(len(y)))
    X = np.column_stack(Xd)
    beta, res_, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - rank
    sigma2 = (resid @ resid) / dof
    covb = sigma2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(covb))
    t = beta[0] / se[0]
    p = 2 * stats.t.sf(abs(t), dof)
    return dict(score=scorecol, covariates="+".join(covars),
                beta_group=round(beta[0], 4), se=round(se[0], 4), t=round(t, 3), pval_adj=p)

arows = []
for sc in ["FERROPTOSIS", "HMGB1_AXIS", "WNT"]:
    arows.append(adj_test(sc, []))
    arows.append(adj_test(sc, ["Epithelium"]))
    arows.append(adj_test(sc, ["Epithelium", "Stroma", "Myeloid"]))
adj = pd.DataFrame(arows)
adj.to_csv(RES + "/P2_bulk_deconv_adjusted.csv", index=False)
print(adj.to_string(index=False))

# ---- B. GSEA on GSE311899 ----
import gseapy
gd = pd.read_csv(RES + "/P2_GSE311899_gene_direction.csv")
gsea_out = []
for comp_name in ["Stroma", "Epithelium"]:
    sub = gd[gd.compartment == comp_name].dropna(subset=["mean_delta"])
    rnk = sub[["gene", "mean_delta"]].sort_values("mean_delta", ascending=False)
    rnk = rnk.drop_duplicates("gene").set_index("gene")
    for lib in ["MSigDB_Hallmark_2020", "KEGG_2021_Human"]:
        try:
            pre = gseapy.prerank(rnk=rnk, terms=None, gene_sets=lib, outdir=None,
                                 min_size=8, max_size=500, threads=4)
            pr = pre.res2d
            pr["compartment"] = comp_name; pr["library"] = lib
            gsea_out.append(pr)
        except Exception as e:
            print("gsea fail", comp_name, lib, e)
gsea = pd.concat(gsea_out, ignore_index=True)
num = pd.to_numeric(gsea["NES"], errors="coerce")
gsea["FDR q-val"] = pd.to_numeric(gsea["FDR q-val"], errors="coerce")
gsea = gsea[num.notna()]
gsea.to_csv(RES + "/P2_GSE311899_GSEA_stroma.csv", index=False)
show = gsea[(gsea.compartment == "Stroma") & (gsea["FDR q-val"] < 0.25)].sort_values("NES", ascending=False)
print("Stroma significant sets (FDR<0.25):")
print(show[["Term", "NES", "FDR q-val"]].head(15).to_string(index=False))

# ---- figS4: composition + adjusted scores (GSEA not feasible: panel-limited rank) ----
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
ax = axes[0]
ctab.plot(x="compartment", y=["median_CTL", "median_IUA"], kind="bar", ax=ax,
          color=["#2f7d4f", "#c0392b"], legend=True)
ax.set_title("Marker-based composition scores (bulk)"); ax.set_xlabel("")
ax.tick_params(axis="x", rotation=30)
ax2 = axes[1]
lab = {"FERROPTOSIS": "Ferroptosis", "HMGB1_AXIS": "HMGB1 axis", "WNT": "Wnt"}
xx = np.arange(3)
for i, sc in enumerate(["FERROPTOSIS", "HMGB1_AXIS", "WNT"]):
    r0 = adj[(adj.score == sc) & (adj.covariates == "")].iloc[0]
    r1 = adj[(adj.score == sc) & (adj.covariates == "Epithelium")].iloc[0]
    ax2.scatter([i - 0.12], r0.pval_adj, c="#555555", s=50, label="unadjusted" if i == 0 else None)
    ax2.scatter([i + 0.12], r1.pval_adj, c="#c0392b", s=50, label="+ epithelium adj." if i == 0 else None)
    ax2.text(i + 0.12, r1.pval_adj, f"p={r1.pval_adj:.3f}", fontsize=7, ha="left")
ax2.axhline(0.05, ls="--", lw=.8, c="#888888")
ax2.set_yscale("log"); ax2.set_xticks(xx)
ax2.set_xticklabels([lab[s] for s in ["FERROPTOSIS", "HMGB1_AXIS", "WNT"]])
ax2.set_ylabel("adjusted p (group effect)")
ax2.set_title("Group effect persists after composition adjustment")
ax2.legend(fontsize=7)
fig.tight_layout()
fig.savefig(FIG + "/figS4_deconv_gsea.png", dpi=200, bbox_inches="tight"); fig.savefig(FIG + "/figS4_deconv_gsea.pdf", bbox_inches="tight")
plt.close(fig)
print("figS4 saved; DONE")
