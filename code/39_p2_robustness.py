# -*- coding: utf-8 -*-
"""P2 robustness reinforcement (no wet-lab): 
  A. Leave-one-out gene-set sensitivity  -> results/P2_LOO_*.csv, figures/figS3_loo_sensitivity.png
  B. External ferroptosis gene-set cross-validation (Enrichr/KEGG/GO) -> results/P2_external_set_validation.csv
  C. Per-celltype Wnt-ferroptosis coupling (sc layer) -> results/P2_wnt_ferro_coupling_celltype.csv, figures/figS5_wnt_ferro_celltype.png
  D. Myeloid subtyping M1/M2 vs iron-load/HMGB1 -> results/P2_myeloid_subtyping.csv, figures/figS6_myeloid.png
Reuse of 03_sc_ferroptosis.py block-scan; z-scores computed once, LOO only changes column means (cheap).
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

# ---- gene sets (sc layer modules, identical to 03) ----
PRO_FERRO = ["ACSL4","LPCAT3","PTGS2","ALOX15","ALOX5","ALOXE3","TFRC","SAT1",
             "NCOA4","SLC39A14","NOX1","NOX4","PLA2G6"]
DEFENSE   = ["GPX4","SLC7A11","SLC3A2","FTH1","FTL","DHODH","GCLC","GCLM",
             "NQO1","SOD2","PRDX1","TXNRD1"]
IRON      = ["TFRC","FTH1","FTL","HMOX1","NCOA4","SLC40A1","SLC39A14"]
INFLAM    = ["IL1B","IL6","CCL2","CCL3","CCL4","CXCL8","TNF","NLRP3","SPP1","CD83"]
WNT = ["WNT1","WNT2","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT7A","WNT7B","WNT10A",
       "WNT11","WNT16","FZD1","FZD2","FZD3","FZD4","FZD5","FZD6","FZD7","FZD8","FZD9",
       "LEF1","MYC","CCND1","DKK1","SFRP1","WIF1","ROR1","ROR2"]
HMGB1_AXIS = ["HMGB1","TLR2","TLR4","TLR9","AGER","NLRP3","AIM2","PYCARD","CASP1","CASP4",
              "CASP5","IL1B","IL18","MYD88","TICAM1","RELA","NFKB1","NFKBIA","CXCL8","CCL2"]
M1 = ["IL1B","TNF","CXCL9","CXCL10","CCL5","NLRP3","CD80","CCR7","SOCS3","CXCL8",
      "IL6","STAT1","IRF5","CD86"]
M2 = ["MRC1","CD163","MSR1","ARG1","CD200R1","IL10","TGFB1","MERTK","CCL18","STAB1",
      "ALOX15","MARCO","FOLR2","TREM2"]

MODULES = {"PRO_FERRO": PRO_FERRO, "DEFENSE": DEFENSE, "IRON": IRON, "INFLAM": INFLAM}

# key conclusions: (module, celltype, expected direction in AS vs CTL)
KEY = [("DEFENSE", "Epithelium", "down"),
       ("PRO_FERRO", "Epithelium", "up"),
       ("IRON", "Macrophages", "up"),
       ("INFLAM", "Macrophages", "up")]

ALL = list(dict.fromkeys(PRO_FERRO + DEFENSE + IRON + INFLAM + WNT + HMGB1_AXIS + M1 + M2))

# ---- external ferroptosis gene sets via Enrichr (internet; graceful fallback) ----
def fetch_external_sets():
    sets = {}
    try:
        import gseapy
        lib = gseapy.get_library("KEGG_2021_Human")
        for k, genes in lib.items():
            if "ferroptosis" in k.lower():
                sets["KEGG_2021_" + k] = genes
        try:
            lib2 = gseapy.get_library("GO_Biological_Process_2023")
            for k, genes in lib2.items():
                kl = k.lower()
                if "ferroptosis" in kl and "negative" not in kl and "positive" not in kl and "regulation" not in kl:
                    sets["GO_BP_2023_" + k] = genes
        except Exception as e2:
            print("GO lib fetch failed:", e2)
    except Exception as e:
        print("Enrichr fetch failed:", e)
    # keep sets with 10-200 genes, min overlap with panel universe later
    return {k: sorted(set(v)) for k, v in sets.items() if 10 <= len(set(v)) <= 200}

EXT_SETS = fetch_external_sets()
print("external sets:", {k: len(v) for k, v in EXT_SETS.items()})

# ---- metadata + block-scan extraction ----
f = h5py.File(H5AD, "r")
n_cells = f["X/indptr"].shape[0] - 1
group = np.array([x.decode() for x in f["obs/Group"][:]])
ctype = np.array([x.decode() for x in f["obs/principal_cell_types"][:]])
sampl = np.array([x.decode() for x in f["obs/orig.ident"][:]])
genes = np.array([x.decode() if isinstance(x, bytes) else x
                  for x in f["var/features"][:]])
gset = set(genes.tolist())
wanted = [g for g in dict.fromkeys(ALL) if g in gset]
col_idx = np.array([int(np.where(genes == g)[0][0]) for g in wanted])
print("cells", n_cells, "| genes wanted:", len(ALL), "found:", len(wanted),
      "| missing:", [g for g in ALL if g not in gset])

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
        pos = np.fromiter((posmap[int(c)] for c in h_idx), dtype=np.int64, count=len(h_idx))
        M[h_row, pos] = h_val
f.close()
print("extracted:", M.shape)

# ---- z-score once ----
Z = np.empty_like(M)
for j in range(M.shape[1]):
    col = M[:, j]
    sd = col.std()
    Z[:, j] = (col - col.mean()) / sd if sd > 1e-8 else 0.0
gi = {g: i for i, g in enumerate(wanted)}
obs = pd.DataFrame(dict(group=group, celltype=ctype, sample=sampl))
as_mask = (obs.group == "AS").values
ctl_mask = (obs.group == "WOI Control").values

def score(names):
    ii = [gi[g] for g in names if g in gi]
    return Z[:, ii].mean(axis=1)

def mwu_pair(v, sel):
    u, p = stats.mannwhitneyu(v[sel & as_mask], v[sel & ctl_mask], alternative="two-sided")
    return float(np.median(v[sel & as_mask]) - np.median(v[sel & ctl_mask])), float(p)

# ---- A. Leave-one-out ----
loo_rows = []
for mod, genes_mod in MODULES.items():
    present = [g for g in genes_mod if g in gi]
    base_v = score(present)
    for (m2, ct2, direction) in KEY:
        if m2 != mod:
            continue
        sel = (ctype == ct2)
        for drop in [None] + present:
            names = [g for g in present if g != drop]
            v = base_v if drop is None else score(names)
            delta, p = mwu_pair(v, sel)
            loo_rows.append(dict(module=mod, celltype=ct2, expected=direction,
                                 dropped_gene="NONE(full)" if drop is None else drop,
                                 n_genes=len(names),
                                 delta_AS_minus_CTL=round(delta, 4), pval=p))
loo = pd.DataFrame(loo_rows)
loo["direction_ok"] = np.where(loo.expected == "down", loo.delta_AS_minus_CTL < 0,
                               loo.delta_AS_minus_CTL > 0)
loo["sig_p05"] = loo.pval < 0.05
loo["robust"] = loo.direction_ok & loo.sig_p05
loo.to_csv(RES + "/P2_LOO_full.csv", index=False)

summ = (loo[loo.dropped_gene != "NONE(full)"]
        .groupby(["module", "celltype", "expected"])
        .agg(n_loo=("dropped_gene", "count"),
             n_dir_ok=("direction_ok", "sum"),
             n_robust=("robust", "sum")).reset_index())
summ["frac_dir_ok"] = (summ.n_dir_ok / summ.n_loo).round(3)
summ["frac_robust"] = (summ.n_robust / summ.n_loo).round(3)
base = loo[loo.dropped_gene == "NONE(full)"][["module", "celltype", "delta_AS_minus_CTL", "pval"]]
base.columns = ["module", "celltype", "baseline_delta", "baseline_pval"]
summ = summ.merge(base, on=["module", "celltype"])
summ.to_csv(RES + "/P2_LOO_summary.csv", index=False)
print(summ.to_string(index=False))

# figS3: LOO delta range (dot plot with baseline line)
fig, axes = plt.subplots(1, len(summ), figsize=(3.1 * len(summ), 4.0), sharey=False)
if len(summ) == 1: axes = [axes]
for ax, (_, r) in zip(axes, summ.iterrows()):
    sub = loo[(loo.module == r.module) & (loo.celltype == r.celltype) & (loo.dropped_gene != "NONE(full)")]
    col = "#c0392b" if r.expected == "up" else "#2f7d4f"
    ax.scatter(range(len(sub)), sub.delta_AS_minus_CTL, s=22, c=col, alpha=.8)
    ax.axhline(r.baseline_delta, ls="--", lw=1, c="k", label="full panel")
    ax.axhline(0, lw=.8, c="#888888")
    ax.set_title(f"{r.module}\n{r.celltype} (AS−CTL)\nrobust {r.frac_robust*100:.0f}% (dir {r.frac_dir_ok*100:.0f}%)",
                 fontsize=9)
    ax.set_xlabel("one gene removed"); ax.tick_params(labelsize=7)
axes[0].set_ylabel("Δ median module score")
axes[0].legend(fontsize=7)
fig.suptitle("Leave-one-out sensitivity, GSE215968 single-cell modules", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(FIG + "/figS3_loo_sensitivity.png", dpi=200, bbox_inches="tight"); fig.savefig(FIG + "/figS3_loo_sensitivity.pdf", bbox_inches="tight")
plt.close(fig)
print("figS3 saved")

# ---- B. external gene-set validation ----
ext_rows = []
for k, gl in EXT_SETS.items():
    gl = [g for g in gl if g in gi]
    if len(gl) < 8:
        continue
    v = score(gl)
    for ct, direction in [("Epithelium", "up"), ("Macrophages", "up"),
                          ("Stromal", None), ("DC", "up")]:
        sel = (ctype == ct)
        if sel.sum() < 15:
            continue
        delta, p = mwu_pair(v, sel)
        ext_rows.append(dict(set_name=k, n_genes_in_data=len(gl), celltype=ct,
                             delta_AS_minus_CTL=round(delta, 4), pval=p))
ext = pd.DataFrame(ext_rows)
if len(ext):
    ext["padj"] = multipletests(ext.pval, method="fdr_bh")[1]
    ext.to_csv(RES + "/P2_external_set_validation.csv", index=False)
    print(ext.to_string(index=False))
else:
    print("no external sets usable")

# ---- C. per-celltype Wnt-ferroptosis coupling ----
obs["ferro_activity"] = score(PRO_FERRO)
obs["wnt"] = score(WNT)
obs["hmgb1"] = score(HMGB1_AXIS)
obs["M1"] = score(M1)
obs["M2"] = score(M2)
obs["iron_load"] = score(IRON)
obs["inflam_score"] = score(INFLAM)
obs.to_csv(RES + "/GSE215968_per_cell_scores_extended.csv.gz", index=False, compression="gzip")

crows = []
for ct, g in obs.groupby("celltype"):
    if len(g) < 30:
        continue
    r, p = stats.spearmanr(g.ferro_activity, g.wnt)
    crows.append(dict(celltype=ct, n=len(g), rho=round(r, 3), pval=p))
cpl = pd.DataFrame(crows)
cpl["padj"] = multipletests(cpl.pval, method="fdr_bh")[1]
cpl = cpl.sort_values("rho", ascending=False)
cpl.to_csv(RES + "/P2_wnt_ferro_coupling_celltype.csv", index=False)
print(cpl.to_string(index=False))

fig, ax = plt.subplots(figsize=(6.4, 4.6))
d = cpl.sort_values("rho")
cols = ["#c0392b" if v > 0 else "#2f7d4f" for v in d.rho]
ax.barh(d.celltype, d.rho, color=cols, alpha=.85)
for y, (rho, p) in enumerate(zip(d.rho, d.padj)):
    star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    ax.text(rho + (0.004 if rho > 0 else -0.004), y, star, va="center",
            ha="left" if rho > 0 else "right", fontsize=8)
ax.axvline(0, lw=.8, c="#555555")
ax.set_xlabel("Spearman ρ (ferroptosis activity vs Wnt module), per cell")
ax.set_title("Wnt–ferroptosis coupling by cell type (GSE215968)", fontsize=11)
fig.tight_layout()
fig.savefig(FIG + "/figS5_wnt_ferro_celltype.png", dpi=200, bbox_inches="tight"); fig.savefig(FIG + "/figS5_wnt_ferro_celltype.pdf", bbox_inches="tight")
plt.close(fig)
print("figS5 saved")

# ---- D. myeloid subtyping ----
mye = obs[obs.celltype.isin(["Macrophages", "DC"])].copy()
mrows = []
for ct, g in mye.groupby("celltype"):
    for sig in ["M1", "M2"]:
        for tgt in ["iron_load", "hmgb1", "inflam_score", "ferro_activity"]:
            r, p = stats.spearmanr(g[sig], g[tgt])
            mrows.append(dict(celltype=ct, signature=sig, target=tgt, n=len(g),
                              rho=round(r, 3), pval=p))
    # AS vs CTL of M1/M2
    for sig in ["M1", "M2"]:
        a_, c_ = g[g.group == "AS"], g[g.group == "WOI Control"]
        u, p = stats.mannwhitneyu(a_[sig], c_[sig], alternative="two-sided")
        mrows.append(dict(celltype=ct, signature=sig + "_ASvsCTL", target="median_delta",
                          n=len(g), rho=round(a_[sig].median() - c_[sig].median(), 3), pval=p))
mysig = pd.DataFrame(mrows)
mysig["padj"] = multipletests(mysig.pval, method="fdr_bh")[1]
mysig.to_csv(RES + "/P2_myeloid_subtyping.csv", index=False)
print(mysig.to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
sub1 = mysig[mysig.target.isin(["iron_load", "hmgb1", "inflam_score", "ferro_activity"])]
for ax, ct in zip(axes, ["Macrophages", "DC"]):
    d = sub1[sub1.celltype == ct]
    piv = d.pivot(index="signature", columns="target", values="rho")
    piv = piv[["iron_load", "hmgb1", "inflam_score", "ferro_activity"]]
    sns.heatmap(piv, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                vmin=-1, vmax=1, ax=ax, cbar_kws={"shrink": .8})
    ax.set_title(f"{ct} (n per rho col)", fontsize=10)
    ax.set_xlabel("")
fig.suptitle("M1/M2 polarization scores vs iron-load / HMGB1 axis (myeloid)", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(FIG + "/figS6_myeloid.png", dpi=200, bbox_inches="tight"); fig.savefig(FIG + "/figS6_myeloid.pdf", bbox_inches="tight")
plt.close(fig)
print("figS6 saved")
print("ALL DONE")
