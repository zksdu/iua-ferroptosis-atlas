# -*- coding: utf-8 -*-
"""23_p2_pubfigs.py — P2 出版规格主文图重排 (补强④, Fig2-6 + FigS1 + Fig7b)
统一风格: Arial, 300dpi, 双栏 180mm / 单栏 85mm, Frontiers 风格.
不依赖 GSE311899 (Fig7a 占位待补).
"""
import warnings
warnings.filterwarnings("ignore")
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo"
RES, FIG = ROOT + "/results", ROOT + "/figures"
MM = 1 / 25.4

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 8,
    "axes.labelsize": 7.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "legend.fontsize": 6.5, "pdf.fonttype": 42, "ps.fonttype": 42,
})
C_IUA, C_CTL = "#c2504d", "#4878a8"   # 病变红 / 对照蓝 (国标方向无关, 疾病标注色)
C_DEF, C_ACT = "#3a7d5c", "#b0653a"

def stars(p):
    return "**" if p < 0.01 else ("*" if p < 0.05 else "ns")

def save(fig, name):
    fig.savefig(f"{FIG}/{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{FIG}/{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print("saved", name)

# ============ Fig2: bulk 层① (A 通路箱线 B 防御/Wnt 基因方向 C 相关散点) ============
ps = pd.read_csv(f"{RES}/GSE224093_pathway_scores.csv")
stats = pd.read_csv(f"{RES}/GSE224093_pathway_stats.csv")
deg = pd.read_csv(f"{RES}/GSE224093_DEG.csv")

fig, axes = plt.subplots(1, 3, figsize=(180*MM, 60*MM))
# A
mods = ["FERROPTOSIS", "HMGB1_AXIS", "NOTCH", "WNT"]
lab = {"FERROPTOSIS": "Ferroptosis", "HMGB1_AXIS": "HMGB1 axis", "NOTCH": "Notch", "WNT": "Wnt"}
ax = axes[0]
for i, m in enumerate(mods):
    a = ps[ps.group == "IUA"][m]; b = ps[ps.group == "Control"][m]
    bp = ax.boxplot([b, a], positions=[i - 0.18, i + 0.18], widths=0.32,
                    patch_artist=True, medianprops=dict(color="k", lw=0.8),
                    whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                    flierprops=dict(ms=1.5, lw=0.4))
    for box, c in zip(bp["boxes"], [C_CTL, C_IUA]):
        box.set_facecolor(c); box.set_alpha(0.75)
    row = stats[stats.iloc[:, 0] == m].iloc[0]
    pcol = [c for c in stats.columns if "pval" in c.lower() or c == "p"][0]
    ax.text(i, max(a.max(), b.max()) + 0.02, stars(row[pcol]), ha="center", fontsize=6.5)
ax.set_xticks(range(4)); ax.set_xticklabels([lab[m] for m in mods], rotation=20, ha="right")
ax.set_ylabel("Pathway score"); ax.set_title("A  GSE224093 bulk (7 vs 7)")
ax.spines[["top", "right"]].set_visible(False)
# B
key = ["GPX4", "DHODH", "KEAP1", "NCOA4", "SLC7A11", "TFRC", "FTH1", "WNT5A", "SFRP1", "DKK1", "HMGB1", "ACTA2"]
kd = deg[deg.GeneName.isin(key)].set_index("GeneName").loc[key]
ax = axes[1]
colors = [C_ACT if v < 0 else "#888" for v in kd.log2FC]  # 下调绿=防御相关
bars = ax.bar(range(len(kd)), kd.log2FC, color=["#3a7d5c" if g in
        ["GPX4", "DHODH", "KEAP1", "NCOA4", "SLC7A11", "TFRC", "WNT5A", "SFRP1"] else "#9a9a9a"
        for g in kd.index], width=0.7)
for i, (g, r) in enumerate(kd.iterrows()):
    if r.pval < 0.05:
        ax.text(i, r.log2FC + (0.04 if r.log2FC >= 0 else -0.09), "*", ha="center", fontsize=7)
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(range(len(kd))); ax.set_xticklabels(kd.index, rotation=45, ha="right")
ax.set_ylabel("log2FC (IUA vs Control)"); ax.set_title("B  defense & Wnt genes (nominal p)")
ax.spines[["top", "right"]].set_visible(False)
# C
corr = pd.read_csv(f"{RES}/GSE224093_score_correlation.csv", index_col=0)
rho = corr.loc["FERROPTOSIS", "WNT"]
ax = axes[2]
for g, c, mk in [("IUA", C_IUA, "o"), ("Control", C_CTL, "s")]:
    d = ps[ps.group == g]
    ax.scatter(d.WNT, d.FERROPTOSIS, s=14, c=c, marker=mk, label=g, alpha=0.85, edgecolors="none")
z = np.polyfit(ps.WNT, ps.FERROPTOSIS, 1)
xs = np.linspace(ps.WNT.min(), ps.WNT.max(), 10)
ax.plot(xs, np.polyval(z, xs), color="k", lw=0.8, ls="--")
ax.set_xlabel("Wnt score"); ax.set_ylabel("Ferroptosis score")
ax.set_title(f"C  coupling  rho={rho:.2f} (p<0.01)")
ax.legend(frameon=False, loc="lower right", handletextpad=0.2)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, "P2fig2_pub_bulk")

# ============ Fig3: 上皮双打击 (per-cell 箱线) ============
pc = pd.read_csv(f"{RES}/GSE215968_per_cell_scores.csv.gz")
focus = ["Epithelium", "AS-Epithelium", "Stromal", "Macrophages", "Endothelium"]
fig, axes = plt.subplots(1, 2, figsize=(180*MM, 65*MM))
ax = axes[0]
ct0 = pc[pc.celltype == "Epithelium"]
for i, m in enumerate(["ferro_activity", "ferro_defense"]):
    a = ct0[ct0.group == "AS"][m]; b = ct0[ct0.group == "WOI Control"][m]
    bp = ax.boxplot([b, a], positions=[i - 0.18, i + 0.18], widths=0.32,
                    patch_artist=True, medianprops=dict(color="k", lw=0.8),
                    whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                    flierprops=dict(ms=0.8, lw=0.3, marker="."))
    for box, c in zip(bp["boxes"], [C_CTL, C_IUA]):
        box.set_facecolor(c); box.set_alpha(0.7)
ax.set_xticks([0, 1]); ax.set_xticklabels(["Activity", "Defense"])
ax.text(0, 1.35, "padj~0", ha="center", fontsize=6.5)
ax.text(1, 1.35, "padj~0", ha="center", fontsize=6.5)
ax.set_ylabel("Module score (per cell)"); ax.set_title("A  Epithelium: dual hit")
ax.set_ylim(-3.2, 1.6)
ax.spines[["top", "right"]].set_visible(False)
ax = axes[1]
width = 0.38
dsub = pc[pc.celltype.isin(focus)]
deltas, labels, ses = [], [], []
for ct in focus:
    d = dsub[dsub.celltype == ct]
    a = d[d.group == "AS"]; b = d[d.group == "WOI Control"]
    deltas.append(a.ferro_defense.median() - b.ferro_defense.median())
    labels.append(ct)
colors = [C_IUA if v < 0 else C_CTL for v in deltas]
ax.barh(range(len(labels))[::-1], deltas, color=colors, height=0.6, alpha=0.85)
ax.set_yticks(range(len(labels))[::-1]); ax.set_yticklabels(labels)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("Median defense delta (AS - Control)")
ax.set_title("B  defense collapse by cell type")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, "P2fig3_pub_epithelium")

# ============ Fig4: 全细胞型模块漂移热图 ============
ct = pd.read_csv(f"{RES}/GSE215968_celltype_ferroptosis.csv")
piv = ct.pivot_table(index="celltype", columns="module", values="delta")
padj = ct.pivot_table(index="celltype", columns="module", values="padj")
order_m = ["ferro_activity", "ferro_defense", "iron_load", "inflam_score"]
order_m = [m for m in order_m if m in piv.columns]
mlab = {"ferro_activity": "Activity", "ferro_defense": "Defense", "iron_load": "Iron load", "inflam_score": "Inflammation"}
piv = piv[order_m]; padj = padj[order_m]
piv = piv.loc[piv.abs().max(axis=1).sort_values().index]
fig, ax = plt.subplots(figsize=(90*MM, 110*MM))
vmax = 0.35
im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(len(order_m))); ax.set_xticklabels([mlab[m] for m in order_m], rotation=30, ha="right")
ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=6)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v, p = piv.values[i, j], padj.values[i, j]
        if not np.isnan(v):
            ax.text(j, i, ("***" if p < 1e-3 else ("**" if p < 1e-2 else ("*" if p < 0.05 else ""))),
                    ha="center", va="center", fontsize=5.5,
                    color="white" if abs(v) > vmax * 0.55 else "k")
cb = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
cb.set_label("delta (AS - Control)", fontsize=6.5)
ax.set_title("Module drift across cell types")
fig.tight_layout()
save(fig, "P2fig4_pub_heatmap")

# ============ Fig5: 耦合 (髓系锚定) ============
cp = pd.read_csv(f"{RES}/GSE215968_ferro_inflam_coupling.csv")
cp = cp[cp.n >= 150].sort_values("rho", ascending=True)
fig, ax = plt.subplots(figsize=(85*MM, 85*MM))
colors = [C_IUA if (r > 0.25 and padj < 0.05) else ("#b8b8b8" if padj >= 0.05 else C_CTL)
          for r, padj in zip(cp.rho, cp.padj)]
ax.barh(range(len(cp)), cp.rho, color=colors, height=0.62)
ax.set_yticks(range(len(cp))); ax.set_yticklabels(cp.celltype, fontsize=6)
for i, (r, pv, n) in enumerate(zip(cp.rho, cp.pval, cp.n)):
    if pv < 0.05:
        ax.text(r + 0.008, i, f"p={pv:.1e}", va="center", fontsize=5)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("Spearman rho (ferroptosis activity x inflammation)")
ax.set_title("Myeloid-anchored coupling (n>=150)")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, "P2fig5_pub_coupling")

# ============ Fig6: HMGB1 重分布 ============
h = pd.read_csv(f"{RES}/P2_hmgb1_celltype.csv").sort_values("delta")
fig, ax = plt.subplots(figsize=(85*MM, 100*MM))
myeloid = {"Macrophages", "DC", "NK", "Cytotoxic NK", "CD8+ T cells", "MAIT cells", "DN T cells"}
colors = [C_IUA if c in myeloid else "#4878a8" for c in h.celltype]
ax.barh(range(len(h)), h.delta, color=colors, height=0.62)
ax.set_yticks(range(len(h))); ax.set_yticklabels(h.celltype, fontsize=6)
for i, (d, p) in enumerate(zip(h.delta, h.padj)):
    if p < 0.05:
        ax.text(d + (0.008 if d > 0 else -0.008), i, stars(p), va="center",
                ha="left" if d > 0 else "right", fontsize=5.5)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("Median HMGB1-axis delta (AS - Control)")
ax.set_title("Compartmentalized HMGB1 redistribution")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(fc=C_IUA, label="myeloid/lymphoid (up)"),
                   Patch(fc="#4878a8", label="parenchymal (down)")],
          frameon=False, loc="lower right", fontsize=5.5)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, "P2fig6_pub_hmgb1")

# ============ Fig7b: GSE160633 方向印证 ============
kg = pd.read_csv(f"{RES}/P2_GSE160633_keygenes.csv")
kg = kg.sort_values("A_over_N")
fig, ax = plt.subplots(figsize=(85*MM, 65*MM))
colors = [C_DEF if v < 1 else "#9a9a9a" for v in kg.A_over_N]
ax.bar(range(len(kg)), kg.A_over_N, color=colors, width=0.65)
ax.axhline(1, color="k", lw=0.8, ls="--")
ax.set_xticks(range(len(kg))); ax.set_xticklabels(kg.gene, rotation=45, ha="right", fontsize=6)
ax.set_ylabel("FPKM ratio (thin / adjacent normal)")
ax.set_title("GSE160633: 8 IUA patients, pooled (direction only)")
ax.text(len(kg) - 0.5, 1.08, "defense module 11/12 genes lower (key genes shown)", ha="right", fontsize=6, color=C_DEF)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, "P2fig7b_pub_gse160633")

# ============ FigS1: 样本级稳健性 ============
sr = pd.read_csv(f"{RES}/P2_sample_robust.csv")
sr["label"] = sr.celltype + " | " + sr["module"].str.replace("ferro_", "").str.replace("iron_", "iron_")
sr = sr.sort_values("hedges_g")
fig, ax = plt.subplots(figsize=(120*MM, 65*MM))
colors = [C_IUA if g < 0 else C_CTL for g in sr.hedges_g]
ax.barh(range(len(sr)), sr.hedges_g, color=colors, height=0.6)
ax.set_yticks(range(len(sr))); ax.set_yticklabels(sr.label, fontsize=6)
for i, (g, wp, np_) in enumerate(zip(sr.hedges_g, sr.welch_p, sr.n_AS_samples)):
    if g < -2.5:  # 超长负条：文字放条内右端白字，防止越出轴与 y 标签重叠
        ax.text(g + 0.06, i, f"p={wp:.3f} (n={np_}v{sr.n_CTL_samples.iloc[i]})",
                va="center", ha="left", fontsize=5, color="white", fontweight="bold")
    else:
        ax.text(g + (0.06 if g > 0 else -0.06), i, f"p={wp:.3f} (n={np_}v{sr.n_CTL_samples.iloc[i]})",
                va="center", ha="left" if g > 0 else "right", fontsize=5)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("Hedges g (sample-level, AS vs Control)")
ax.set_title("Sample-level robustness: 7/7 direction-consistent")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, "P2figS1_pub_sample_robust")

print("ALL PUBFIGS DONE")
