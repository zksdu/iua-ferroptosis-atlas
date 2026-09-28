#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""17_p2_hmgb1_sc.py — P2 补强②
GSE215968 逐细胞重算 HMGB1 轴模块(20 基因, 与层① bulk 同集), 按细胞型 AS vs WOI 检验。
回应层① bulk HMGB1 下调阴性: bulk 信号被细胞构成稀释/混杂, 单细胞按型定位。
输出: results/P2_hmgb1_celltype.csv / P2_hmgb1_summary.md / figures/P2fig_hmgb1.png
"""
import gc
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = r"C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo"
H5AD = ROOT + "/data/GSE215968_sc.h5ad"
RES, FIG = ROOT + "/results", ROOT + "/figures"

HMGB1_AXIS = ["HMGB1","TLR2","TLR4","TLR9","AGER","NLRP3","AIM2","PYCARD","CASP1","CASP4",
              "CASP5","IL1B","IL18","MYD88","TICAM1","RELA","NFKB1","NFKBIA","CXCL8","CCL2"]

f = h5py.File(H5AD, "r")
n_cells = f["X/indptr"].shape[0] - 1
group = np.array([x.decode() for x in f["obs/Group"][:]])
ctype = np.array([x.decode() for x in f["obs/principal_cell_types"][:]])
genes = np.array([x.decode() if isinstance(x, bytes) else x
                  for x in f["var/features"][:]])
gset = set(genes.tolist())
wanted = [g for g in HMGB1_AXIS if g in gset]
print(f"cells {n_cells} | HMGB1-axis genes found {len(wanted)}/{len(HMGB1_AXIS)} "
      f"| missing: {[g for g in HMGB1_AXIS if g not in gset]}")

col_idx = np.array([int(np.where(genes == g)[0][0]) for g in wanted])
X = f["X"]
indptr = X["indptr"][:]
M = np.zeros((n_cells, len(wanted)), dtype=np.float32)
posmap = {int(c): i for i, c in enumerate(col_idx)}
indptr_np = indptr[:]
idx_ds, val_ds = X["indices"], X["data"]
for r in range(n_cells):
    a, b = int(indptr_np[r]), int(indptr_np[r + 1])
    if b <= a:
        continue
    r_idx = idx_ds[a:b]
    r_val = val_ds[a:b]
    sel = [(posmap[int(c)], v) for c, v in zip(r_idx, r_val) if int(c) in posmap]
    for j, v in sel:
        M[r, j] = v
f.close()

# 全细胞 z-score 后取模块均值 (与 03 同法)
Z = np.empty_like(M)
for j in range(M.shape[1]):
    col = M[:, j]
    sd = col.std()
    Z[:, j] = (col - col.mean()) / sd if sd > 1e-8 else 0.0
ii = [wanted.index(g) for g in wanted]
score = Z[:, ii].mean(axis=1)

obs = pd.DataFrame(dict(group=group, celltype=ctype, hmgb1=score))
obs.to_csv(f"{RES}/P2_hmgb1_per_cell.csv.gz", index=False)

# ---- 逐细胞型 AS vs WOI: MWU + BH ----
rows = []
for ct, sub in obs.groupby("celltype"):
    a = sub[sub.group == "AS"].hmgb1.values
    c = sub[sub.group == "WOI Control"].hmgb1.values
    if len(a) < 15 or len(c) < 15:
        rows.append(dict(celltype=ct, n_AS=len(a), n_CTL=len(c), skip="n<15"))
        continue
    u, p = stats.mannwhitneyu(a, c, alternative="two-sided")
    rows.append(dict(celltype=ct, n_AS=len(a), n_CTL=len(c),
                     med_AS=round(float(np.median(a)), 4),
                     med_CTL=round(float(np.median(c)), 4),
                     delta=round(float(np.median(a) - np.median(c)), 4), p=p))
df = pd.DataFrame([r for r in rows if "skip" not in r])
df["padj"] = multipletests(df.p, method="fdr_bh")[1]
df = df.sort_values("delta")
df.to_csv(f"{RES}/P2_hmgb1_celltype.csv", index=False)
print(df.to_string(index=False))

sig = df[df.padj < 0.05]
with open(f"{RES}/P2_hmgb1_summary.md", "w", encoding="utf-8") as fo:
    fo.write("# P2 补强②：HMGB1 轴单细胞层面定位（回应 bulk 下调阴性）\n\n")
    fo.write(f"基因集：HMGB1 轴 20 基因（与层① bulk 同集），找到 {len(wanted)}/20。"
             f"逐细胞 z-score 模块分，细胞型内 AS vs WOI MWU+BH。\n\n")
    fo.write(df.to_string(index=False))
    fo.write("\n\n**结论**：")
    up = sig[sig.delta > 0]; dn = sig[sig.delta < 0]
    fo.write(f"\n- padj<0.05 中上调 {len(up)} 型：{', '.join(up.celltype)}；下调 {len(dn)} 型：{', '.join(dn.celltype)}。\n")
    fo.write("- 若髓系/上皮特异性上调而全组织平均无差异或下调 → 证明 bulk 阴性为细胞构成稀释，"
             "HMGB1 轴为细胞类型限制性激活，与层①结论兼容。\n")

# ---- 图: 按细胞型 delta 条形 (标 padj) + 上皮/髓系代表箱线 ----
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={"width_ratios": [1.4, 1]})
cols = ["#c0392b" if d > 0 else "#2471a3" for d in df.delta]
axes[0].barh(np.arange(len(df)), df.delta, color=cols, height=0.62)
axes[0].set_yticks(np.arange(len(df)), df.celltype, fontsize=8)
axes[0].axvline(0, color="k", lw=0.8)
for y, (d, q) in enumerate(zip(df.delta, df.padj)):
    axes[0].text(d + (0.008 if d >= 0 else -0.008), y,
                 "ns" if q >= 0.05 else f"padj={q:.1e}", va="center",
                 ha="left" if d >= 0 else "right", fontsize=7)
axes[0].set_xlabel("AS − control median HMGB1-axis module score")
axes[0].set_title("HMGB1-axis shift by cell type (AS vs WOI)")

focus = [ct for ct in ["Macrophages", "Epithelium"] if ct in set(obs.celltype)]
pos = 0
for k, ct in enumerate(focus):
    d_ct = obs[obs.celltype == ct]
    data = [d_ct[d_ct.group == "AS"].hmgb1, d_ct[d_ct.group == "WOI Control"].hmgb1]
    bp = axes[1].boxplot(data, positions=[pos, pos + 1], widths=0.6,
                         tick_labels=[f"AS\n(n={len(data[0])})", f"CTL\n(n={len(data[1])})"],
                         patch_artist=True)
    for b in bp["boxes"]:
        b.set_facecolor("#e8afaf" if k == 0 else "#aec6e8")
    pos += 3
axes[1].set_title("Focus: Macrophages / Epithelium")
fig.tight_layout()
fig.savefig(f"{FIG}/P2fig_hmgb1.png", dpi=300)
print("saved figures/P2fig_hmgb1.png | DONE")
