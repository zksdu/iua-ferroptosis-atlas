# -*- coding: utf-8 -*-
"""24_p2_fig1_design.py — Fig1 研究设计示意 (4 数据集 → 双分辨率分析 → 核心发现 → 结论)"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

FIG = r"C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo/figures"
MM = 1 / 25.4
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})

fig, ax = plt.subplots(figsize=(180*MM, 105*MM))
ax.set_xlim(0, 100); ax.set_ylim(0, 58); ax.axis("off")

def box(x, y, w, h, fc, text, fs=6.8, tc="k", lw=0.8, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4",
                                fc=fc, ec="#555", lw=lw))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if bold else "normal", linespacing=1.35)

def arrow(x1, y1, x2, y2, c="#666"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=8, color=c, lw=0.9))

# ---- 列 1: 数据集 ----
box(1, 44, 24, 10, "#eef3fa", "GSE224093\nSevere IUA 7 vs Control 7\nbulk RNA-seq (FPKM)", bold=False)
box(1, 31, 24, 10, "#eef3fa", "GSE215968\nAS 40,336 vs WOI 66,064 cells\nscRNA-seq, 17 cell types")
box(1, 18, 24, 10, "#fdf1ef", "GSE311899\n9 AS patients, Pre/Post\nCD133+ therapy, Epi/Stroma")
box(1, 5, 24, 10, "#f2f7f2", "GSE160633\nThin vs adjacent normal\n8 IUA patients (direction)")
ax.text(13, 56.2, "Public cohorts (GEO)", ha="center", fontsize=7.5, fontweight="bold")

# ---- 列 2: 分析 ----
box(33, 40, 26, 16, "#fafafa", "Bulk scoring\nrank-based pathway scores\nFerroptosis / HMGB1 / Notch / Wnt\nDEG direction + Spearman")
box(33, 18, 26, 18, "#fafafa", "Single-cell modules\nper-cell z-score, 4 modules\nactivity / defense / iron / inflammation\nMWU+BH per cell type\n+ sample-level robustness")
ax.text(46, 56.2, "Analysis", ha="center", fontsize=7.5, fontweight="bold")
arrow(25, 49, 33, 48); arrow(25, 36, 33, 30)
arrow(25, 23, 33, 22.5); arrow(25, 10, 33, 20)

# ---- 列 3: 发现 ----
box(67, 47, 31, 8, "#fff8e6", "Defense collapse, not execution\nGPX4/DHODH/KEAP1 axis down; score lower")
box(67, 38, 31, 8, "#fff8e6", "Epithelium = dual hit\nactivity +0.145 / defense -0.244")
box(67, 29, 31, 8, "#fff8e6", "Macrophages = iron sink\niron load +0.163, highest of all types")
box(67, 20, 31, 8, "#fff8e6", "Myeloid-anchored coupling\nDC 0.31 > Mac 0.30 >> Epi 0.07")
box(67, 11, 31, 8, "#fff8e6", "HMGB1 redistribution\nparenchyma down / myeloid up")
ax.text(82.5, 56.2, "Findings", ha="center", fontsize=7.5, fontweight="bold")
arrow(59, 48, 67, 51); arrow(59, 30, 67, 26); arrow(59, 21, 67, 15)

# ---- 底部结论条 ----
box(20, 0.5, 62, 7, "#e8f0e8", "Therapeutic implication: restore epithelial ferroptosis defense — susceptible state is potentially reversible;\nmyeloid iron niche shared across fibrotic pregnancy disease (IUA / PE)", fs=7, bold=False)
arrow(82.5, 11, 60, 7.5)
arrow(46, 18, 46, 7.5)

fig.tight_layout()
fig.savefig(f"{FIG}/P2fig1_pub_design.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{FIG}/P2fig1_pub_design.pdf", bbox_inches="tight")
print("saved P2fig1_pub_design")
