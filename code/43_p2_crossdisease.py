# -*- coding: utf-8 -*-
"""43_p2_crossdisease.py — Figure S7: IUA vs preeclampsia module deltas, same 4 modules.
IUA side:  results/GSE215968_celltype_ferroptosis.csv   (AS vs mid-secretory CTL)
PE side:   nsfc-direction-B/results/B1_ferroptosis_by_celltype.csv (PE vs CTL placenta)
Matched functional compartments (same scoring scheme, z-mean modules, MWU+BH):
  Parenchyma  : IUA Epithelium   <-> PE Trophoblast
  Stroma      : IUA Stromal      <-> PE Fibro_Stromal
  Macrophage  : IUA Macrophages  <-> PE Hofbauer
  Endothelial : IUA Endothelium  <-> PE Endothelial
  Erythroid   : IUA Erythrocytes <-> PE Erythrocyte
Output: results/P2_crossdisease_module_deltas.csv, figures/figS7_crossdisease.png(+pdf)
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22"
iua = pd.read_csv(BASE + "/iua-exosome-bioinfo/results/GSE215968_celltype_ferroptosis.csv")
pe = pd.read_csv(BASE + "/nsfc-direction-B/results/B1_ferroptosis_by_celltype.csv")

PAIRS = [("Parenchyma", "Epithelium", "Trophoblast"),
         ("Stroma", "Stromal", "Fibro_Stromal"),
         ("Macrophage", "Macrophages", "Hofbauer"),
         ("Endothelial", "Endothelium", "Endothelial")]
MODS = [("ferro_activity", "Pro-ferroptotic activity"),
        ("ferro_defense", "Ferroptosis defense"),
        ("iron_load", "Iron load"),
        ("inflam_score", "Inflammation")]
ORDER = [p[0] for p in PAIRS]

rows = []
for comp, iua_ct, pe_ct in PAIRS:
    for mod, _ in MODS:
        ri = iua[(iua.celltype == iua_ct) & (iua.module == mod)].iloc[0]
        rp = pe[(pe.celltype == pe_ct) & (pe.module == mod)].iloc[0]
        rows.append(dict(compartment=comp, module=mod,
                         IUA_ct=iua_ct, PE_ct=pe_ct,
                         delta_IUA=ri.delta, padj_IUA=ri.padj,
                         delta_PE=rp.delta, padj_PE=rp.padj))
cd = pd.DataFrame(rows)
cd.to_csv(BASE + "/iua-exosome-bioinfo/results/P2_crossdisease_module_deltas.csv", index=False)
print(cd.round(4).to_string(index=False))

# ---- figure: 2x2 panels ----
fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2))
CI, CP = "#c0392b", "#2f6690"
for ax, (mod, title) in zip(axes.flat, MODS):
    d = cd[cd.module == mod].set_index("compartment").loc[ORDER]
    y = np.arange(len(d))
    h = 0.36
    ax.barh(y + h/2, d.delta_IUA, height=h, color=CI, alpha=.85, label="IUA (AS − CTL)")
    ax.barh(y - h/2, d.delta_PE, height=h, color=CP, alpha=.85, label="PE placenta (PE − CTL)")
    for yy, (di, pi, dp, pp) in enumerate(zip(d.delta_IUA, d.padj_IUA, d.delta_PE, d.padj_PE)):
        for val, pv, xoff in [(di, pi, 0.006), (dp, pp, 0.006)]:
            star = "***" if pv < 1e-3 else "**" if pv < 1e-2 else "*" if pv < 5e-2 else "ns"
            ax.text(val + (xoff if val >= 0 else -xoff), yy + (h/2 if (val == di and pv == pi) else -h/2),
                    star, va="center", ha="left" if val >= 0 else "right", fontsize=7.5)
    ax.axvline(0, lw=.8, c="#555555")
    ax.set_yticks(y); ax.set_yticklabels(ORDER, fontsize=9)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=11)
    if mod == "ferro_activity":
        ax.legend(fontsize=8, loc="lower right")
fig.suptitle("Same-module cross-disease comparison: IUA endometrium vs preeclamptic placenta "
             "(per-cell medians, MWU+BH; *p<0.05 **p<0.01 ***p<0.001)", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(BASE + "/iua-exosome-bioinfo/figures/figS7_crossdisease.png", dpi=200, bbox_inches="tight")
fig.savefig(BASE + "/iua-exosome-bioinfo/figures/figS7_crossdisease.pdf", bbox_inches="tight")
plt.close(fig)
print("figS7 saved")
