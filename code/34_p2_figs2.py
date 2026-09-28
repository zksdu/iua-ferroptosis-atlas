# -*- coding: utf-8 -*-
"""34_p2_figs2.py — FigS2: GSE311899 per-gene direction heatmap (Pre -> Post).

Input : results/P2_GSE311899_gene_direction.csv (compartment, gene, n_pairs, mean_delta, frac_up)
Output: figures/P2figS2_gene_direction.png/.pdf
Layout: 2 columns (Epithelium / Stroma) x 37 genes (rows, grouped by module),
        diverging colormap centered at 0 (red=up Post, blue=down Post),
        right strip shows frac_up in Stroma as dot size.
Publication spec: Arial 7pt, 85 mm single-column width, 300 dpi, PNG+PDF.
Run AFTER 28_p2_gse311899_samples.py analyze (works on partial data; rerun when complete).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

PRO_FERRO = ["ACSL4","LPCAT3","PTGS2","ALOX15","ALOX5","ALOXE3","TFRC","SAT1","NCOA4","SLC39A14","NOX1","NOX4","PLA2G6"]
DEFENSE   = ["GPX4","SLC7A11","SLC3A2","FTH1","FTL","DHODH","GCLC","GCLM","NQO1","SOD2","PRDX1","TXNRD1"]
IRON      = ["TFRC","FTH1","FTL","HMOX1","NCOA4","SLC40A1","SLC39A14"]
INFLAM    = ["IL1B","IL6","CCL2","CCL3","CCL4","CXCL8","TNF","NLRP3","SPP1","CD83"]
# order genes: module blocks, dedup keeping first occurrence
ORDER, SEEN = [], set()
for g in PRO_FERRO + DEFENSE + IRON + INFLAM:
    if g not in SEEN:
        ORDER.append(g); SEEN.add(g)
MOD_OF = {}
for name, gs in [("Pro-ferroptosis", PRO_FERRO), ("Defense", DEFENSE), ("Iron", IRON), ("Inflammation", INFLAM)]:
    for g in gs:
        MOD_OF.setdefault(g, name)
MOD_COLORS = {"Pro-ferroptosis": "#C0392B", "Defense": "#2471A3", "Iron": "#7D6608", "Inflammation": "#5B2C6F"}

plt.rcParams.update({"font.family": "Arial", "font.size": 7,
                     "axes.linewidth": 0.5, "pdf.fonttype": 42})
MM = 1 / 25.4

def main():
    gd = pd.read_csv(os.path.join(RES, "P2_GSE311899_gene_direction.csv"))
    comps = ["Epithelium", "Stroma"]
    mat = np.full((len(ORDER), 2), np.nan)
    frac = np.full((len(ORDER), 2), np.nan)
    npair = {}
    for _, r in gd.iterrows():
        if r["gene"] not in SEEN:
            continue
        i = ORDER.index(r["gene"])
        j = comps.index(r["compartment"])
        mat[i, j] = r["mean_delta"]
        frac[i, j] = r["frac_up"]
        npair[r["compartment"]] = int(r["n_pairs"])
    vmax = np.nanmax(np.abs(mat))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)

    h = max(85, 6 * len(ORDER)) * MM
    fig, ax = plt.subplots(figsize=(85 * MM, h))
    im = ax.imshow(mat, cmap="RdBu_r", norm=norm, aspect="auto")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Epithelium\n(n={npair.get('Epithelium','')})",
                        f"Stroma\n(n={npair.get('Stroma','')})"], fontsize=7)
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels(ORDER, fontsize=6)
    ax.tick_params(length=1, width=0.5)
    for s in ax.spines.values():
        s.set_visible(False)
    # module color strip on the left
    for i, g in enumerate(ORDER):
        ax.add_patch(plt.Rectangle((-0.55, i - 0.5), 0.15, 1, color=MOD_COLORS[MOD_OF[g]],
                                   clip_on=False, lw=0))
    ax.set_xlim(-0.55, 1.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("Δ Post − Pre (log1p CPM)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6, length=1, width=0.5)
    # legend for module strip
    handles = [Patch(color=c, label=m) for m, c in MOD_COLORS.items()]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.30, 1.02),
              frameon=False, fontsize=6, handlelength=0.9, handleheight=0.9)
    out = os.path.join(FIG, "P2figS2_gene_direction")
    fig.savefig(out + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print("saved", out + ".png/.pdf")

if __name__ == "__main__":
    main()
