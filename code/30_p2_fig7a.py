# -*- coding: utf-8 -*-
"""30_p2_fig7a.py — Fig7a: GSE311899 paired Pre/Post module scores (Epi/Stroma).

Input : results/P2_GSE311899_module_scores.csv   (patient, comp, stage, 4 modules)
Output: figures/P2fig7a_pub_gse311899.png/.pdf
Layout: 1x4 paired-line panels (PRO_FERRO / DEFENSE / IRON / INFLAM),
        two lines per patient (Epi=red-ish, Stroma=blue-ish), Post-point colored by direction.
Run AFTER 28_p2_gse311899_samples.py analyze.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

MODS = [("PRO_FERRO", "Pro-ferroptosis"), ("DEFENSE", "Ferroptosis defense"),
        ("IRON", "Iron handling"), ("INFLAM", "Inflammation / HMGB1 axis")]
COMP_STYLE = {"Epithelium": ("#C0392B", "Epithelium"), "Stroma": ("#2471A3", "Stroma")}

def main():
    sc = pd.read_csv(os.path.join(RES, "P2_GSE311899_module_scores.csv"))
    st = pd.read_csv(os.path.join(RES, "P2_GSE311899_prepost_stats.csv"))
    st.index = pd.MultiIndex.from_arrays([st["compartment"], st["module"]])

    plt.rcParams.update({"font.family": "Arial", "font.size": 7,
                         "axes.linewidth": 0.6, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 4, figsize=(180/25.4, 55/25.4))

    for ax, (m, label) in zip(axes, MODS):
        d = sc[["patient", "comp", "stage", m]].dropna()
        for comp, (color, clabel) in COMP_STYLE.items():
            dd = d[d["comp"] == comp]
            if dd.empty:
                continue
            for pat, g in dd.groupby("patient"):
                g = g.set_index("stage")
                if {"Pre", "Post"} <= set(g.index):
                    ax.plot([0, 1], [g.loc["Pre", m], g.loc["Post", m]],
                            color=color, lw=0.7, alpha=0.45, zorder=1)
                    ax.scatter([0], [g.loc["Pre", m]], color=color, s=8, zorder=2)
                    post = g.loc["Post", m]
                    ax.scatter([1], [post], color=color, s=8, zorder=2)
            # mean line
            piv = dd.pivot(index="patient", columns="stage", values=m).dropna()
            if len(piv):
                ax.plot([0, 1], [piv["Pre"].mean(), piv["Post"].mean()],
                        color=color, lw=1.8, zorder=3)
        # annotation from stats table
        yslots = {"Epithelium": 0.03, "Stroma": 0.11}
        for comp in COMP_STYLE:
            key = (comp, m)
            if key in st.index:
                r = st.loc[key]
                star = ""
                p = r["wilcox_p"]
                if p < 0.001: star = "***"
                elif p < 0.01: star = "**"
                elif p < 0.05: star = "*"
                else: star = "ns"
                y = yslots[comp]
                ax.text(0.5, y, f"{COMP_STYLE[comp][1][0]}: {star} (p={p:.3g}, n={int(r['n_pairs'])})",
                        transform=ax.transAxes, fontsize=5.2, ha="center",
                        color=COMP_STYLE[comp][0])
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Pre", "Post"])
        ax.set_title(label, fontsize=7, fontweight="bold")
        if m == "PRO_FERRO":
            ax.set_ylabel("Module score (mean log1p CPM)")
        ax.tick_params(length=2, width=0.6)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"P2fig7a_pub_gse311899.{ext}"), dpi=300)
    print("saved P2fig7a_pub_gse311899.png/.pdf")

if __name__ == "__main__":
    main()
