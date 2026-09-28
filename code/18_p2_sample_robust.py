#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""18_p2_sample_robust.py — P2 补强③：样本级稳健性
把细胞级核心结论(上皮防御塌陷/巨噬铁负荷/上皮活性升)落到样本水平:
每样本×细胞型中位数(既有 sample_medians.csv) -> 样本级 Welch t + Mann-Whitney + Hedges g。
AS 样本数与 CTL 样本数一目了然, 响应"AS 侧样本少"的审稿质疑。
输出: results/P2_sample_robust.csv + 汇总追加至 P2_hmgb1_summary.md 同级新 md
"""
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo"
d = pd.read_csv(f"{ROOT}/results/GSE215968_sample_medians.csv")

TESTS = [  # (celltype, module, 预期方向)
    ("Epithelium",      "ferro_defense", "down"),
    ("AS-Epithelium",   "ferro_defense", "down"),
    ("Epithelium",      "ferro_activity", "up"),
    ("Macrophages",     "iron_load", "up"),
    ("Macrophages",     "ferro_activity", "up"),
    ("Stromal",         "ferro_defense", "down"),
    ("Endothelium",     "ferro_defense", "down"),
]

def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.std(ddof=1) ** 2 + (nb - 1) * b.std(ddof=1) ** 2) / (na + nb - 2))
    return float((a.mean() - b.mean()) / sp) if sp > 0 else np.nan

rows = []
for ct, mod, direction in TESTS:
    sub = d[d.celltype == ct]
    a = sub[sub.group == "AS"][mod].values
    b = sub[sub.group == "WOI Control"][mod].values
    if len(a) < 2 or len(b) < 2:
        rows.append(dict(celltype=ct, module=mod, n_AS=len(a), n_CTL=len(b), note="insufficient"))
        continue
    t, tp = stats.ttest_ind(a, b, equal_var=False)
    u, up = stats.mannwhitneyu(a, b, alternative="two-sided")
    g = hedges_g(a, b)
    ok = (direction == "up" and a.mean() > b.mean()) or (direction == "down" and a.mean() < b.mean())
    rows.append(dict(celltype=ct, module=mod, n_AS_samples=len(a), n_CTL_samples=len(b),
                     mean_AS=round(float(a.mean()), 4), mean_CTL=round(float(b.mean()), 4),
                     welch_p=round(float(tp), 5), mwu_p=round(float(up), 5),
                     hedges_g=round(g, 2), direction_consistent=bool(ok)))
df = pd.DataFrame(rows)
df.to_csv(f"{ROOT}/results/P2_sample_robust.csv", index=False)
print(df.to_string(index=False))

n_samples = d.groupby("group")["sample"].nunique()
with open(f"{ROOT}/results/P2_sample_robust_summary.md", "w", encoding="utf-8") as fo:
    fo.write("# P2 补强③：样本级稳健性（sample-level validation）\n\n")
    fo.write(f"样本数（按 orig.ident）：{n_samples.to_dict()}\n\n")
    fo.write("细胞级核心结论在样本水平（每样本×细胞型中位数为一个观察）的复验：\n\n")
    fo.write(df.to_string(index=False))
    fo.write("\n\n**解读**：若各核心结论方向一致（direction_consistent=True）且 Welch/MWU p 可接受，"
             "说明细胞级结果非 AS 侧少数样本的聚类伪影。AS 侧样本量有限时以效应量 Hedges g 报告。\n")
print("saved P2_sample_robust.csv / P2_sample_robust_summary.md")
