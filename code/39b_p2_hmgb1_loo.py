# -*- coding: utf-8 -*-
"""39b: LOO sensitivity for the 20-gene HMGB1 axis in macrophages (adds to P2_LOO_full).
Uses same z-once scheme; only the HMGB1_AXIS block is extracted.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import h5py
from scipy import stats

ROOT = r"C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo"
H5AD = ROOT + "/data/GSE215968_sc.h5ad"
HMGB1_AXIS = ["HMGB1","TLR2","TLR4","TLR9","AGER","NLRP3","AIM2","PYCARD","CASP1","CASP4",
              "CASP5","IL1B","IL18","MYD88","TICAM1","RELA","NFKB1","NFKBIA","CXCL8","CCL2"]

f = h5py.File(H5AD, "r")
n_cells = f["X/indptr"].shape[0] - 1
group = np.array([x.decode() for x in f["obs/Group"][:]])
ctype = np.array([x.decode() for x in f["obs/principal_cell_types"][:]])
genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var/features"][:]])
gset = set(genes.tolist())
wanted = [g for g in HMGB1_AXIS if g in gset]
col_idx = np.array([int(np.where(genes == g)[0][0]) for g in wanted])
X = f["X"]; indptr = X["indptr"][:]
M = np.zeros((n_cells, len(wanted)), dtype=np.float32)
posmap = {int(c): i for i, c in enumerate(col_idx)}
for s in range(0, n_cells, 5000):
    e = min(s + 5000, n_cells)
    i0, i1 = int(indptr[s]), int(indptr[e])
    idx, val = X["indices"][i0:i1], X["data"][i0:i1]
    lp = indptr[s:e + 1].astype(np.int64) - i0
    rows = np.repeat(np.arange(s, e, dtype=np.int64), np.diff(lp))
    mask = np.isin(idx, col_idx)
    if mask.any():
        h_idx, h_val, h_row = idx[mask], val[mask], rows[mask]
        pos = np.fromiter((posmap[int(c)] for c in h_idx), dtype=np.int64, count=len(h_idx))
        M[h_row, pos] = h_val
f.close()

Z = np.empty_like(M)
for j in range(M.shape[1]):
    col = M[:, j]; sd = col.std()
    Z[:, j] = (col - col.mean()) / sd if sd > 1e-8 else 0.0
gi = {g: i for i, g in enumerate(wanted)}
as_mask = group == "AS"; ctl_mask = group == "WOI Control"
sel = ctype == "Macrophages"
present = [g for g in HMGB1_AXIS if g in gi]

rows = []
for drop in [None] + present:
    names = [g for g in present if g != drop]
    v = Z[:, [gi[g] for g in names]].mean(axis=1)
    u, p = stats.mannwhitneyu(v[sel & as_mask], v[sel & ctl_mask], alternative="two-sided")
    delta = float(np.median(v[sel & as_mask]) - np.median(v[sel & ctl_mask]))
    rows.append(dict(module="HMGB1_AXIS", celltype="Macrophages", expected="up",
                     dropped_gene="NONE(full)" if drop is None else drop,
                     n_genes=len(names), delta_AS_minus_CTL=round(delta, 4), pval=p))
loo = pd.DataFrame(rows)
loo["direction_ok"] = loo.delta_AS_minus_CTL > 0
loo["sig_p05"] = loo.pval < 0.05
loo["robust"] = loo.direction_ok & loo.sig_p05
old = pd.read_csv(ROOT + "/results/P2_LOO_full.csv")
allloo = pd.concat([old, loo], ignore_index=True)
allloo.to_csv(ROOT + "/results/P2_LOO_full.csv", index=False)
sub = loo[loo.dropped_gene != "NONE(full)"]
print("HMGB1 macrophage LOO: dir_ok %d/%d, robust %d/%d" %
      (sub.direction_ok.sum(), len(sub), sub.robust.sum(), len(sub)))
print(loo.to_string(index=False))
