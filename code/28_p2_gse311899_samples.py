# -*- coding: utf-8 -*-
"""28_p2_gse311899_samples.py — GSE311899 per-sample pipeline (Plan B, memory-light)
Phases:
  python 28_..._samples.py download   # 35 tar.gz, 4 workers
  python 28_..._samples.py process    # stream each tar.gz -> 37-gene pseudobulk counts
  python 28_..._samples.py analyze    # log1p CPM -> module scores -> paired Wilcoxon
"""
import sys, os, gzip, tarfile, time, csv, io, subprocess, re
import numpy as np
import pandas as pd
from scipy import stats
from concurrent.futures import ThreadPoolExecutor

ROOT = "C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo"
DATA = os.path.join(ROOT, "data", "gse311899")
RES = os.path.join(ROOT, "results")
os.makedirs(DATA, exist_ok=True)
SOFT = os.path.join(ROOT, "data", "GSE311899_family.soft.gz")
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9334nnn/{gsm}/suppl/{fn}"

PRO_FERRO = ["ACSL4","LPCAT3","PTGS2","ALOX15","ALOX5","ALOXE3","TFRC","SAT1","NCOA4","SLC39A14","NOX1","NOX4","PLA2G6"]
DEFENSE   = ["GPX4","SLC7A11","SLC3A2","FTH1","FTL","DHODH","GCLC","GCLM","NQO1","SOD2","PRDX1","TXNRD1"]
IRON      = ["TFRC","FTH1","FTL","HMOX1","NCOA4","SLC40A1","SLC39A14"]
INFLAM    = ["IL1B","IL6","CCL2","CCL3","CCL4","CXCL8","TNF","NLRP3","SPP1","CD83"]
ALLG = sorted(set(PRO_FERRO + DEFENSE + IRON + INFLAM))

def filelist_map():
    """gsm -> (real filename, size) from series filelist.txt."""
    fl = os.path.join(ROOT, "data", "GSE311899_filelist.txt")
    if not os.path.exists(fl):
        subprocess.run(["curl", "-s", "--ssl-no-revoke", "-o", fl,
                        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311899/suppl/filelist.txt"],
                       check=True, timeout=120)
    m = {}
    with open(fl, encoding="utf-8", errors="replace") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 4 and p[0] == "File" and p[1].startswith("GSM"):
                m[p[1][:10]] = (p[1], int(p[3]))
    return m

def parse_soft():
    with gzip.open(SOFT, "rt", encoding="utf-8", errors="replace") as f:
        txt = f.read()
    samples = []
    cur = {}
    for line in txt.splitlines():
        if line.startswith("^SAMPLE"):
            if cur: samples.append(cur)
            cur = {}
        elif line.startswith("!Sample_geo_accession"):
            cur["gsm"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1"):
            v = line.split("=", 1)[1].strip()
            if v.startswith("biopsy origin:"):
                cur["comp"] = v.split(":", 1)[1].strip()
            elif v.startswith("individual:"):
                cur["patient"] = v.split(":", 1)[1].strip().replace("Patient ", "P")
            elif v.startswith("treatment stage:"):
                cur["stage"] = "Post" if "Post" in v else ("Pre" if "Pre" in v else v)
    if cur: samples.append(cur)
    df = pd.DataFrame(samples)
    fm = filelist_map()
    df["file"] = df.gsm.map(lambda g: fm[g][0])
    df["size"] = df.gsm.map(lambda g: fm[g][1])
    assert len(df) == 35, f"expected 35 samples, got {len(df)}"
    return df

def download(nw=4):
    df = parse_soft()
    df.to_csv(os.path.join(RES, "P2_GSE311899_sample_table.csv"), index=False)
    print(df.groupby(["comp", "stage"]).size().to_string(), flush=True)
    jobs = []
    for _, r in df.iterrows():
        dest = os.path.join(DATA, r["file"])
        if not (os.path.exists(dest) and os.path.getsize(dest) == r["size"]):
            jobs.append((r["gsm"], r["file"]))
    print(f"to download: {len(jobs)}", flush=True)
    def fetch(job):
        gsm, fn = job
        url = BASE.format(gsm=gsm, fn=fn)
        dest = os.path.join(DATA, fn)
        want = int(df.loc[df.gsm == gsm, "size"].iloc[0])
        for attempt in range(6):
            try:
                subprocess.run(["curl", "-sf", "--ssl-no-revoke", "--retry", "3", "-o", dest, url],
                               check=True, timeout=1800)
                if os.path.getsize(dest) != want:
                    print(f"sizemismatch {fn}: {os.path.getsize(dest)} != {want}", flush=True)
                    time.sleep(10); continue
                print(f"OK {fn} {want/1e6:.1f}MB", flush=True)
                return True
            except Exception as e:
                print(f"retry{attempt} {fn}: {type(e).__name__}", flush=True)
                time.sleep(20)
        return False
    ok = 0
    with ThreadPoolExecutor(max_workers=nw) as ex:
        for r in ex.map(fetch, jobs):
            ok += int(bool(r))
    print(f"downloaded {ok}/{len(jobs)}", flush=True)

def find_member(tf, pred):
    for m in tf.getmembers():
        if pred(m.name):
            return m
    return None

def process_sample(fn):
    """Stream tar.gz -> 37-gene sum vector. Returns dict or raises."""
    path = os.path.join(DATA, fn)
    with tarfile.open(path, "r:gz") as tf:
        names = tf.getnames()
        feats = [n for n in names if re.search(r"features\.tsv(\.gz)?$|genes\.tsv(\.gz)?$", n)]
        mtx = [n for n in names if n.endswith("matrix.mtx") or n.endswith("matrix.mtx.gz")]
        h5 = [n for n in names if n.endswith(".h5")]
        counts = None
        if h5:
            import h5py, tempfile, shutil
            tmp = os.path.join(DATA, "_tmp.h5")
            with tf.extractfile(h5[0]) as src, open(tmp, "wb") as out:
                shutil.copyfileobj(src, out)
            with h5py.File(tmp, "r") as f:
                g = f["matrix"]
                names_h5 = [x.decode() for x in (g["features/name"][:] if "name" in g["features"] else g["features/gene_names"][:])]
                shape = g["shape"][:]
                data = g["data"][:]; indices = g["indices"][:]; indptr = g["indptr"][:]
            os.remove(tmp)
            gi = {n: i for i, n in enumerate(names_h5)}
            sel = np.array([gi.get(x, -1) for x in ALLG])
            have = sel >= 0
            out = np.zeros(len(ALLG))
            # CSC: iterate columns, accumulate entries in selected rows
            ncol = shape[1]
            colsum_cells = np.zeros(ncol)  # not needed for pseudobulk sum
            # total counts per selected gene
            for j in range(ncol):
                s, e = indptr[j], indptr[j+1]
                rows = indices[s:e]; vals = data[s:e]
                m = have[rows] if False else None
            # vectorized: map all entries once
            mask = have[indices]
            np.add.at(out, sel[indices[mask]], vals[mask])
            counts = out
            ncells = ncol
        elif mtx:
            fm = find_member(tf, lambda n: n.endswith("matrix.mtx") or n.endswith("matrix.mtx.gz"))
            ff = feats[0] if feats else None
            if ff is None:
                raise RuntimeError("no features file")
            with tf.extractfile(ff) as fh:
                op = gzip.open if ff.endswith(".gz") else open
                with op(fh, "rt", encoding="utf-8", errors="replace") as fr:
                    gnames = [ln.split("\t")[1].strip() for ln in fr if ln.strip()]
            gi = {n: i for i, n in enumerate(gnames)}
            want = {x: k for k, x in enumerate(ALLG)}
            out = np.zeros(len(ALLG))
            with tf.extractfile(fm) as fh:
                op = gzip.open if fm.name.endswith(".gz") else open
                with op(fh, "rt", encoding="utf-8", errors="replace") as fr:
                    ncells = None
                    nrep = 0
                    for ln in fr:
                        if ln.startswith("%"):
                            continue
                        p = ln.split()
                        if ncells is None:
                            ncells = int(p[1])  # dims line: rows cols nnz
                            continue
                        if len(p) < 3: continue
                        r = int(p[0]) - 1; v = float(p[2])
                        g = gnames[r] if r < len(gnames) else None
                        if g in want:
                            out[want[g]] += v
            counts = out
        else:
            raise RuntimeError(f"unknown format: {names[:6]}")
    return counts, ncells

def process():
    df = parse_soft()
    cache_path = os.path.join(RES, "P2_GSE311899_pseudobulk_raw.csv")
    done = {}
    if os.path.exists(cache_path):
        try:
            prev = pd.read_csv(cache_path)
            done = {r["file"]: r for _, r in prev.iterrows()}
            print(f"cache: {len(done)} rows reused", flush=True)
        except Exception:
            done = {}
    rows = []
    fm = filelist_map()
    for _, r in df.iterrows():
        if r["file"] in done:
            rows.append(done[r["file"]].to_dict())
            continue
        fpath = os.path.join(DATA, r["file"])
        fm_entry = fm.get(r["gsm"])
        want = fm_entry[1] if fm_entry else None
        if want is not None and (not os.path.exists(fpath) or os.path.getsize(fpath) != want):
            print(f"skip (incomplete): {r['file']}", flush=True)
            continue
        t0 = time.time()
        counts, ncells = process_sample(r["file"])
        rows.append(dict(gsm=r["gsm"], title=r["title"], patient=r["patient"],
                         comp=r["comp"], stage=r["stage"], file=r["file"],
                         n_cells=ncells, **{g: counts[k] for k, g in enumerate(ALLG)}))
        # incremental save after every sample
        pd.DataFrame(rows).to_csv(cache_path, index=False)
        print(f"{r['file']}: {ncells} cells, {time.time()-t0:.0f}s [saved]", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(cache_path, index=False)
    print("saved P2_GSE311899_pseudobulk_raw.csv", flush=True)

def analyze():
    df = pd.read_csv(os.path.join(RES, "P2_GSE311899_pseudobulk_raw.csv"))
    G = [c for c in df.columns if c in ALLG]
    pb = df[G].to_numpy(float).T  # genes x samples
    lib = pb.sum(0); lib[lib == 0] = 1
    logcpm = np.log1p(pb / lib * 1e6)
    L = pd.DataFrame(logcpm, index=G, columns=df["title"])
    mods = dict(PRO_FERRO=PRO_FERRO, DEFENSE=DEFENSE, IRON=IRON, INFLAM=INFLAM)
    sc = pd.DataFrame({m: L.loc[[g for g in gs if g in L.index]].mean(0) for m, gs in mods.items()})
    sc.insert(0, "stage", df["stage"].values); sc.insert(0, "comp", df["comp"].values)
    sc.insert(0, "patient", df["patient"].values)
    sc.to_csv(os.path.join(RES, "P2_GSE311899_module_scores.csv"), index=False)

    # paired Pre vs Post per compartment
    rows = []
    for cp in sc["comp"].unique():
        d = sc[sc["comp"] == cp]
        piv = d.pivot(index="patient", columns="stage", values=list(mods.keys()))
        piv = piv.dropna()
        for m in mods:
            pre, post = piv[(m, "Pre")], piv[(m, "Post")]
            wt = stats.wilcoxon(post, pre)
            tt = stats.ttest_rel(post, pre)
            rows.append(dict(compartment=cp, module=m, n_pairs=len(piv),
                             mean_Pre=pre.mean(), mean_Post=post.mean(),
                             delta_Post_minus_Pre=post.mean() - pre.mean(),
                             wilcox_p=wt.pvalue, paired_t_p=tt.pvalue))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(RES, "P2_GSE311899_prepost_stats.csv"), index=False)
    print(res.to_string(index=False), flush=True)

    # per-gene direction (full 37-gene panel; FigS2 consumes all modules)
    grows = []
    Ldf = L.T
    Ldf.insert(0, "stage", df["stage"].values); Ldf.insert(0, "comp", df["comp"].values)
    Ldf.insert(0, "patient", df["patient"].values)
    for cp in sc["comp"].unique():
        d = Ldf[Ldf["comp"] == cp]
        for g in ALLG:
            if g not in d.columns: continue
            piv = d.pivot(index="patient", columns="stage", values=g).dropna()
            delta = piv["Post"] - piv["Pre"]
            grows.append(dict(compartment=cp, gene=g, n_pairs=len(piv),
                              mean_delta=delta.mean(), frac_up=float((delta > 0).mean())))
    gdf = pd.DataFrame(grows)
    gdf.to_csv(os.path.join(RES, "P2_GSE311899_gene_direction.csv"), index=False)
    print(gdf.to_string(index=False), flush=True)
    print("ANALYZE_DONE", flush=True)

if __name__ == "__main__":
    ph = sys.argv[1] if len(sys.argv) > 1 else "download"
    dict(download=download, process=process, analyze=analyze)[ph]()
