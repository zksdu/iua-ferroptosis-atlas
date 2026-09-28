#!/usr/bin/env Rscript
# 22_p2_gse311899_analysis.R — GSE311899 治疗拯救验证 (自适应对象结构)
# 设计: 9 例中重度 AS, CD133+ BMDSC 治疗 Pre vs Post 配对, Epithelium/Stroma 分室
# 分析: 37 基因铁死亡 4 模块(42 计 5 跨模块重复) -> 每样本拟批量(log1p CPM) -> 成对 Wilcoxon + 逐基因方向
# 依赖: 先跑 21_probe.R 确认结构; 本脚本含 Seurat/SCE 双分支
suppressWarnings(suppressMessages({
  library(methods)
}))

ROOT <- "C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo"
RDS  <- file.path(ROOT, "data/GSE311899_sc_AS_biop_pre_post.rds")
RES  <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

PRO_FERRO <- c("ACSL4","LPCAT3","PTGS2","ALOX15","ALOX5","ALOXE3","TFRC","SAT1","NCOA4","SLC39A14","NOX1","NOX4","PLA2G6")
DEFENSE   <- c("GPX4","SLC7A11","SLC3A2","FTH1","FTL","DHODH","GCLC","GCLM","NQO1","SOD2","PRDX1","TXNRD1")
IRON      <- c("TFRC","FTH1","FTL","HMOX1","NCOA4","SLC40A1","SLC39A14")
INFLAM    <- c("IL1B","IL6","CCL2","CCL3","CCL4","CXCL8","TNF","NLRP3","SPP1","CD83")
ALLG <- unique(c(PRO_FERRO, DEFENSE, IRON, INFLAM))

cat("== readRDS ==\n"); t0 <- Sys.time()
obj <- readRDS(RDS)
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="mins")),1), "min\n")
cat("class:", paste(class(obj), collapse=","), "\n")

# ---- 提取 counts + metadata (自适应) ----
if (inherits(obj, "Seurat")) {
  md  <- obj@meta.data
  assay <- if ("RNA" %in% Assays(obj)) "RNA" else Assays(obj)[1]
  cat("assay:", assay, "\n")
  counts <- GetAssayData(obj, assay = assay, slot = "counts")
  genes <- rownames(counts)
} else if (inherits(obj, "SingleCellExperiment")) {
  suppressMessages(library(SingleCellExperiment))
  md <- as.data.frame(colData(obj)); md$cell_id <- rownames(md)
  counts <- counts(obj)
  genes <- rownames(counts)
} else {
  stop("unexpected object class: ", paste(class(obj), collapse=","))
}

# ---- 识别元数据列 (自适应: patient / stage / compartment) ----
pick_col <- function(pats, md) {
  for (p in pats) {
    hit <- grep(p, colnames(md), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  NA_character_
}
col_pat  <- pick_col(c("individual","patient","sample_id","orig.ident"), md)
col_stg  <- pick_col(c("stage","treatment","timepoint","pre|post"), md)
col_cmp  <- pick_col(c("origin","compartment","biopsy","layer"), md)
# 兜底: 找不到的列用常量填充, 保证 data.frame 不崩
if (is.na(col_pat)) { md$.patient_fb <- "P1"; col_pat <- ".patient_fb" }
if (is.na(col_stg)) { md$.stage_fb   <- "Pre"; col_stg <- ".stage_fb" }
if (is.na(col_cmp)) { md$.comp_fb    <- "all";  col_cmp <- ".comp_fb" }
cat("meta cols -> patient:", col_pat, "| stage:", col_stg, "| compartment:", col_cmp, "\n")
for (cn in c(col_pat, col_stg, col_cmp)) {
  if (!is.na(cn)) { cat("\n--", cn, "--\n"); print(table(md[[cn]], useNA="ifany")) }
}

getv <- function(cn) {
  v <- md[[cn]]
  if (is.factor(v)) v <- as.character(v)
  as.character(v)
}

# ---- 基因子集 ----
gset <- intersect(ALLG, genes)
miss <- setdiff(ALLG, gset)
cat("genes found:", length(gset), "/", length(ALLG), "| missing:", paste(miss, collapse=","), "\n")

# ---- 每样本拟批量: sum counts per (patient x stage x compartment) ----
grp <- apply(data.frame(p=getv(col_pat), s=getv(col_stg), c=getv(col_cmp)),
             1, paste, collapse="||")
ug <- unique(grp)
cat("groups:", length(ug), "\n")
gi <- match(gset, genes)
sub <- counts[gi, , drop = FALSE]  # 稀疏矩阵行子集
# 逐组列求和 (Matrix colSums on index; 用稀疏列切片避免大临时)
pb <- matrix(0, nrow = length(gset), ncol = length(ug),
             dimnames = list(gset, ug))
grpF <- factor(grp, levels = ug)
cs <- table(grpF)
for (g in ug) {
  idx <- which(grp == g)
  if (length(idx) == 1) pb[, g] <- as.numeric(sub[, idx])
  else pb[, g] <- as.numeric(Matrix::rowSums(sub[, idx, drop = FALSE]))
}
write.csv(t(pb), file.path(RES, "P2_GSE311899_pseudobulk_raw.csv"))

# ---- log1p CPM ----
lib <- colSums(pb); lib[lib == 0] <- 1
cpm <- sweep(pb, 2, lib, "/") * 1e6
logcpm <- log1p(cpm)

# ---- 模块得分: 模块内基因 log-CPM 均值 ----
mod_list <- list(PRO_FERRO = PRO_FERRO, DEFENSE = DEFENSE, IRON = IRON, INFLAM = INFLAM)
score <- sapply(names(mod_list), function(mn) {
  gs <- intersect(mod_list[[mn]], gset)
  if (length(gs) == 0) return(rep(NA, ncol(logcpm)))
  colMeans(logcpm[gs, , drop = FALSE])
})
info <- data.frame(group = ug)
info$patient <- sapply(info$group, function(x) strsplit(x, "\\|\\|")[[1]][1])
info$stage   <- sapply(info$group, function(x) strsplit(x, "\\|\\|")[[1]][2])
info$comp    <- sapply(info$group, function(x) strsplit(x, "||", fixed = TRUE)[[1]][3])
out <- cbind(info, score)
write.csv(out, file.path(RES, "P2_GSE311899_module_scores.csv"), row.names = FALSE)

# ---- 归一化 stage 值: Pre/Post ----
std_stage <- function(x) {
  x <- tolower(trimws(x))
  ifelse(grepl("pre", x), "Pre", ifelse(grepl("post", x), "Post", x))
}
out$stage_std <- std_stage(out$stage)

# ---- 成对分析: 每 compartment, Pre vs Post 成对 Wilcoxon ----
res_rows <- list()
for (cp in unique(out$comp)) {
  d <- out[out$comp == cp & out$stage_std %in% c("Pre","Post"), ]
  piv <- reshape(d[, c("patient","stage_std", names(mod_list))],
                 idvar = "patient", timevar = "stage_std", direction = "wide")
  piv <- piv[complete.cases(piv), ]
  if (nrow(piv) < 3) { cat("comp", cp, "pairs:", nrow(piv), "skip\n"); next }
  for (mn in names(mod_list)) {
    pre  <- piv[[paste0(mn, ".Pre")]]
    post <- piv[[paste0(mn, ".Post")]]
    wt <- suppressWarnings(wilcox.test(pre, post, paired = TRUE))
    tt <- suppressWarnings(t.test(pre, post, paired = TRUE))
    res_rows[[length(res_rows)+1]] <- data.frame(
      compartment = cp, module = mn, n_pairs = nrow(piv),
      mean_Pre = mean(pre), mean_Post = mean(post),
      delta_Post_minus_Pre = mean(post) - mean(pre),
      wilcox_p = wt$p.value, paired_t_p = tt$p.value)
  }
}
res_df <- do.call(rbind, res_rows)
write.csv(res_df, file.path(RES, "P2_GSE311899_prepost_stats.csv"), row.names = FALSE)
cat("\n== Pre vs Post paired results ==\n"); print(res_df)

# ---- 逐基因方向一致性 (DEFENSE 重点) ----
gene_rows <- list()
for (cp in unique(out$comp)) {
  d <- out$group[grepl(cp, out$group, fixed = TRUE)]
  for (g in intersect(DEFENSE, gset)) {
    dd <- out[out$comp == cp, c("patient","stage_std", g)]
    piv <- reshape(dd, idvar="patient", timevar="stage_std", direction="wide")
    piv <- piv[complete.cases(piv), ]
    if (nrow(piv) < 3) next
    delta <- piv[[paste0(g, ".Post")]] - piv[[paste0(g, ".Pre")]]
    gene_rows[[length(gene_rows)+1]] <- data.frame(
      compartment = cp, gene = g, n_pairs = nrow(piv),
      mean_delta = mean(delta), frac_up = mean(delta > 0))
  }
}
gdf <- do.call(rbind, gene_rows)
write.csv(gdf, file.path(RES, "P2_GSE311899_gene_direction.csv"), row.names = FALSE)
cat("\n== defense gene direction (Post-Pre) ==\n"); print(gdf)

cat("\nANALYSIS_DONE\n")
