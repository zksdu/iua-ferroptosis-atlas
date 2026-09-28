#!/usr/bin/env Rscript
# 21_p2_gse311899_probe.R — 探查 GSE311899 RDS 对象结构 (类/维度/元数据列)
RDS <- "C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo/data/GSE311899_sc_AS_biop_pre_post.rds"
RES <- "C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo/results"

cat("== loading RDS (may take minutes) ==\n")
obj <- readRDS(RDS)
cat("class:", paste(class(obj), collapse=", "), "\n")
if (inherits(obj, "Seurat")) {
  cat("assays:", paste(Assays(obj), collapse=","), "\n")
  cat("dims:", dim(obj), "\n")
  md <- obj@meta.data
  cat("meta columns:", paste(colnames(md), collapse=" | "), "\n")
  write.csv(md, paste0(RES, "/P2_GSE311899_metadata_full.csv"), row.names=TRUE)
  # 打印每个字符列的取值分布 (前8个)
  for (cn in colnames(md)) {
    if (is.character(md[[cn]]) || is.factor(md[[cn]])) {
      tb <- table(md[[cn]], useNA="ifany")
      if (length(tb) <= 60) {
        cat("\n--", cn, "--\n"); print(tb)
      }
    }
  }
} else if (inherits(obj, "SingleCellExperiment")) {
  cat("SCE dims:", dim(obj), "\n")
  md <- as.data.frame(colData(obj))
  cat("meta columns:", paste(colnames(md), collapse=" | "), "\n")
  write.csv(md, paste0(RES, "/P2_GSE311899_metadata_full.csv"), row.names=TRUE)
} else if (is.list(obj)) {
  cat("list names:", paste(names(obj), collapse=" | "), "\n")
  str(obj, max.level=2, list.len=15)
} else {
  cat("dim:", paste(dim(obj), collapse=" x "), "\n")
  str(obj, max.level=1, list.len=10)
}
cat("\nPROBE_DONE\n")
