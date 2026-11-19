#!/usr/bin/env Rscript
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
#OUT <- "/work/zhouqing1/cfRNA_project/output"
args <- commandArgs(trailingOnly = TRUE)
c <- args[1]

f <- file.path(OUT, paste0("DEG_All_", c, ".csv"))
stopifnot(file.exists(f))
d <- read.csv(f)
d <- d[!is.na(d$logFC) & !is.na(d$gene) & d$gene != "", ]
d <- d[!duplicated(d$gene), ]
ranklist <- setNames(d$logFC, d$gene)
ranklist <- sort(ranklist, decreasing = TRUE)
eg <- mapIds(org.Hs.eg.db, keys = names(ranklist), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
names(ranklist) <- ifelse(is.na(eg), names(ranklist), eg)
ranklist <- ranklist[!is.na(ranklist) & !duplicated(names(ranklist))]

cat(sprintf("===== All %s:  %d =====\n", c, length(ranklist)))
gsea <- tryCatch(
  gseGO(geneList = ranklist, OrgDb = org.Hs.eg.db, ont = "BP", minGSSize = 15,
        maxGSSize = 300, pvalueCutoff = 0.25, seed = TRUE, nPermSimple = 500, verbose = FALSE),
  error = function(e) { cat("GSEA :", conditionMessage(e), "\n"); NULL }
)
out_csv <- file.path(OUT, paste0("GSEA_longitudinal_All_", c, "_GOBP.csv"))
if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
  gs <- as.data.frame(gsea)
  write.csv(gs, out_csv, row.names = FALSE)
  cat(sprintf("(q<0.25): %d\n", sum(gs$qvalue < 0.25)))
  top <- head(gs[order(gs$qvalue), c("ID", "Description", "setSize", "NES", "pvalue", "qvalue")], 10)
  print(top)
} else {
  cat("\n")
}
cat("\n")
