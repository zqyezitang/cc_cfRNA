#!/usr/bin/env Rscript

suppressMessages({ library(edgeR) })
args <- commandArgs(trailingOnly = TRUE)
ARM <- args[1]
DATA <- "/work/zhouqing1/cfRNA_project/data"
OUT  <- "/work/zhouqing1/cfRNA_project/output"

counts <- read.csv(file.path(DATA, "counts_filtered.csv"), row.names = 1, check.names = FALSE)
meta   <- read.csv(file.path(DATA, "metadata_subgroup_bmi.csv"), stringsAsFactors = FALSE)

sams <- colnames(counts)
meta <- meta[match(sams, meta$core), ]
stopifnot(all(meta$core == sams))

# ）
meta <- meta[!is.na(meta$age), ]
counts <- counts[, meta$core]

# 
sel <- meta$tx_arm %in% c(ARM, "Control")
meta2 <- meta[sel, ]
counts2 <- counts[, meta2$core]

# 
meta2$grp <- ifelse(meta2$tx_arm == "Control", "Control",
              paste0(ARM, "_", meta2$group))
meta2$grp <- factor(meta2$grp, levels = c("Control", paste0(ARM, "_", c("Pre", "Post3M", "Post3Y"))))

print(table(meta2$grp))

y <- DGEList(counts = counts2, group = meta2$grp, samples = meta2[, c("core", "grp", "treatment", "age", "BMI")])
y <- calcNormFactors(y)
keep <- filterByExpr(y, min.count = 5, min.prop = 0.2)
y <- y[keep, , keep.lib.sizes = FALSE]


# 
design <- model.matrix(~ 0 + grp + age + BMI, data = meta2)
colnames(design) <- c(levels(meta2$grp), "age", "BMI")
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# 
contrasts <- list(
  "Pre_vs_Control"     = c(-1, 1, 0, 0, 0, 0),
  "Post3M_vs_Pre"      = c(0, -1, 1, 0, 0, 0),
  "Post3Y_vs_Pre"      = c(0, -1, 0, 1, 0, 0)
)
res <- list()
FC_TH <- log2(1.2)
for (nm in names(contrasts)) {
  ql <- glmQLFTest(fit, contrast = contrasts[[nm]])
  tt <- topTags(ql, n = Inf)
  df <- as.data.frame(tt$table); df$gene <- rownames(df)
  df <- df[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")]
  write.csv(df, file.path(OUT, paste0("DEG_", ARM, "_", nm, "_cov.csv")), row.names = FALSE)
  n_sig <- sum(df$FDR < 0.1 & abs(df$logFC) > FC_TH)
  n_up <- sum(df$FDR < 0.1 & df$logFC > FC_TH); n_dn <- sum(df$FDR < 0.1 & df$logFC < -FC_TH)
  # 
  old_f <- file.path(OUT, paste0("DEG_", ARM, "_", nm, ".csv"))
  if (file.exists(old_f)) {
    old <- read.csv(old_f)
    o_sig <- sum(old$FDR < 0.1 & abs(old$logFC) > FC_TH)
    
  } else {
    cat(sprintf("[%s] %s:  DEG=%d (up %d / down %d)\n", ARM, nm, n_sig, n_up, n_dn))
  }
  res[[nm]] <- df
}
saveRDS(res, file.path(OUT, paste0("DEG_", ARM, "_all_cov.rds")))
cat( ARM, "\n")


###################GO enrichment
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
OUT <- "/work/zhouqing1/cfRNA_project/output"
FC_TH <- log2(1.2)

sym2entrez <- function(genes) {
  mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
}

run_ora <- function(genes, label, out_prefix) {
  eg <- unique(unname(sym2entrez(genes)))
  eg <- eg[!is.na(eg)]
  cat(sprintf("\n[%s] symbols=%d, mapped Entrez=%d\n", label, length(genes), length(eg)))
  if (length(eg) < 10) { cat("  too few genes, skip\n"); return(invisible(NULL)) }
  go <- enrichGO(gene = eg, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
    write.csv(as.data.frame(go), file.path(OUT, paste0(out_prefix, "_GOBP.csv")), row.names = FALSE)
    
  } else cat("  GO BP: \n")
  kegg <- tryCatch(
    enrichKEGG(gene = eg, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) { cat("  KEGG 失败(可能无网络):", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    write.csv(as.data.frame(kegg), file.path(OUT, paste0(out_prefix, "_KEGG.csv")), row.names = FALSE)
    cat("  KEGG :", nrow(as.data.frame(kegg)), "\n")
  } else cat("  KEGG: \n")
}

for (arm in c("Surgery", "Chemo")) {
  for (nm in c("Pre_vs_Control", "Post3M_vs_Pre", "Post3Y_vs_Pre")) {
    f <- file.path(OUT, paste0("DEG_", arm, "_", nm, "_cov.csv"))
    if (!file.exists(f)) next
    d <- read.csv(f)
    sig <- d$gene[d$FDR < 0.1 & abs(d$logFC) > FC_TH]
    cat(sprintf("\n===== %s %s (cov):  DEG=%d =====", arm, nm, length(sig)))
    if (length(sig) >= 10) run_ora(sig, paste0(arm, "_", nm, "_cov"), paste0("enrich_", arm, "_", nm, "_cov"))
    else cat(sprintf("\n[%s_%s]  DEG=%d，不足 10 跳过\n", arm, nm, length(sig)))
  }
}
cat("\n enrich_*_cov_GOBP.csv / _KEGG.csv\n")