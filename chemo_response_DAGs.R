#!/usr/bin/env Rscript
# 
suppressMessages({ library(edgeR) })
DATA <- "/work/zhouqing1/cfRNA_project/data"
OUT  <- "/work/zhouqing1/cfRNA_project/output"

counts <- read.csv(file.path(DATA, "counts_filtered.csv"), row.names = 1, check.names = FALSE)
meta   <- read.csv(file.path(DATA, "metadata_subgroup_bmi.csv"), stringsAsFactors = FALSE)
raw    <- readxl::read_excel("/work/zhouqing1/uploads/Samples_info.xlsx", sheet = 1)

clean_response <- function(v) {
  if (is.na(v) || trimws(as.character(v)) == "") return("NA")
  s <- trimws(as.character(v))
  m <- regmatches(s, regexpr("\\b(CR|PR|PD|SD)\\b", s))
  if (length(m) > 0) return(m[1])
  if (grepl("", s)) return("NA")
  substr(s, 1, 10)
}
meta$response <- sapply(meta$response, clean_response)

# 
sel <- meta$tx_arm == "Chemo" & meta$group != "Pre" & meta$response %in% c("CR", "PR")
meta2 <- meta[sel, ]
meta2 <- meta2[!is.na(meta2$age) & !is.na(meta2$BMI), ]
counts2 <- counts[, meta2$core]
stopifnot(all(colnames(counts2) == meta2$core))


meta2$response <- factor(meta2$response, levels = c("CR", "PR"))
meta2$group    <- factor(meta2$group, levels = c("Post3M", "Post3Y"))
meta2$age_c    <- meta2$age - mean(meta2$age, na.rm = TRUE)
meta2$bmi_c    <- meta2$BMI - mean(meta2$BMI, na.rm = TRUE)

y <- DGEList(counts = counts2, group = meta2$response, samples = meta2[, c("core", "response", "group", "age_c", "bmi_c")])
y <- calcNormFactors(y)
keep <- filterByExpr(y, min.count = 5, min.prop = 0.2)
y <- y[keep, , keep.lib.sizes = FALSE]
cat("", nrow(y), "\n")

design <- model.matrix(~ response + age_c + bmi_c + group, data = meta2)
colnames(design)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# 
ql <- glmQLFTest(fit, coef = "responsePR")
tt <- topTags(ql, n = Inf)
df <- as.data.frame(tt$table); df$gene <- rownames(df)
df <- df[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")]
# 
df$direction <- ifelse(df$logFC > 0, "PR_up", "CR_up")
write.csv(df, file.path(OUT, "DEG_Chemo_CR_vs_PR.csv"), row.names = FALSE)

# 
FC_TH <- log2(1.2)
sig  <- df$FDR < 0.1 & abs(df$logFC) > FC_TH
cat(sprintf("FDR<0.1 & fold-change>1.2: DEG = %d（PR up %d / CR up %d）\n",
            sum(sig), sum(sig & df$logFC > 0), sum(sig & df$logFC < 0)))
sig5 <- df$FDR < 0.05 & abs(df$logFC) > FC_TH
cat(sprintf("FDR<0.05 : %d\n", sum(sig5)))
saveRDS(list(res = df, n_keep = nrow(y), n_cr = sum(meta2$response == "CR"), n_pr = sum(meta2$response == "PR")),
        file.path(OUT, "DEG_Chemo_CR_vs_PR.rds"))
cat("\n")
