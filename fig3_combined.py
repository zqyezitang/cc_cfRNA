#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combined chemotherapy-arm figure script (main figure panels + Fig S1).

Main figure (chemotherapy arm):
  - clinical information heatmap (transposed)
  - DEG count barchart across three contrasts
  - GO BP combined dotplot across three periods
  - CR vs PR boxplot of the 5 response DEGs

Fig S1 (chemotherapy response details):
  - volcano plot (PR vs CR)
  - 5-DEG heatmap (pheatmap)
  - GO BP GSEA bar (PR vs CR)

Pipeline order (as referenced):
  10 arm-stratified GO/KEGG enrichment (non-cov) -> R script (embedded)
   1 clinical heatmap                                   -> Python
   2 DEG count barchart                                 -> Python
   3 cov-adjusted enrichment (cross-check)              -> R script (embedded)
   4 combined GO dotplot                                -> Python
   5 volcano plot                                       -> Python
   6 5-DEG pheatmap                                     -> R script (embedded)
   7 GSEA bar                                           -> Python
   8 CR vs PR boxplot                                   -> Python

Inputs : cfRNA_project/data tables, DEG/enrichment tables, Paper_TPM.csv, Samples_info.xlsx
Outputs: cfRNA_project/output/ (PNG/PDF/SVG)

Usage: python3 fig3_combined.py
"""
import pandas as pd, numpy as np, os, re, subprocess, tempfile
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
font_manager.fontManager.addfont("/public/skills/skill_tools/font/fonts/NotoSansCJKsc-Regular.otf")
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK SC"]
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams["font.size"] = 10
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import Normalize
import matplotlib.cm as cm

DATA = "/work/zhouqing1/cfRNA_project/data"
OUT  = "/work/zhouqing1/cfRNA_project/output"
TPM  = "/work/zhouqing1/uploads/Paper_TPM.csv"
FC_TH = np.log2(1.2)

# ---------- embedded R sources ----------

R_ARM_ENRICH = r'''
#!/usr/bin/env Rscript
# Step 1: arm-stratified GO/KEGG enrichment (non-cov)
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
DATA <- "/work/zhouqing1/cfRNA_project/data"
OUT  <- "/work/zhouqing1/cfRNA_project/output"
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
    cat("  GO BP significant terms:", nrow(as.data.frame(go)), "\n")
  } else cat("  GO BP: no significant terms\n")
  kegg <- tryCatch(
    enrichKEGG(gene = eg, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) { cat("  KEGG failed (network?):", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    write.csv(as.data.frame(kegg), file.path(OUT, paste0(out_prefix, "_KEGG.csv")), row.names = FALSE)
    cat("  KEGG significant terms:", nrow(as.data.frame(kegg)), "\n")
  } else cat("  KEGG: no significant terms or unavailable\n")
}

for (arm in c("Surgery", "Chemo")) {
  for (nm in c("Pre_vs_Control", "Post3M_vs_Pre", "Post3Y_vs_Pre")) {
    f <- file.path(OUT, paste0("DEG_", arm, "_", nm, ".csv"))
    if (!file.exists(f)) next
    d <- read.csv(f)
    sig <- d$gene[d$FDR < 0.1 & abs(d$logFC) > FC_TH]
    cat(sprintf("\n===== %s %s: significant DEGs=%d =====\n", arm, nm, length(sig)))
    if (length(sig) >= 10) run_ora(sig, paste0(arm, "_", nm), paste0("enrich_", arm, "_", nm))
    else cat(sprintf("\n[%s_%s] DEGs=%d, skipped (<10)\n", arm, nm, length(sig)))
  }
}
cat("\nArm-stratified enrichment done. Output: enrich_Surgery_*/enrich_Chemo_*_GOBP.csv and _KEGG.csv\n")
'''

R_COV_ENRICH = r'''
#!/usr/bin/env Rscript
# Cov-adjusted arm-stratified GO/KEGG enrichment (reads *_cov.csv)
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
    cat("  GO BP significant terms:", nrow(as.data.frame(go)), "\n")
  } else cat("  GO BP: no significant terms\n")
  kegg <- tryCatch(
    enrichKEGG(gene = eg, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) { cat("  KEGG failed (network?):", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    write.csv(as.data.frame(kegg), file.path(OUT, paste0(out_prefix, "_KEGG.csv")), row.names = FALSE)
    cat("  KEGG significant terms:", nrow(as.data.frame(kegg)), "\n")
  } else cat("  KEGG: no significant terms or unavailable\n")
}

for (arm in c("Surgery", "Chemo")) {
  for (nm in c("Pre_vs_Control", "Post3M_vs_Pre", "Post3Y_vs_Pre")) {
    f <- file.path(OUT, paste0("DEG_", arm, "_", nm, "_cov.csv"))
    if (!file.exists(f)) next
    d <- read.csv(f)
    sig <- d$gene[d$FDR < 0.1 & abs(d$logFC) > FC_TH]
    cat(sprintf("\n===== %s %s (cov): significant DEGs=%d =====\n", arm, nm, length(sig)))
    if (length(sig) >= 10) run_ora(sig, paste0(arm, "_", nm, "_cov"), paste0("enrich_", arm, "_", nm, "_cov"))
    else cat(sprintf("\n[%s_%s] DEGs=%d, skipped (<10)\n", arm, nm, length(sig)))
  }
}
cat("\nCov-adjusted enrichment done. Output: enrich_*_cov_GOBP.csv / _KEGG.csv\n")
'''

R_DEG5_HEATMAP = r'''
#!/usr/bin/env Rscript
# 5 response DEGs heatmap (Paper_TPM, log2(TPM+1) row z-score), pheatmap
# Columns: chemo post-treatment samples CR -> PR, top annotation aligned
.libPaths(c("/home/zhouqing1/software/R_libs", .libPaths()))
suppressMessages({ library(pheatmap); library(RColorBrewer); library(readxl) })

DATA <- "/work/zhouqing1/cfRNA_project/data"
OUT  <- "/work/zhouqing1/cfRNA_project/output"
TPM  <- "/work/zhouqing1/uploads/Paper_TPM.csv"

genes <- c("CYTB", "ND1", "LINC02018", "LOC124902198", "COX1")

meta <- read.csv(file.path(DATA, "metadata_subgroup_bmi.csv"), stringsAsFactors = FALSE)
raw  <- read_excel("/work/zhouqing1/uploads/Samples_info.xlsx", sheet = 1)
raw2 <- raw[, c("sample", "采血时评价")]
names(raw2)[2] <- "eval_raw"
meta <- merge(meta, raw2, by = "sample", all.x = TRUE)
if ("采血时评价" %in% names(meta)) {
  meta$`采血时评价` <- ifelse(is.na(meta$eval_raw), meta$`采血时评价`, meta$eval_raw)
  meta$eval_raw <- NULL
} else {
  names(meta)[names(meta) == "eval_raw"] <- "采血时评价"
}
clean_response <- function(v) {
  if (is.na(v) || trimws(as.character(v)) == "") return("NA")
  s <- trimws(as.character(v))
  m <- regmatches(s, regexpr("\\b(CR|PR|PD|SD)\\b", s))
  if (length(m) > 0) return(m[1])
  if (grepl("未", s)) return("NA")
  substr(s, 1, 10)
}
meta$response <- sapply(meta$`采血时评价`, clean_response)
sel <- meta$tx_arm == "Chemo" & meta$group != "Pre" & meta$response %in% c("CR", "PR")
meta2 <- meta[sel, , drop = FALSE]
tpm <- read.csv(TPM, row.names = 1, check.names = FALSE)
meta2 <- meta2[meta2$expr_col %in% colnames(tpm), , drop = FALSE]
cat("CR:", sum(meta2$response == "CR"), "PR:", sum(meta2$response == "PR"), "\n")

ord <- c("CR", "PR")
meta2 <- meta2[order(factor(meta2$response, levels = ord), meta2$core), , drop = FALSE]
samples <- meta2$expr_col
mat <- log2(tpm[genes, samples] + 1)

annot <- data.frame(Response = factor(meta2$response, levels = ord), row.names = samples)
ann_colors <- list(Response = c(CR = "#52BE80", PR = "#E67E22"))
gaps <- cumsum(as.integer(table(factor(meta2$response, levels = ord))))[1]

png(file.path(OUT, "Heatmap_Chemo_5DEG_CRvsPR_pheatmap.png"), width = 12, height = 5, units = "cm", res = 300)
pheatmap(mat,
         scale = "row",
         annotation_col = annot,
         annotation_colors = ann_colors,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE,
         fontsize_row = 8,
         fontsize = 7,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
         breaks = seq(-2.5, 2.5, length.out = 101),
         legend = TRUE,
         annotation_legend = TRUE,
         border_color = NA,
         gaps_col = gaps,
         main = NA)
dev.off()

pdf(file.path(OUT, "Heatmap_Chemo_5DEG_CRvsPR_pheatmap.pdf"), width = 12/2.54, height = 5/2.54)
pheatmap(mat,
         scale = "row",
         annotation_col = annot,
         annotation_colors = ann_colors,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE,
         fontsize_row = 8,
         fontsize = 7,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
         breaks = seq(-2.5, 2.5, length.out = 101),
         legend = TRUE,
         annotation_legend = TRUE,
         border_color = NA,
         gaps_col = gaps,
         main = NA)
dev.off()

cat("5DEG pheatmap saved (PNG/PDF)\n")
'''

def run_r(source):
    RSCRIPT = "/public/tools/r-packages/R-4.5.3/Rscript"
    fd, path = tempfile.mkstemp(suffix=".R")
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write(source)
    proc = subprocess.run([RSCRIPT, path], capture_output=True, text=True)
    if proc.stdout:
        print(proc.stdout, end="")
    os.unlink(path)
    if proc.returncode != 0:
        print("R ERROR:", proc.stderr, end="")
        raise RuntimeError("embedded R step failed")

# ---------- shared helpers ----------

def clean_response(v):
    if pd.isna(v) or str(v).strip() == "":
        return "NA"
    s = str(v).strip()
    m = re.search(r'\b(CR|PR|PD|SD)\b', s)
    return m.group(1) if m else ("NA" if "未" in s else s[:10])

def hpv_status(v):
    if pd.isna(v) or str(v).strip() == "" or str(v).strip().upper() in ("NA", "NAN", "未查", "/"):
        return "NA"
    if str(v).strip().lower() in ("negative", "阴性", "-"):
        return "Negative"
    return "Positive"

def load_meta_with_response():
    raw = pd.read_excel("/work/zhouqing1/uploads/Samples_info.xlsx", sheet_name=0)
    meta = pd.read_csv(f"{DATA}/metadata_subgroup_bmi.csv")
    raw_eval = raw[["sample", "采血时评价"]].rename(columns={"采血时评价": "eval_raw"})
    meta = meta.merge(raw_eval, on="sample", how="left")
    if "采血时评价" in meta.columns:
        meta["采血时评价"] = meta["eval_raw"].fillna(meta["采血时评价"])
    else:
        meta = meta.rename(columns={"eval_raw": "采血时评价"})
    meta = meta.drop(columns=["eval_raw"])
    return meta

# ---------- Fig (main) : clinical information heatmap ----------

def fig_clinical_heatmap():
    meta = load_meta_with_response()
    FIGO_ORDER = {"I": 0, "II": 1, "III": 2, "NA": 3}
    PERIOD_ORDER = {"Pre": 0, "Post3M": 1, "Post3Y": 2}

    def draw(arm, out_png, out_svg, out_pdf):
        s = meta[meta["tx_arm"] == arm].copy()
        s["hpv_status"] = s["HPV型别（多型别以逗号分割），阴性写negative，未测写NA"].apply(hpv_status)
        s["lymph"] = s["初治时淋巴结是否转移"].replace({"无": "否", "/": "NA", "未查": "NA"}).fillna("NA").astype(str)
        s["recur"] = s["是否复发（是/否）"].fillna("NA")
        s["figo"] = s["figo_roman"].fillna("NA")
        s["response"] = s["采血时评价"].apply(clean_response)
        s["_pk"] = s["group"].map(PERIOD_ORDER) * 10 + s["figo"].map(FIGO_ORDER) * 2 + np.arange(len(s)) * 0.001
        s = s.sort_values("_pk").reset_index(drop=True)
        s = s.drop(columns="_pk")
        n = len(s)

        period_cats = ["Pre", "Post3M", "Post3Y"]
        figo_cats = ["I", "II", "III", "NA"]
        lymph_cats = ["否", "是", "NA"]
        hpv_cats = ["Positive", "Negative", "NA"]
        resp_cats = ["CR", "PR", "PD", "NA"]
        recur_cats = ["否", "是", "NA"]

        def encode_cat(series, cats):
            m = {c: i for i, c in enumerate(cats)}
            return series.map(m).fillna(len(cats)).astype(int)

        M = np.zeros((n, 8))
        M[:, 0] = encode_cat(s["group"].astype(str), period_cats)
        M[:, 1] = (s["age"] - s["age"].min()) / (s["age"].max() - s["age"].min())
        M[:, 2] = (s["BMI"] - s["BMI"].min()) / (s["BMI"].max() - s["BMI"].min())
        M[:, 3] = encode_cat(s["figo"], figo_cats)
        M[:, 4] = encode_cat(s["lymph"], lymph_cats)
        M[:, 5] = encode_cat(s["hpv_status"], hpv_cats)
        M[:, 6] = encode_cat(s["response"], resp_cats)
        M[:, 7] = encode_cat(s["recur"], recur_cats)

        MT = M.T
        col_names = ["Period", "Age", "BMI", "FIGO", "Lymph", "HPV", "Response", "Recur"]
        palette = {
            "Pre": "#F1948A", "Post3M": "#F7DC6F", "Post3Y": "#85C1E9",
            "I": "#2ECC71", "II": "#F39C12", "III": "#E74C3C",
            "否": "#AED6F1", "是": "#C0392B",
            "Positive": "#E74C3C", "Negative": "#7FB3D5",
            "CR": "#52BE80", "PR": "#F7DC6F", "PD": "#E74C3C", "NA": "#BDC3C7",
        }
        cats_list = [period_cats, None, None, figo_cats, lymph_cats, hpv_cats, resp_cats, recur_cats]

        w_in = 20 / 2.54
        h_in = 10 / 2.54
        fig, ax = plt.subplots(figsize=(w_in, h_in))
        for j in range(8):
            for i in range(n):
                v = MT[j, i]
                if j in (1, 2):
                    c = cm.RdBu_r(v)
                else:
                    c = palette[cats_list[j][int(v)]]
                ax.add_patch(plt.Rectangle((i - 0.5, j - 0.5), 1, 1, facecolor=c, edgecolor="white", lw=0.3))

        bounds = [0]
        for g in period_cats:
            bounds.append(bounds[-1] + int((s["group"].astype(str) == g).sum()))
        for b in bounds[1:-1]:
            ax.axvline(b - 0.5, color="black", lw=1.2)

        ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(-0.5, 7.5)
        ax.set_yticks(range(8)); ax.set_yticklabels(col_names, fontsize=10)
        ax.set_xticks([])
        ax.tick_params(length=0)
        ax.set_title(f"{arm} arm (n={n})", fontsize=11)

        handles = []
        for label, key in [("Period: Pre", "Pre"), ("Period: Post3M", "Post3M"), ("Period: Post3Y", "Post3Y"),
                           ("FIGO: I", "I"), ("FIGO: II", "II"), ("FIGO: III", "III"),
                           ("Lymph: No", "否"), ("Lymph: Yes", "是"),
                           ("HPV: Positive", "Positive"), ("HPV: Negative", "Negative"),
                           ("Response: CR", "CR"), ("Response: PR", "PR"), ("Response: PD", "PD"),
                           ("Recur: No", "否"), ("Recur: Yes", "是"), ("NA", "NA")]:
            handles.append(Patch(facecolor=palette[key], label=label))
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                  fontsize=9, frameon=False, ncol=1, title="Categories", title_fontsize=10)
        sm = cm.ScalarMappable(cmap="RdBu_r", norm=Normalize(0, 1))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.15, location="bottom")
        cbar.set_label("Age / BMI (min-max normalized)", fontsize=9)
        cbar.ax.tick_params(labelsize=9)
        plt.tight_layout(rect=[0, 0.05, 0.82, 1])
        plt.savefig(out_png, dpi=300)
        plt.savefig(out_svg, dpi=150)
        plt.savefig(out_pdf)
        plt.close()
        print(f"Heatmap saved: {out_png} (20cm x 10cm transposed, 300dpi, n={n})")

    draw("Surgery", f"{OUT}/Surgery_clinical_heatmap_transpose_20x10cm.png", f"{OUT}/Surgery_clinical_heatmap_transpose_20x10cm.svg", f"{OUT}/Surgery_clinical_heatmap_transpose_20x10cm.pdf")
    draw("Chemo", f"{OUT}/Chemo_clinical_heatmap_transpose_20x10cm.png", f"{OUT}/Chemo_clinical_heatmap_transpose_20x10cm.svg", f"{OUT}/Chemo_clinical_heatmap_transpose_20x10cm.pdf")

# ---------- Fig (main) : DEG count barchart ----------

def fig_deg_bar():
    arm = "Chemo"
    contrasts = ["Pre_vs_Control", "Post3M_vs_Pre", "Post3Y_vs_Pre"]
    labels = ["Pre vs Control", "Post3M vs Pre", "Post3Y vs Pre"]
    ups, dns = [], []
    for nm in contrasts:
        d = pd.read_csv(f"{OUT}/DEG_{arm}_{nm}_cov.csv")
        sig = d["FDR"] < 0.1
        ups.append(int((sig & (d["logFC"] > FC_TH)).sum()))
        dns.append(int((sig & (d["logFC"] < -FC_TH)).sum()))
    print(f"Chemo: up={ups}, down={dns}")

    x = np.arange(3)
    ymax = max(ups) * 1.25
    ymin = -max(dns) * 1.25
    gap = max(max(ups), max(dns), 1) * 0.015
    fig, ax = plt.subplots(figsize=(9/2.54, 9/2.54))
    ax.bar(x, ups, width=0.55, color="#C0392B", label="Upregulated")
    ax.bar(x, [-v for v in dns], width=0.55, color="#2980B9", label="Downregulated")
    for i in range(3):
        ax.text(i, ups[i] + gap, str(ups[i]), ha="center", va="bottom", fontsize=12, color="#C0392B")
        ax.text(i, -dns[i] - gap, str(dns[i]), ha="center", va="top", fontsize=12, color="#2980B9")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=10, rotation=45, ha="right")
    ax.set_ylabel("Number of DEGs", fontsize=10)
    ax.set_ylim(-40, 120)
    ax.set_title("")
    ax.legend(frameon=False, fontsize=10, loc="upper left")
    plt.tight_layout()
    plt.savefig(f"{OUT}/Chemo_DEG_bar_cov_v3.png", dpi=300)
    plt.savefig(f"{OUT}/Chemo_DEG_bar_cov_v3.pdf")
    plt.close()
    print("DEG bar saved: Chemo_DEG_bar_cov_v3.png/pdf")

# ---------- Fig (main) : combined GO dotplot ----------

def parse_gene_ratio(s):
    try:
        a, b = str(s).split("/")
        return float(a) / float(b)
    except Exception:
        return np.nan

def fig_go_dotplot():
    CONTRAST_LABEL = {
        "Pre_vs_Control": "Pre vs Control",
        "Post3M_vs_Pre": "Post3M vs Pre",
        "Post3Y_vs_Pre": "Post3Y vs Pre",
    }

    def load_block(nm, topn=8):
        f = os.path.join(OUT, f"enrich_Chemo_{nm}_GOBP.csv")
        d = pd.read_csv(f).sort_values("pvalue").head(topn).copy()
        d["ratio"] = d["GeneRatio"].apply(parse_gene_ratio)
        d["neglogp"] = -np.log10(d["pvalue"])
        return d.sort_values("ratio")

    entries = []
    d1 = load_block("Pre_vs_Control", 8)
    entries.append((True, CONTRAST_LABEL["Pre_vs_Control"], np.nan, 0, 0))
    for _, r in d1.iterrows():
        entries.append((False, str(r["Description"])[:58] + ("..." if len(str(r["Description"])) > 58 else ""),
                        r["ratio"], r["Count"], r["neglogp"]))
    entries.append((True, "Post3M vs Pre (no significant DEGs)", np.nan, 0, 0))
    d2 = load_block("Post3Y_vs_Pre", 8)
    entries.append((True, CONTRAST_LABEL["Post3Y_vs_Pre"], np.nan, 0, 0))
    for _, r in d2.iterrows():
        entries.append((False, str(r["Description"])[:58] + ("..." if len(str(r["Description"])) > 58 else ""),
                        r["ratio"], r["Count"], r["neglogp"]))

    rows = [(i, e) for i, e in enumerate(entries)]
    y = np.array([i for i, e in rows])
    ratio = np.array([e[2] if not e[0] else np.nan for i, e in rows])
    count = np.array([e[3] for i, e in rows])
    neglogp = np.array([e[4] if not e[0] else 0 for i, e in rows])
    mask = ~np.isnan(ratio)

    fig, ax = plt.subplots(figsize=(12, 8.5))
    sc = ax.scatter(ratio[mask], y[mask], s=count[mask] * 6 + 20, c=neglogp[mask],
                    cmap="YlOrRd", edgecolor="grey", linewidth=0.5, alpha=0.92)
    ax.set_yticks(y)
    ax.set_yticklabels([e[1] for i, e in rows], fontsize=9)
    for i, e in rows:
        if e[0]:
            ax.axhspan(i - 0.5, i + 0.5, color="lightgrey", alpha=0.35, zorder=0)
            ax.get_yticklabels()[i].set_fontweight("bold")
            ax.get_yticklabels()[i].set_fontsize(10)
            ax.get_yticklabels()[i].set_color("#1F618D")
    ax.set_xlabel("GeneRatio", fontsize=10)
    ax.set_xlim(0, np.nanmax(ratio) * 1.22 + 0.02)
    ax.grid(axis="x", ls="--", alpha=0.35)
    cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("-log10(pvalue)", fontsize=10)
    cb.ax.tick_params(labelsize=9)
    ax.set_title("Chemo (CRT) arm: GO BP across three periods", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{OUT}/Chemo_GO_combined.png", dpi=300)
    plt.close()
    print(f"GO dotplot saved: Chemo_GO_combined.png ({len(entries)} entries)")

# ---------- Fig S1A : volcano plot ----------

def figs1_volcano():
    df = pd.read_csv(f"{OUT}/DEG_Chemo_CR_vs_PR.csv")
    sig = (df.FDR < 0.1) & (df.logFC.abs() > FC_TH)

    df["neg_log10fdr"] = -np.log10(df.FDR.clip(lower=1e-300))
    colors = ["#E74C3C" if (s and r > 0) else "#2E86C1" if (s and r < 0) else "#BDC3C7"
              for s, r in zip(sig, df.logFC)]

    fig, ax = plt.subplots(figsize=(12/2.54, 9/2.54))
    ax.scatter(df.logFC, df.neg_log10fdr, s=8, c=colors, alpha=0.7, linewidths=0)
    ax.axvline(FC_TH, color="#999999", lw=0.8, ls="--")
    ax.axvline(-FC_TH, color="#999999", lw=0.8, ls="--")
    ax.axhline(-np.log10(0.1), color="#999999", lw=0.8, ls="--")

    ax.set_xlabel("log2 FC (PR vs CR)")
    ax.set_ylabel("-log10 FDR")
    ax.set_title("Chemo arm: PR vs CR cfRNA DEGs", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(f"{OUT}/Volcano_Chemo_PR_vs_CR.png", dpi=300)
    plt.savefig(f"{OUT}/Volcano_Chemo_PR_vs_CR.pdf")
    plt.close()
    print(f"Volcano saved: Volcano_Chemo_PR_vs_CR.png/pdf (significant DEGs: {sig.sum()})")

# ---------- Fig S1C : GSEA bar ----------

def figs1_gsea_bar():
    gs = pd.read_csv(f"{OUT}/GSEA_Chemo_PR_vs_CR_GOBP.csv")
    s = gs[gs.qvalue < 0.15].copy()
    s = s.reindex(s.NES.abs().sort_values(ascending=False).index).head(18).iloc[::-1]
    labels = list(s.Description)

    n_row = len(s)
    fig, ax = plt.subplots(figsize=(20/2.54, max(10, 0.48 * n_row + 2)/2.54))
    colors = ["#E74C3C" if v > 0 else "#2E86C1" for v in s.NES]
    ax.barh(range(n_row), s.NES, color=colors, edgecolor="none")
    ax.set_yticks(range(n_row)); ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Normalized Enrichment Score (NES)")
    ax.set_title("Chemo arm GO BP GSEA: PR vs CR (q<0.15)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(f"{OUT}/GSEA_Chemo_PR_vs_CR_bar.png", dpi=300)
    plt.savefig(f"{OUT}/GSEA_Chemo_PR_vs_CR_bar.pdf")
    plt.close()
    print(f"GSEA bar saved (PNG+PDF); entries: {n_row}")

# ---------- Fig (main, bottom) : CR vs PR boxplot ----------

def figs1_boxplot():
    meta = load_meta_with_response()
    meta["response"] = meta["采血时评价"].apply(clean_response)
    sel = (meta.tx_arm == "Chemo") & (meta.group != "Pre") & (meta.response.isin(["CR", "PR"]))
    meta2 = meta[sel].copy().reset_index(drop=True)
    print("samples: CR", (meta2.response == "CR").sum(), "PR", (meta2.response == "PR").sum())

    log2cpm = pd.read_csv(f"{DATA}/log2cpm.csv", index_col=0)
    deg = pd.read_csv(f"{OUT}/DEG_Chemo_CR_vs_PR.csv")
    genes = deg[(deg.FDR < 0.1) & (deg.logFC.abs() > 0.263)]["gene"].tolist()
    print("genes:", genes)

    order = ["CR", "PR"]
    colors = {"CR": "#52BE80", "PR": "#E67E22"}

    fig, axes = plt.subplots(2, 3, figsize=(16/2.54, 10/2.54))
    axes = axes.flatten()
    for ax, g in zip(axes, genes):
        vals = {r: log2cpm.loc[g, meta2[meta2.response == r]["core"]].astype(float).values for r in order}
        bp = ax.boxplot([vals["CR"], vals["PR"]], tick_labels=order, widths=0.55,
                        patch_artist=True, showfliers=False, medianprops=dict(color="black", lw=1.2))
        for patch, r in zip(bp["boxes"], order):
            patch.set_facecolor(colors[r]); patch.set_alpha(0.75)
        for i, r in enumerate(order, start=1):
            jit = np.random.default_rng(42).normal(0, 0.06, len(vals[r]))
            ax.scatter(np.full(len(vals[r]), i) + jit, vals[r], s=12, color=colors[r],
                       edgecolor="white", lw=0.3, alpha=0.85, zorder=3)
        ax.set_title(g, fontsize=10)
        ax.set_ylabel("log2 CPM" if g == genes[0] else "", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=9)

    axes[5].axis("off")
    plt.tight_layout()
    plt.savefig(f"{OUT}/Boxplot_Chemo_DEG_CR_vs_PR_nopval.png", dpi=300)
    plt.savefig(f"{OUT}/Boxplot_Chemo_DEG_CR_vs_PR_nopval.svg", dpi=150)
    plt.close()
    print("CR vs PR boxplot saved (no p-value)")

# ---------- main ----------

if __name__ == "__main__":
    # 10: arm-stratified GO/KEGG enrichment (non-cov) -> supplies GO dotplot data
    run_r(R_ARM_ENRICH)
    # 1: clinical heatmap
    fig_clinical_heatmap()
    # 2: DEG count barchart
    fig_deg_bar()
    # 3: cov-adjusted enrichment (cross-check; Chemo cov has no significant GO)
    run_r(R_COV_ENRICH)
    # 4: combined GO dotplot
    fig_go_dotplot()
    # 5: volcano (Fig S1A)
    figs1_volcano()
    # 6: 5-DEG heatmap (Fig S1B)
    run_r(R_DEG5_HEATMAP)
    # 7: GSEA bar (Fig S1C)
    figs1_gsea_bar()
    # 8: CR vs PR boxplot (main figure bottom)
    figs1_boxplot()
    print("All chemotherapy-arm figures saved.")
