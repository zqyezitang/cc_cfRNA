#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combined Figure 4 generation script.

Fig 4A: Sankey diagram of the 28 chemo-arm Pre-up genes (Y3 outcomes)
Fig 4B: 13-gene 7-group mean-expression heatmap (R/pheatmap, embedded)
Fig 4C/4D: UpSet plots of time-course DEG sets (R/UpSetR, embedded)
Fig 4E/4G: four-time-point boxplots of MEMO1 and MMP11 (parameterised)

Pipeline order (as referenced):
  1 Sankey                                    -> Python
  2 13-gene heatmap                           -> R script (embedded)
  3 UpSet                                     -> R script (embedded)
  4 four-time-point gene boxplot              -> Python

Inputs : cfRNA_project/data tables, DEG tables, Paper_TPM.csv
Outputs: cfRNA_project/output/ (PNG/PDF/SVG)

Usage: python3 fig4_combined.py
"""
import pandas as pd, numpy as np, os, subprocess, tempfile
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
font_manager.fontManager.addfont("/public/skills/skill_tools/font/fonts/NotoSansCJKsc-Regular.otf")
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK SC"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as mpatches
from matplotlib.patches import Patch

DATA = "/work/zhouqing1/cfRNA_project/data"
OUT  = "/work/zhouqing1/cfRNA_project/output"
TPM  = "/work/zhouqing1/uploads/Paper_TPM.csv"

# ---------- embedded R sources ----------

R_13GENES_HEATMAP = r'''
#!/usr/bin/env Rscript
# 13 overlapping-gene expression heatmap: 7 groups (Control / Pre-M3-Y3 x Surgery-Chemo),
# each group summarised by the mean of log2(TPM+1)
.libPaths(c("/home/zhouqing1/software/R_libs", .libPaths()))
suppressMessages({ library(pheatmap); library(RColorBrewer) })

DATA <- "/work/zhouqing1/cfRNA_project/data"
OUT  <- "/work/zhouqing1/cfRNA_project/output"
TPM  <- "/work/zhouqing1/uploads/Paper_TPM.csv"

genes <- c("ANKRD28","AVPR1A","CLEC1B","F2R","IGF2BP3","LOC105371661","LOC105377384",
           "SNORA104","TRND","TRNK","TRNQ","TRNR","TRNV")

meta <- read.csv(file.path(DATA, "metadata_subgroup_bmi.csv"), stringsAsFactors = FALSE)
tpm  <- read.csv(TPM, row.names = 1, check.names = FALSE)
meta <- meta[meta$expr_col %in% colnames(tpm), , drop = FALSE]

# 7-group definition (Control single; others group x tx_arm, drop Other)
meta$grp <- ifelse(meta$group == "Control", "Control",
            ifelse(meta$tx_arm == "Surgery", paste0(meta$group, "_Surgery"),
            ifelse(meta$tx_arm == "Chemo",    paste0(meta$group, "_Chemo"), NA)))
ord <- c("Control", "Pre_Surgery", "Pre_Chemo", "Post3M_Surgery", "Post3M_Chemo",
         "Post3Y_Surgery", "Post3Y_Chemo")
meta <- meta[!is.na(meta$grp) & meta$grp %in% ord, , drop = FALSE]

# per-group mean of log2(TPM+1) (13 genes x 7 groups)
mat <- sapply(ord, function(g) {
  cols <- meta$expr_col[meta$grp == g]
  rowMeans(log2(tpm[genes, cols] + 1), na.rm = TRUE)
})
rownames(mat) <- genes
colnames(mat) <- ord
cat("group sample counts:", table(factor(meta$grp, levels = ord)), "\n")

ann_col <- data.frame(Group = factor(colnames(mat), levels = ord), row.names = colnames(mat))
ann_colors <- list(Group = c(Control = "#95A5A6",
                             Pre_Surgery = "#85C1E9", Pre_Chemo = "#2E86C1",
                             Post3M_Surgery = "#82E0AA", Post3M_Chemo = "#229954",
                             Post3Y_Surgery = "#F0B27A", Post3Y_Chemo = "#CA6F1E"))

png(file.path(OUT, "Heatmap_13genes_7groups_mean_pheatmap.png"), width = 12, height = 8, units = "cm", res = 300)
pheatmap(mat,
         scale = "row",
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         cluster_rows = FALSE, cluster_cols = FALSE,
         fontsize_row = 7, fontsize = 7,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
         breaks = seq(-2.5, 2.5, length.out = 101),
         legend = TRUE, annotation_legend = TRUE,
         border_color = NA, main = NA)
dev.off()

pdf(file.path(OUT, "Heatmap_13genes_7groups_mean_pheatmap.pdf"), width = 12/2.54, height = 8/2.54)
pheatmap(mat,
         scale = "row",
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         cluster_rows = FALSE, cluster_cols = FALSE,
         fontsize_row = 7, fontsize = 7,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
         breaks = seq(-2.5, 2.5, length.out = 101),
         legend = TRUE, annotation_legend = TRUE,
         border_color = NA, main = NA)
dev.off()

cat("13-gene 7-group mean heatmap saved (PNG/PDF)\n")
'''

R_UPSET = r'''
#!/usr/bin/env Rscript
# UpSet of time-course DEG sets: Surgery M3/Y3 vs Chemo Y3 (cov-adjusted)
.libPaths(c("/home/zhouqing1/software/R_libs", .libPaths()))
suppressMessages(library(UpSetR))

out  <- "/work/zhouqing1/cfRNA_project/output"
fc   <- log2(1.2)

read_sets <- function(f) {
  d <- read.csv(f, check.names = FALSE)
  up   <- d$gene[d$FDR < 0.1 & d$logFC >  fc]
  down <- d$gene[d$FDR < 0.1 & d$logFC < -fc]
  list(up = up, down = down)
}

sm <- read_sets(file.path(out, "DEG_Surgery_Post3M_vs_Pre_cov.csv"))
sy <- read_sets(file.path(out, "DEG_Surgery_Post3Y_vs_Pre_cov.csv"))
cy <- read_sets(file.path(out, "DEG_Chemo_Post3Y_vs_Pre_cov.csv"))

cat("Surgery M3: up", length(sm$up), "down", length(sm$down), "\n")
cat("Surgery Y3: up", length(sy$up), "down", length(sy$down), "\n")
cat("Chemo   Y3: up", length(cy$up), "down", length(cy$down), "\n")

plot_upset <- function(sets_list, sets_names, main_color, file_base) {
  lst <- setNames(lapply(sets_names, function(n) sets_list[[n]]), sets_names)
  names(lst) <- paste0(names(lst), " (", sapply(lst, length), ")")
  png(paste0(file_base, ".png"), width = 10, height = 7.5, units = "cm", res = 300)
  print(upset(fromList(lst),
        nsets = 3, nintersects = 20, order.by = "freq", decreasing = TRUE,
        mb.ratio = c(0.6, 0.4), point.size = 3.2, line.size = 1.4,
        sets.bar.color = "#9E9E9E", main.bar.color = main_color,
        matrix.color = "#444444",
        set_size.show = FALSE,
        mainbar.y.label = "Intersection size",
        sets.x.label = "Set size",
        text.scale = c(1.3, 1.3, 1.15, 1.15, 1.35, 1.45)))
  dev.off()
  pdf(paste0(file_base, ".pdf"), width = 10/2.54, height = 7.5/2.54)
  print(upset(fromList(lst),
        nsets = 3, nintersects = 20, order.by = "freq", decreasing = TRUE,
        mb.ratio = c(0.6, 0.4), point.size = 3.2, line.size = 1.4,
        sets.bar.color = "#9E9E9E", main.bar.color = main_color,
        matrix.color = "#444444",
        set_size.show = FALSE,
        mainbar.y.label = "Intersection size",
        sets.x.label = "Set size",
        text.scale = c(1.3, 1.3, 1.15, 1.15, 1.35, 1.45)))
  dev.off()
  cat("saved:", file_base, ".png/.pdf\n")
}

# down-regulated: Surgery M3 / Surgery Y3 / Chemo Y3 (blue)
plot_upset(list(M3 = sm$down, Y3 = sy$down, Cy3 = cy$down),
           c("M3", "Y3", "Cy3"), "#5DADE2",
           file.path(out, "UpSetR_Surgery_M3_Y3_vs_Chemo_down"))
# up-regulated: Surgery M3 / Surgery Y3 / Chemo Y3 (red)
plot_upset(list(M3 = sm$up, Y3 = sy$up, Cy3 = cy$up),
           c("M3", "Y3", "Cy3"), "#E74C3C",
           file.path(out, "UpSetR_Surgery_M3_Y3_vs_Chemo_up"))
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

# ---------- Fig 4A: Sankey of the 28 genes ----------

def fig4a_sankey():
    surg_only_n, trnq_n, chemo_only_n, non_n = 2, 1, 10, 15
    Y0 = 2.0
    surg_seg  = (Y0, Y0+surg_only_n)
    trnq_seg  = (Y0+surg_only_n, Y0+surg_only_n+trnq_n)
    chemo_seg = (Y0+surg_only_n+trnq_n, Y0+surg_only_n+trnq_n+chemo_only_n)
    non_seg   = (Y0+surg_only_n+trnq_n+chemo_only_n, Y0+28)

    non_nd   = (Y0, Y0+non_n)
    chemo_nd = (Y0+non_n+1.5, Y0+non_n+1.5+chemo_only_n+trnq_n)
    surg_nd  = (Y0+non_n+1.5+11+1.5, Y0+non_n+1.5+11+1.5+surg_only_n+trnq_n)
    chemo_part_only = (chemo_nd[0], chemo_nd[1]-trnq_n)
    chemo_part_trnq = (chemo_nd[1]-trnq_n, chemo_nd[1])
    surg_part_only  = (surg_nd[0], surg_nd[1]-trnq_n)
    surg_part_trnq  = (surg_nd[1]-trnq_n, surg_nd[1])

    C_RED, C_SURG, C_CHEMO, C_BOTH, C_NON = "#E74C3C", "#E67E22", "#2E86C1", "#C0392B", "#95A5A6"

    X0, XFLOW0, XR, XR_W = 0.3, 1.5, 8.7, 1.2
    XR_END = XR + XR_W

    def flow(ax, y0l, y1l, y0r, y1r, color):
        verts = [(XFLOW0, y0l), (4.4, y0l), (6.9, y0r), (XR, y0r),
                 (XR, y1r), (6.9, y1r), (4.4, y1l), (XFLOW0, y1l)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        ax.add_patch(mpatches.PathPatch(Path(verts, codes), fc=color, ec="none", alpha=0.5, lw=0))

    fig, ax = plt.subplots(figsize=(8/2.54, 7/2.54))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.22)
    ax.set_xlim(0, XR_END+0.4); ax.set_ylim(0, 42); ax.axis("off")

    ax.add_patch(mpatches.FancyBboxPatch((X0, Y0), 1.0, 28, boxstyle="round,pad=0.15", fc=C_RED, ec="#444"))
    def node(ax, y0, y1, fc):
        ax.add_patch(mpatches.FancyBboxPatch((XR, y0), XR_W, y1-y0, boxstyle="round,pad=0.15", fc=fc, ec="#444"))
    node(ax, *surg_nd, C_SURG)
    node(ax, *chemo_nd, C_CHEMO)
    node(ax, *non_nd, C_NON)

    flow(ax, *surg_seg, *surg_part_only, C_SURG)
    flow(ax, *trnq_seg, *surg_part_trnq, C_BOTH)
    flow(ax, *trnq_seg, *chemo_part_trnq, C_BOTH)
    flow(ax, *chemo_seg, *chemo_part_only, C_CHEMO)
    flow(ax, *non_seg, *non_nd, C_NON)

    handles = [
        mpatches.Patch(fc=C_RED,  label="Chemo-arm Pre-up (28)"),
        mpatches.Patch(fc=C_SURG, label="Surgery Y3 down (3)"),
        mpatches.Patch(fc=C_CHEMO,label="Chemo Y3 down (11)"),
        mpatches.Patch(fc=C_BOTH, label="TRNQ both (1)"),
        mpatches.Patch(fc=C_NON,  label="Non-DEG (15)"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.04), ncol=2,
               frameon=False, fontsize=6.5, handlelength=1.2, handleheight=1.2, columnspacing=0.8)

    plt.savefig(f"{OUT}/Sankey_28genes_Y3_outcome_v9.png", dpi=300)
    plt.savefig(f"{OUT}/Sankey_28genes_Y3_outcome_v9.pdf")
    plt.close()
    print("Fig4a saved: Sankey_28genes_Y3_outcome_v9.png/pdf")

# ---------- Fig 4B: 13-gene heatmap (R) ----------

def fig4b_heatmap():
    run_r(R_13GENES_HEATMAP)

# ---------- Fig 4C/4D: UpSet (R) ----------

def fig4c_upset():
    run_r(R_UPSET)

# ---------- Fig 4E/4G: four-time-point boxplots of MEMO1 / MMP11 ----------

def fig4e_gene_4tp(genes=("MEMO1", "MMP11")):
    meta = pd.read_csv(f"{DATA}/metadata_subgroup_bmi.csv")
    tpm = pd.read_csv(TPM, index_col=0)
    meta = meta[meta["expr_col"].isin(tpm.columns)].reset_index(drop=True)
    meta["recur"] = meta["是否复发（是/否）"].fillna("NA").astype(str)

    def series(g, mask):
        return pd.Series(np.log2(tpm.loc[g, meta.loc[mask, "expr_col"]].astype(float).values + 1),
                         index=meta.loc[mask, "core"]).dropna()

    def point_data(g, mask):
        sub = meta[mask]
        return np.log2(tpm.loc[g, sub["expr_col"]].astype(float).values + 1), sub["recur"].values

    for gene in genes:
        if gene not in tpm.index:
            print(f"!! {gene} not found in TPM matrix, skipped")
            continue
        ctrl = series(gene, meta["group"]=="Control")
        pre = meta[meta["group"]=="Pre"].copy(); pre["figo_roman"] = pre["figo_roman"].astype(str)
        pre_early = series(gene, (meta["group"]=="Pre") & (meta["expr_col"].isin(pre.loc[pre["figo_roman"].isin(["I","II"]), "expr_col"])))
        pre_late  = series(gene, (meta["group"]=="Pre") & (meta["expr_col"].isin(pre.loc[pre["figo_roman"].isin(["III","IV"]), "expr_col"])))
        m3s = series(gene, (meta["group"]=="Post3M") & (meta["tx_arm"]=="Surgery"))
        m3c = series(gene, (meta["group"]=="Post3M") & (meta["tx_arm"]=="Chemo"))
        y3s = series(gene, (meta["group"]=="Post3Y") & (meta["tx_arm"]=="Surgery"))
        y3c = series(gene, (meta["group"]=="Post3Y") & (meta["tx_arm"]=="Chemo"))

        blocks = [
            ("Control", [ctrl], ["#95A5A6"], 1.0),
            ("Pre", [pre_early, pre_late], ["#85C1E9", "#2E86C1"], 2.3),
            ("Post3M", [m3s, m3c], ["#82E0AA", "#229954"], 3.6),
            ("Post3Y", [y3s, y3c], ["#F0B27A", "#CA6F1E"], 4.9),
        ]
        legend_items = [("Control", "#95A5A6"), ("Pre Early (I-II)", "#85C1E9"), ("Pre Late (III-IV)", "#2E86C1"),
                        ("M3 Surgery", "#82E0AA"), ("M3 Chemo", "#229954"),
                        ("Y3 Surgery", "#F0B27A"), ("Y3 Chemo", "#CA6F1E")]

        pts = [
            (1.0, point_data(gene, meta["group"]=="Control")),
            (2.14, point_data(gene, (meta["group"]=="Pre") & (meta["expr_col"].isin(pre.loc[pre["figo_roman"].isin(["I","II"]), "expr_col"])))),
            (2.46, point_data(gene, (meta["group"]=="Pre") & (meta["expr_col"].isin(pre.loc[pre["figo_roman"].isin(["III","IV"]), "expr_col"])))),
            (3.44, point_data(gene, (meta["group"]=="Post3M") & (meta["tx_arm"]=="Surgery"))),
            (3.76, point_data(gene, (meta["group"]=="Post3M") & (meta["tx_arm"]=="Chemo"))),
            (4.74, point_data(gene, (meta["group"]=="Post3Y") & (meta["tx_arm"]=="Surgery"))),
            (5.06, point_data(gene, (meta["group"]=="Post3Y") & (meta["tx_arm"]=="Chemo"))),
        ]

        fig, ax = plt.subplots(figsize=(12/2.54, 6.96/2.54))
        fig.subplots_adjust(left=0.13, right=0.68, top=0.86, bottom=0.16)
        all_pos, all_vals, all_cols, x_centers = [], [], [], []
        for tp_name, gs, cs, x0 in blocks:
            n = len(gs)
            pos = [x0 + 0.32*k - 0.16*(n-1) for k in range(n)]
            for p, g, c in zip(pos, gs, cs):
                all_pos.append(p); all_vals.append(g.values); all_cols.append(c)
            x_centers.append(x0)
            for g, nm in zip(gs, [tp_name]*n):
                print(f"{gene} {nm} (n={len(g)}): median {g.median():.3f}")

        bp = ax.boxplot(all_vals, positions=all_pos, widths=0.26, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black", lw=1.2))
        for patch, c in zip(bp["boxes"], all_cols):
            patch.set_facecolor(c); patch.set_alpha(0.85); patch.set_edgecolor("#444")
        pts_colors = ["#95A5A6", "#85C1E9", "#2E86C1", "#82E0AA", "#229954", "#F0B27A", "#CA6F1E"]
        n_rec = 0
        for p, (vals, recs), c in zip([x[0] for x in pts], [x[1] for x in pts], pts_colors):
            jit = np.random.default_rng(11).normal(0, 0.012, len(vals))
            for v, rc, j in zip(vals, recs, jit):
                if rc == "是":
                    ax.scatter(p + j, v, marker="*", s=35, color="#E74C3C", edgecolor="none", zorder=5)
                    n_rec += 1
                else:
                    ax.scatter(p + j, v, s=10, color=c, edgecolor="white", lw=0.3, alpha=0.85, zorder=3)
        print(f"{gene} recurrent markers: {n_rec}")

        ax.set_xticks(x_centers)
        ax.set_xticklabels(["Control", "Pre", "Post3M", "Post3Y"], fontsize=8.5)
        ax.set_ylabel("log2(TPM+1)", fontsize=9.5)
        ax.set_title(f"{gene} expression across time points", fontsize=10)
        ax.spines[["top","right"]].set_visible(False)
        ymin = min(v.min() for v in all_vals); ymax = max(v.max() for v in all_vals)
        ax.set_ylim(ymin - 0.02*(ymax-ymin), ymax + 0.04*(ymax-ymin+1e-9))
        handles = [Patch(facecolor=c, edgecolor="#444", alpha=0.85, label=l) for l, c in legend_items]
        handles.append(plt.Line2D([0],[0], marker="*", color="w", markerfacecolor="#E74C3C", markersize=7, label="recurrent patient"))
        ax.legend(handles=handles, bbox_to_anchor=(1.02, 0.5), loc="center left", fontsize=6.5, frameon=False,
                  handlelength=1.0, handleheight=1.0, borderaxespad=0)
        plt.savefig(f"{OUT}/{gene}_4timepoints_boxplot.png", dpi=300)
        plt.savefig(f"{OUT}/{gene}_4timepoints_boxplot.pdf")
        plt.close()
        print(f"saved: {gene}_4timepoints_boxplot.png/pdf")

# ---------- main ----------

if __name__ == "__main__":
    # 1: Sankey
    fig4a_sankey()
    # 2: 13-gene heatmap
    fig4b_heatmap()
    # 3: UpSet
    fig4c_upset()
    # 4: MEMO1 / MMP11 four-time-point boxplots
    fig4e_gene_4tp(["MEMO1", "MMP11"])
    print("All Figure 4 panels saved.")
