#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combined Figure 2 (panels a-e) generation script for the surgery arm.

Fig 2a: clinical information heatmap (transposed, 8 variables x samples)
Fig 2b: DEG count barchart (up/down) across three contrasts
Fig 2c: GO BP combined dotplot across three periods (cov-adjusted)
Fig 2d: three-way Venn (Pre up / Post3M down / Post3Y down)
Fig 2e: mean z-score trend of the overlapping genes

Inputs : cfRNA_project/data tables, DEG and enrichment tables, Paper_TPM.csv, Samples_info.xlsx
Outputs: cfRNA_project/output/ (PNG/PDF/SVG)

Usage: python3 fig2_combined.py
"""
import pandas as pd, numpy as np, os, re
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
from matplotlib_venn import venn3

DATA = "/work/zhouqing1/cfRNA_project/data"
OUT  = "/work/zhouqing1/cfRNA_project/output"
TPM  = "/work/zhouqing1/uploads/Paper_TPM.csv"
FC_TH = np.log2(1.2)

# ---------- shared helpers ----------

def sig_genes(arm, nm, direction):
    d = pd.read_csv(os.path.join(OUT, f"DEG_{arm}_{nm}_cov.csv"))
    m = (d["FDR"] < 0.1) & (d["logFC"] > FC_TH if direction == "up" else d["logFC"] < -FC_TH)
    return set(d.loc[m, "gene"])

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

# ---------- Fig 2a: clinical information heatmap (transposed) ----------

def fig2a_clinical_heatmap():
    raw = pd.read_excel("/work/zhouqing1/uploads/Samples_info.xlsx", sheet_name=0)
    meta = pd.read_csv(f"{DATA}/metadata_subgroup_bmi.csv")
    raw_eval = raw[["sample", "采血时评价"]].rename(columns={"采血时评价": "eval_raw"})
    meta = meta.merge(raw_eval, on="sample", how="left")
    if "采血时评价" in meta.columns:
        meta["采血时评价"] = meta["eval_raw"].fillna(meta["采血时评价"])
    else:
        meta = meta.rename(columns={"eval_raw": "采血时评价"})
    meta = meta.drop(columns=["eval_raw"])

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
        print(f"Fig2a saved: {out_png} (20cm x 10cm transposed, 300dpi, n={n})")

    draw("Surgery", f"{OUT}/Surgery_clinical_heatmap_transpose_20x10cm.png", f"{OUT}/Surgery_clinical_heatmap_transpose_20x10cm.svg", f"{OUT}/Surgery_clinical_heatmap_transpose_20x10cm.pdf")
    draw("Chemo", f"{OUT}/Chemo_clinical_heatmap_transpose_20x10cm.png", f"{OUT}/Chemo_clinical_heatmap_transpose_20x10cm.svg", f"{OUT}/Chemo_clinical_heatmap_transpose_20x10cm.pdf")

# ---------- Fig 2b: DEG count barchart ----------

def fig2b_deg_bar():
    arm = "Surgery"
    contrasts = ["Pre_vs_Control", "Post3M_vs_Pre", "Post3Y_vs_Pre"]
    labels = ["Pre vs Control", "Post3M vs Pre", "Post3Y vs Pre"]
    ups, dns = [], []
    for nm in contrasts:
        d = pd.read_csv(f"{OUT}/DEG_{arm}_{nm}_cov.csv")
        sig = d["FDR"] < 0.1
        ups.append(int((sig & (d["logFC"] > FC_TH)).sum()))
        dns.append(int((sig & (d["logFC"] < -FC_TH)).sum()))
    print(f"Surgery: up={ups}, down={dns}")

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
    ax.set_ylim(ymin, ymax)
    ax.set_title("")
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    plt.tight_layout()
    plt.savefig(f"{OUT}/Surgery_DEG_bar_cov_v3.png", dpi=300)
    plt.savefig(f"{OUT}/Surgery_DEG_bar_cov_v3.pdf")
    plt.close()
    print("Fig2b saved: Surgery_DEG_bar_cov_v3.png/pdf")

# ---------- Fig 2c: GO BP combined dotplot ----------

def parse_gene_ratio(s):
    try:
        a, b = str(s).split("/")
        return float(a) / float(b)
    except Exception:
        return np.nan

def fig2c_go_dotplot():
    CONTRAST_LABEL = {
        "Pre_vs_Control": "Pre vs Control",
        "Post3M_vs_Pre": "Post3M vs Pre",
        "Post3Y_vs_Pre": "Post3Y vs Pre",
    }

    def load_panel(arm, contrasts, topn=8):
        entries = []
        for nm in contrasts:
            f = os.path.join(OUT, f"enrich_{arm}_{nm}_cov_GOBP.csv")
            if not os.path.exists(f):
                continue
            d = pd.read_csv(f)
            if d.empty:
                continue
            d = d.sort_values("pvalue").head(topn).copy()
            d["ratio"] = d["GeneRatio"].apply(parse_gene_ratio)
            d["neglogp"] = -np.log10(d["pvalue"])
            d = d.sort_values("ratio")
            entries.append((True, CONTRAST_LABEL[nm], np.nan, 0, 0))
            for _, r in d.iterrows():
                entries.append((False, str(r["Description"])[:58] + ("..." if len(str(r["Description"])) > 58 else ""),
                                r["ratio"], r["Count"], r["neglogp"]))
        return entries

    entries = load_panel("Surgery", ["Pre_vs_Control", "Post3M_vs_Pre", "Post3Y_vs_Pre"], topn=8)
    print("Surgery GO entries:", len(entries))
    rows = [(i, e) for i, e in enumerate(entries)]
    y = np.array([i for i, e in rows])
    ratio = np.array([e[2] if not e[0] else np.nan for i, e in rows])
    count = np.array([e[3] for i, e in rows])
    neglogp = np.array([e[4] if not e[0] else 0 for i, e in rows])
    mask = ~np.isnan(ratio)

    fig, ax = plt.subplots(figsize=(12, 9))
    sc = ax.scatter(ratio[mask], y[mask], s=count[mask] * 5 + 15, c=neglogp[mask],
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
    ax.set_title("Surgery arm (cov-adjusted): GO BP across three periods", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{OUT}/Surgery_GO_combined_cov.png", dpi=300)
    plt.close()
    print("Fig2c saved: Surgery_GO_combined_cov.png")

# ---------- Fig 2d: three-way Venn ----------

def fig2d_venn():
    pre_up = sig_genes("Surgery", "Pre_vs_Control", "up")
    m_down = sig_genes("Surgery", "Post3M_vs_Pre", "down")
    y_down = sig_genes("Surgery", "Post3Y_vs_Pre", "down")
    print(f"cov: pre_up={len(pre_up)}, m_down={len(m_down)}, y_down={len(y_down)}")

    fig, ax = plt.subplots(figsize=(7/2.54, 7/2.54))
    v = venn3([pre_up, m_down, y_down],
              set_labels=(f"Pre up\n(n={len(pre_up)})", f"Post3M down\n(n={len(m_down)})", f"Post3Y down\n(n={len(y_down)})"))
    for t in v.set_labels:
        t.set_fontsize(10)
    for t in v.subset_labels:
        if t is not None:
            t.set_fontsize(10)
    ax.set_title("")
    plt.tight_layout()
    plt.savefig(f"{OUT}/Surgery_Venn_cov_v9.png", dpi=300)
    plt.savefig(f"{OUT}/Surgery_Venn_cov_v9.pdf")
    plt.savefig(f"{OUT}/Surgery_Venn_cov_v9.svg", dpi=300)
    plt.close()
    print("Fig2d saved: Surgery_Venn_cov_v9.png/pdf/svg")

    om = pre_up & m_down; oy = pre_up & y_down; omy = m_down & y_down; o3 = pre_up & m_down & y_down
    reg = {"100": len(pre_up)-len(om)-len(oy)+len(o3), "010": len(m_down)-len(om)-len(omy)+len(o3),
           "001": len(y_down)-len(oy)-len(omy)+len(o3), "110": len(om)-len(o3), "101": len(oy)-len(o3),
           "011": len(omy)-len(o3), "111": len(o3)}
    print("regions:", reg, "| union:", sum(reg.values()))
    print("Pre up & M3 down:", len(om), "Pre up & Y3 down:", len(oy), "M3 & Y3:", len(omy))

# ---------- Fig 2e: overlap gene trend ----------

def fig2e_trend():
    arm = "Surgery"
    pre_up = sig_genes(arm, "Pre_vs_Control", "up")
    m_down = sig_genes(arm, "Post3M_vs_Pre", "down")
    y_down = sig_genes(arm, "Post3Y_vs_Pre", "down")
    o_all = pre_up & (m_down | y_down)
    print(f"overlap genes n={len(o_all)}")

    tpm = pd.read_csv(TPM, index_col=0)
    meta = pd.read_csv(f"{DATA}/metadata_subgroup_bmi.csv")
    meta = meta[meta["expr_col"].isin(tpm.columns) & (meta["tx_arm"] == arm)]
    genes = sorted(o_all & set(tpm.index))
    print(f"overlap genes available in TPM: {len(genes)}")
    expr = np.log2(tpm.loc[genes] + 1)

    means = {}
    for gname in ["Pre", "Post3M", "Post3Y"]:
        cols = meta.loc[meta["group"] == gname, "expr_col"]
        means[gname] = expr[cols].mean(axis=1)
    prof = pd.DataFrame(means)
    zs = prof.apply(lambda r: (r - r.mean()) / (r.std() + 1e-9), axis=1)

    fig, ax = plt.subplots(figsize=(8/2.54, 8/2.54))
    t = np.arange(3)
    tick = ["Pre", "Post3M", "Post3Y"]
    for _, row in zs.iterrows():
        ax.plot(t, row.values, color="#CCCCCC", lw=0.6, alpha=0.5, zorder=1)
    m = zs.mean(axis=0).values
    se = zs.std(axis=0).values / np.sqrt(len(zs))
    ax.plot(t, m, color="#C0392B", lw=2.8, zorder=3, label=f"mean (n={len(zs)} genes)")
    ax.fill_between(t, m - se, m + se, color="#C0392B", alpha=0.2, zorder=2)
    ax.axhline(0, color="grey", lw=0.8, ls="--")
    ax.set_xticks(t); ax.set_xticklabels(tick, fontsize=10)
    ax.set_ylabel("Mean z-score of log2(TPM+1)", fontsize=10)
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    ax.set_ylim(-2.4, 2.4)
    plt.tight_layout()
    plt.savefig(f"{OUT}/Surgery_overlap_trend_cov_v3.png", dpi=300)
    plt.savefig(f"{OUT}/Surgery_overlap_trend_cov_v3.pdf")
    plt.close()
    print("Fig2e saved: Surgery_overlap_trend_cov_v3.png/pdf")

# ---------- main ----------

if __name__ == "__main__":
    fig2a_clinical_heatmap()
    fig2b_deg_bar()
    fig2c_go_dotplot()
    fig2d_venn()
    fig2e_trend()
    print("All Figure 2 panels (a-e) saved.")
