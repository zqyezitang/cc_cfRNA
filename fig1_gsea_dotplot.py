#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import pandas as pd
import numpy as np
import textwrap
from matplotlib import font_manager
import matplotlib
matplotlib.use("Agg")
font_manager.fontManager.addfont("/public/skills/skill_tools/font/fonts/NotoSansCJKsc-Regular.otf")
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK SC"]
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams["font.size"] = 8
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

######Figure 1D
TOP = 20

gs = pd.read_csv(f"{OUT}/GSEA_GOBP_Pre_vs_Control.csv")
gs = gs[gs["p.adjust"] < 0.05].copy()
print("FDR<0.05 :", len(gs), " NES<0 :", (gs.NES < 0).sum())
neg = gs[gs.NES < 0]
print("NES<0 :")
print(neg[["Description", "setSize", "NES", "p.adjust"]].round(3).to_string(index=False))

# 
s = gs.reindex(gs.NES.abs().sort_values(ascending=False).index).head(TOP)
s = s.reindex(s.NES.abs().sort_values(ascending=False).index).iloc[::-1]  # 
print("\n:", len(s), "  NES<0:", (s.NES < 0).sum())
print(s[["Description", "NES", "setSize", "p.adjust"]].round(3).to_string(index=False))

labels_txt = ["\n".join(textwrap.wrap(d, 100)) for d in s.Description]  # 
sizes = 40 + 140 * (s["setSize"] - s["setSize"].min()) / (s["setSize"].max() - s["setSize"].min() + 1e-9)
vmin, vmax = 1, 30
colors = np.clip(-np.log10(s["p.adjust"].clip(lower=1e-300)), vmin, vmax)

#
fig = plt.figure(figsize=(14/2.54, 7/2.54))
ax = fig.add_axes([0.60, 0.12, 0.32, 0.78])
ys = np.arange(len(s)) * 1.3  
sc = ax.scatter(s["NES"], ys, s=sizes, c=colors, cmap="RdYlBu_r",
                norm=Normalize(vmin, vmax), edgecolor="gray", lw=0.3, alpha=0.95)
ax.set_yticks(ys)
ax.set_yticklabels(labels_txt, fontsize=6.5)
ax.tick_params(axis="y", pad=2)
ax.axvline(0, color="black", lw=0.7)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(labelsize=7)
ax.set_ylim(-0.8, len(s) * 1.3 - 0.4)
# 
ax.set_xlim(s["NES"].min() - 0.15, s["NES"].max() + 0.15)
ax.set_xlabel("NES", fontsize=8)
ax.set_title("GSEA GO BP: Pre vs Control\n(|NES| top20, FDR<0.05)", fontsize=9)
cbar = fig.colorbar(cm.ScalarMappable(norm=Normalize(vmin, vmax), cmap="RdYlBu_r"),
                    ax=ax, fraction=0.04, pad=0.03)
cbar.set_label("-log10(adjusted P)", fontsize=6)
cbar.ax.tick_params(labelsize=6)

plt.savefig(f"{OUT}/GSEA_Pre_vs_Control_dotplot.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{OUT}/GSEA_Pre_vs_Control_dotplot.pdf", bbox_inches="tight")
plt.close()
print(" GSEA_Pre_vs_Control_dotplot.png/pdf")

######Figure1E-F
OUT = "/work/zhouqing1/cfRNA_project/output"
prefixes = [("AllAdj")]
contrasts = [("Post3M_vs_Pre", "T1 vs T0"), ("Post3Y_vs_Pre", "T2 vs T0")]
TOP = 12

fig, axes = plt.subplots(1, 2, figsize=(20/2.54, 9/2.54))
for ai, (pfx, pname) in enumerate(prefixes):
    for ci, (c, cname) in enumerate(contrasts):
        ax = axes[ci]
        gs = pd.read_csv(f"{OUT}/GSEA_longitudinal_{pfx}_{c}_GOBP.csv")
        gs["neg_log10q"] = -np.log10(gs["qvalue"].clip(lower=1e-300))
        s = gs.reindex(gs.NES.abs().sort_values(ascending=False).index).head(TOP).iloc[::-1]
        labels_txt = ["\n".join(textwrap.wrap(d, 34)) for d in s.Description]
        sizes = 60 + 260 * (s["setSize"] - s["setSize"].min()) / (s["setSize"].max() - s["setSize"].min() + 1e-9)
        vmin, vmax = 1, 30  # -log10(qvalue) 
        colors = np.clip(s["neg_log10q"], vmin, vmax)
        sc = ax.scatter(s["NES"], range(len(s)), s=sizes, c=colors, cmap="RdYlBu_r",
                        norm=Normalize(vmin, vmax), edgecolor="gray", lw=0.3, alpha=0.95)
        ax.set_yticks(range(len(s))); ax.set_yticklabels(labels_txt, fontsize=6.5)
        ax.axvline(0, color="black", lw=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=7)
        ax.set_xlim(s["NES"].min() - 0.5, s["NES"].max() + 0.5)
        if ci == 0:
            ax.set_xlabel("NES", fontsize=8)
        # 
        cbar = fig.colorbar(cm.ScalarMappable(norm=Normalize(vmin, vmax), cmap="RdYlBu_r"),
                            ax=ax, fraction=0.03, pad=0.03)
        cbar.set_label("-log10(qvalue)", fontsize=6)
        cbar.ax.tick_params(labelsize=6)

plt.tight_layout(h_pad=1.2, w_pad=0.8)
plt.savefig(f"{OUT}/GSEA_All_T1T2_dotplot.png", dpi=300)
plt.savefig(f"{OUT}/GSEA_All_T1T2_dotplot.pdf")
plt.close()

