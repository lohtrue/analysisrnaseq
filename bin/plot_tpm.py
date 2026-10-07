#!/usr/bin/env python3
"""
Create some overview plots from the TPM matrix (genes x samples).  

Usage: plot_tpm.py tpm_matrix.tsv 

Output 
  - tpm_top_genes_mqc.png       top 20 genes by mean TPM
  - tpm_distribution_mqc.png    distribution of log10(TPM + 1) per sample
  - tpm_correlation_mqc.png     Spearman correlation between samples (only if >= 2 samples) 
"""
import os
import sys

os.environ["MPLCONFIGDIR"] = os.getcwd()  # matplotlib needs a writable cache directory inside the container

import matplotlib
matplotlib.use("Agg")  # no screen in the container: render directly to files
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TOP_N = 20

tpm = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
names = tpm.pop("gene_name")          # keep readable gene names separately
samples = list(tpm.columns)
log_tpm = np.log10(tpm + 1)           # log scale: TPM values span several orders of magnitude

# 1) Top expressed genes
top = tpm.loc[tpm.mean(axis=1).sort_values(ascending=False).index[:TOP_N]]
ax = top.set_index(names[top.index]).plot.barh(figsize=(8, 7))
ax.invert_yaxis()                     # highest expressed gene on top
ax.set_xlabel("TPM")
ax.set_ylabel("")
ax.set_title(f"Top {TOP_N} genes by mean TPM")
plt.tight_layout()
plt.savefig("tpm_top_genes_mqc.png", dpi=150)
plt.close()

# 2) Distribution of expression values per sample (only genes with TPM > 0)
fig, ax = plt.subplots(figsize=(1.5 + 1.5 * len(samples), 5))
ax.boxplot([log_tpm[s][tpm[s] > 0] for s in samples], tick_labels=samples)
ax.set_ylabel("log10(TPM + 1)")
ax.set_title("Expression distribution of detected genes")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig("tpm_distribution_mqc.png", dpi=150)
plt.close()

# 3) Sample-to-sample correlation (needs at least two samples)
if len(samples) >= 2:
    corr = log_tpm.corr(method="spearman")
    fig, ax = plt.subplots(figsize=(1.5 + 1.2 * len(samples), 1 + 1.2 * len(samples)))
    im = ax.imshow(corr, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(samples)), samples, rotation=30, ha="right")
    ax.set_yticks(range(len(samples)), samples)
    for i in range(len(samples)):
        for j in range(len(samples)):
            ax.text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center", color="white")
    fig.colorbar(im, ax=ax, label="Spearman correlation")
    ax.set_title("Sample correlation (log10 TPM)")
    plt.tight_layout()
    plt.savefig("tpm_correlation_mqc.png", dpi=150)
    plt.close()