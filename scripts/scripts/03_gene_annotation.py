#!/usr/bin/env python3
"""
03_gene_annotation.py
Annotates a list of significant genes with information from MyGene.info
"""

import pandas as pd
from mygene import MyGeneInfo
import os

# --------------------------------------------------
# Create folders
# --------------------------------------------------
if not os.path.exists("results"):
    os.makedirs("results")

# --------------------------------------------------
# Load DESeq2 results
# --------------------------------------------------
results_file = "results/deseq2_results.csv"
if not os.path.isfile(results_file):
    raise FileNotFoundError(f"Cannot find {results_file}. Run 02_DE_analysis.R first.")

df = pd.read_csv(results_file)

# --------------------------------------------------
# Filter significant genes (p-value < 0.05)
# --------------------------------------------------
sig_genes = df[df["pvalue"] < 0.05]["Gene"].tolist()

if not sig_genes:
    print("No significant genes found. Exiting.")
    exit()

# --------------------------------------------------
# Annotate using MyGene.info
# --------------------------------------------------
mg = MyGeneInfo()
annotations = mg.querymany(sig_genes, scopes="symbol", fields="name,entrezgene,summary", species="human")

# Convert to DataFrame
annot_df = pd.DataFrame(annotations)
