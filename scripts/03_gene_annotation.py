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
os.makedirs("results", exist_ok=True)

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
    print("No significant genes found. Creating empty annotation file.")
    pd.DataFrame(columns=["query", "name", "entrezgene", "summary"]).to_csv(
        "results/gene_annotations.csv", index=False
    )
    exit()

# --------------------------------------------------
# Annotate using MyGene.info
# --------------------------------------------------
mg = MyGeneInfo()

# Query MyGene.info safely
annotations = mg.querymany(
    sig_genes,
    scopes="symbol",
    fields="name,entrezgene,summary",
    species="human",
    as_dataframe=True
)

# Fill missing columns with Na if necessary
for col in ["name", "entrezgene", "summary"]:
    if col not in annotations.columns:
        annotations[col] = pd.NA

# Save annotations
annotations.reset_index(inplace=True)
annotations.rename(columns={"query": "Gene"}, inplace=True)
annotations.to_csv("results/gene_annotations.csv", index=False)

print(f"Gene annotation complete. Results saved to 'results/gene_annotations.csv'.")
