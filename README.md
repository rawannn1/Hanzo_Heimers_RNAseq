# Hanzo Heimers RNA-seq Pipeline

A simple bioinformatics workflow for RNA-seq analysis.

## Overview
This project shows how to:
- Download and QC RNA-seq FASTQ files
- Perform differential expression analysis using DESeq2
- Annotate significant genes with biological information

## Folder Structure
```
Hanzo_Heimers_RNAseq/
│── README.md
│── scripts/
│   ├── 01_download_and_qc.sh
│   ├── 02_DE_analysis.R
│   └── 03_gene_annotation.py
│── environment.yml
│── data/            <-- raw & processed data (not uploaded here)
│── results/         <-- analysis outputs
```

## How to Run

### 1. Install dependencies
Make sure you have [Conda](https://docs.conda.io/en/latest/miniconda.html) installed.

```bash
conda env create -f environment.yml
conda activate hanzo_env
```

### 2. Download data & run QC
```bash
bash scripts/01_download_and_qc.sh
```

### 3. Run differential expression analysis
```bash
Rscript scripts/02_DE_analysis.R
```

### 4. Annotate significant genes
```bash
python scripts/03_gene_annotation.py
```

## Requirements
- Bash  
- R (with DESeq2)  
- Python 3 (with pandas, mygene)  
- FastQC, MultiQC (via bioconda)  
