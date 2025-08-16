#!/bin/bash
# 01_download_and_qc.sh
# Downloads RNA-seq data from SRA, runs FastQC, and summarizes with MultiQC

# Exit if any command fails
set -e

# Create folders
mkdir -p data/raw data/qc

# SRA accession numbers 
SRA_IDS=("SRR5944937" "SRR5944938" "SRR5944939" "SRR5944940")



# Download FASTQ files
for SRA_ID in "${SRA_IDS[@]}"; do
    echo "Downloading $SRA_ID..."
    fasterq-dump $SRA_ID -O data/raw
done

# Run FastQC
echo "Running FastQC..."
fastqc data/raw/*.fastq -o data/qc

# Summarize with MultiQC
echo "Running MultiQC..."
multiqc data/qc -o data/qc

echo "QC pipeline complete!"
