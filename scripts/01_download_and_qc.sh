#!/bin/bash
# 01_download_and_qc.sh
# Downloads RNA-seq data from SRA, runs FastQC, and summarizes with MultiQC

# Exit if any command fails
set -e

# Create folders
mkdir -p data/raw data/qc

# Example SRA accession numbers (replace with your own if needed)
SRA_IDS=("SRR12345678" "SRR12345679")

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
