#!/bin/bash
# 01_download_and_qc.sh
# Downloads RNA-seq data from SRA, runs FastQC, and summarizes with MultiQC

# Exit if any command fails
set -e

# Create folders
mkdir -p data/raw data/qc

# Example SRA accession numbers 
SRA_IDS=("SRR1039508" "SRR1039509" "SRR1039510" "SRR1039511")


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
