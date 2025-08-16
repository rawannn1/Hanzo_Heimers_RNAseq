#!/usr/bin/env Rscript
# 02_DE_analysis.R
# Performs differential expression analysis using DESeq2 

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(tibble)
})

# Create folders
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# ----------------------------
# Load count data
# ----------------------------
# Adjust path to your featureCounts output
counts <- read.csv("data/raw/featureCounts_counts.csv", row.names = 1, check.names = FALSE)

# Automatically detect numeric columns (counts)
count_cols <- sapply(counts, is.numeric)
counts <- counts[, count_cols]

# ----------------------------
# Load sample metadata
# ----------------------------
colData <- data.frame(
  condition = factor(c("Control", "Control", "Treatment", "Treatment"))
)
rownames(colData) <- colnames(counts)

# ----------------------------
# Run DESeq2
# ----------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# ----------------------------
# Save results table
# ----------------------------
res_tbl <- as_tibble(res, rownames = "Gene")
write.csv(res_tbl, file = "results/deseq2_results.csv", row.names = FALSE)

# ----------------------------
# Volcano plot
# ----------------------------
res_tbl$significant <- ifelse(res_tbl$pvalue < 0.05, "yes", "no")
p <- ggplot(res_tbl, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Differential Expression (Volcano Plot)",
       x = "log2 Fold Change",
       y = "-log10 p-value")

ggsave("results/volcano_plot.png", plot = p, width = 6, height = 4)

cat("Differential expression analysis complete. Results saved in 'results/' folder.\n")
