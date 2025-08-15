#!/usr/bin/env Rscript
# 02_DE_analysis.R
# Performs a simple differential expression analysis using DESeq2

# Load libraries
suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(tibble)
})

# Create folders
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# --------------------------------------------------
# Example input: pretend we already have a count matrix
# Replace this with real data if running in practice
# --------------------------------------------------
set.seed(42)
counts <- matrix(
  rpois(1000, lambda = 20), 
  ncol = 4
)
colnames(counts) <- c("Ctrl_1", "Ctrl_2", "Treat_1", "Treat_2")
rownames(counts) <- paste0("Gene_", 1:nrow(counts))

# Metadata (experimental design)
colData <- data.frame(
  condition = factor(c("Control", "Control", "Treatment", "Treatment"))
)
rownames(colData) <- colnames(counts)

# --------------------------------------------------
# Run DESeq2
# --------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# --------------------------------------------------
# Save results table
# --------------------------------------------------
res_tbl <- as_tibble(res, rownames = "Gene")
write.csv(res_tbl, file = "results/deseq2_results.csv", row.names = FALSE)

# --------------------------------------------------
# Volcano plot
# --------------------------------------------------
res_tbl$significant <- ifelse(res_tbl$pvalue < 0.05, "yes", "no")
p <- ggplot(res_tbl, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Differential Expression (Volcano Plot)",
       x = "log2 Fold Change",
       y = "-log10 p-value")

ggsave("results/volcano_plot.png", plot = p, width = 6, height = 4)

cat("Differential expression analysis complete. Results saved in 'results/' folder.\n")
