# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)

# Load data
mlh1_res <- readRDS("MLH1/results/4_deseq2/MLH1_res.Rds")
msh2_res <- readRDS("MSH2/results/4_deseq2/MSH2_res.Rds")
combined_results <- read.csv("combined_DE_results.csv")

# Function to create volcano plot
create_volcano_plot <- function(res, title) {
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = title,
                  pCutoff = 0.05,
                  FCcutoff = 1.0,
                  pointSize = 2.0,
                  labSize = 3.0)
}

# Function to create DE plot
create_de_plot <- function(res, title) {
  res$significant <- ifelse(res$padj < 0.05, "Significant", "Not Significant")
  ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
    geom_point() +
    theme_minimal() +
    labs(title=title, x="Log2 Fold Change", y="-Log10 Adjusted P-value")
}

# Generate plots for MLH1
mlh1_volcano <- create_volcano_plot(mlh1_res, "Volcano Plot - MLH1")
mlh1_de <- create_de_plot(mlh1_res, "DE Plot - MLH1")

# Generate plots for MSH2
msh2_volcano <- create_volcano_plot(msh2_res, "Volcano Plot - MSH2")
msh2_de <- create_de_plot(msh2_res, "DE Plot - MSH2")

# Display plots
print(mlh1_volcano)
print(mlh1_de)
print(msh2_volcano)
print(msh2_de)

#HEATMAP
# Get top 20 DE genes for each analysis
top_genes_mlh1 <- head(rownames(mlh1_res[order(mlh1_res$padj), ]), 20)
top_genes_msh2 <- head(rownames(msh2_res[order(msh2_res$padj), ]), 20)

# Filter for genes DE in only one condition
only_mlh1 <- combined_results[combined_results$Result_KOvsWT == "DE (down) in MLH1" | combined_results$Result_KOvsWT == "DE (up) in MLH1", ]
only_msh2 <- combined_results[combined_results$Result_KOvsWT == "DE (down) in MSH2" | combined_results$Result_KOvsWT == "DE (up) in MSH2", ]

# Ensure we only have top 20 DE genes unique to each condition
top_genes_mlh1 <- intersect(top_genes_mlh1, only_mlh1$gene)
top_genes_msh2 <- intersect(top_genes_msh2, only_msh2$gene)

# Function to create heatmap
create_heatmap <- function(dds, genes, title) {
  data <- assay(dds)[genes, ]
  pheatmap(data, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames=TRUE, main=title)
}

# Create heatmaps
mlh1_dds <- readRDS("MLH1/results/4_deseq2/MLH1_dds.Rds")
msh2_dds <- readRDS("MSH2/results/4_deseq2/MSH2_dds.Rds")

create_heatmap(mlh1_dds, top_genes_mlh1, "Top 20 DE Genes - MLH1")
create_heatmap(msh2_dds, top_genes_msh2, "Top 20 DE Genes - MSH2")

#VIOLIN PLOT
# Function to create violin plots
create_violin_plot <- function(dds, genes, condition, title) {
  data <- assay(dds)[genes, ]
  data <- as.data.frame(t(data))
  data$condition <- condition
  
  melted_data <- melt(data, id.vars = "condition")
  
  ggplot(melted_data, aes(x=condition, y=value, fill=condition)) +
    geom_violin(trim=FALSE) +
    facet_wrap(~ variable, scales="free_y") +
    theme_minimal() +
    labs(title=title, x="Condition", y="Counts")
}

# Extract condition labels
mlh1_condition <- colData(mlh1_dds)$condition
msh2_condition <- colData(msh2_dds)$condition

# Create violin plots for MLH1
create_violin_plot(mlh1_dds, top_genes_mlh1, mlh1_condition, "Violin Plots - MLH1")

# Create violin plots for MSH2
create_violin_plot(msh2_dds, top_genes_msh2, msh2_condition, "Violin Plots - MSH2")

