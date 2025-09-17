# Load the libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(msigdbr)
library(BiocParallel)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(grid)
#load data
MSH2_res <- readRDS("../results/4_deseq2/MSH2_res.Rds")
MSH2_dds <- readRDS("../results/4_deseq2/MSH2_dds.Rds")

# Convert DESeqMSH2_results to a data frame
MSH2_res_df <- as.data.frame(MSH2_res)
MSH2_res_df$external_gene_name <- MSH2_res_df$geneSymbol

#read in gene mapping df
#gene_mapping <- readRDS("../results/4_deseq2/MSH2_gene_mapping.Rds")

# Merge with the gene mapping to retain ENSEMBL IDs and gene names
#MSH2_res_df <- merge(MSH2_res_df, gene_mapping, by = "external_gene_name")

# Create a named vector of log2FoldChange values with ENSEMBL IDs
geneList <- MSH2_res_df$log2FoldChange
names(geneList) <- MSH2_res_df$external_gene_name

# Sort the gene list in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# Download MSigDB gene sets
msigdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
msigdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
msigdb_h <- msigdbr(species = "Homo sapiens", category = "H")

# Combine C2 and H sets for further use
msigdb_combined <- rbind(msigdb_c2, msigdb_h)

# Create gene sets for immune-related pathways (if using C5, you'll need to filter specifically)
immune_gene_sets <- msigdb_c5[grepl("immune", msigdb_c5$gs_name, ignore.case = TRUE),]

# Convert to list format required by GSEA function
immune_gene_sets_list <- split(immune_gene_sets$gene_symbol, immune_gene_sets$gs_name)
c2_gene_sets_list <- split(msigdb_c2$gene_symbol, msigdb_c2$gs_name)
h_gene_sets_list <- split(msigdb_h$gene_symbol, msigdb_h$gs_name)
combined_gene_sets_list <- split(msigdb_combined$gene_symbol, msigdb_combined$gs_name)

# Ensure your gene sets are in the correct format: a data frame with two columns
immune_gene_sets_df <- do.call(rbind, lapply(names(immune_gene_sets_list), function(gs_name) {
  data.frame(term = gs_name, gene = immune_gene_sets_list[[gs_name]])
}))

c2_gene_sets_df <- do.call(rbind, lapply(names(c2_gene_sets_list), function(gs_name) {
  data.frame(term = gs_name, gene = c2_gene_sets_list[[gs_name]])
}))

h_gene_sets_df <- do.call(rbind, lapply(names(h_gene_sets_list), function(gs_name) {
  data.frame(term = gs_name, gene = h_gene_sets_list[[gs_name]])
}))

combined_gene_sets_df <- do.call(rbind, lapply(names(combined_gene_sets_list), function(gs_name) {
  data.frame(term = gs_name, gene = combined_gene_sets_list[[gs_name]])
}))

# Register a parallel backend
register(MulticoreParam(4))  # Use 4 cores; adjust based on your system

# Assuming 'geneList' is already prepared as shown previously
gsea_immune <- GSEA(geneList, TERM2GENE = immune_gene_sets_df, pvalueCutoff = 0.05, verbose = TRUE, BPPARAM = MulticoreParam(4))
gsea_c2 <- GSEA(geneList, TERM2GENE = c2_gene_sets_df, pvalueCutoff = 0.05, verbose = TRUE, BPPARAM = MulticoreParam(4))
gsea_h <- GSEA(geneList, TERM2GENE = h_gene_sets_df, pvalueCutoff = 0.05, verbose = TRUE, BPPARAM = MulticoreParam(4))
gsea_combined <- GSEA(geneList, TERM2GENE = combined_gene_sets_df, pvalueCutoff = 0.05, verbose = TRUE, BPPARAM = MulticoreParam(4))

# Print and save results
# Function to create GSEA plots for each pathway
create_gsea_plots <- function(gsea_result) {
  plots <- list()
  n_pathways <- length(gsea_result@result$ID)  # Determine the number of pathways
  for (i in seq_len(n_pathways)) {
    pathway_name <- gsea_result@result$Description[i]
    # Create one plot per pathway
    p <- gseaplot(gsea_result, geneSetID = gsea_result@result$ID[i], title = pathway_name)
    
    # If the result is a gglist, extract the ggplot objects
    if (inherits(p, "gglist")) {
      plots <- c(plots, p)  # Concatenate the list of ggplot objects into the main list
    } else if (inherits(p, "ggplot")) {
      plots[[length(plots) + 1]] <- p  # If it's a single ggplot object, add it to the list
    }
  }
  return(plots)
}

# Save plots to a multi-page PDF with 2 plots per page
save_plots_to_pdf <- function(plots, file_name, ncol = 1, nrow = 2) {
  pdf(file_name, width = 8.5, height = 11)  # Adjust width and height as needed
  for (i in seq(1, length(plots), by = ncol * nrow)) {
    gridExtra::grid.arrange(grobs = plots[i:min(i + ncol * nrow - 1, length(plots))], ncol = ncol, nrow = nrow)
  }
  dev.off()
}
# Create GSEA plots for each pathway
c2_plots <- create_gsea_plots(gsea_c2)
#h_plots <- create_gsea_plots(gsea_h)
combined_plots <- create_gsea_plots(gsea_combined)

# Save the plots to PDF
save_plots_to_pdf(c2_plots, "../results/7_GSEA/MSH2_gsea_c2.pdf")
#save_plots_to_pdf(h_plots, "../results/7_GSEA/MSH2_gsea_h.pdf")
save_plots_to_pdf(combined_plots, "../results/7_GSEA/MSH2_gsea_combined.pdf")
