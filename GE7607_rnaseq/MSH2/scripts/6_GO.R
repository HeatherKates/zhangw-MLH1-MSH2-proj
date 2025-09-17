# Load the libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(BiocParallel)

#Which dataset
dataset="MSH2"
# Load data
MSH2_res <- readRDS("../results/4_deseq2/MSH2_res.Rds")
MSH2_dds <- readRDS("../results/4_deseq2/MSH2_dds.Rds")

# Convert DESeq results to a data frame
MSH2_res_df <- as.data.frame(MSH2_res)
MSH2_res_df$external_gene_name <- rownames(MSH2_res_df) #Don't know why this is different than MLH1, should be edited upstream but for now OK

# Filter genes by adjusted p-value < 0.05
significant_genes <- na.omit(MSH2_res_df[MSH2_res_df$padj < 0.05, ])

# Create lists for upregulated and downregulated genes
up_genes <- significant_genes[significant_genes$log2FoldChange > 0.58, ]
down_genes <- significant_genes[significant_genes$log2FoldChange < -0.58, ]

# Perform GO enrichment analysis for upregulated genes
go_enrichment_up <- enrichGO(
  gene          = up_genes$external_gene_name,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "ALL",  # Perform enrichment for all GO categories (BP, MF, CC)
  pAdjustMethod = "BH",
  #pvalueCutoff  = 0.05,
  #qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Perform GO enrichment analysis for downregulated genes
go_enrichment_down <- enrichGO(
  gene          = down_genes$external_gene_name,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "ALL",  # Perform enrichment for all GO categories (BP, MF, CC)
  pAdjustMethod = "BH",
  #pvalueCutoff  = 0.05,
  #qvalueCutoff  = 0.2,
  readable      = TRUE
)

# Function to create a combined dot plot with faceting by ontology using tryCatch
create_combined_go_plot <- function(go_result, dataset, direction) {
  tryCatch({
    if (nrow(go_result@result) > 0) {
      go_result_df <- as.data.frame(go_result) %>% arrange(p.adjust) %>% slice_head(n=20)
      
      # Create the dot plot with faceting by ontology
      p <- ggplot(go_result_df, aes(x = GeneRatio, y = Description)) +
        geom_point(aes(size = Count, color = -log10(p.adjust))) +
        facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
        scale_color_gradient(low = "blue", high = "red") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 10),
              strip.text.y = element_text(size = 12, face = "bold")) +
        labs(color = "-log10(p.adjust)", size = "Gene Count", x = "Gene Ratio", y = "GO Term") +
        ggtitle(paste0("Up to top 20 GO results for genes ", direction, "-regulated in ", dataset, " KO vs. HA (all ont categories)"))
      
      return(p)
    } else {
      stop("No enrichment results found.")
    }
  }, error = function(e) {
    message("No enrichment results found or another error occurred: ", e$message)
    return(NULL)
  })
}
go_combined_plot_down=NULL
go_combined_plot_down <- create_combined_go_plot(go_enrichment_down,dataset,direction="down")

go_combined_plot_up=NULL
go_combined_plot_up <- create_combined_go_plot(go_enrichment_up,dataset,direction="up")
# Save the plots as single-page PDFs

# Save the plots as single-page PDFs
if(!is.null(go_combined_plot_down)){
ggsave("../results/8_GO/MSH2_combined_go_enrichment_down.pdf",
       plot = go_combined_plot_down, width = 12, height = 16) } else {
         print("no down plot to print")
}

if(!is.null(go_combined_plot_up)){
ggsave("../results/8_GO/MSH2_combined_go_enrichment_up.pdf",
       plot = go_combined_plot_up, width = 12, height = 16) } else {
         print("no down plot to print")
       }

#Write the results to a file
write.csv(file="../results/8_GO/MSH2_combined_go_enrichment_up.csv",go_enrichment_up@result)
write.csv(file="../results/8_GO/MSH2_combined_go_enrichment_down.csv",go_enrichment_down@result)
