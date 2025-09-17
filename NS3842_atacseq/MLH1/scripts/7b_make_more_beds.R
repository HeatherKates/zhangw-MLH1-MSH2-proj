# Read in the annotated peak data
library(readxl)
df <- read_excel("results/7_diffbind/MLH1_Differential_binding_results.xlsx",sheet = "peaks_with_anno")

# Make a peak df of just MLH1KO unique peaks that are also DE (p adj < 0.05) and log2FC > 0 in MLH1KO vs MLH1R4
MLH1KO_unique_peaks_DE_genes <- df %>% filter(Unique_to_group == "MLH1KO") %>% filter (DE_result.padj < 0.0509) %>% filter (DE_result.log2FoldChange_KOvsR4 > 0)
MLH1R4_unique_peaks_DE_genes <- df %>% filter(Unique_to_group == "MLH1R4") %>% filter (DE_result.padj < 0.0509) %>% filter (DE_result.log2FoldChange_KOvsR4 < 0)


# 1. Select and rename relevant columns
# Select the necessary columns and reverse the DE log2FoldChange values
MLH1KO_df_for_bed <- MLH1KO_unique_peaks_DE_genes %>%
  dplyr::select(seqnames = geneChr, 
                start, 
                end, 
                name = peak_ID, 
                score = Fold,  # Will rename later in write function
                gene_expression_log2FC = DE_result.log2FoldChange_KOvsR4) %>%
  dplyr::mutate(gene_expression_log2FC = -gene_expression_log2FC)  # Reverse the log2FoldChange values

MLH1R4_df_for_bed <- MLH1R4_unique_peaks_DE_genes %>%
  dplyr::select(seqnames = geneChr, 
                start, 
                end, 
                name = peak_ID, 
                score = Fold,  # Will rename later in write function
                gene_expression_log2FC = DE_result.log2FoldChange_KOvsR4) %>%
  dplyr::mutate(gene_expression_log2FC = -gene_expression_log2FC)  # Reverse the log2FoldChange values


write_gene_expression_log2FC_to_bedgraph <- function(peaks_df, output_file) {
  
  # Ensure required columns are present
  if (!all(c("seqnames", "start", "end", "gene_expression_log2FC") %in% colnames(peaks_df))) {
    stop("The data frame does not contain the necessary columns: 'seqnames', 'start', 'end', 'gene_expression_log2FC'")
  }
  
  # Create the BEDGraph format data frame
  bedgraph <- data.frame(
    chr = peaks_df$seqnames,              # Chromosome column
    start = peaks_df$start,               # Start position
    end = peaks_df$end,                   # End position
    value = peaks_df$gene_expression_log2FC  # Gene expression log2FC as the value
  )
  
  # Write the BEDGraph file
  write.table(bedgraph, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  message("BEDGraph file written to: ", output_file)
}

# Example usage

write_binding_log2FC_to_bedgraph <- function(peaks_df, output_file) {
  
  # Ensure required columns are present
  if (!all(c("seqnames", "start", "end", "score") %in% colnames(peaks_df))) {
    stop("The data frame does not contain the necessary columns: 'seqnames', 'start', 'end', 'score'")
  }
  
  # Create the BEDGraph format data frame
  bedgraph <- data.frame(
    chr = peaks_df$seqnames,      # Chromosome column
    start = peaks_df$start,       # Start position
    end = peaks_df$end,           # End position
    value = peaks_df$score        # Binding log2FC as the value
  )
  
  # Write the BEDGraph file
  write.table(bedgraph, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  message("BEDGraph file written to: ", output_file)
}

write_gene_expression_log2FC_to_bedgraph(MLH1KO_df_for_bed, "results/7_diffbind/MLH1K0_peaks_gene_expression_log2FC.bedgraph")
write_binding_log2FC_to_bedgraph(MLH1KO_df_for_bed, "results/7_diffbind/MLH1KO_peaks_binding_log2FC.bedgraph")

write_gene_expression_log2FC_to_bedgraph(MLH1R4_df_for_bed, "results/7_diffbind/MLH1R4_gene_expression_log2FC.bedgraph")
write_binding_log2FC_to_bedgraph(MLH1R4_df_for_bed, "results/7_diffbind/MLH1R4_binding_log2FC.bedgraph")

