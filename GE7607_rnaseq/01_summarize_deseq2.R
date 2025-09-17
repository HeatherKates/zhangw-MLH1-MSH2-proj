library(dplyr)
library(DESeq2)
MLH1_dds <- readRDS("MLH1/results/4_deseq2/MLH1_dds.Rds")
MSH2_dds <- readRDS("MSH2/results/4_deseq2/MSH2_dds.Rds")
MLH1_res <- readRDS("MLH1/results/4_deseq2/MLH1_res.Rds")
MSH2_res <- readRDS("MSH2/results/4_deseq2/MSH2_res.Rds")

# Function to label DE genes
label_DE <- function(res, antibody) {
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  res$antibody <- antibody
  res$DE <- ifelse(res$padj < 0.05 & res$log2FoldChange > 0, "DE (up)",
                   ifelse(res$padj < 0.05 & res$log2FoldChange < 0, "DE (down)", "Not DE"))
  return(res)
}

# Label DE genes for each result set
MLH1_res_labeled <- label_DE(MLH1_res, "MLH1")
MSH2_res_labeled <- label_DE(MSH2_res, "MSH2")

# Combine the results into one data frame
combined_res <- full_join(MLH1_res_labeled, MSH2_res_labeled, by = "gene", suffix = c("_MLH1", "_MSH2"))

# Determine the result status
combined_res <- combined_res %>%
  mutate(Result_KOvsHA = case_when(
    DE_MLH1 == "DE (up)" & DE_MSH2 == "Not DE" ~ "DE (up) in MLH1",
    DE_MLH1 == "DE (down)" & DE_MSH2 == "Not DE" ~ "DE (down) in MLH1",
    DE_MSH2 == "DE (up)" & DE_MLH1 == "Not DE" ~ "DE (up) in MSH2",
    DE_MSH2 == "DE (down)" & DE_MLH1 == "Not DE" ~ "DE (down) in MSH2",
    DE_MLH1 == "DE (up)" & DE_MSH2 == "DE (up)" ~ "DE (up) in both",
    DE_MLH1 == "DE (down)" & DE_MSH2 == "DE (down)" ~ "DE (down) in both",
    TRUE ~ "Not DE"
  ))

# Select relevant columns
final_res <- combined_res %>%
  dplyr::select(gene, padj_MLH1, log2FoldChange_MLH1, padj_MSH2, log2FoldChange_MSH2, Result_KOvsHA)

# Save to CSV
write.csv(final_res, file = "combined_DE_results.csv", row.names = FALSE)
