library(readxl)
MLH1 <- readxl::read_excel("MLH1/results/01_diffbind/MLH1R4_peaks.xlsx",sheet = "peaks_with_anno")
MSH2 <- readxl::read_excel("MSH2/results/01_diffbind/MSH2R4_peaks.xlsx",sheet="peaks_with_anno")
library(dplyr)

# Identify geneSymbols unique to each group
mlh1_genes <- setdiff(MLH1$geneSymbol, MSH2$geneSymbol)  # Genes only in MLH1
msh2_genes <- setdiff(MSH2$geneSymbol, MLH1$geneSymbol)  # Genes only in MSH2

# Filter rows based on unique geneSymbols
mlh1_unique <- MLH1 %>%
  filter(geneSymbol %in% mlh1_genes) %>%
  mutate(unique_to_group = "MLH1")

msh2_unique <- MSH2 %>%
  filter(geneSymbol %in% msh2_genes) %>%
  mutate(unique_to_group = "MSH2")

# Combine the results
combined_df <- bind_rows(mlh1_unique, msh2_unique)


README <- readxl::read_excel("MSH2/results/01_diffbind/MSH2R4_peaks.xlsx")
lastrow <- data.frame(cbind("unique_to_group","The analysis in which this peak was uniquely identified (i.e. this peak is not present in the other group","setdiff (R)"))
colnames(lastrow) <- colnames(README)
README <- rbind(README,lastrow)
data <- list(README,combined_df)
write_xlsx(data,"MLH1_MSH2_unique_peaks.xlsx")
