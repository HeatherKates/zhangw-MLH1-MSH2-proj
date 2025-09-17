# Load necessary libraries
library(tximport)
library(DESeq2)
library(biomaRt)


# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Please provide both 'dir' and 'out_dir' as arguments.")
}

# Set the variables from the command line
dir <- args[1]
out_dir <- args[2]

# Define the directory containing the quantification results
dir <- "results/3_salmon/"
#out_dir="../results/4_deseq2/"
#out_dir= "test_4a_out/"

# Get the list of quantification directories starting with "MLH1"
dirs <- list.files(dir, pattern="^MLH1.*_quant$", full.names=TRUE)

# Create a named vector of files for tximport
files <- file.path(dirs, "quant.sf")
names(files) <- basename(dirs)

# Define the conditions based on filenames
sample_names <- names(files)
conditions <- ifelse(grepl("KO", sample_names), "KO", "HA")

# Create colData dataframe
colData <- data.frame(
  row.names = sample_names,
  condition = factor(conditions, levels = c("HA", "KO"))
)

# Use biomaRt to get the transcript-to-gene mapping
options(timeout = 600)  # 10 minutes
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
tx2gene <- getBM(attributes=c("ensembl_transcript_id_version", "ensembl_gene_id"),
                 mart=mart)

# Import quantification data using tximport and summarize to gene level
txi <- tximport(files, type="salmon", tx2gene=tx2gene, txOut=FALSE)

# Create DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi, colData, design=~condition)

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# replace with gene names
# Get the gene names for the Ensembl gene IDs
gene_ids <- rownames(dds)
genes <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
               filters="ensembl_gene_id", 
               values=gene_ids, 
               mart=mart)

# Ensure that all gene IDs are mapped
genes <- genes[match(gene_ids, genes$ensembl_gene_id),]

# Store the mapping in a separate data frame
gene_mapping <- data.frame(ensembl_gene_id = gene_ids, 
                           external_gene_name = genes$external_gene_name)


# Add gene mapping to dds@elementMetadata
dds@elementMetadata <- DataFrame(gene_mapping)
# Add gene symbols to the results data frame
res$geneSymbol <- gene_mapping$external_gene_name

# Save results
MLH1_dds <- dds
MLH1_res <- res
MLH1_gene_mapping <- gene_mapping

write.csv(as.data.frame(MLH1_res), file=paste0(out_dir,"MLH1_DESeq2_gene_results.csv")
save.image(file=paste0(out_dir,"MLH1_DESeq2_result.RDATA")

saveRDS(MLH1_dds,paste0(out_dir,"MLH1_dds.Rds")
saveRDS(MLH1_res,paste0(out_dir,"MLH1_res.Rds")
saveRDS(MLH1_gene_mapping,paste0(out_dir,"MLH1_gene_mapping.Rds")

#Write more results (counts, etc.)
# Extract the counts data frame
counts_df <- as.data.frame(MLH1_dds@assays@data@listData[["counts"]])
counts_df$geneSymbol <- res$geneSymbol
counts_df$ensembl_gene_id <- rownames(counts_df)
counts_df <- counts_df[,c("geneSymbol","ensembl_gene_id",colnames(counts_df)[1:6])]
avgTxLength_df <- as.data.frame(MLH1_dds@assays@data@listData[["avgTxLength"]])
normalizationFactors_df <- as.data.frame(MLH1_dds@assays@data@listData[["normalizationFactors"]])
mu_df <- as.data.frame(MLH1_dds@assays@data@listData[["mu"]])
H_df <- as.data.frame(MLH1_dds@assays@data@listData[["H"]])
cooks_df <- as.data.frame(MLH1_dds@assays@data@listData[["cooks"]])

# Extract the geneSymbol and ensembl_gene_id columns from counts_df
cols_to_add <- counts_df[, c("geneSymbol", "ensembl_gene_id")]
avgTxLength_df <- cbind(cols_to_add,avgTxLength_df)
normalizationFactors_df <- cbind(cols_to_add,normalizationFactors_df )
mu_df <- cbind(cols_to_add,mu_df)
H_df <- cbind(cols_to_add,H_df)
cooks_df <- cbind(cols_to_add,cooks_df)

# Create a list of the data frames
list_of_dfs <- list(
  counts = counts_df,
  avgTxLength = avgTxLength_df,
  normalizationFactors = normalizationFactors_df,
  mu = mu_df,
  H = H_df,
  cooks = cooks_df
)
# Set the column names based on counts_df (including the added columns)
new_colnames <-colnames(counts_df)

# Apply the new column names to each data frame
list_of_dfs <- lapply(list_of_dfs, function(df) {
  colnames(df) <- new_colnames
  return(df)
})
# Create the README data frame
README <- data.frame(
  sheet = names(list_of_dfs),
  description = c(
    "Raw gene counts per sample. Rows correspond to genes and columns correspond to samples.",
    "Average transcript length for each gene in each sample. Rows correspond to genes and columns correspond to samples.",
    "Normalization factors used to scale the counts for each sample. Rows correspond to genes and columns correspond to samples.",
    "Fitted mean values (mu) for each gene in each sample. Rows correspond to genes and columns correspond to samples.",
    "H values, which are intermediate values used in DESeq2’s dispersion estimation. Rows correspond to genes and columns correspond to samples.",
    "Cook's distance values, indicating the influence of each sample on the estimation of the model parameters for each gene. Rows correspond to genes and columns correspond to samples."
  )
)
# Add the README df to the beginning of the list_of_dfs
list_of_dfs <- c(list(README = as.data.frame(README)), list_of_dfs)

#Write to excel
library(openxlsx)
# Define the output file path
output_file <- paste0(out_dir,"MLH1_DESeq2_analysis_dataframes.xlsx")

# Create a new Excel workbook
wb <- createWorkbook()

# Add each data frame as a separate sheet in the workbook
for (sheet_name in names(list_of_dfs)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = list_of_dfs[[sheet_name]])
}

# Save the workbook
saveWorkbook(wb, output_file, overwrite = TRUE)

