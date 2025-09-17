library(edgeR)
library(BiocParallel)
register(SerialParam())
library(dplyr)
# Set up parallelization to use 4 cores
bp_param <- MulticoreParam(workers = 6)
library(DiffBind)

samples <- data.frame(
  SampleID = c("MLH1KO-1", "MLH1KO-2", "MLH1KO-3", "MLH1R4-1", "MLH1R4-2", "MLH1R4-3"),
  Tissue = c(rep("MLH1KO", 3), rep("MLH1R4", 3)),  #  'Tissue' can represent conditions or leave this column out.
  Factor = c(rep("MLH1KO", 3), rep("MLH1R4", 3)),
  Condition = c(rep("MLH1KO", 3), rep("MLH1R4", 3)),  # The experimental conditions
  Treatment = rep("None", 6),  # Since there's no treatment in this setup
  Replicate = c(1, 2, 3, 1, 2, 3),  # Indicating replicates
  bamReads = c("results/6_bigwig/MLH1KO-1_sort_n.bam", "results/6_bigwig/MLH1KO-2_sort_n.bam", "results/6_bigwig/MLH1KO-3_sort_n.bam",
               "results/6_bigwig/MLH1R4-1_sort_n.bam", "results/6_bigwig/MLH1R4-2_sort_n.bam", "results/6_bigwig/MLH1R4-3_sort_n.bam"),
  Peaks = c("results/5_genrich/MLH1KO-1.narrowPeak", "results/5_genrich/MLH1KO-2.narrowPeak", "results/5_genrich/MLH1KO-3.narrowPeak",
            "results/5_genrich/MLH1R4-1.narrowPeak", "results/5_genrich/MLH1R4-2.narrowPeak", "results/5_genrich/MLH1R4-3.narrowPeak"),
  PeakCaller = rep("narrow", 6),  # For Genrich, 'narrow' for narrowPeak format
  stringsAsFactors = FALSE
)


# Initialize DiffBind object
dbaObj <- NULL
dbaObj <- dba(sampleSheet=samples)
# Count reads in peaks

# Count reads in consensus peaks using 6 cores
dbaObj <- dba.count(dbaObj, summits=250, bParallel = TRUE,minOverlap = 2)

# Define contrasts (KO vs HLA)
dbaObj <- dba.contrast(dbaObj, categories = DBA_CONDITION)

# Perform differential binding analysis
dbaObj <- dba.analyze(dbaObj)
saveRDS(dbaObj,file="results/7_diffbind/MLH1_dbaObj.RDATA")
# Generate and inspect the report of differentially bound regions
diffPeaks <- dba.report(dbaObj)
write.csv(diffPeaks, file = "results/7_diffbind/MLH1_differential_binding_results.csv")
diffPeaksIn <- read.csv("results/7_diffbind/MLH1_differential_binding_results.csv")

# Function to write peaks data frame to a BED file
write_peaks_to_bed <- function(peaks_df, output_file) {
  
  # Ensure required columns are present
  if (!all(c("seqnames", "start", "end", "X", "Fold") %in% colnames(peaks_df))) {
    stop("The data frame does not contain the necessary columns: 'seqnames', 'start', 'end', 'X', 'Fold'")
  }
  
  # Create the BED format data frame
  bed <- data.frame(
    chr = peaks_df$seqnames,      # Chromosome column
    start = peaks_df$start,       # Start position
    end = peaks_df$end,           # End position
    name = peaks_df$X,            # Peak ID or use `.` if no name is required
    score = peaks_df$Fold         # Log2 fold-change as score
  )
  
  # Write the BED file
  write.table(bed, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  message("BED file written to: ", output_file)
}

# Example usage with your dataframe 'diffPeaksIn'
write_peaks_to_bed(diffPeaksIn, "results/7_diffbind/MLH1_differential_binding_results.bed")

# MLH1 unique peaks

# Set thresholds
min_Conc_MLH1R4 <- 2  # Adjust this value based on data and expectations
max_Conc_MLH1KO <- 0.2  # Adjust based on definition of "low" in KO

# Filter peaks for positive fold change, low concentration in KO, and reasonably high concentration in MLH1
MLH1R4_filtered_peaks <- diffPeaksIn %>%
  filter(Fold > 0 & Conc_MLH1KO < max_Conc_MLH1KO & Conc_MLH1R4 > min_Conc_MLH1R4)

# Check the filtered results
write.csv(MLH1R4_filtered_peaks,"results/7_diffbind/MLH1R4_filtered_differential_binding_results.csv")
write_peaks_to_bed(MLH1R4_filtered_peaks, "results/7_diffbind/MLH1R4_filtered_differential_binding_results.bed")

# MLH1KO unique peaks

# Set thresholds
min_Conc_MLH1KO <- 2  # Adjust this value based on data and expectations
max_Conc_MLH1R4 <- 0.2  # Adjust based on definition of "low" in MLH1

# Filter peaks for positive fold change, low concentration in KO, and reasonably high concentration in MLH1
MLH1KO_filtered_peaks <- diffPeaksIn %>%
  filter(Fold < 0 & Conc_MLH1R4 < max_Conc_MLH1R4 & Conc_MLH1KO > min_Conc_MLH1KO)

# Check the filtered results
write.csv(MLH1KO_filtered_peaks,"results/7_diffbind/MLH1KO_filtered_differential_binding_results.csv")
write_peaks_to_bed(MLH1KO_filtered_peaks, "results/7_diffbind/MLH1KO_filtered_differential_binding_results.bed")

#From this point, you can view filtered peaks and per-sample reads in IGV by loading:
#"results/7_diffbind/filtered_differential_binding_results.bed"
#"results/6_bigwig/MLH1KO-1.bw","results/6_bigwig/MLH1KO-2.bw","results/6_bigwig/MLH1KO-3.bw","results/6_bigwig/MLH1R4-1.bw","results/6_bigwig/MLH1R4-2.bw","results/6_bigwig/MLH1R4-3.bw")
#Write all the results to an excel

##Annotate peaks with chipseeker
# Load necessary libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # You can change the genome version accordingly
library(clusterProfiler)

# Load the peaks (ensure peaks are in BED format)
peaks <- readPeakFile("results/7_diffbind/MLH1_differential_binding_results.bed")

library(GenomicRanges)

# Ensure peaks are in the GRanges format, then add the "chr" prefix
seqlevelsStyle(peaks) <- "UCSC"  # This will add "chr" to the chromosome names

# Annotate the peaks to genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Now, run the annotation
peak_annotation <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=txdb)

# Add gene symbols

# Load the necessary library for gene symbol mapping
library(org.Hs.eg.db)  # For human genes, change if using another species

# Extract the Entrez IDs from the annotation result
entrez_ids <- as.data.frame(peak_annotation)$geneId

# Map Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# Save the annotated results
peak_annotation_df <- as.data.frame(peak_annotation)

# Add the gene symbols to the peak_annotation result
peak_annotation_df$geneSymbol <- gene_symbols

#Add peak annotation to diffbind results
# Add a "X" column
peak_annotation_df$X=rownames(peak_annotation_df)

# merge peaks and annotation
peaks_with_anno <- merge(diffPeaksIn,peak_annotation_df,by="X")

# Remove redundant .y columns (from ChIPseekresults/diffbind for seqnames)
peaks_with_anno <- peaks_with_anno %>%
  dplyr::select(-seqnames.y, -start.y, -end.y, -width.y, -strand.y,-seqnames.x)

# Add group specific significance results
# Create a new column to indicate in which group the peak is significant
peaks_with_anno$Unique_to_group <- ifelse(
  peaks_with_anno$X %in% MLH1KO_filtered_peaks$X, "MLH1KO",
  ifelse(peaks_with_anno$X %in% MLH1R4_filtered_peaks$X, "MLH1R4", "None")
)

# Add DE results

#Read in DE results
DE_results <- read.csv("../RNAseq/MLH1/results/4_deseq2/MLH1_DESeq2_gene_results.csv")

# Specify the direction of the log2FoldChange
DE_results <- DE_results %>% dplyr::rename(log2FoldChange_KOvsR4=log2FoldChange)

# Rename ensemble ID col
DE_results <- DE_results %>% dplyr::rename(ensembl_gene_id=X)

# Add entrez ID for merge with peak annotation
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Extract the Ensembl Gene IDs from DE_results dataframe
ensembl_ids <- DE_results$ensembl_gene_id

# Use biomaRt to retrieve the corresponding Entrez Gene IDs
gene_info <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                   filters = "ensembl_gene_id", values = ensembl_ids, mart = ensembl)

# Merge the Entrez Gene IDs with your DE_results dataframe based on the Ensembl Gene IDs
DE_results <- merge(DE_results,gene_info,by="ensembl_gene_id",all.x=TRUE)

# Add DE_result to all colnames except the last 2 col which is geneSymbol and entrez ID
colnames(DE_results)[1:7] <- paste0("DE_result.",colnames(DE_results)[1:7])

# Merge DE result with peaks
peaks_with_anno <- merge(peaks_with_anno,DE_results,by.y="entrezgene_id",by.x="geneId",all.x=TRUE)

# Add identifier for gene symbol
peaks_with_anno <- peaks_with_anno %>% dplyr::rename(DE_result.geneSymbol=geneSymbol.y)

# Rename "X" to be more informative
peaks_with_anno <- peaks_with_anno %>% dplyr::rename(peak_ID=X)

# Remove suffix x from merges
colnames(peaks_with_anno) <- gsub(".x$","",colnames(peaks_with_anno))

# Remove chipseeker cols redundant with genrich cols
peaks_with_anno <- peaks_with_anno %>% dplyr::select(c(-V4,-V5))
## Add README

# Create a data frame with the column names, descriptions, and the source of the data
README_df <- data.frame(
  Column = colnames(peaks_with_anno),
  Description = c(
    "Gene ID from annotation",
    "Peak ID from Genrich peak calling",
    "Start position of the peak (Genrich)",
    "End position of the peak (Genrich)",
    "Width of the peak (Genrich)",
    "Strand of the peak (Genrich)",
    "Concentration of signal across all conditions (DiffBind)",
    "Concentration of signal in MLH1R4 condition (DiffBind)",
    "Concentration of signal in MLH1KO condition (DiffBind)",
    "Fold change between MLH1KO and MLH1R4 (DiffBind)",
    "P-value for the difference between conditions (DiffBind)",
    "FDR-adjusted p-value (DiffBind)",
    "Peak annotation (e.g., promoter, exon, intron) from ChIPseeker",
    "Chromosome of the annotated gene (ChIPseeker)",
    "Start position of the annotated gene (ChIPseeker)",
    "End position of the annotated gene (ChIPseeker)",
    "Length of the annotated gene (ChIPseeker)",
    "Strand of the annotated gene (ChIPseeker)",
    "Transcript ID from ChIPseeker annotation",
    "Distance from the peak to the nearest TSS (ChIPseeker)",
    "Gene symbol from ChIPseeker annotation",
    "Indicates whether the peak is unique to a specific group (DiffBind analysis)",
    "Ensembl gene ID from DESeq2 results",
    "Base mean expression level (DESeq2)",
    "Log2 fold change between MLH1KO and MLH1R4 (DESeq2)",
    "Standard error of the log2 fold change (DESeq2)",
    "Statistical test statistic (DESeq2)",
    "P-value for differential expression (DESeq2)",
    "FDR-adjusted p-value for differential expression (DESeq2)",
    "Gene symbol from DESeq2 results"
  ),
  `Analysis that produced field` = c(
    "ChIPseeker",
    "Genrich",
    "Genrich",
    "Genrich",
    "Genrich",
    "Genrich",
    "DiffBind",
    "DiffBind",
    "DiffBind",
    "DiffBind",
    "DiffBind",
    "DiffBind",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "DiffBind",
    "DESeq2",
    "DESeq2",
    "DESeq2",
    "DESeq2",
    "DESeq2",
    "DESeq2",
    "DESeq2",
    "DESeq2"
  ),
  stringsAsFactors = FALSE
)


library(writexl)

dfs <- list(
  README = README_df,
  peaks_with_anno = peaks_with_anno
)

#Add the gene symbol to the annotation
# Write the list of data frames to an Excel file using writexl
write_xlsx(dfs, "results/7_diffbind/MLH1_Differential_binding_results.xlsx")

