#!/usr/bin/env Rscript

library(readr)
library(openxlsx)
library(tools)

# ===============================
# 1. Paths and Setup
# ===============================
input_dir <- "../OUTPUT/variant_tables"  
output_dir <- "/orange/cancercenter-dept/web/public/BCB-SR/Zhangw/TANZIA/MLH1MSH2/Jiao_rnaseq_var/"
output_file <- file.path(output_dir, "variant_summary.xlsx")
public_link <- "https://data.rc.ufl.edu/pub/cancercenter-dept/BCB-SR/Zhangw/TANZIA/MLH1MSH2/Jiao_rnaseq_var/variant_summary.xlsx"

# ===============================
# 2. Column Descriptions
# ===============================
readme_df <- data.frame(
  Column = c(
    "Sample", "Chrom", "Pos", "Ref", "Alt", "Qual", "Allele", "Effect", "Impact",
    "GeneName", "GeneID", "FeatureType", "FeatureID", "TranscriptBiotype", "Rank",
    "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance", "ErrorsWarnings"
  ),
  Description = c(
    "Sample name or ID",
    "Chromosome of the variant",
    "Genomic position (1-based)",
    "Reference allele",
    "Alternate (variant) allele",
    "Phred-scaled variant call quality score",
    "Allele used in annotation (typically same as Alt)",
    "Predicted effect on gene or transcript",
    "Impact of the variant (HIGH, MODERATE, LOW, MODIFIER)",
    "Gene symbol",
    "Gene ID (e.g., Ensembl ID)",
    "Type of affected feature (e.g., transcript)",
    "ID of the affected feature (e.g., transcript ID)",
    "Biotype of the transcript (e.g., protein_coding)",
    "Exon or intron rank in transcript (e.g., 2/7)",
    "HGVS notation at cDNA level (e.g., c.215C>G)",
    "HGVS notation at protein level (e.g., p.Pro72Arg)",
    "Position in cDNA",
    "Position in coding sequence",
    "Position of affected amino acid",
    "Distance to gene if intergenic",
    "Errors or warnings from SnpEff"
  ),
  stringsAsFactors = FALSE
)

# ===============================
# 3. Load CSVs (ignore all_samples_*)
# ===============================
csv_files <- list.files(input_dir, pattern = "^[^aA].*_variants\\.csv$", full.names = TRUE)

variant_list <- lapply(csv_files, read_csv)
names(variant_list) <- file_path_sans_ext(basename(csv_files))

# ===============================
# 3. Filter CSVs
# ===============================

# Define which effects to keep
keep_effects <- c("missense_variant", "stop_gained", "splice_acceptor_variant", "splice_donor_variant")

# Filter each data frame in the list
filtered_variant_list <- lapply(variant_list, function(df) {
  df <- df[df$Impact %in% c("MODERATE", "HIGH") &
             df$TranscriptBiotype == "protein_coding" &
             (is.na(df$Distance) | df$Distance == 0) &
             df$Effect %in% keep_effects, ]
  return(df)
})


# ===============================
# 5. Create Workbook
# ===============================
wb <- createWorkbook()

# Add README sheet first
addWorksheet(wb, "README")
writeData(wb, "README", readme_df)

# Add each sample sheet
for (sample_name in names(filtered_variant_list)) {
  sheet_name <- substr(sample_name, 1, 31)  # Excel sheet name limit
  message("Writing sheet: ", sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, filtered_variant_list[[sample_name]])
}

# ===============================
# 6. Save Workbook
# ===============================
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveWorkbook(wb, output_file, overwrite = TRUE)

# Set open file permissions (rw-rw-r--)
Sys.chmod(output_file, mode = "0664")

cat("✅ Excel file created and accessible at:\n", public_link, "\n")
