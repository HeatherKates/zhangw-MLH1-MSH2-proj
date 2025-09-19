#!/usr/bin/env Rscript

# ---------------------- #
# Annotate + Extract RNAseq Variant Effects
# ---------------------- #

suppressMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(VariantAnnotation)
  library(glue)
})

# ---- Paths ----
vcf_dir <- "../../nf-core_pipeline/OUTPUT/variant_calling"
out_dir <- "../OUTPUT/variant_tables"
bash_script <- "./run_snpeff.sh"
snpEff_db <- "mm10"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Locate filtered VCFs ----
vcf_files <- list.files(
  path = vcf_dir,
  pattern = "haplotypecaller.filtered.vcf.gz$",
  recursive = TRUE,
  full.names = TRUE
)

# ---- Annotate using snpEff (via shell wrapper) ----
annotated_vcfs <- vcf_files %>%
  map_chr(function(vcf) {
    sample_id <- str_replace(basename(vcf), ".haplotypecaller.filtered.vcf.gz", "")
    out_vcf <- file.path(out_dir, paste0(sample_id, ".annotated.vcf"))
    
    if (!file.exists(out_vcf)) {
      cmd <- glue("{bash_script} {vcf} {out_vcf} {snpEff_db}")
      message("Annotating ", sample_id)
      system(cmd)
    } else {
      message("✔ Skipping existing: ", sample_id)
    }
    
    return(out_vcf)
  })

# ---- Extract annotation tables ----
extract_annotations <- function(vcf_file) {
  sample_id <- str_replace(basename(vcf_file), ".annotated.vcf", "")
  
  vcf <- readVcf(vcf_file, genome = "GRCm38")
  
  if (!"ANN" %in% names(info(vcf))) {
    warning("No ANN field in ", sample_id)
    return(NULL)
  }

  # Extract ANN field as character vector and flatten
  ann_raw <- unlist(info(vcf)$ANN)
  ann_split <- str_split(ann_raw, "\\|", simplify = TRUE)

  ann_df <- as.data.frame(ann_split, stringsAsFactors = FALSE)
  colnames(ann_df) <- c(
    "Allele", "Effect", "Impact", "GeneName", "GeneID",
    "FeatureType", "FeatureID", "TranscriptBiotype",
    "Rank", "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos",
    "AA.pos", "Distance", "ErrorsWarnings"
  )[1:ncol(ann_df)]

  # Get variant metadata (one row per variant)
  variant_ranges <- rowRanges(vcf)
  variant_df <- data.frame(
    Chrom = as.character(seqnames(variant_ranges)),
    Pos   = start(variant_ranges),
    Ref   = as.character(ref(vcf)),
    Alt   = sapply(alt(vcf), function(x) paste(as.character(x), collapse = ",")),
    Qual  = qual(vcf)
  )

  # Repeat variant rows to match number of ANN rows
  variant_df_expanded <- variant_df[rep(seq_len(nrow(variant_df)), elementNROWS(info(vcf)$ANN)), ]

  final_df <- bind_cols(
    Sample = sample_id,
    variant_df_expanded,
    ann_df
  )

  return(final_df)
}

all_tables <- map_dfr(annotated_vcfs, extract_annotations)

# ---- Save results ----
write_csv(all_tables, file.path(out_dir, "all_samples_annotated_variants.csv"))

all_tables %>%
  group_split(Sample) %>%
  walk(function(df) {
    out_path <- file.path(out_dir, paste0(unique(df$Sample), "_variants.csv"))
    write_csv(df, out_path)
  })

message("Done. Results written to: ", out_dir)
