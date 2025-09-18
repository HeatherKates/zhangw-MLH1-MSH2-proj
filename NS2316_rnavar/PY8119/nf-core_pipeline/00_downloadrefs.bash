#!/bin/bash
set -euo pipefail

# Define target directory
TARGET_DIR="/orange/zhangw/GENOME_REFS/Mus_musculus/GenomicVariants/GRCm38"

# Create directory if it doesn't exist
mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

# List of URLs to download
URLS=(
  "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
  "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz.tbi"
  "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz"
  "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz.tbi"
)

# Download each file if it doesn't exist
for url in "${URLS[@]}"; do
  # Extract filename from URL
  filename=$(basename "$url")
  
  if [[ -f "$filename" ]]; then
    echo "File $filename already exists, skipping download."
  else
    echo "Downloading $filename..."
    wget "$url"
  fi
done

echo "All files processed in $TARGET_DIR"
