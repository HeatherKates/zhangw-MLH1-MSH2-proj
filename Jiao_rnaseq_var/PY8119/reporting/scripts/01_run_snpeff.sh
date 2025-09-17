#!/bin/bash

# run_snpeff.sh
# Usage: ./run_snpeff.sh <input_vcf> <output_vcf> <snpEff_db>

SNPEFF_JAR="/apps/snpeff/5.1d/jar/snpEff.jar"
SNPEFF_CONFIG="/apps/snpeff/5.1d/conf/snpEff.config"
SNPEFF_DATA="/apps/snpeff/5.1d/conf/data"

java -Xmx4g -jar "$SNPEFF_JAR" \
  -c "$SNPEFF_CONFIG" \
  -dataDir "$SNPEFF_DATA" \
  "$3" "$1" > "$2"
