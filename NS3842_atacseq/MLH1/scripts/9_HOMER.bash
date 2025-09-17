#!/bin/bash
module load homer

# Add "chr" prefix to the chromosome column in the BED file
sed 's/^/chr/' results/7_diffbind/MLH1R4_filtered_differential_binding_results.bed > results/9_HOMER/MLH1R4_filtered_differential_binding_results_chr.bed
# Add "chr" prefix to the chromosome column in the BED file
sed 's/^/chr/' results/7_diffbind/MLH1KO_filtered_differential_binding_results.bed > results/9_HOMER/MLH1KO_filtered_differential_binding_results_chr.bed

findMotifsGenome.pl results/9_HOMER/MLH1KO_filtered_differential_binding_results_chr.bed hg38 results/9_HOMER/MLH1KO_peaks_motifs -size 200
findMotifsGenome.pl results/9_HOMER/MLH1R4_filtered_differential_binding_results_chr.bed hg38 results/9_HOMER/MLH1R4_peaks_motifs -size 200

#Make csv files for R
sed 's:/Homer:/Homer,:g' < results/9_HOMER/MLH1KO_peaks_motifs/knownResults.txt | tr '\t' ',' > results/9_HOMER/MLH1KO_peaks_motifs/knownResults.csv
sed 's:/Homer:/Homer,:g' < results/9_HOMER/MLH1R4_peaks_motifs/knownResults.txt | tr '\t' ',' > results/9_HOMER/MLH1R4_peaks_motifs/knownResults.csv
