#!/bin/bash
module load multiqc
# Run multiqc to summarize the fastqc results
multiqc results/2_fastQC/ -o results/2_fastQC/
