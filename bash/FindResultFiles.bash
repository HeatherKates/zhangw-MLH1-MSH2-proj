#!/bin/bash

# Find and print full file paths for the specified file types
find "$(pwd)" \( \
    -name "*.bw" -o \
    -name "*.bam" -o \
    -name "*.bai" -o \
    -name "*.csv" -o \
    -name "*.xlsx" -o \
    -name "*.RDATA" -o \
    -iname "*peaks*" \
\) -type f

