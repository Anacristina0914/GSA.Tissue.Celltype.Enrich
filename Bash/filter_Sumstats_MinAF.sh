#!/bin/bash
#Script takes four arguments: 1)The name of the summary statistics file; 2)The column number of the minAlleleFreq; 3) The column number of the BP SNP position; and 4)The lower threshold (TH) for MinAlleleFreq filtering (e.g. 0.01 for MinAlleleFreq > 0.01). 

#Gets the file name without extension and creates a new variable containing the name of the new file.
file=$(ls $1 |cut -d. -f 1)"_filt"$3".txt"

#Skips the first line (header) and filters from the file the SNPs whose MinAlleleFreq >TH and <1-TH
awk -v T=$3 -v I=$(echo "1-$3"|bc -l) -v MAF=$2 'NR==1 || $MAF>T && $MAF<I' $1 > $file

#Filter for MHC from BP
