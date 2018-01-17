#!/usr/bin/env bash


# Insert a column taken from one file into another
# Run as e.g.:
#bash insert_col.sh 1 AIRWAVE_1DNMR_BatchCorrected_log_SampleInfo.tsv AIRWAVE_1DNMR_BatchCorrected_log_Data_Var.tsv AIRWAVE_1DNMR_BatchCorrected_log_Data_Var_Sample.tsv 


# Set bash script options:
# https://kvz.io/blog/2013/11/21/bash-best-practices/
set -o errexit
set -o pipefail
set -o nounset


# Variables to substitute:
# Column to copy. Pass an integer:
index_col=$1
# File with column to copy:
infile1=$2
# File that will have the column inserted:
infile2=$3
# Name of new file:
outfile=$4

# Run command:
cut -f${index_col} ${infile1} | paste - ${infile2} > ${outfile}
echo 'Done'
