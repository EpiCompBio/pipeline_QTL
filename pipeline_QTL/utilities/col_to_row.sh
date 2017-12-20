#!/usr/bin/env bash


# Transposes a column from a tsv file and inserts it as a row into a second
# file.
# Run as e.g.:
#bash col_to_row.sh AIRWAVE_1DNMR_BatchCorrected_log_VarInfo.tsv 2 AIRWAVE_1DNMR_BatchCorrected_log_Data.tsv AIRWAVE_1DNMR_BatchCorrected_log_Data_labelled.tsv 


# Set bash script options:
# https://kvz.io/blog/2013/11/21/bash-best-practices/
set -o errexit
set -o pipefail
set -o nounset


# Variables to substitute:
# File to transpose:
infile1=$1
# Column to transpose. Pass an integer (e.g. 1 or 2 if header from first file is to be kept)
index_col=$2
# File that will be concatenated below infile1:
infile2=$3
# Name of new file:
outfile=$4

# Run command:
cut -f1 ${infile1} | paste -s - | cut -f${index_col}- | cat - ${infile2} > ${outfile}
echo 'Done'

