#!/usr/bin/env bash

# Set bash script options:
# https://kvz.io/blog/2013/11/21/bash-best-practices/
set -o errexit
set -o pipefail
set -o nounset

# Transposes the first column from a tsv file and inserts it as a row into a second
# tsv file.
# Run as e.g.:
#bash col_to_row.sh AIRWAVE_1DNMR_BatchCorrected_log_VarInfo.tsv \
             #      2 \
             #      AIRWAVE_1DNMR_BatchCorrected_log_Data.tsv \
             #      AIRWAVE_1DNMR_BatchCorrected_log_Data_labelled.tsv 

# Variables to substitute:
# File to transpose:
infile1=$1
# Column to transpose. Pass an integer (e.g. 1 or 2 if header from first file
# is to be kept, no header present, etc.)
index_col1=$2
# File that will be concatenated below infile1 (if number of rows and columns
# between files don't match your outfile will also be mismatched):
infile2=$3
# Index of where the pasting will start in the second file
# (i.e. position/column number of the
# second file  at which the first element
# of the column from the first file will be pasted across):
#index_col2=$4
# Name of new file:
outfile=$4

# Run command:
cut -f1 ${infile1} | paste -s - | cut -f${index_col1}- | cat - ${infile2} > ${outfile}

echo 'Done'
