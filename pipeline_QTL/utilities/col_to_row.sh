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
             #      1 \
             #      2 \
             #      AIRWAVE_1DNMR_BatchCorrected_log_Data.tsv \
             #      AIRWAVE_1DNMR_BatchCorrected_log_Data_labelled.tsv 


# Variables to substitute

# File to transpose:
infile1=$1

# Column to transpose:
col_to_transpose=$2

# Pass an integer for the index of where to cut from (e.g. 1 to keep the header
# of the column, 1 if there is no header; or 2 to skip it header):
# is to be kept, no header present, etc.)
index_col1=$3

# File that will be concatenated below infile1 (if number of rows and columns
# between files don't match your outfile will also be mismatched):
infile2=$4

# Index of where the pasting will start in the second file
# (i.e. position/column number of the
# second file  at which the first element
# of the column from the first file will be pasted across):
#index_col2=$4
# Name of new file:
outfile=$5

# Run command:
cut -f${col_to_transpose} ${infile1} | paste -s - | cut -f${index_col1}- | cat - ${infile2} > ${outfile}

# TO DO:
# Add row counts and column counts for each file and out file as sanity check:

echo 'Done'
