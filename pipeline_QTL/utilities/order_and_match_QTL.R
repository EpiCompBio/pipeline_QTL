#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
order_and_match_QTL.R
======================

Author: |author_names| 
Release: |version|
Date: |today|


Purpose
=======

|description|

Order and match columns between two files, usually a genotype and a phenotype file.
The output is meant to be passed to MatrixeQTL.


Usage and options
=================

To run, type:
Rscript order_and_match_QTL.R --file1 <FILE> --file2 <FILE> [options]

Usage: order_and_match_QTL.R (--file1 <FILE>) (--file2 <FILE>)
       order_and_match_QTL.R [options]

Options:
--file1 <FILE>                Usually a genotype input file name, columns are samples, rows are features
--file2 <FILE>                usually a phenotype input file name, columns are samples, rows are features
-O <OUTPUT_FILE>              Output file name
--session <R_SESSION_NAME>    R session name if to be saved
-h --help                     Show this screen

Input:

Two tab separated files. These are read with data.table and stringsAsFactors = FALSE
Rows must be features (phenotypes, variables, etc.) and columns must be samples (individuals)
The first row and first column must be the ID labels.

Output:

Checks that the labels and order match perfectly. If not, it orders them and outputs files to disk.

Requirements:

library(docopt)
library(data.table)


Documentation
=============

For more information see:

|url|
' -> doc
# Load docopt:
library(docopt, quietly = TRUE)
# Retrieve the command-line arguments:
args <- docopt(doc)

# Save arguments needed:
# num_PCs <- as.integer(args[['--num_PCs']])

# Print all arguments to screen:
str(args)
######################

######################
# Logging
# This can be taken care of by CGAT Experiment.py if running as a pipeline.
# Otherwise there seem to be few good alternatives. A workaround is the code in:
# XXXX/project_quickstart/templates/script_templates/logging.R
# It does not run on its own though, needs copy/pasting for now.
######################

######################
# Load a previous R session, data and objects:
#load('R_session_saved_image_order_and_match.RData', verbose=T)
######################

######################
# This function allows other R scripts to obtain the path to a script directory
# (ie where this script lives). Useful when using source('some_script.R')
# without having to pre-specify the location of where that script is.
# This is taken directly from:
# How to source another_file.R from within your R script molgenis/molgenis-pipelines Wiki
# https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
# Couldn't find a licence at the time (12 June 2018)
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
    this.file = NULL
    # This file may be 'sourced' 
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }

    if (!is.null(this.file)) return(dirname(this.file))

    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))

    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}
Rscripts_dir <- LocationOfThisScript()
print('Location where this script lives:')
Rscripts_dir
# R scripts sourced with source() have to be in the same directory as this one
# (or the path constructed appropriately with file.path)
######################

######################
# Import libraries
# source('http://bioconductor.org/biocLite.R')
# biocLite
library(data.table)
source(file.path(Rscripts_dir, 'moveme.R')) #, chdir = TRUE)
######################

######################
##########
# Read files, this is with data.table:
if (is.null(args[['--file1']]) == FALSE) {
  input_name_1 <- as.character(args[['--file1']])
  # For tests:
  # setwd('~/Documents/quickstart_projects/chronic_inflammation_Airwave.p_q/results/QTL_all_files_04Feb2018/')
  # input_name_1 <- 'airwave-illumina_exome-all_chrs.A-transpose.matrixQTL.geno'
  input_data_1 <- fread(input_name_1, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
} else {
  # Stop if arguments not given:
  print('You need to provide an input file. This has to be tab separated with headers.')
  stopifnot(!is.null(args[['--file1']]) == TRUE)
}
print('First file being used: ')
print(input_name_1)
##########

##########
# Read file2:
if (is.null(args[['--file2']]) == FALSE) {
  input_name_2 <- as.character(args[['--file2']])
  # For tests:
  # input_name_2 <- 'AIRWAVE-CPMG_BatchCorrected_log_Var_Data_Sample-plasma.transposed.tsv'
  input_data_2 <- fread(input_name_2, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
} else {
  # Stop if arguments not given:
  print('You need to provide an input file. This has to be tab separated with headers.')
  stopifnot(!is.null(args[['--file2']]) == TRUE)
}
print('Second file file being used (expecting molecular variables): ')
print(input_name_2)
##########

##########
# Set output file name prefix:
# TO DO: sort out as two input_name s
if (is.null(args[['-O']])) {
  stopifnot(!is.null(args[['--file1']]), !is.null(args[['--file2']]))
  print('Output file name prefix not given. Using:')
  output_file_name <- 'matched_'
  print('Output file names will contain: ')
  print(output_file_name)
} else {
  output_file_name <- as.character(args[['-O']])
  # output_file_name <- 'testing'
  output_file_name <- sprintf('%s', output_file_name)
  print(sprintf('Output file names will contain %s', output_file_name))
}
##########
######################

######################
##########
# Check file1 and order
# setkey(input_data_1, )
input_data_1[1:5, 1:5, with = F]
# TO DO for data.table ordering if large files:
# col_order <- order(colnames(input_data_1))
# col_order
# setcolorder(input_data_1, c("SNP", master_IDs_read))
input_data_1 <- as.data.frame(input_data_1[, order(colnames(input_data_1)), with = F])
class(input_data_1)
colnames(input_data_1)[ncol(input_data_1)] <- 'FID'
# names(input_data_1)
input_data_1 <- input_data_1[, moveme(names(input_data_1), 'FID first')]
dim(input_data_1)
input_data_1[1:5, 1:5]
##########

##########
# Check file2 and order
input_data_2[1:5, 1:5, with = F]
input_data_2 <- as.data.frame(input_data_2[, order(colnames(input_data_2)), with = F])
class(input_data_2)
colnames(input_data_2)[ncol(input_data_2)] <- 'FID'
# names(input_data_2)
input_data_2 <- input_data_2[, moveme(names(input_data_2), 'FID first')]
dim(input_data_2)
input_data_2[1:5, 1:5]
##########
######################

######################
# Match each
# file1 and file2:
print('Number of column names present in both file 1 and file 2:')
length(which((colnames(input_data_2) %in% colnames(input_data_1))))
print('Shape of file 1:')
dim(input_data_1)
print('Shape of file 2:')
dim(input_data_2)
# input_data_2_temp <- input_data_2
input_data_1 <- input_data_1[, which(colnames(input_data_1) %in% colnames(input_data_2))]
input_data_2 <- input_data_2[, which(colnames(input_data_2) %in% colnames(input_data_1))]



# Check and stop if not the same:
sanity_check <- identical(as.character(colnames(input_data_2)),
                          as.character(colnames(input_data_1)))
if (sanity_check == TRUE) {
  print('Column names in both files match, continuing.')
} else {
  stop('Exiting. The order of the column names in file 1 and file 2 do not match after 
       ordering and keeping only those which match in both. Are column names the same in
       both files? They must be identical in both files.')
}
input_data_2[1:5, 1:5]
input_data_1[1:5, 1:5]
##########
######################


######################
# Write to file each, this saves with an empty first header and row names as first column, cut when read next:
# col.names = NA makes headers match but can't be used with row.names = F
# Save file:
fwrite(input_data_1,
       sprintf('%s%s', output_file_name, input_name_1),
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)

fwrite(input_data_2,
       sprintf('%s%s', output_file_name, input_name_2),
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)
######################


######################
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])
#rm(list=ls(xxx))
#objects_to_save <- (c('xxx_var'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# Filename to save current R session, data and objects at the end:
if (is.null(args[['--session']]) == FALSE) {
  save_session <- as.character(args[['--session']]) #args $ `--session`
  R_session_saved_image <- sprintf('R_session_saved_image_%s.RData', save_session)
  print(sprintf('Saving an R session image as: %s', R_session_saved_image))
  save.image(file = R_session_saved_image, compress = 'gzip')
} else {
  print('Not saving an R session image, this is the default. Specify the --session option otherwise')
}

# If using Rscript and creating plots, Rscript will create the file Rplots.pdf 
# by default, it doesn't look like there is an easy way to suppress it, so deleting here:
print('Deleting the file Rplots.pdf...')
system('rm -f Rplots.pdf')
print('Finished successfully.')
sessionInfo()
q()

# Next: run the script for xxx
######################
