#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
merge_dataframes.R
====================

Author: |author_names| 
Release: |version|
Date: |today|


Purpose
=======

|description|

Merge two data frames. This is specifically for principal component files derived from
genotype and phenotype data. File 2 will get a label added to make headers unique.

Usage and options
=================

To run, type:
Rscript order_and_match_QTL.R --file1 <FILE> --file2 <FILE> [options]

Usage: merge_dataframes.R (--file1 <FILE>) (--file2 <FILE>) [options]
       merge_dataframes.R [options]

Options:
--file1 <FILE>                Usually a genotype input file name, columns are samples, rows are features
--file2 <FILE>                usually a phenotype input file name, columns are samples, rows are features
--file1-PCs                   Number of PCs from file 1 to keep. [default: 10]
--file2-PCs                   Number of PCs from file 1 to keep. [default: 30]
-O <OUTPUT_FILE>              Output file name
--session <R_SESSION_NAME>    R session name if to be saved
-h --help                     Show this screen

Input:

Two tab separated files. These are read with data.table and stringsAsFactors = FALSE
Rows must be features (phenotypes, variables, etc.) and columns must be samples (individuals)
The first row and first column must be the ID labels.

Output:

Merges two dataframes and writes output to disk.

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
file1_PCs <- as.integer(args[['--file1-PCs']])
file2_PCs <- as.integer(args[['--file2-PCs']])

# Print to screen:
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
# Import libraries
# source('http://bioconductor.org/biocLite.R')
# biocLite
library(data.table)
# TO DO: sort paths out so they are read from utilities folder after installation:
source('moveme.R')
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
# TO DO: sort out as two input_names:
if (is.null(args[['-O']])) {
  stopifnot(!is.null(args[['--file1']]), !is.null(args[['--file2']]))
  print('Output file name prefix not given. Using:')
  output_file_name <- sprintf('merged_%s_%s.tsv', input_name_1, input_name_2)
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
# Check file dimensions and merge
print('Number of column names present in both file 1 and file 2:')
length(which((colnames(input_data_2) %in% colnames(input_data_1))))
print('Shape of file 1:')
dim(input_data_1)
print('Shape of file 2:')
dim(input_data_2)

# Add a string to have unique names:
# input_data_2$string <- '_pheno'
# input_data_2$ID <- with(input_data_2, paste0(FID, string))
# input_data_2 <- input_data_2[, moveme(names(input_data_2), 'ID first')]
# input_data_2[, 'FID'] <- NULL
# input_data_2[, 'string'] <- NULL
# colnames(input_data_2)[1] <- 'FID'
# input_data_2[1:5, 1:5]

# Merge into one file, choice of PCs is arbitrary:
all_covar_PCs <- rbind(input_data_1[1:file1_PCs, ], input_data_2[1:file2_PCs, ])
fwrite(all_covar_PCs,
       sprintf('%s', output_file_name),
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
