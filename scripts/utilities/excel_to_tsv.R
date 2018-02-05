#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
excel_to_tsv.R
================

Author: |author_names| 
Release: |version|
Date: |today|


Purpose
=======

|description|

Convert the first spreadsheet of an Excel file into a tab separated plain text file.

Usage and options
=================

To run, type:
Rscript excel_to_tsv.R -I <INPUT_FILE> [options]

Usage: excel_to_tsv.R (-I <INPUT_FILE>)
       excel_to_tsv.R [options]

Options:
  -I <INPUT_FILE>                 Input file name in xlsx format
  -O <OUTPUT_FILE>                Output file name
  --session <R_SESSION_NAME>      R session name if to be saved
  -h --help                       Show this screen
  --spreadsheet <integer>         Number of spreadsheet to convert [default: 1]
  --header <boolean>              Specify whether file has column names [default: TRUE]

Input:
 An Excel file that will be read with the R gdata package.
  Specify sheet number and presence of column names if needed.

Output:
  Tab separated file.

Requirements:

library(docopt)
library(data.table)
library(gdata)

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
sheet <- as.integer(args[['--spreadsheet']])
header <- as.logical(args[['--header']])

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
# Import libraries
# source('http://bioconductor.org/biocLite.R')
# biocLite
library(data.table)
library(gdata)
######################

######################
##########
# Read files, this is with data.table:
if (is.null(args[['-I']]) == FALSE) {
  input_name <- as.character(args[['-I']])
  # For tests:
  # setwd('~/Documents/quickstart_projects/chronic_inflammation_Airwave.p_q/data/processed/metabolomics_data.dir/')
  # input_name <- 'AIRWAVE_1DNMR_BatchCorrected_log_VarInfo.xlsx'
  input_data <- read.xls(input_name, sheet = sheet, header = header)
} else {
  # Stop if arguments not given:
  print('You need to provide an input file. This has to be from Excel (xlsx suffix).')
  stopifnot(!is.null(args[['-I']]) == TRUE)
}
print('File being used: ')
print(input_name)
##########

##########
# Set output file name:
if (is.null(args[['-O']])) {
  stopifnot(!is.null(args[['-I']]))
  print('Output file name prefix not given. Using:')
  # Split infile name at the last '.':
  input_name <- strsplit(input_name, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
  output_file_name <- sprintf('%s.tsv', input_name)
  print('Name being used to save output files: ')
  print(output_file_name)
} else {
  output_file_name <- as.character(args[['-O']])
  # output_file_name <- 'testing'
  #output_file_name <- sprintf('%s.tsv', output_file_name)
  print(sprintf('Output file names will contain %s', output_file_name))
}
##########
######################

######################
##########
# Explore input data:
class(input_data)
dim(input_data) # nrow(), ncol()
str(input_data)
head(input_data)
colnames(input_data)

# Save file:
fwrite(input_data,
       output_file_name,
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)
##########
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
