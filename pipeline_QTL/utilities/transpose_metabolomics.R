#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
transpose_metabolomics.R
==========================

Author: |author_names| 
Release: |version|
Date: |today|


Purpose
=======

|description|

Transposes metabolomics file. 

Usage and options
=================

To run, type:
Rscript transpose_metabolomics.R -I <INPUT_FILE> [options]

Usage: transpose_metabolomics.R (-I <INPUT_FILE>)
       transpose_metabolomics.R [options]

Options:
  -I <INPUT_FILE>                 Input file name, columns are samples, rows are features
  -O <OUTPUT_FILE>                Output file name
  --session <R_SESSION_NAME>      R session name if to be saved
  -h --help                       Show this screen

Input:

Tab separated file

Output:

A tranposed tab separated file.

Requirements:

library(docopt)
library(data.table)
library(ggplot2)
library(cowplot)

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
# Import libraries
# source('http://bioconductor.org/biocLite.R')
# biocLite
library(data.table)
# TO DO: sort paths out so they are read from utilities folder after installation:
# source('~/Documents/github.dir/EpiCompBio/pipeline_QTL/pipeline_QTL/utilities/moveme.R')
######################

######################
##########
# Read files, this is with data.table:
if (is.null(args[['-I']]) == FALSE) {
  input_name <- as.character(args[['-I']])
  # For tests:
  # setwd('~/Documents/quickstart_projects/chronic_inflammation_Airwave.p_q/results/QTL_core_illumina')
  # input_name <- 'AIRWAVE_1DNMR_BatchCorrected_log_Data_Var_Sample.tsv_test'
  input_data <- fread(input_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
} else {
  # Stop if arguments not given:
  print('You need to provide an input file. This has to be tab separated with headers.')
  stopifnot(!is.null(args[['-I']]) == TRUE)
}
print('File being used: ')
print(input_name)
##########

##########
# Set output file names:
if (is.null(args[['-O']])) {
  stopifnot(!is.null(args[['-I']]))
  print('Output file name prefix not given. Using:')
  # Split infile name at the last '.':
  input_name <- strsplit(input_name, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
  output_file_name <- sprintf('%s.transposed.tsv', input_name)
  print('Name being used to save output files: ')
  print(output_file_name)
} else {
  output_file_name <- as.character(args[['-O']])
  # output_file_name <- 'testing'
  #output_file_name <- sprintf('%s.transposed.tsv', output_file_name)
  print(sprintf('Output file names will contain %s', output_file_name))
}
##########
######################

######################
##########
# Explore input data:
class(input_data)
dim(input_data) # nrow(), ncol()
# str(input_data)
# input_data # equivalent to head() and tail()
# setkey(input_data) # memory efficient and fast
# key(input_data)
# tables()
# colnames(input_data)
input_data[1:5, 1:5]

# Save original IDs from first column:
rows <- as.character(unlist(input_data[, 1]))
# Save original IDs from first row:
cols <- as.character(colnames(input_data))
# Transpose file without first column (containing IDs):
input_data_t <- transpose(input_data[, -1])
input_data_t[1:5, 1:5]
dim(input_data)
dim(input_data_t)
# Insert original IDs as new colnames in transposed:
colnames(input_data_t) <- rows
# Insert original IDs as first column into transposed:
input_data_t <- cbind(as.character(cols[-1]), input_data_t) # Exclude first label
input_data_t[1:5, 1:5]
dim(input_data_t)

# Sanity:
input_data[1:5, 1:5]
input_data_t[1:5, 1:5]
# Second column should be the first row:
stopifnot(identical(as.character(unlist(input_data[, 2])),
                    as.character(unlist(input_data_t[1, -1])) # Minus first column which is ID
                    )
          )
# print('There was a problem, the file does not match after being transposed. Are there IDs in the first row and column?')

# Random sample range:
random_n <- sample(ncol(input_data), 1)
# random_n
stopifnot(identical(as.character(unlist(input_data[, c(random_n + 1), with = F])), # Plus one for offset in IDs
                    as.character(unlist(input_data_t[c(random_n), -1]))
                    )
          )
# print('There was a problem, the file does not match after being transposed. Are there IDs in the first row and column?')  

# Save file:
fwrite(input_data_t,
       output_file_name,
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)
##########
######################

######################
## Save some text:
# Methods
# Legend
# Interpretation
# cat(file <- output_file, some_var, '\t', another_var, '\n', append = TRUE)
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
