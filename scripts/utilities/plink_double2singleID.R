#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
plink_double2singleID.R
========================

Author: |author_names| 
Release: |version|
Date: |today|


Purpose
=======

|description|

Plink outputs double IDs (FID_IID) after using e.g. --recode A-transpose, which generates a .traw file.
Double IDs will then fail to match against phenotype and other files.
This script does a simple conversion, keeping everything after "_", i.e. the FID.
This works if FID and IIDs are identical, otherwise you may get repeated IDs.


Usage and options
=================

To run, type:
Rscript plink_double2singleID.R -I <FILE> [options]

Usage: plink_double2singleID.R (-I <FILE>)
       plink_double2singleID.R [options]

Options:
-I <FILE>                     A file with double IDs (FID_IID) as column names
-O <OUTPUT_FILE>              Output file name
--session <R_SESSION_NAME>    R session name if to be saved
-h --help                     Show this screen

Input:

A tab separate file with double IDs in column names.
Files are read with data.table and stringsAsFactors = FALSE

Output:

Outputs a file with single ID column names (i.e. keeps the IID from "FID_IID")

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
# Import libraries
# source('http://bioconductor.org/biocLite.R')
# biocLite
library(data.table)
######################

######################
##########
# Read files, this is with data.table:
if (is.null(args[['-I']]) == FALSE) {
  input_name_1 <- as.character(args[['-I']])
  # For tests:
  # setwd('~/Documents/quickstart_projects/chronic_inflammation_Airwave.p_q/results/QTL_core_illumina/')
  # input_name_1 <- 'all.clean-base.A-transpose.matrixQTL.geno'
  input_data_1 <- fread(input_name_1, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
} else {
  # Stop if arguments not given:
  print('You need to provide an input file. This has to be tab separated with headers as double IDs.')
  stopifnot(!is.null(args[['-I']]) == TRUE)
}
print('Input file is: ')
print(input_name_1)
##########

##########
# Set output file name prefix:
# TO DO: sort out as two input_name s
if (is.null(args[['-O']])) {
  stopifnot(!is.null(args[['-I']]))
  print('Output file name prefix not given. Using:')
  output_file_name <- 'IID'
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
# Check input file:
# setkey(input_data_1, )
input_data_1[1:5, 1:5, with = F]

# Names don't match between genotypes and metabolomics data (geno data has plink double IDs):
sample_IDs <- colnames(input_data_1)[-1]
head(sample_IDs)
# # TO DO: check QC:
# # Some have IDs which get split as "11421" "rep"   "11421" "rep":
# grep_rep <- grep('rep', sample_IDs)
# grep_rep
# sample_IDs[c(grep_rep[1])] <- "11421_11421"
# sample_IDs[c(grep_rep[2])] <- "46750_46750"
# Split plink double IDs:
sample_IDs_split <- strsplit(sample_IDs, split = '[_]')
# Create data frame:
sample_IDs_split_df <- data.frame(matrix(unlist(sample_IDs_split), 
                                         nrow = length(sample_IDs_split), 
                                         byrow = T),
                                  stringsAsFactors = FALSE)
# Add 'FID' back in:
sample_IDs_split_df <- rbind(c('FID', 'FID'), sample_IDs_split_df)
# head(sample_IDs_split_df)
dim(sample_IDs_split_df)
# str(sample_IDs_split_df)

# Insert single IDs into file1 data:
head(colnames(input_data_1))
head(sample_IDs_split_df$X1)
# input_data_1_temp <- input_data_1
colnames(input_data_1) <- sample_IDs_split_df$X1
head(colnames(input_data_1))
input_data_1[1:5, 1:5]
##########
######################


######################
# Write to file each, this saves with an empty first header and row names as first column, cut when read next:
# col.names = NA makes headers match but can't be used with row.names = F
# Save file:
# TO DO: sort out naming of files:
fwrite(input_data_1,
       sprintf('%s_%s', output_file_name, input_name_1),
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
