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
# Import libraries
# source('http://bioconductor.org/biocLite.R')
# biocLite
library(data.table)
# TO DO: sort paths out so they are read from utilities folder after installation:
source('~/Documents/github.dir/EpiCompBio/pipeline_QTL/pipeline_QTL/utilities/moveme.R')
# source('~/Documents/github.dir/EpiCompBio/pipeline_QTL/pipeline_QTL/utilities/functions_for_MatrixeQTL.R')
######################

######################
##########
# Read files, this is with data.table:
if (is.null(args[['--file1']]) == FALSE) {
  input_name_1 <- as.character(args[['--file1']])
  # For tests:
  # setwd('~/Documents/quickstart_projects/chronic_inflammation_Airwave.p_q/results/QTL_core_illumina/')
  # input_name_1 <- 'all.clean-base.A-transpose.matrixQTL.geno'
  input_data_1 <- fread(input_name_1, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
} else {
  # Stop if arguments not given:
  print('You need to provide an input file. This has to be tab separated with headers.')
  stopifnot(!is.null(args[['--file1']]) == TRUE)
}
print('First file being used (expecting genotypes here): ')
print(input_name_1)
##########

##########
# Read file2:
if (is.null(args[['--file2']]) == FALSE) {
  input_name_2 <- as.character(args[['--file2']])
  # For tests:
  # input_name_2 <- 'AIRWAVE_1DNMR_BatchCorrected_log_Data_Var_Sample.transposed.tsv'
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
  output_file_name <- 'matched'
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
# TO DO: this is specific to set of files:
# Names don't match between genotypes and metabolomics data (geno data has plink double IDs):
input_data_1[1:5, 1:5]
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
length(which((colnames(input_data_2) %in% colnames(input_data_1))))
# input_data_2_temp <- input_data_2
input_data_2 <- input_data_2[, which(colnames(input_data_2) %in% colnames(input_data_1))]
input_data_1 <- input_data_1[, which(colnames(input_data_1) %in% colnames(input_data_2))]
# # Delete 'rep' samples:
# to_del <- grep('46750.1', colnames(input_data_1))
# input_data_1 <- input_data_1[, -to_del]
stopifnot(identical(as.character(colnames(input_data_2)),
          as.character(colnames(input_data_1))))
input_data_2[1:5, 1:5]
input_data_1[1:5, 1:5]
##########

##########
# TO DO: separate or move this file to project specific
# Covar PCs to file2:
covar_PCs_1_name <- 'pcs.all.clean-base.pruned.flashpca.transposed.tsv'
covar_PCs_1 <- fread(covar_PCs_1, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# This file after using transposing script has double IDs in rows:
covar_PCs_1 <- covar_PCs_1[-1, ]
# Change 'V1' to 'FID':
setnames(x = covar_PCs_1, old = 'V1', new = 'FID')
covar_PCs_1 <- as.data.frame(covar_PCs_1[, order(colnames(covar_PCs_1)), with = F])
covar_PCs_1 <- covar_PCs_1[, moveme(names(covar_PCs_1), 'FID first')]
covar_PCs_1[1:5, 1:5]


covar_PCs_2_name <- 'AIRWAVE_1DNMR_BatchCorrected_log_Data_Var_Sample.pca.transposed.tsv'
covar_PCs_2 <- fread(covar_PCs_2, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# Change 'V1' to 'FID':
setnames(x = covar_PCs_2, old = 'V1', new = 'FID')
covar_PCs_2 <- as.data.frame(covar_PCs_2[, order(colnames(covar_PCs_2)), with = F])
covar_PCs_2 <- covar_PCs_2[, moveme(names(covar_PCs_2), 'FID first')]
covar_PCs_2[1:5, 1:5]
dim(covar_PCs_2)

# Match PCs: 
# length(which(as.character(unlist(covar_PCs_1[, 'FID'])) %in% as.character(colnames(input_data_1))))
length(which(colnames(covar_PCs_1) %in% colnames(input_data_1)))
covar_PCs_1 <- covar_PCs_1[, which(colnames(covar_PCs_1) %in% colnames(input_data_1))]
covar_PCs_1[1:5, 1:5]
dim(covar_PCs_1)

length(which(colnames(covar_PCs_2) %in% colnames(input_data_2)))
covar_PCs_2 <- covar_PCs_2[, which(colnames(covar_PCs_2) %in% colnames(input_data_2))]
covar_PCs_2[1:5, 1:5]
dim(covar_PCs_2)

# Stop if not 'TRUE':
stopifnot(identical(colnames(covar_PCs_1), colnames(covar_PCs_2)))
stopifnot(identical(colnames(covar_PCs_1), colnames(input_data_2)))
stopifnot(identical(colnames(covar_PCs_2), colnames(input_data_1)))

# TO DO: this is analysis specific, clean up:
covar_PCs_1[1:5, 1:5]
covar_PCs_2[30:35, 1:5]
# Add a string to have unique names:
covar_PCs_2$string <- '_pheno'
covar_PCs_2$ID <- with(covar_PCs_2, paste0(FID, string))
covar_PCs_2 <- covar_PCs_2[, moveme(names(covar_PCs_2), 'ID first')]
covar_PCs_2[, 'FID'] <- NULL
covar_PCs_2[, 'string'] <- NULL
colnames(covar_PCs_2)[1] <- 'FID'
covar_PCs_2[1:5, 1:5]

# Merge into one file, choice of PCs is arbitrary:
all_covar_PCs <- rbind(covar_PCs_1, covar_PCs_2[1:35, ])
fwrite(all_covar_PCs,
       sprintf('%s_all_covar_PCs.tsv', output_file_name),
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)
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

fwrite(input_data_2,
       sprintf('%s_%s', output_file_name, input_name_2),
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)

fwrite(covar_PCs_1,
       sprintf('%s_%s', output_file_name, covar_PCs_1_name),
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)

fwrite(covar_PCs_2,
       sprintf('%s_%s', output_file_name, covar_PCs_2_name),
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)
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