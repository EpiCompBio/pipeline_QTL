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

Order and match columns between a genotype and a phenotype file.
Output is meant to be passed to MatrixeQTL.


Usage and options
=================

To run, type:
Rscript order_and_match_QTL.R --geno <FILE> --pheno <FILE> [options]

Usage: order_and_match_QTL.R (--geno <FILE>) (--pheno <FILE>)
       order_and_match_QTL.R [options]

Options:
--geno <FILE>                 Genotype input file name, columns are samples, rows are features
--pheno <FILE>                Phenotype input file name, columns are samples, rows are features
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
######################

######################
##########
# Read files, this is with data.table:
if (is.null(args[['--geno']]) == FALSE) {
  input_name <- as.character(args[['--geno']])
  # For tests:
  # setwd('~/Documents/quickstart_projects/chronic_inflammation_Airwave.p_q/results/QTL_core_illumina/')
  # input_name <- 'all.clean-base.bed'
  input_data_geno <- fread(input_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
} else {
  # Stop if arguments not given:
  print('You need to provide an input file. This has to be tab separated with headers.')
  stopifnot(!is.null(args[['--pcs']]) == TRUE)
}
print('PC file being used: ')
print(input_name)
##########

##########
# pve is a single column with cumulative proportion explained
# the proportion of total variance explained by each of the top k eigenvectors
# To get the cumulative variance explained, simply do the cumulative sum of the variances (cumsum in R)
if (is.null(args[['--pve']]) == FALSE) {
  input_name <- as.character(args[['--pve']])
  # For tests:
  # input_name <- 'pve.all.clean-base.pruned.flashpca.tsv'
  input_data_pve <- fread(input_name, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
} else {
  # Stop if arguments not given:
  print('You need to provide an input file. This has to be tab separated with headers.')
  stopifnot(!is.null(args[['--pve']]) == TRUE)
}
print('pve file being used: ')
print(input_name)
##########

##########
# Set output file name prefix:
# TO DO: sort out as two input_name from pcs and pve
if (is.null(args[['-O']])) {
  stopifnot(!is.null(args[['--pcs']]), !is.null(args[['--pve']]))
  print('Output file name prefix not given. Using:')
  # Split infile name at the last '.':
  input_name <- strsplit(input_name, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
  output_file_name <- sprintf('%s', input_name)
  print('Name being used to save output files: ')
  print(output_file_name)
} else {
  output_file_name <- as.character(args[['-O']])
  # output_file_name <- 'testing'
  output_file_name <- sprintf('%s', output_file_name)
  print(sprintf('Output file names will contain %s', output_file_name))
}
##########
######################


########################
library(data.table)
source('moveme.R')
source('functions_for_MatrixeQTL.R')
########################


########################
# # TO DO for other files: Master file ordered:
# phenotype_file <- 'BEST-D_phenotype_file_final.tsv'
# phenotype_data <- fread(phenotype_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# phenotype_data[1:5, 1:5, with = F]
# master_IDs <- phenotype_data[, c('pt_id', 'kit_id_randomisation', 'kit_id_finalVisit'), with = F]
# master_IDs <- master_IDs[order(kit_id_randomisation)] # Watch the lack of comma as DT not DF
# # master_IDs <- transpose(master_IDs) # Looses headers
# # master_IDs <- as.list(master_IDs)
# master_IDs <- as.data.frame(master_IDs)
# master_IDs <- as.data.frame(t(master_IDs))
# class(master_IDs)
# head(master_IDs)
# row.names(master_IDs)[2] <- 'FID'
# row.names(master_IDs)

#########################


#############################################
# Run with command line arguments:
options(echo=TRUE) # to see commands in output file. TO DO: check how it works with sink() above.
args <- commandArgs(trailingOnly = TRUE)

# TO DO: pass to configuration file

# geno_file <- 'genotype_data_all_treated_baseline.tsv'
# expr_file <- 'GEx_baseline_4000_and_2000.tsv'
# covar_PCs_file <- 'principal_components_normalised_filtered_PC20.tsv'

geno_file <- as.character(args[1])
#geno_file <- 'genotype_data_all_treated_baseline.tsv'
#geno_file <- 'genotype_data_all_treated_final.tsv'
# geno_file <- 'chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome.A-transpose.matrixQTL.geno'

expr_file <- as.character(args[2])
#expr_file <- 'GEx_baseline_4000_and_2000.tsv'
#expr_file <- 'GEx_treated_4000_and_2000.tsv'
# expr_file <- 'Airwave_CPMG_Plasma.txt'

# covar_PCs_file <- as.character(args[3])
#covar_PCs_file <- 'principal_components_normalised_filtered_PC20.tsv'

print(args)
########################

########################
# Read geno baseline and order:
geno_data <- fread(geno_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)#, drop = 1) # Column doesn't need to be dropped, may have been from double IDs before
# setkey(geno_data, )
geno_data[1:5, 1:5, with = F]
# TO DO for data.table ordering if large files:
# col_order <- order(colnames(geno_data))
# col_order
# setcolorder(geno_data, c("SNP", master_IDs_read))
geno_data <- as.data.frame(geno_data[, order(colnames(geno_data)), with = F])
class(geno_data)
colnames(geno_data)[ncol(geno_data)] <- 'FID'
names(geno_data)
geno_data <- geno_data[, moveme(names(geno_data), 'FID first')]
head(geno_data)
dim(geno_data)
geno_data[1:5, 1:5]
########################

########################
# Read expr baseline and order:
expr_data <- fread(expr_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
expr_data[1:5, 1:5, with = F]
# expr_data <- as.data.frame(expr_data[, order(colnames(expr_data)), with = F])
# class(expr_data)
# colnames(expr_data)[ncol(expr_data)] <- 'FID'
# names(expr_data)
# expr_data <- expr_data[, moveme(names(expr_data), 'FID first')]
# head(expr_data)
# dim(expr_data)
# expr_data[1:5, 1:5]

# CPMG needs transposing:
expr_data_t <- transpose_file(expr_data, 1)
class(expr_data_t)
expr_data_t[1:5, 1:5, with = F]
expr_data_t[1:5, (ncol(expr_data_t)-5):ncol(expr_data_t), with = F]
expr_data <- expr_data_t
# Rename datapoints column and re-order:
expr_data <- as.data.frame(expr_data[, order(colnames(expr_data)), with = F])
class(expr_data)
colnames(expr_data)[ncol(expr_data)] <- 'FID'
names(expr_data)
expr_data <- expr_data[, moveme(names(expr_data), 'FID first')]
head(expr_data)
dim(expr_data)
expr_data[1:5, 1:5]
########################

########################
# Names don't match between genotypes and metabolomics data (geno data has plink double IDs):
geno_data[1:5, 1:5]
sample_IDs <- colnames(geno_data)[-1]
head(sample_IDs)
# Split plink double IDs:
sample_IDs_split <- strsplit(sample_IDs, split = '[_]')
# Create data frame:
sample_IDs_split_df <- data.frame(matrix(unlist(sample_IDs_split), 
                                         nrow = length(sample_IDs_split), 
                                         byrow = T),
                                  stringsAsFactors=FALSE)
# Add 'FID' back in:
sample_IDs_split_df <- rbind(c('FID', 'FID'), sample_IDs_split_df)
head(sample_IDs_split_df)
dim(sample_IDs_split_df)
str(sample_IDs_split_df)

# Insert single IDs into geno data:
head(colnames(geno_data))
head(sample_IDs_split_df$X1)
colnames(geno_data) <- sample_IDs_split_df$X1
head(colnames(geno_data))
########################

# ########################
# # Read PCs from expression data:
# covar_PCs <- fread(covar_PCs_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# covar_PCs <- data.table::transpose(covar_PCs) # Tranpose looses headers but these were in order already
# # TO DO: transpose() gives odd error (not in namespace...), runs OK outside of Rscript but don't know where the error is.
# # See:
# # http://stackoverflow.com/questions/23252231/r-data-table-breaks-in-exported-functions
# # http://stackoverflow.com/questions/29549690/function-on-data-table-environment-errors
# # http://stackoverflow.com/questions/15223367/wrapping-data-table-using-an-evaluated-call-in-a-package
# 
# covar_PCs <- as.data.frame(covar_PCs)
# names(covar_PCs) <- covar_PCs[1,]
# covar_PCs <- covar_PCs[-1,]
# covar_PCs[1:5, 1:5]
# class(covar_PCs)
# covar_PCs <- covar_PCs[, order(colnames(covar_PCs))]
# covar_PCs[1:5, 1:5]
# dim(covar_PCs)
# ########################

########################
# Match each
# Geno and expr:
length(which((colnames(expr_data) %in% colnames(geno_data))))
expr_data <- expr_data[, which(colnames(expr_data) %in% colnames(geno_data))]
geno_data <- geno_data[, which(colnames(geno_data) %in% colnames(expr_data))]
identical(colnames(expr_data), colnames(geno_data))
expr_data[1:5, 1:5]
geno_data[1:5, 1:5]

# Covar PCs to expr:
length(which((colnames(covar_PCs) %in% colnames(expr_data)[-1])))
covar_PCs <- covar_PCs[, which(colnames(covar_PCs) %in% colnames(expr_data)[-1])]

# TO DO: stop if not 'TRUE':
identical(colnames(covar_PCs), colnames(expr_data)[-1])
expr_data[1:5, 1:5]
covar_PCs[1:5, 1:5]
########################

########################
# # TO DO:
# # Match gene expression probes measured to location file:
# snp_pos_data <- fread(snp_pos_file, sep = ' ', header = TRUE, stringsAsFactors = FALSE)
# head(snp_pos_data)
# dim(snp_pos_data)
# probe_pos_data <- fread(probe_pos_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# head(probe_pos_data)
# dim(probe_pos_data)
# 
# # Switch classes as these are read by data.table and cause problems later:
# snp_pos_data <- as.data.frame(snp_pos_data)
# probePos <- as.data.frame(probe_pos_data)
# 
# length(which(as.character(probe_pos_data[, 1]) %in% as.character(expr_data[, 1])))
# probe_pos_data <- probe_pos_data[which(as.character(probe_pos_data[, 1]) %in% as.character(expr_data[, 1])), ]
# dim(probe_pos_data)
# head(probe_pos_data)
# 
# length(which(as.character(snp_pos_data[, 1]) %in% as.character(geno_data[, 1])))
# snp_pos_data <- snp_pos_data[which(as.character(snp_pos_data[, 1]) %in% as.character(geno_data[, 1])), ]
# dim(snp_pos_data)
# head(snp_pos_data)
########################

########################
# Write to file each, this saves with an empty first header and row names as first column, cut when read next:
# col.names = NA makes headers match but can't be used with row.names = F
write.table(expr_data, paste(expr_file, '_matched.tsv', sep = ''), sep='\t', quote = FALSE, col.names = NA)
write.table(geno_data, paste(geno_file, '_matched.tsv', sep = ''), sep='\t', quote = FALSE, col.names = NA)

covar_PCs_file_substr <- substring(as.character(covar_PCs_file), 1, 21)
expr_file_substr <- substring(as.character(expr_file), 1, 12)
write.table(covar_PCs, paste(covar_PCs_file_substr, expr_file_substr, '_matched.tsv', sep = ''), 
            sep='\t', quote = FALSE, col.names = NA)

# Cut first column (row numbers):
expr_written <- paste(expr_file, '_matched.tsv', sep = '')
cmd_cut <- sprintf('cat %s | cut -f2- > cut_%s', expr_written, expr_written)
system(as.character(cmd_cut))
geno_written <- paste(geno_file, '_matched.tsv', sep = '')
cmd_cut <- sprintf('cat %s | cut -f2- > cut_%s', geno_written, geno_written)
system(as.character(cmd_cut))
covar_written <- paste(covar_PCs_file_substr, expr_file_substr, '_matched.tsv', sep = '')
cmd_cut <- sprintf('cat %s | cut -f2- > cut_%s', covar_written, covar_written)
system(as.character(cmd_cut))
########################



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