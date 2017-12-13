#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
run_SVA.R
===============

Author: |author_names| 
Release: |version|
Date: |today|


Purpose
=======

|description|

Run surrogate variable analysis (SVA) using the SmartSVA R package.


From example in manual:
# https://cran.r-project.org/web/packages/SmartSVA/SmartSVA.pdf
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3808-1

Also see some tutorials and recent papers:
# http://genomicsclass.github.io/book/pages/pca_svd.html
# https://www.biorxiv.org/content/early/2017/03/26/120899

Also see CMS paper and tool:
# https://github.com/haschard/CMS/blob/master/CMS_v1.0.R
# https://www.nature.com/articles/ng.3975.pdf


Usage and options
=================

To run, type:
    Rscript run_SVAR.R -I <INPUT_FILE> [options]

Usage: run_SVAR.R (-I <INPUT_FILE>)
       run_SVAR.R [options]
       run_SVAR.R [-h | --help]

Options:
  -I <INPUT_FILE>                 Input file name, columns are samples, rows are features
  -O <OUTPUT_FILE>                Output file name
  --session <R_SESSION_NAME>      R session name if to be saved
  -h --help                       Show this screen
  -mod0	<null_model_file>         File name with null model matrix
  -n.sv <integer>                 Number of surrogate variables to estimate [default: 50].
  -B <numeric>                    The maximum iteration number [default: 100].
  -alpha <numeric>                Initial point for optimization for convergence rate [default: 0.25].
  -VERBOSE <boolean>              SmartSVA prints additional details if TRUE [default: TRUE].
  -epsilon <numeric>              Convergence threshold [default: 0.001].

Input:

    A tab separated file with headers. This is read with data.table and stringsAsFactors = FALSE
    Rows must be features (phenotypes) and columns must be samples (individuals)
    The first row and first column must be labels.
    Input file is converted to a matrix, scaled (scale()) and passed to smartsva.cpp()
    If providing a null model, give the name of a tab separated file with headers.

    The model matrix used to fit the data is calculated as
    model.matrix(~1+X, data = data.frame(input_data))


Output:

SmartSVA returns a list containing the surrogate variables and some meta data about the convergence criterion.
Save the session if you need the SV object. This wrapper saves a file with the significant surrogate variables
as a tab separated file with headers.

Requirements:

    library(docopt)
    library(data.table)
    library(SmartSVA)

Documentation
=============

    For more information see:

    |url|
' -> doc
# Load docopt:
library(docopt, quietly = TRUE)
# Retrieve the command-line arguments:
args <- docopt(doc)
# See:
# https://cran.r-project.org/web/packages/docopt/docopt.pdf
# docopt(doc, args = commandArgs(TRUE), name = NULL, help = TRUE,
# version = NULL, strict = FALSE, strip_names = !strict,
# quoted_args = !strict)

# Within the script specify options as:
# args[['--session']]
# args $ `-I` == TRUE

# Save arguments needed to run SmartSVA:
# mod0 <- as.character(args[['-mod0']])
B  <- as.numeric(args[['-B']])
alpha <- as.numeric(args[['-alpha']])
VERBOSE <- as.character(args[['-VERBOSE']])
n.sv <- as.numeric(args[['-n.sv']])
epsilon <- as.numeric(args[['-epsilon']])

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
# biocLite()
library(data.table)
library(SmartSVA)
######################

######################
##########
# Read files, this is with data.table:
if (is.null(args[['-I']]) == FALSE) {
  input_name <- as.character(args[['-I']])
  # For tests:
  # setwd('~/Desktop/Downloads_to_delete/miscellaneous_tests/pipe_QTL_tests/results/tests3')
  # input_name <- 'airwave-NMR-blood.pheno'
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
  output_file_name <- sprintf('%s.sva', input_name)
  print('Name being used to save output files: ')
  print(output_file_name)
} else {
  output_file_name <- as.character(args[['-O']])
  # output_file_name <- 'testing'
  output_file_name <- sprintf('%s.sva', output_file_name)
  print(sprintf('Output file names will contain %s', output_file_name))
}
##########
######################

######################
# Run SVA on input data
# Explore data:
class(input_data)
dim(input_data) # nrow(), ncol()
str(input_data)
input_data # equivalent to head() and tail()
setkey(input_data) # memory efficient and fast
key(input_data)
tables()
colnames(input_data)
# First column with feature labels:
input_data[, 1, with = FALSE] # by column position, preferable by column name to avoid silent bugs

# Convert to matrix, exclude the feature IDs
# Save feature names as rownames:
input_data <- as.data.frame(input_data)
rownames(input_data) <- input_data[, 1]
rownames(input_data)
input_data
input_data <- data.matrix(input_data[1:nrow(input_data), 2:ncol(input_data)])
class(input_data)
str(input_data)
dim(input_data)
head(input_data)
input_data <- scale(input_data)
str(input_data)
class(input_data)
head(input_data) # Columns must be individuals and rows phenotypes

# TO DO: From the manual:
# Y <- matrix(rnorm(100*2700), 2700, 100)
# input_data2 <- matrix(rnorm(100*50), 100, 50)
# df <- data.frame(pred=gl(2, ncol(input_data) / 2))
# Determine the number of SVs
# Y.r <- t(resid(lm(t(Y) ~ pred, data=df)))
# mod <- model.matrix( ~ pred, df)
# For testing:
# dim(input_data)
# df_cov <- data.frame(pred = gl(2, ncol(input_data) / 2))
# df_cov
# mod  <- model.matrix(~1 + pred, data = df_cov)
mod  <- model.matrix(~1, data = data.frame(input_data))
mod
# What I want is like PC loop where PCs are obtained per sample for each sample present.
# From Andy's example in VD:
# X			<- cbind( time, time*vit_D )
# mod0  <- model.matrix(~1, data=data.frame(expr_sv) )
# mod   <- model.matrix(~1+X, data=data.frame(expr_sv) )

# Provide a null model matrix:
# if (is.null(args[['-mod0']]) == FALSE) {
if (is.null(args[['-mod0']])) {
  # Leave mod0 as NULL, this is the default:
  mod0 <- NULL # Passed to variable as.character above, clean up
  print('Leaving mod0 as null.')
} else {
  # Reading null model matrix provided:
  mod0 <- fread(mod0, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
}

# Run SmartSVA:
svs  <- smartsva.cpp(dat = input_data,
                     n.sv = n.sv,
                     mod = mod,
                     mod0 = mod0,
                     B = B,
                     alpha = alpha,
                     VERBOSE = VERBOSE,
                     epsilon = epsilon)

# Explore the file for sanity (will be saved to log file):
str(svs)
class(svs)
dim(svs)
head(svs$sv)
tail(svs$sv)
# Save surrogate variables:
svs_sv <- svs$sv
######################

###################### 
# Save as a file, give a name if -O option not provided:
if (is.null(args[['-O']]) == FALSE) {
  output_name <- as.character(args[['-O']])
  input_data <- fread(input_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
} else {
  # Warn if arguments not given:
  outfile <- sprintf('desc_stats_%s.tsv', input_name)
  warning(sprintf('Outfile name not given, using a default: %s', outfile))
}

# Save file:
fwrite(desc_stats, outfile, 
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)
######################

######################
# Some basic exploratory plots
# Plot here or use a separate plot_template.R script (preferable if processing
# large datasets, process first, save, plot separately):
plot_name <- svg(sprintf('%s_%s_%s.svg', input_name, x_var_name, y_var_name))
# par(mfrow = c(1, 3)) # rows, cols
boxplot(input_data_df[,  y_var_name] ~ input_data_df[,  x_var_name])
dev.off()
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
