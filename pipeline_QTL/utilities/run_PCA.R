#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
run_bigPCA.R
===============

Author: |author_names| 
Release: |version|
Date: |today|


Purpose
=======

|description|

Run principal component analysis (PCA) on input data using the bigPCA R package.


From example in manual:
# https://github.com/nicholasjcooper/bigpca

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
# See:
# https://cran.r-project.org/web/packages/docopt/docopt.pdf
# docopt(doc, args = commandArgs(TRUE), name = NULL, help = TRUE,
# version = NULL, strict = FALSE, strip_names = !strict,
# quoted_args = !strict)

# Within the script specify options as:
# args[['--session']]
# args $ `-I` == TRUE

# Save arguments needed:
# B  <- as.numeric(args[['-B']])
# VERBOSE <- as.character(args[['-VERBOSE']])

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
# library(bigpca)
library(ggplot2)
library(cowplot)
# TO DO: sort paths out so they are read from utilities folder after installation:
source('~/Documents/github.dir/EpiCompBio/pipeline_QTL/pipeline_QTL/utilities/moveme.R')
source('~/Documents/github.dir/EpiCompBio/pipeline_QTL/pipeline_QTL/utilities/ggtheme.R')
source('~/Documents/github.dir/EpiCompBio/pipeline_QTL/pipeline_QTL/utilities/ggtheme_bestd.R')
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
  output_file_name <- sprintf('%s.pca', input_name)
  print('Name being used to save output files: ')
  print(output_file_name)
} else {
  output_file_name <- as.character(args[['-O']])
  # output_file_name <- 'testing'
  output_file_name <- sprintf('%s.pca', output_file_name)
  print(sprintf('Output file names will contain %s', output_file_name))
}
##########
######################

######################
# Run PCA on input data

##########
# Explore input data:
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
##########

##########
# Compute the PCs:
input_data_prcomp <- prcomp(input_data, center = TRUE, scale = TRUE)
input_data_prcomp
# Obtain values for all PCs output:
pc <- data.frame(round(input_data_prcomp$x, 4))
pc$sample_id <- rownames(pc) 
pc <- pc[, moveme(names(pc), 'sample_id first')]
names(pc)[1:10]
class(pc)
# Save file:
fwrite(pc,
       sprintf('%s.tsv', output_file_name),
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)
##########

##########
# Explore dimensions and plot first 10 or so components:
dim(pc) # Extra column for sampleID
dim(input_data)
str(input_data_prcomp)
head(pc)

# Get cumulative proportion:
sum_pca <- summary(input_data_prcomp)
sum_pca$importance[, 1:10]
sum_pca
sum_pca_df <- as.data.frame(sum_pca$importance)
sum_pca_df <- t(sum_pca_df)
sum_pca_df <- as.data.frame(sum_pca_df)
# View(sum_pca_df)
sum_pca_df$percent_var <- round(100 * (sum_pca_df$`Proportion of Variance`), 3)
sum_pca_df$PC <- factor(row.names(sum_pca_df), levels = row.names(sum_pca_df),
                        labels = row.names(sum_pca_df))
head(sum_pca_df)
tail(sum_pca_df)
sum_pca_df[1:20, ]
names(sum_pca_df)
str(sum_pca_df)

# Plot proportion of variance of first x PCs:
plot_prop_vars <- ggplot(sum_pca_df[1:100, ], aes(y = percent_var, x = PC)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Principal component', y = 'Proportion of variance (%)') +
  theme_Publication() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid = element_blank(), panel.border = element_blank())
plot_name <- sprintf('prop_var_%s.svg', output_file_name)
ggsave(plot_name, plot_prop_vars)

# Plot first 10 PCs:
head(pc)
# Function for plotting any PC pair:
plot_PCs <- function(data = pc, PCa, PCb){
  xlab <- sprintf('%s', PCa)
  ylab <- sprintf('%s', PCb)
  ggplot(data = data, aes(x = data[, PCa], y = data[, PCb])) + geom_point(alpha = 0.7) +
    theme_Publication() +
    xlab(label = xlab) +
    ylab(label = ylab)
}

for (i in 1:10){
  PCa <- sprintf('PC%s', i)
  PCb <- sprintf('PC%s', i + 1)
  plot_name <- sprintf('plot_%s_%s', PCa, PCb)
  print(plot_name)
  # Name variable on the fly with assign:
  assign(plot_name, plot_PCs(pc, as.character(PCa), as.character(PCb)))
}

# Put all plots together in one figure:
cow_grid <- plot_grid(plot_PC1_PC2,
                      plot_PC2_PC3,
                      plot_PC3_PC4,
                      plot_PC4_PC5,
                      plot_PC5_PC6,
                      plot_PC6_PC7,
                      plot_PC7_PC8,
                      plot_PC8_PC9,
                      plot_PC9_PC10,
                      plot_PC10_PC11,
                      # align = 'vh',
                      # rel_widths = c(1, 0.05, 1),
                      # labels = c("A", "B"),
                      ncol=2)
# Save figure to disk as svg:
plot_name <- sprintf('top_10_PCs_%s.svg', output_file_name)
# A4 paper measures 210 × 297 millimeters or 8.27 × 11.69 inches
svg(plot_name, width = 12, height = 10)
cow_grid
dev.off()
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