#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
run_PCA.R
===============

Author: |author_names| 
Release: |version|
Date: |today|


Purpose
=======

|description|

Run principal component analysis (PCA) on input data using prcomp in R.

If this is too slow look at the bigPCA package:
# https://github.com/nicholasjcooper/bigpca

Usage and options
=================

To run, type:
Rscript run_PCA.R -I <INPUT_FILE> [options]

Usage: run_PCA.R (-I <INPUT_FILE>)
       run_PCA.R [options]

Options:
  -I <INPUT_FILE>                 Input file name, columns are samples, rows are features
  -O <OUTPUT_FILE>                Output file name
  --session <R_SESSION_NAME>      R session name if to be saved
  -h --help                       Show this screen
  --num_PCs	<PCs_to_plot>         Number of PCs to plot, max 10 [default: 10]

Input:

A tab separated file. This is read with data.table and stringsAsFactors = FALSE
Rows must be features (phenotypes) and columns must be samples (individuals)
The first row and first column must be the ID labels.
PCA is run using prcomp, data is scaled and centred.

Output:

A tab separated file with the principal components and a plot of the top 10 PCs and variance explained.
Save the session if you want all the details from prcomp.

Requirements:

library(docopt)
library(data.table)
library(ggplot2)
library(cowplot)
library(mice)

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
num_PCs <- as.integer(args[['--num_PCs']])

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
# library(mice)
# TO DO: sort paths out so they are read from utilities folder after installation:
#source('~/Documents/github.dir/EpiCompBio/pipeline_QTL/pipeline_QTL/utilities/moveme.R')
#source('~/Documents/github.dir/EpiCompBio/pipeline_QTL/pipeline_QTL/utilities/ggtheme.R')
#source('~/Documents/github.dir/EpiCompBio/pipeline_QTL/pipeline_QTL/utilities/ggtheme_bestd.R')
source('moveme.R')
source('ggtheme.R')
source('ggtheme_bestd.R')
######################

######################
##########
# Read files, this is with data.table:
if (is.null(args[['-I']]) == FALSE) {
  input_name <- as.character(args[['-I']])
  # For tests:
  # setwd('~/Documents/quickstart_projects/chronic_inflammation_Airwave.p_q/results/QTL_core_illumina')
  # input_name <- 'AIRWAVE_1DNMR_BatchCorrected_log_Data_Var_Sample.tsv'
  input_data <- fread(input_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE,
                      na.strings = c('', ' ', 'NA', 'NaN'))
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
  #output_file_name <- sprintf('%s.pca', output_file_name)
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
# str(input_data)
# input_data # equivalent to head() and tail()
# setkey(input_data) # memory efficient and fast
# key(input_data)
# tables()
# colnames(input_data)
# First column with feature labels:
# input_data[, 1, with = FALSE] # by column position, preferable by column name to avoid silent bugs
##########

##########
# Check whether there are missing values:
NAs <- which(complete.cases(input_data) == FALSE)
# input_data[c(NAs), ]
length(NAs)

if(length(NAs) == 0) {
  print('Dataset set is complete, no imputation done.')
} else {
  # Impute
  # TO DO: Check this is a sane approach
  # Run imputation
  # Roughly one imputation per percent of incomplete data (White et al.,2011),
  # but the more the better, 100 can easily be run on small datasets on a laptop
  # Roughly 20-30 iterations should be enough, use plot() to check convergence:
  # http://stats.stackexchange.com/questions/219013/how-do-the-number-of-imputations-the-maximum-iterations-affect-accuracy-in-mul
  # Fill in with median for now:
  impute.median <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
  input_data <- apply(input_data, 2, impute.median)
  print(sprintf('Imputed %s missing values', length(NAs)))
  # # TO DO: run proper imputation
  # input_data <- mice(input_data,
  #                      m = 2, # Number of imputed datasets, 5 is default
  #                      maxit = 3,
  #                      meth = 'pmm', # predictive mean matching, leave empty for
  #                      # auto selection depending on variable type
  #                      diagnostics = T,
  #                      seed = 500)
  # input_data <- complete(input_data, 1)
  # summary(input_data)
}
##########

##########
# Convert to a matrix, exclude the first column (which is expected to have the feature names) 
# and save feature names as rownames:
input_data <- as.data.frame(input_data)
rownames(input_data) <- input_data[, 1]
# rownames(input_data)
# input_data
input_data[1:5, 1:5]
input_data <- data.matrix(input_data[1:nrow(input_data), 2:ncol(input_data)])
class(input_data)
# str(input_data)
dim(input_data)
# head(input_data)
# str(input_data)
class(input_data)
# head(input_data) # Columns must be individuals and rows phenotypes
input_data[1:5, 1:5]
dim(input_data)
##########

##########
# Compute the PCs
input_data_prcomp <- prcomp(input_data, center = TRUE, scale = TRUE)

# Obtain values for all PCs output:
pc <- data.frame(round(input_data_prcomp$x, 4))
pc$sample_id <- rownames(pc)
pc <- pc[, moveme(names(pc), 'sample_id first')]
names(pc)[1:num_PCs]
class(pc)
pc[1:5, 1:5]
# Save file:
fwrite(pc,
       sprintf('%s', output_file_name),
       sep = '\t', na = 'NA',
       col.names = TRUE, row.names = FALSE,
       quote = FALSE)
##########

##########
# Explore dimensions and plot first 10 or so components:
dim(pc) # Extra column for sampleID
dim(input_data)
# str(input_data_prcomp)
# head(pc)

# Get cumulative proportion:
sum_pca <- summary(input_data_prcomp)
sum_pca$importance[, 1:num_PCs]
# sum_pca
sum_pca_df <- as.data.frame(sum_pca$importance)
sum_pca_df <- t(sum_pca_df)
sum_pca_df <- as.data.frame(sum_pca_df)
# View(sum_pca_df)
sum_pca_df$percent_var <- round(100 * (sum_pca_df$`Proportion of Variance`), 3)
sum_pca_df$PC <- factor(row.names(sum_pca_df), levels = row.names(sum_pca_df),
                        labels = row.names(sum_pca_df))
# head(sum_pca_df)
# tail(sum_pca_df)
# sum_pca_df[1:num_PCs, ]
# names(sum_pca_df)
# str(sum_pca_df)

# Plot proportion of variance of first x PCs:
plot_prop_vars <- ggplot(sum_pca_df[1:num_PCs, ], aes(y = percent_var, x = PC)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Principal component', y = 'Proportion of variance (%)') +
  theme_Publication() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid = element_blank(), panel.border = element_blank())
# plot_name <- sprintf('prop_var_%s.svg', output_file_name)
# ggsave(plot_name, plot_prop_vars)

# Plot top PCs:
# head(pc)
# Function for plotting any PC pair:
plot_PCs <- function(data = pc, PCa, PCb){
  xlab <- sprintf('%s', PCa)
  ylab <- sprintf('%s', PCb)
  ggplot(data = data, aes(x = data[, PCa], y = data[, PCb])) + geom_point(alpha = 0.7) +
    theme_Publication() +
    xlab(label = xlab) +
    ylab(label = ylab)
}

for (i in 1:num_PCs){
  PCa <- sprintf('PC%s', i)
  PCb <- sprintf('PC%s', i + 1)
  plot_name <- sprintf('plot_%s_%s', PCa, PCb)
  print(plot_name)
  # Name variable on the fly with assign:
  assign(plot_name, plot_PCs(pc, as.character(PCa), as.character(PCb)))
}

# Put all plots together in one figure:
cow_grid <- plot_grid(plot_prop_vars,
                      plot_PC1_PC2,
                      plot_PC2_PC3,
                      plot_PC3_PC4,
                      plot_PC4_PC5,
                      plot_PC5_PC6,
                      plot_PC6_PC7,
                      plot_PC7_PC8,
                      plot_PC8_PC9,
                      plot_PC9_PC10,
                      # align = 'vh',
                      # rel_widths = c(1, 0.05, 1),
                      # labels = c("A", "B"),
                      ncol=2)
# Save figure to disk as svg:
plot_name <- sprintf('top_%s_PCs_%s.svg', num_PCs, output_file_name)
# A4 paper measures 210 × 297 millimeters or 8.27 × 11.69 inches
svg(plot_name, width = 12, height = 12)
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
