#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
plot_flashpca.R
==================

Author: |author_names| 
Release: |version|
Date: |today|


Purpose
=======

|description|

Plots the results of flashpca output.


See flashpca:
# https://github.com/gabraham/flashpca#R


Usage and options
=================

To run, type:
Rscript plot_flashpca.R (--pcs <FILE>) (--pve <FILE>) [options]

Usage: plot_flashpca.R (--pcs <FILE>) (--pve <FILE>)
       plot_flashpca.R [options]

Options:
  --pcs <FILE>                     Input file name for pcs
  --pve <FILE>                     Input file name for pve
  -O <OUTPUT_FILE>                Output file name
  --session <R_SESSION_NAME>      R session name if to be saved
  -h --help                       Show this screen
  --num_PCs	<PCs_to_plot>         Number of PCs to plot, max 10 [default: 10]


Input:

pcs and pve files as output by flashpca.
Files must be named as:
pcs.*.flashpca.tsv
pve.*.flashpca.tsv

* meaning any character.

Output:

A multi-plot figure of the top PCs and cumulative variance explained.

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
# TO DO: sort paths out so they are read from utilities folder after installation:
source('moveme.R')
source('ggtheme.R')
source('ggtheme_bestd.R')
######################

######################
##########
# Read files, this is with data.table:
# pcs have FID and IID as first two columns
if (is.null(args[['--pcs']]) == FALSE) {
  input_name <- as.character(args[['--pcs']])
  # For tests:
  # setwd('~/Documents/quickstart_projects/chronic_inflammation_Airwave.p_q/results/QTL_core_illumina/')
  # setwd('~/Desktop/Downloads_to_delete/miscellaneous_tests/pipe_QTL_tests/results/tests_full_pipe/')
  # input_name <- 'pcs.all.clean-base.pruned.flashpca.tsv'
  # input_name <- 'simulated-dummy_binary-QC.flashpca.pcs.tsv'
  input_data_pcs <- fread(input_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
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
  # input_name <- 'simulated-dummy_binary-QC.flashpca.pve.tsv'
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

######################
##########
# Explore input data for pcs:
class(input_data_pcs)
dim(input_data_pcs) # nrow(), ncol()
str(input_data_pcs)
input_data_pcs # equivalent to head() and tail()
setkey(input_data_pcs) # memory efficient and fast
key(input_data_pcs)
tables()
colnames(input_data_pcs)
# First column with feature labels:
input_data_pcs[, 1, with = FALSE] # by column position, preferable by column name to avoid silent bugs
input_data_pcs <- as.data.frame(input_data_pcs)

# Explore input data for proportions:
class(input_data_pve)
dim(input_data_pve) # nrow(), ncol()
str(input_data_pve)
input_data_pve # equivalent to head() and tail()
colnames(input_data_pve)
# First column with feature labels:
input_data_pve[, 1, with = FALSE] # by column position, preferable by column name to avoid silent bugs
##########

##########
# Get proportion explained:
input_data_pve <- as.data.frame(input_data_pve)
input_data_pve$percent_var <- round(100 * (input_data_pve$V1), 3)
input_data_pve$PC <- factor(row.names(input_data_pve), levels = row.names(input_data_pve),
                        labels = row.names(input_data_pve))
head(input_data_pve)
tail(input_data_pve)
names(input_data_pve)
str(input_data_pve)
##########

##########
# Plot proportion of variance of first x PCs:
plot_prop_vars <- ggplot(input_data_pve[c(1:num_PCs), ], aes(y = percent_var, x = PC)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Principal component', y = 'Proportion of variance (%)') +
  theme_Publication() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid = element_blank(), panel.border = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "mm"))
# plot_name <- sprintf('prop_var_%s.svg', output_file_name)
# ggsave(plot_name, plot_prop_vars)

# Plot top PCs:
pc <- input_data_pcs[, -c(1)]
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
# plot_PCs(pc, 'PC1', 'PC2')

for (i in 1:num_PCs){
  PCa <- sprintf('PC%s', i)
  PCb <- sprintf('PC%s', i + 1)
  plot_name <- sprintf('plot_%s_%s', PCa, PCb)
  print(plot_name)
  # Name variable on the fly with assign:
  assign(plot_name, plot_PCs(pc, as.character(PCa), as.character(PCb)))
}

# TO DO: num_PCs has to be 10 for now, will fail otherwise.
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
plot_name <- sprintf('top_10_PCs_%s.svg', output_file_name)
# A4 paper measures 210 ?? 297 millimeters or 8.27 ?? 11.69 inches
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
