#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
plot_template.R
===============

Author: |author_names| 
Release: |version|
Date: |today|


Purpose
=======

|description|


Usage and options
=================

These are based on docopt_ for R:

https://github.com/docopt/docopt.R
https://cran.r-project.org/web/packages/docopt/index.html

To run, type:
    Rscript plot_template.R -I <INPUT_FILE> [options]

Usage: plot_template.R [-I <INPUT_FILE>] [--session=<R_SESSION_NAME>]
       plot_template.R [-h | --help]

Options:
  -I <INPUT_FILE>                           Input file name
  --session=<R_SESSION_NAME>                R session name if to be saved
  -h --help                                 Show this screen

Input:

    A tab separated file with headers. This is read with data.table and stringsAsFactors = FALSE

Output:

    Plots as svg files.

Requirements:

    library(docopt)
    library(ggplot2)
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
# See:
# https://cran.r-project.org/web/packages/docopt/docopt.pdf
# docopt(doc, args = commandArgs(TRUE), name = NULL, help = TRUE,
# version = NULL, strict = FALSE, strip_names = !strict,
# quoted_args = !strict)

# Print to screen:
str(args)
# Within the script specify options as e.g.:
# args[['--session']]
# args $ `-I` == TRUE
######################

######################
# Logging
# This can be taken care of by CGAT Experiment.py if running as a pipeline.
# Otherwise there seem to be few good alternatives. A workaround is this code:
# logging.R
# In the script_templates dir of project_quickstart.
# It does not run on it own though, needs copy/pasting for now.
######################

######################
# Load a previous R session, data and objects:
#load('R_session_saved_image_order_and_match.RData', verbose=T)
######################


######################
# Import libraries
# source('http://bioconductor.org/biocLite.R')
library(ggplot2)
# library(ggthemes)
library(data.table)
######################

######################
# Read files:
if (is.null(args[['-I']]) == FALSE) {
  # args[['-I']] <- as.character('pandas_DF.tsv') # For testing
  input_name <- as.character(args[['-I']])#(args $ `-I`)
  input_data <- fread(input_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
} else {
  # Stop if arguments not given:
  print('You need to provide an input file. This has to be tab separated with headers.')
  stopifnot(!is.null(args[['-I']]) == TRUE)
}

print('File being used: ')
print(input_name)
# Split at the last '.':
input_name <- strsplit(input_name, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
print('Name being used to save output files: ')
print(input_name)
######################

######################
# Explore data:
class(input_data)
dim(input_data) # nrow(), ncol()
str(input_data)
input_data # equivalent to head() and tail()
setkey(input_data) # memory efficient and fast
key(input_data)
tables()
colnames(input_data)
input_data[, 2, with = FALSE] # by column position, preferable by column name to avoid silent bugs
###################### 

######################
# Convert data.table to data.frame as easier to handle for plotting:
setDF(x = input_data)
# input_data <- as.data.frame(input_data)
class(input_data)
######################

######################
####
# Histogram overlaid with kernel density curve
# http://www.cookbook-r.com/Graphs/Plotting_distributions_(ggplot2)/

# Setup:
plot_name <- sprintf('%s_%s_histogram.svg', input_name, x_var)
x_var_label <- x_var
# Plot:
ggplot(input_data, aes(x = input_data[, x_var])) +
       geom_histogram(aes( y = ..density..), # Histogram with density instead of count on y-axis
                 binwidth = 0.5,
                 colour = "black", fill = "white") +
       geom_density(alpha = 0.2, fill = "#FF6666") + # Overlay with transparent density plot
       ylab('density') +
       xlab(x_var_label)
# Save to file:
ggsave(plot_name)
# Prevent Rplots.pdf from being generated. ggsave() without weight/height opens a device.
# Rscript also saves Rplots.pdf by default, these are deleted at the end of this script.
dev.off()
####

####
# A boxplot
# Setup:
var3_label <- var3
var3_factor <- factor(input_data[, var3])
y_var_label <- y_var
plot_name <- sprintf('%s_%s_%s_boxplot.svg', input_name, var3, y_var)
# Plot:
ggplot(input_data, aes(x = var3_factor, y = input_data[, y_var], fill = var3_factor)) +
       geom_boxplot() +
       ylab(y_var_label) +
       xlab(var3_label) +
       theme_classic() +
       theme(legend.position = 'none')
# Save to file:
ggsave(plot_name)
# Prevent Rplots.pdf from being generated. ggsave() without weight/height opens a device.
# Rscript also saves Rplots.pdf by default, these are deleted at the end of this script.
dev.off()
####

####
# Scatterplot and legend:
# http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
# Setup:
x_var_label <- x_var
y_var_label <- y_var
plot_name <- sprintf('%s_%s_%s_scatterplot.svg', input_name, x_var, y_var)
# Plot:
ggplot(input_data, aes(x = input_data[, x_var], y = input_data[, y_var], colour = var3_factor)) +
       geom_point() +
       geom_smooth(method = lm) +
       ylab(y_var_label) +
       xlab(x_var_label) +
       labs(colour = var3) +
       theme_classic()
# Save:
ggsave(plot_name)
# Prevent Rplots.pdf from being generated. ggsave() without weight/height opens a device.
# Rscript also saves Rplots.pdf by default, these are deleted at the end of this script.
dev.off()
####
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
sessionInfo()
print('Finished successfully')
q()

# Next: run the script for xxx
######################
