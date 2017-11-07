#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
script_name
===============

Author: |author_names| 
Release: |version|
Date: 20 May 2016


Purpose
=======

|description|

Script to run the mePipe package by Peter Humburg. See:

https://github.com/jknightlab/mePipe

Usage and options
=================

Usage: script_name (--geno <SNPs file>)
       script_name (--gex <phenotype file>)
       script_name [--cov <covariates_file>]
       script_name [--PCs <start and end of principal components to adjust for>]
       script_name [--model <MatrixEQTLs model>]
       script_name [--threshold <p-value>]
       script_name [--cisThreshold <cis p-value>]
       script_name [--cis <cis distance definition>]1e+06
       script_name [-O <output file name>]
       script_name [-h | --help]
       script_name [--session <R_SESSION_NAME>] 

Options:
--geno <SNPs file>          File name with genetic data   
--gex <phenotype file>      File name with phenotypes (e.g. gene expression, metabolites, etc.)
--PCs <range>               Start and end of principal components to adjust for, comma separated. Step is always 5 [default: 0,50].
--model <MatrixEQTLs model> One of linear, modelANOVA or modelLINEAR or modelLINEAR_CROSS [default: modelLINEAR].
--threshold <p-value>       p-value threshold at which pheno-SNP associations are saved [default: 1e-08].
--cisThreshold <cis p-value> p-value threshold at which cis associations are saved [default: 1e-05].
--cis <cis distance definition> Distance for defining cis vs trans [default: 1e+06].


-O <OUTPUT_FILE>           Output file name [default: <infile>_QTL.txt]
--session <R_SESSION_NAME> R session name if to be saved
-h --help                  Show this screen

Input:

Requires gene expression and genotype data

The same as MatrixEQTL requires. See:

http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

These are tab separated files without headers. Files here are read with data.table and stringsAsFactors = FALSE

Output:

Namely a qqlot and tables of genotype molecular phenotype associations. These are saved in the working directory.

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
# See:
# https://cran.r-project.org/web/packages/docopt/docopt.pdf
# docopt(doc, args = commandArgs(TRUE), name = NULL, help = TRUE,
# version = NULL, strict = FALSE, strip_names = !strict,
# quoted_args = !strict)

# Within the script specify options as:
# args[['--session']]
# args $ `-I` == TRUE

# Set up MatrixEQTL options:
geno_file <- as.character(args[['--geno']])
# geno_file <- 'SNP.txt'

# Gene expression levels:
Gex_file <- as.character(args[['--gex']])
# Gex_file <- 'GE.txt'

PC_seq_to_test  <- seq(args[['--PCs']], by = 5)
# PC_seq_to_test = seq(0, 50, by = 5)

useModel <- as.character(args[['--model']]) 
threshold  <- as.numeric(args[['--threshold']])
cisThreshold <- as.numeric(args[['--cisThreshold']])
cis <- as.numeric(args[['--cis']])

# Set output file names:
output_file_trans = sprintf('trans_%s_%s.txt', threshold, geno_file)
output_file_cis = sprintf('cis_%s_%s_%s.txt', cis, cisThreshold, geno_file)

# Error covariance matrix, pass as numeric() but not needed for mePipe function runME(), leave empty.

# Get PCs for gene expression data:
covariates_file = sprintf('PCs_%s', Gex_file)
print('Getting principal components from expression/metabolite file')
allCovariates(Gex_file, output = covariates_file, sep = '\t')

# Identify covariates that are significantly associated with genotype:
covOut_filename = sprintf('PCs_to_adjust_for_%s.txt', geno_file)
covAssoc(genotype = geno_file, covariate = covariates_file, output = sprintf('covars_assoc_to_genotype_%s', geno_file),
         covOut = covOut_filename)

# Determine the optimal number of covariates after excluding those associated with genotypes and runME 
# at the same time:
output_dir_covSelect = sprintf('covSelect_%s_%s_%s_%s', threshold, cisThreshold, cis, geno_file)


# Print to screen:
str(args)

##################
# TO DO: generate covariate PCs from GEx file first, then assoc cov to genotype, then filter these from run_mePipe run.
# chooseCov or runMe don't have these options (assocCov and filterCov) so run covariates.R, assocCov, filterCov and pass these to
# chooseCov
# default FDR threshold for PC exclusion is 0.001
#(command line does though, check this again to make running straightforward) 
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
# https://github.com/Rdatatable/data.table/wiki/Getting-started

# For mePipe, see:
# biocLite(c("MatrixEQTL", "trio", "optparse", "XML", "snow"))
# library(devtools)
# biocLite('optparse')
# # devtools::install_github("humburg/Rsge")
# system("wget https://github.com/jknightlab/mePipe/releases/download/mePipe_v1.3.5/mePipe_1.3.5.tar.gz")
# install.packages("mePipe_1.3.5.tar.gz", repos=NULL)
library(mePipe)

# Switch Rsge off for mePipe:
sge.options(sge.use.cluster = FALSE) #sge.user.options=opt$sgeoptions)
# Also check the following as errors in cgat150
# make_option("--sgeoptions", default="-S /bin/bash -V",
# help="Options to pass to qsub (ignored unless `--cluster` is used. [default: %default])"),
# mePipe::
mePipe::getOptions()

library(ggplot2)
######################

######################
# Read files, this is with data.table:
if (is.null(args[['-I']]) == FALSE) {
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


#--------------------------------------------------------------------------------------------

#################
# Run mePipe
choose_covariates = chooseCov(expression = Gex_file, genotype = geno_file, covariate = covOut_filename, 
                              output = output_dir_covSelect, # output directory with results
                              # snpsPos = snpspos, genePos = genepos, # runME options here
                              # output = output_file_trans, cisOutput = output_file_cis, # Leave this line out
                              # as chooseCov() and runME() output options collide. output will be saved to covSelect
                              # directory. trans files are .eQTl and cis files are .eQTL.cis
                              doCis = F, doTrans = F, candidates = PC_seq_to_test,
                              threshold = threshold, cisThreshold = cisThreshold, model = useModel, 
                              cis = cis, bins = 1000, qqplot = TRUE, 
                              cluster = sge.getOption("sge.use.cluster"), verbose = TRUE)
print(choose_covariates)

# Save as dataframe to disk to read for next script:
choose_covariates_df <- as.data.frame(choose_covariates)
head(choose_covariates_df)
dim(choose_covariates_df)

write.table(choose_covariates_df, sprintf('choose_covariates_df_%s.txt', geno_file), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

# Create soft links to cis and trans files and qqplot:
base_path <- sprintf('%s/', getwd())

# File with maximum number of eQTL results based on covariate loop:
trans_file <- choose_covariates$best[[1]]
trans_file
trans_file_path <- sprintf('%s%s', base_path, trans_file)
cmd_soft_link <- sprintf('ln -s %s QTL_%s', trans_file_path, dirname(trans_file))
cmd_soft_link
system(cmd_soft_link)
# File with maximum number of eQTL results based on covariate loop:
#cis_file <- choose_covariates$best[[1]]
#cis_file
#cis_file_path <- sprintf('%s%s', base_path, cis_file)
#cmd_soft_link <- sprintf('ln -s %s cis_%s', cis_file_path, dirname(cis_file))
#cmd_soft_link
#system(cmd_soft_link)
# qqplots are saved as PDFs in covSelect directory, save both trans and cis:
qqplot_file <- paste(choose_covariates$best[[1]], '.pdf', sep = '')
qqplot_file
qqplot_file_path <- sprintf('%s%s', base_path, qqplot_file)
cmd_soft_link <- sprintf('ln -s %s qqplot_%s_QTL.pdf', qqplot_file_path, dirname(qqplot_file))
cmd_soft_link
system(cmd_soft_link)

# qqplots all (cis and trans are plotted together) end as xxx.eQTL.pdf:
#qqplot_file <- paste(choose_covariates$best[[1]], '.pdf', sep = '')
#qqplot_file
#qqplot_file_path <- sprintf('%s%s', base_path, qqplot_file)
#cmd_soft_link <- sprintf('ln -s %s qqplot_%s_cis.pdf', qqplot_file_path, dirname(qqplot_file))
#cmd_soft_link
#system(cmd_soft_link)

# Plot change of identified eQTLs for different number of covariates:

# Errors:
# pdf(sprintf('choose_covariates_%s.pdf', geno_file))
# plotCov(candidates = choose_covariates$covariates, eqtls = choose_covariates$eqtls$significant, 
#         selected = choose_covariates$selected)
# dev.off()

#ggplot(data = as.data.frame(choose_covariates), aes(y = eqtls.significant.cis, x = covariates)) +
#  geom_point() + geom_line() + 
#  xlab('Gene expression principal components') + ylab('Unique probes with associated cis SNPs')
#ggsave(sprintf('chooseCov_%s_%s_%s.pdf', cisThreshold, cis, geno_file))


ggplot(data = as.data.frame(choose_covariates), aes(y = eqtls.significant.all, x = covariates)) +
  geom_point() + geom_line() + 
  xlab('Metabolome principal components') + ylab('Unique data-points with associated SNPs')
ggsave(sprintf('chooseCov_%s_%s.pdf', threshold, geno_file))
#################


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