#!/usr/bin/env Rscript
# R script to run with docopt for command line options:
'
run_matrixEQTL.R
=================

Author: |author_names| 
Release: |version|
Date: 20 May 2016


Purpose
=======

|description|

This is a wrapper for matrixEQTL

http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

The script is intended to run as part of a Ruffus (python based) pipeline. Naming files in the expected way is important (see below).


Usage and options
=================

Usage: run_matrixEQTL.R (--geno <SNPs_file> --gex <phenotype_file>)
       run_matrixEQTL.R [options]
       run_matrixEQTL.R [-h | --help]

Options:
   --geno SNPs_file                     File name with genetic data
   --gex phenotype_file                 File name with phenotypes
   -O OUTPUT_FILE                       Output file name, ".MxEQTL" is added to outputs
   --cov covariates_file                File name with covariates to adjust for
   --condition tissue_or_context        Specify name of condition (if running multiple files/groups/tissues/etc.)
   --snpspos SNP_position_file          File name with SNP position information
   --genepos phenotype_position_file    File name with molecular phenotype annotation information
   --model MatrixEQTLs_model            One of modelANOVA, modelLINEAR or modelLINEAR_CROSS [default: modelLINEAR].
   --pvOutputThreshold p-value          trans associations p-value threshold at which pheno-SNP associations are saved [default: 1e-08].
   --pvOutputThreshold.cis cis_p-value  cis p-value threshold at which associations are saved if files on positions are given [default: 1e-05].
   --cisDist cis_distance_definition    Distance for defining cis vs trans [default: 1e+06].
   -h --help                            Show this screen
   --session R_SESSION_NAME             R session name if to be saved



Input
=====

Quality controlled (molecular) phenotype (e.g. gene expression) and genotyping data. 

Optionally covariates and error covariance matrix.

Annotation files are also needed for cis vs trans analysis (SNPs positions and probe positions and any associated annotations).

Files need to be in the same format as MatrixEQTL requires:

SNPs/genes/covariates in rows, individuals in columns with (dummy) headers and row names (the first column and first row are skipped when read).

Missing values must be set as "NA".

SNP and probe position files, e.g. 
   - snp146Common_MatrixEQTL_snp_pos.txt
   - biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt
must be tab separated. 

See:

http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/


Output
======

Namely a qqlot and tables of genotype molecular phenotype associations. These are saved in the working directory.


Naming convention for input files
=================================

Output files get named based on the input files. The script assumes "cohort" is the same for input files (but only takes it from the genotype file).

Please rename your files in the following way (use soft links to rename for example):

Infile: cohort-platform-other_descriptor.suffix

Outfile: cohort-platform_infile1-descriptor1-platform_infile2-descriptor2.new_suffix

For example:

genotype file: airwave-illumina_exome-all_chrs.geno

phenotype file: airwave-NMR-blood.pheno

covariates file: airwave-NMR-blood.cov

and depending on the input and arguments you might get:

airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.cis.tst
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.trans
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.qqplot.svg
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.degrees_condition.txt
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.log

File names can get long so use abbreviations or short versions.

You can also override this and simply choose your outfile prefix.

If you do not use an outfile name and your files do not follow the naming above you might get something like:

"SNP.txt-NA-NA-NA-NA.MxEQTL"

This script will take files with any suffix. If running in a pipeline .geno, .pheno and .cov are required.

Requirements
============

docopt
data.table
ggplot2


Documentation
=============

For more information see:

|url|
' -> doc

##################
# TO DO: generate covariate PCs from GEx file first, then assoc cov to genotype, then filter these from run_mePipe run.
# chooseCov or runMe don't have these options (assocCov and filterCov) so run covariates.R, assocCov, filterCov and pass these to
# chooseCov
# default FDR threshold for PC exclusion is 0.001
#(command line does though, check this again to make running straightforward) 
# PC_seq_to_test  <- seq(args[['--PCs']], by = 5)
# PC_seq_to_test = seq(0, 50, by = 5)
######################


######################
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
library(ggplot2)
library(MatrixEQTL)
######################


######################
##########
# Set-up file names:
SNP_file <- as.character(args[['--geno']])
expression_file <- as.character(args[['--gex']])
covariates_file <- as.character(args[['--cov']])
condition <- as.character(args[['--condition']])
# For testing:
# SNP_file <- 'airwave-NMR-blood.geno'
# expression_file <- 'airwave-NMR-blood.pheno'
# covariates_file <- 'airwave-NMR-blood.cov'
# SNP_file <- 'SNP.txt'
# expression_file <- 'GE.txt'
# covariates_file <- 'Covariates.txt'

# Load annotations for SNP location, probe annotation and probe location:
snpspos <- as.character(args[['--snpspos']])
genepos <- as.character(args[['--genepos']])
# snpspos <- 'snp146Common_MatrixEQTL_snp_pos.txt'
# snpspos <- 'snpsloc.txt'
# genepos <- 'biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt'
# genepos <- 'geneloc.txt'

# Set up MatrixEQTL options:
useModel <- as.character(args[['--model']])
if (useModel == 'modelLINEAR') {
  useModel <- modelLINEAR
} else if (useModel == 'modelANOVA') {
    useModel = modelANOVA
    } else if (useModel == 'modelLINEAR_CROSS') {
      useModel = modelLINEAR_CROSS
      } else (useModel <- modelLINEAR)
# These will appear as:
# modelLINEAR = 117348L;
# modelANOVA  = 47074L;
# modelLINEAR_CROSS = 1113461L;

# The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. 
# Note that for larger datasets the threshold should be lower. 
# Setting the threshold to a high value for a large dataset may cause excessively large output files.
pvOutputThreshold  <- as.numeric(args[['--pvOutputThreshold']])
pvOutputThreshold.cis <- as.numeric(args[['--pvOutputThreshold.cis']])
cisDist <- as.numeric(args[['--cisDist']])

# Define the covariance matrix for the error term. 
# Consider an error covariance matrix if necessary (correlated variables or errors)
# This parameter is rarely used. If the covariance matrix is a multiple of identity, set it to numeric().
errorCovariance <- numeric()

# Print all arguments to screen:
str(args)
##########

##########
# Set output file names:
if (is.null(args[['-O']])) {
  stopifnot(!is.null(args[['--geno']]), !is.null(args[['--gex']]))
  print('Output file name prefix not given. Using:')
  # Separate name components based on convention for input in this script:
  infile_1 <- strsplit(SNP_file, '-')
  cohort <- infile_1[[1]][1]
  platform1 <- infile_1[[1]][2]
  descriptor_1 <- strsplit(infile_1[[1]][3], '\\.')[[1]][1]
  # Do the same for the phenotype file, take cohort only from genotype file:
  infile_2 <- strsplit(expression_file, '-')
  platform2 <- infile_2[[1]][2]
  descriptor_2 <- strsplit(infile_2[[1]][3], '\\.')[[1]][1]
  output_file_name <- sprintf('%s-%s-%s-%s-%s.MxEQTL', cohort, platform1, descriptor_1, platform2, descriptor_2)
  print(output_file_name)
} else {
  output_file_name <- as.character(args[['-O']])
  # output_file_name <- 'testing'
  output_file_name <- sprintf('%s.MxEQTL', output_file_name)
  print(sprintf('Output file names will contain %s', output_file_name))
}
##########
##################


######################
# Read covariates and position/annotation files if given:
# For covariates file use matrixEQTL functions to read in:
if (is.null(args[['--cov']])) {
  covariates_file <- character()
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values; This is from the plink encoding.
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
  # cvrt$LoadFile(covariates_file);
  cvrt$CreateFromMatrix(data.matrix(covariates_file)); # use this if covariates_file = character()
  print('Covariates file not provided, running analysis without.')
  cvrt
} else {
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values; This is from the plink encoding.
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
  cvrt$LoadFile(covariates_file);
  # cvrt$CreateFromMatrix(data.matrix(covariates_file)); # use this if covariates_file = character()
  cvrt
}

# Values for annotation files:
if (is.null(args[['--snpspos']]) | is.null(args[['--genepos']])) {
  print('SNP and/or probe (gene expression, metabolites, etc.) position files not provided, running analysis without (no cis vs trans).')
} else {
  snpPos <- fread(snpspos, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  head(snpPos)
  dim(snpPos)
  probePos <- fread(genepos, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  head(probePos)
  dim(probePos)
  # Switch classes as these are read by data.table and cause problems later:
  snpPos <- as.data.frame(snpPos)
  class(snpPos)
  probePos <- as.data.frame(probePos)
  class(probePos)
  # Set cis and trans output file names:
  output_file_name.cis <- sprintf('%s.cis.tsv', output_file_name)
  output_file_name.trans <- sprintf('%s.trans.tsv', output_file_name)
}

# Provide a name for the condition being run if none given:
# If neither --condition or the naming convention appear:
time_run <- as.character(gsub(' ', '_', as.character(Sys.time())))
if (is.null(args[['--condition']]) & is.na(descriptor_2)) {
    print(sprintf("You did not provide a name for the condition (i.e. tissue) and your input files
                   don't follow this script's convention, using a timestamp instead: %s.", time_run))
    condition <- time_run
  } else {
    # If --condition isn't given but the naming convention is followed so can pick it up:
    if (is.null(args[['--condition']]) & !is.na(descriptor_2)) {
  print(sprintf('You did not provide a name for the condition file, using a descriptor from the file names instead: %s', descriptor_2))
  condition <- sprintf('%s', descriptor_2)
    } else {
      print('Could not find a name for the condition variable, using a timestamp instead')
      condition <- time_run
    }
    }
# If --condition is given that takes priority (through docopt without specifying more)
# MatrixEQTL reads in the main files.
######################


######################
##########
# Run MatrixeQTL
# Load the files:
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file)
snps
# snps$CreateFromMatrix(data.matrix(SNP_file)); # Data need to be passed as a matrix

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile(expression_file);
gene
# gene$CreateFromMatrix(data.matrix(expression_file));
##########

##########
# Sanity check:
# MatriEQTL will do this anyway:
# Stop if number of column numbers are not the same:
# if (ncol(snps) != ncol(gene)) {
#   # Stop:
#   print('Number of columns between SNPs and phenotype files differ, these need to be the same (and in the same order, but this is not checked)')
#   stopifnot(ncol(snps) == ncol(gene)) 
# }
##########

##########
# Call the main Matrix eQTL function
# If cis information available run:
if (!is.null(args[['--snpspos']]) | !is.null(args[['--genepos']])) {
  me <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name.trans,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = 'qqplot',
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    pvOutputThreshold = pvOutputThreshold,
    output_file_name.cis = output_file_name.cis,
    pvOutputThreshold.cis = pvOutputThreshold.cis,
    snpspos = snpPos,
    genepos = probePos,
    cisDist = cisDist)
} else {
  # Make the main results table end in tsv if no cis or trans being run:
  output_file_name_no_cis <- sprintf('%s.tsv', output_file_name)
  me <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name_no_cis,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = 'qqplot',
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    pvOutputThreshold = pvOutputThreshold
    # output_file_name.cis = output_file_name.cis,
    # pvOutputThreshold.cis = pvOutputThreshold.cis,
    # snpspos = snpPos,
    # genepos = probePos,
    # cisDist = cisDist
    )
}
##########

##########
# # Inspect results:
# # Each significant gene-SNP association is recorded in a separate line in the output file and in the returned object me. 
# # In case of cis/trans eQTL analysis described below, two output files are produced, one with cis-eQTLs, another only with trans. 
# # Every record contains a SNP name, a transcript name, estimate of the effect size, t- or F-statistic, p-value, and FDR.
# show(me$all$eqtls)
# head(me$all$eqtls)
# plot(me)
# qqnorm(me$all$eqtls[, 4])
# head(me$cis$eqtls)
# nrow(me$cis$eqtls)
# head(me$cis$ntests)
# head(me$trans$eqtls)
# nrow(me$trans$eqtls)
# head(me$trans$ntests)
##########
######################

######################
# Modify plot for publication (SF8):
# setwd('/Users/antoniob/Documents/quickstart_projects/BEST_D_molecular.p_q/results/repro_re_runs/')
# load('R_session_saved_image_MatrixQTL_analysis.RData')
# source('../../code/BEST_D_molecular/utilities/ggtheme.R')
# library(ggplot2)
# plot:
svg(sprintf('%s.qqplot.svg', output_file_name))
plot(me,
     main = '',
     # cex.main = 1.25,
     cex.lab = 1.25,
     cex.axis = 0.75)
dev.off()


## Save degrees of freedom in order to be able to run multi-tissue Matrix EQTL:
degrees_filename <- sprintf('%s.degrees_condition.txt', output_file_name)
cat(file = degrees_filename, condition, "\t", me$param$dfFull, '\n', append = TRUE)
######################

######################
## Save some text:
# Methods
# Legend
# Interpretation
# cat(file <- output_file, some_var, '\t', another_var, '\n', append = TRUE)
######################


######################
##########
# Save a minimal file with parameters called and basic counts (counts get saved as a separate dataframe though):
log_file <- sprintf('%s.log', output_file_name)
sink(log_file, append = TRUE, split = TRUE, type = c("output", "message"))
print(paste('Minimal log file for parameters called with:',
            'run_matrixEQTL.R'))
Sys.time()
print(paste('Working directory :', getwd()))
print('Arguments called from the command line:')
str(args)
print('Arguments saved in matrixEQTL object:')
me$param
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
# Total QTLs:
cat('QTLs tested (cis plus trans):', '\t')
me$all$ntests
cat(sprintf('Significant QTLs below raw p-value %s:', pvOutputThreshold), '\t')
me$all$neqtls
# Cis counts:
cat('Local (cis) QTLs tested:', '\t');
me$cis$ntests
cat(sprintf('Significant local (cis) QTLs below raw p-value %s:', pvOutputThreshold.cis), '\t')
me$cis$neqtls
# Trans counts:
cat('Distant (trans) QTLs tested:', '\t');
me$trans$ntests
cat(sprintf('Significant distant (trans) QTLs below raw p-value %s:', pvOutputThreshold), '\t')
me$trans$neqtls
# Close:
sink()
##########

##########
# Save a minimal counts file with results:
counts_df <- list(as.character(output_file_name),
                        me$all$ntests,
                        me$all$neqtls,
                        me$cis$ntests,
                        me$cis$neqtls,
                        me$trans$ntests,
                        me$trans$neqtls)
# If there are NULL values convert them to strings:
counts_df <- sapply(counts_df, function(x) ifelse(is.null(x), 'NULL', x))
counts_df <- as.data.frame(counts_df, stringsAsFactors = FALSE)
# Variable names:
df_rows <- c('dataset',
             'total_QTLs_tested',
             sprintf('QTLs_below_p_%s', pvOutputThreshold),
             'local_QTLs_tested',
             sprintf('local_QTLs_below_p_%s', pvOutputThreshold.cis),
             'distant_QTLs_tested',
             sprintf('distant_QTLs_below_p_%s', pvOutputThreshold)
             )
counts_df <- cbind(df_rows, counts_df)
# This is extra but gets delted when saving to file:
rownames(counts_df) <- df_rows
counts_file <- sprintf('%s.counts.tsv', output_file_name)
write.table(file = counts_file,
            x = counts_df,
            sep = '\t',
            na = 'NA',
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
##########
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