########################
# Script for plotting eQTLs following MatrixEQTL analysis
# Antonio J Berlanga-Taylor
# 30 Sept 2015

########################

#############################################
##Set working directory and file locations and names of required inputs:
options(echo = TRUE, digits = 3)

# Working directory:
#setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_eQTL_plotting",Sys.Date(),".txt", sep=""))
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))
getwd()

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
#load('R_session_saved_image_order_and_match.RData', verbose=T)
#load('R_session_saved_image_eQTL_responseQTLs.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_eQTL_plotting.RData', sep='')
R_session_saved_image
####################

######################
# This function allows other R scripts to obtain the path to a script directory
# (ie where this script lives). Useful when using source('some_script.R')
# without having to pre-specify the location of where that script is.
# This is taken directly from:
# How to source another_file.R from within your R script · molgenis/molgenis-pipelines Wiki
# https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
# Couldn't find a licence at the time (12 June 2018)
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
    this.file = NULL
    # This file may be 'sourced'
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }

    if (!is.null(this.file)) return(dirname(this.file))

    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))

    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}
Rscripts_dir <- LocationOfThisScript()
print('Location where this script lives:')
Rscripts_dir
# R scripts sourced with source() have to be in the same directory as this one
# (or the path constructed appropriately with file.path)
######################

####################
# Load packages:
library(ggplot2)
library(data.table)
library(car) # Run correlations and plots for lm
library(gvlma) # Compare models and fit

# TO DO: sort out paths

source(file.path(Rscripts_dir, 'functions_for_MatrixeQTL.R'))
source(file.path(Rscripts_dir, 'moveme.R')) #, chdir = TRUE)
########################


####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)
# setwd('~/Documents/quickstart_projects/chronic_inflammation_Airwave.p_q/results/QTL_core_illumina')
# getwd()

# Get files:
# TO DO: automate to plot top ten SNPs/probes:
eQTL_file_name <- as.character(args[1])
# eQTL_file_name <- "matched_all.clean-base.A-transpose-NA-NA.MxEQTL.tsv"

snp_file_name <- as.character(args[2])
# snp_file_name <- 'matched_all.clean-base.A-transpose.matrixQTL.geno'

probe_file_name <- as.character(args[3])
# probe_file_name <- 'matched_AIRWAVE_1DNMR_BatchCorrected_log_Data_Var_Sample.transposed.tsv'

snp <- as.character(args[4])
# snp <- 'exm847632' # 'rs7072216'
snp_col <- sprintf('X%s', snp)

probe <- as.character(args[5])
# probe <- '2.711437729'
probe_col <- sprintf('X%s', probe)

PC_file <- as.character(args[6])
# PC_file <- 'matched_all_covar_PCs.tsv'

PCs_to_correct <- as.numeric(args[7])
# PCs_to_correct <- '45'

print(args)
####
# Final PCs to correct for were:
# See: /ifs/projects/proj043/analysis.dir/mePipe_runs_2.dir
# trans final = me_25_covariates.eQTL
# trans baseline = me_35_covariates.eQTL
# cis final = me_35_covariates.eQTL.cis
# cis baseline = me_40_covariates.eQTL.cis
####

# eQTL_file_name="reQTLs_FDR5_cis_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_baseline.tsv_matched.tsv.txt"
# snp_file_name="cut_genotype_data_all_treated_final.tsv_matched.tsv"
# probe_file_name="cut_GEx_treated_4000_and_2000.tsv_matched.tsv"
# snp="rs17385885"
# probe="ILMN_1674985"
########################

# ####################
# TO DO # Set up file naming for outputs:
# # to create: 2000+4000-baseline-1.trans_in_cis
# eQTL_file1_base <- strsplit(eQTL_file1, '[.]')
# eQTL_file1_base <- eQTL_file1_base[[1]][1]
# eQTL_file1_base
# eQTL_file1_ext <- strsplit(eQTL_file1, '_')
# eQTL_file1_ext <- eQTL_file1_ext[[1]][2]
# eQTL_file1_ext
# 
# eQTL_file2_base <- strsplit(eQTL_file2, '[.]')
# eQTL_file2_base <- eQTL_file2_base[[1]][1]
# eQTL_file2_base
# eQTL_file2_ext <- strsplit(eQTL_file2, '_')
# eQTL_file2_ext <- eQTL_file2_ext[[1]][2]
# eQTL_file2_ext
# 
# output_file_name <- sprintf('%s_in_%s_%s.eQTL', eQTL_file1_ext, eQTL_file2_ext, eQTL_file1_base)
# output_file_name
# ####################

########################
# Read files:
eQTL_data <- fread(eQTL_file_name, sep = '\t', header = TRUE, stringsAsFactors = F)
eQTL_data

snp_file <- fread(snp_file_name, sep = '\t', header = TRUE, stringsAsFactors = F)
setkey(snp_file, FID)
snp_file[1:5, 1:5]
dim(snp_file)
# TO DO: Eliminate non-unique SNPs before eQTL run:
length(unique(snp_file[['FID']]))
snp_file <- unique(snp_file, by = 'FID')
dim(snp_file)
# str(snp_file)
snp_file[1:5, 1:5, with = F]

probe_file <- fread(probe_file_name, sep = '\t', header = T, stringsAsFactors = F)
setkey(probe_file, FID)
dim(probe_file)
# probe_file[, 1] <- as.character(probe_file[, 1])
probe_file[1:5, 1:5, with = F]
str(probe_file[, 1])

# Subset files so that only the SNP and probe of interest are kept, columns are samples and rows probes/SNPs:
which(snp_file[snp, ] == snp)
snp_data <- as.data.frame(snp_file[snp, ])
# head(snp_data)
dim(snp_data)

probe_data <- as.data.frame(probe_file[as.numeric(probe), ])
# head(probe_data)
dim(probe_data)
probe_data[, 1:5]

# Read in PC data to adjust for, columns are samples and rows PCs:
PC_data <- fread(PC_file, sep = '\t', header = T, stringsAsFactors = F)
# setkey(PC_data, FID) # Setting a key re-orders (unix like) so be careful
dim(PC_data)
# PC_data
PC_data[1:5, 1:5, with = F]
# str(PC_data)

# Subset and transpose principal components to those which were corrected for in the genome-wide analysis:
dim(PC_data[1:PCs_to_correct, , ])
PC_data <- PC_data[1:PCs_to_correct, , ]
PC_data[1:5, 1:5, with = F]
# head(PC_data)
dim(PC_data)
PC_data <- as.data.frame(PC_data) # Convert to data.frame otherwise function doesn't work (with data.table)
PC_data_t <- transpose_file(PC_data, 1)
setcolorder(PC_data_t, moveme(names(PC_data_t), "rownames first"))
class(PC_data_t)
# head(PC_data_t)
dim(PC_data_t)
rownames(PC_data_t)
# colnames(PC_data_t)
PC_data_t[1:5, 1:5, with = F]
# Check variable definitions for linear model and plotting:
PC_data_t <- PC_data_t[, rownames:=as.character(rownames)]
# str(PC_data_t)
# Change to dataframe and rename rownames for merging later:
PC_data_t <- as.data.frame(PC_data_t)
row.names(PC_data_t) <- PC_data_t[, 'rownames']
PC_data_t <- PC_data_t[, -1]
PC_data_t[1:5, 1:5]
dim(PC_data_t)

# Sanity checks, all should be true:
# TO DO: stop if false:
identical(colnames(snp_data), colnames(probe_data))
identical(colnames(snp_data), colnames(PC_data))
# identical(colnames(probe_data)[-1], PC_data_t[['rownames']]) # Rownames for tranposed PC data
identical(colnames(probe_data)[-1], rownames(PC_data_t)) # Rownames for tranposed PC data
########################

########################
# Get SNP and probe data:
snp_probe_data <- rbind(snp_data, probe_data)
snp_probe_data[, 1:5]

# Process for ggplot2:
snp_probe_data <- t(snp_probe_data) #, check.names = FALSE)
snp_probe_data <- as.data.frame(snp_probe_data)
head(snp_probe_data)
class(snp_probe_data)
str(snp_probe_data)

# Set headers: 
# colnames(snp_probe_data)[1] <- as.character(snp)
colnames(snp_probe_data)[1] <- snp_col
colnames(snp_probe_data)[2] <- probe_col
colnames(snp_probe_data)
head(snp_probe_data)

# Remove first line:
snp_probe_data <- snp_probe_data[-1, ]
head(snp_probe_data)

# Remove NAs:
length(which(is.na(snp_probe_data[, 1])))
length(which(is.na(snp_probe_data[, 2])))
snp_probe_data <- na.omit(snp_probe_data)
length(which(is.na(snp_probe_data[, 1])))
length(which(is.na(snp_probe_data[, 2])))
head(snp_probe_data)
colnames(snp_probe_data)
dim(snp_probe_data)

# Define variables for uncorrected data:
colnames(snp_probe_data)
str(snp_probe_data)
snp_probe_data[, snp_col] <- as.numeric(as.character(snp_probe_data[, snp_col]))
# TO DO: check whether to run as additive (numeric encoding) or logistic (factor encoding):
# snp_probe_data[, snp] <- as.factor(as.numeric(as.character(snp_probe_data[, snp])))
snp
probe
snp_probe_data[1:5, snp_col]
snp_probe_data[1:5, probe_col]
class(snp_probe_data)
snp_probe_data[, probe_col] <- as.numeric(as.character(snp_probe_data[, probe_col]))
str(snp_probe_data)
sapply(snp_probe_data, class)

# Add PCs to adjust for (corrected data):
snp_probe_PC_data <- merge(snp_probe_data, PC_data_t, by = 'row.names')
row.names(snp_probe_PC_data) <- snp_probe_PC_data[, 1]
snp_probe_PC_data <- snp_probe_PC_data[, -1]
snp_probe_PC_data[1:5, 1:5]
dim(snp_probe_PC_data)
# str(snp_probe_PC_data)

# Define variables for corrected data:
colnames(snp_probe_PC_data)
# str(snp_probe_PC_data)
snp_probe_PC_data[, snp_col] <- as.numeric(as.character(snp_probe_PC_data[, snp_col]))
# TO DO: check whether to run as additive (numeric encoding) or logistic (factor encoding):
# snp_probe_PC_data[, snp] <- as.factor(as.numeric(as.character(snp_probe_PC_data[, snp])))
snp_probe_PC_data[, probe_col] <- as.numeric(as.character(snp_probe_PC_data[, probe_col]))
# str(snp_probe_PC_data)
sapply(snp_probe_PC_data, class)
snp_probe_PC_data[1:5, 1:5]
########################

########################
if (grepl('[:]', as.character(snp))) {
  snp_col <- sprintf('x%s', gsub(":", "_", as.character(snp)))
} else {
  snp_col <- snp_col
}
snp_col

# Run linear regression model adjusting for PCs for individual SNP and probe:
# Uncorrected:
names(snp_probe_data)
names(snp_probe_data)[1] <- snp_col # rename if SNP name was eg 22:18906839 and now x22_18906839
names(snp_probe_data)
sapply(snp_probe_data, class)
snp_probe_data[1:5, ]
summary(snp_probe_data)
pass_formula <- as.formula(sprintf('%s ~ %s', probe_col, snp_col))
pass_formula
lm_fit <- lm(formula = pass_formula, data = snp_probe_data)
# lm_fit <- lm(formula = snp_probe_data$Plasma_CPMG_NMR_4_127146 ~ snp_probe_data$`x22:18906839`, data = snp_probe_data)
summary.lm(lm_fit)
# Test assumptions:
gvmodel <- gvlma(lm_fit)
summary(gvmodel)

# Corrected for PCs:
names(snp_probe_PC_data)
names(snp_probe_PC_data)[1] <- snp_col # rename if SNP name was eg 22:18906839 and now x22_18906839
sapply(snp_probe_PC_data, class)
pass_formula <- as.formula(sprintf('%s ~ .', probe_col, snp_col))
pass_formula
lm_fit_PCs <- lm(formula = pass_formula, data = snp_probe_PC_data)
summary.lm(lm_fit_PCs)
# Test assumptions:
gvmodel <- gvlma(lm_fit_PCs)
summary(gvmodel)

# Model comparison:
AIC(lm_fit, lm_fit_PCs)
anova(lm_fit_PCs, lm_fit)

# Comparison of MatrixEQTL results to lm results:
coefficients(lm_fit) # lm_fit_PCs$coefficients
str(lm_fit)
# Get coefficients (estimate, t-stat, p-value) from lm:
lm_summary <- summary.lm(lm_fit)
lm_summary$coefficients
lm_summary_coef <- as.data.frame(lm_summary$coefficients)
str(lm_summary_coef)
# Get for covariates lm:
lm_summary_PCs <- summary.lm(lm_fit_PCs)
lm_summary_PCs_coef <- as.data.frame(lm_summary_PCs$coefficients)
# Compare lm and MatrixEQTL:
lm_summary_coef[snp_col, ]
lm_summary_PCs_coef[snp_col, ]
# eQTL_data[which(eQTL_data[, SNP] == snp_col), ]
eQTL_data[which(eQTL_data[, SNP] == snp), ]

# TO DO: save to file
# file_name <- sprintf('%s_%s_%s_%s_%s.svg', snp, probe, gene_name, plot_file_name, plot_file_name_group)
# 
# cat(file = sprintf('counts_%s.txt', output_file_name),
#     append = FALSE)
########################


########################
# Labels:
names(eQTL_data)
# Passing reQTLs and eQTL files so check for column names:
if (names(eQTL_data)[1] == 'Probe_ID') {
  gene_name <- eQTL_data[which(eQTL_data[, 'Probe_ID', with = F] == probe), 'gene', with = F]
  } else {
    gene_name <- eQTL_data[which(eQTL_data[, 'gene', with = F] == probe), 'gene', with = F]
    }
gene_name <- as.character(gene_name[1])
gene_name

# TO DO: changing file names will break this:
eQTL_file_name
snp_file_name

plot_file_name <- eQTL_file_name
plot_file_name

plot_file_name_group <- strsplit(snp_file_name, split = '_')
# plot_file_name_group <- strsplit(plot_file_name_group[[1]][6], split = '[.]')[[1]][1]
plot_file_name_group <- plot_file_name_group[[1]][4]
plot_file_name_group

# plot_file_name <- strsplit(eQTL_file_name, split = '_')
# plot_file_name <- paste(plot_file_name[[1]][1], plot_file_name[[1]][2], plot_file_name[[1]][3], sep = '_')
# plot_file_name

plot_title <- sprintf('Effect of %s on levels of %s', snp, probe)
plot_title

# TO DO: add line plot for fitted values against genotype (average for each genotype):
# lines(lm_fit$fitted ~ snp, type = 'b', pch = 15, col = 'darkgrey') 

# Plot unadjusted values:
snp_probe_data[1:5, ]
ggplot(snp_probe_data, aes(as.factor(snp_probe_data[, snp_col]), snp_probe_data[, probe_col])) + 
  geom_jitter(colour='darkgrey', position=position_jitter(width=0.25)) + 
  geom_boxplot(outlier.size=0, alpha=0.6, fill='grey') + 
  ylab('Average unadjusted metabolite (datapoint, uncorrected) abundance levels') + 
  xlab('Genotype') + ggtitle(plot_title) + theme_minimal() +
  theme(text = element_text(size = 8))
file_name <- sprintf('uncorrected_%s_%s_%s.svg', snp, probe, plot_file_name)
file_name
ggsave(plot = last_plot(), filename = file_name)

# Plot SNP and probe associated adjusted with PCs:
snp_probe_PC_data[1:5, 1:5]
# snp_probe_PC_data[1:5, 1:15]
# TO DO: Check the adjustment is correct (check whether running SNP as numeric or factor above and here adjust lm coef accordingly):
# gex_corrected <- snp_probe_PC_data[, probe] - rowSums(coef(lm_fit_PCs)[-(1:4)] * snp_probe_PC_data[, 3:ncol(snp_probe_PC_data)])
# If SNP numeric:
gex_corrected <- snp_probe_PC_data[, probe_col] - rowSums(coef(lm_fit_PCs)[-(1:2)] * snp_probe_PC_data[, 3:ncol(snp_probe_PC_data)])
# Collect as data frame to pass to ggplot2:
gex_corrected_snp_df <- data.frame(expression = gex_corrected, genotype = factor(snp_probe_PC_data[, snp_col]))
ggplot(gex_corrected_snp_df, aes(genotype, expression)) +
  geom_jitter(colour="darkgrey", position=position_jitter(width=0.25)) +
  geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") + 
  ylab('Average adjusted metabolite (datapoint) abundance levels') + 
  xlab('Genotype') + ggtitle(plot_title) + theme_minimal() +
  theme(text = element_text(size = 8))
file_name <- sprintf('corrected_%s_%s_%s.svg', snp, probe, plot_file_name)
file_name
ggsave(plot = last_plot(), filename = file_name)

########################


########################
# TO DO: add legend, interpretation, method

########################

####################
# The end:
# Remove objects that are not necessary to save:

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
# save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next: 
####################
