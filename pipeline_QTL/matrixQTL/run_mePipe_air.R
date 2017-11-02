#################
# Script to run mePipe
# Antonio J Berlanga-Taylor
# 20 May 2016
# Requires gene expression and genotype data
# Outputs plots and files to disk
#################

#################
# See:
# https://github.com/jknightlab/mePipe
#################


#################
##Set working directory and file locations and names of required inputs:
options(echo = TRUE)

# Working directory:
# setwd('/Users/antoniob/Desktop/Airwave_inflammation/results_3.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_",Sys.Date(),".txt", sep=""))
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))
getwd()

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
#load('R_session_saved_image_mePipe.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_mePipe','.RData', sep='')
R_session_saved_image
####################

#################
# Load libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("MatrixEQTL", "trio", "optparse", "XML", "snow"))
# library(devtools)
# biocLite('optparse')
# # devtools::install_github("humburg/Rsge")
# system("wget https://github.com/jknightlab/mePipe/releases/download/mePipe_v1.3.5/mePipe_1.3.5.tar.gz")
# install.packages("mePipe_1.3.5.tar.gz", repos=NULL)
library(mePipe)
library(ggplot2)
# library(data.table)
#################

#################
# Switch Rsge off:
sge.options(sge.use.cluster = FALSE) #sge.user.options=opt$sgeoptions)
# Also check the following as errors in cgat150
# make_option("--sgeoptions", default="-S /bin/bash -V",
# help="Options to pass to qsub (ignored unless `--cluster` is used. [default: %default])"),
# mePipe::
mePipe::getOptions()

#################


####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

geno_file <- as.character(args[1])
# geno_file <- 'cut_chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome.geno_matched.tsv'
# geno_file <- 'head_cut_chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome.geno_matched.tsv'

# Gene expression levels:
Gex_file <- as.character(args[2])
# Gex_file <- 'cut_Airwave_CPMG_Plasma.txt_matched.tsv'
# Gex_file <- 'head_cut_Airwave_CPMG_Plasma.txt_matched.tsv'

#PC_seq_to_test = as.character(args[5])
PC_seq_to_test = seq(0, 50, by = 5)

print(args)
##################
# TO DO: generate covariate PCs from GEx file first, then assoc cov to genotype, then filter these from run_mePipe run.
# chooseCov or runMe don't have these options (assocCov and filterCov) so run covariates.R, assocCov, filterCov and pass these to
# chooseCov
# default FDR threshold for PC exclusion is 0.001
#(command line does though, check this again to make running straightforward) 

#################
# Run mePipe

# Set up MatrixEQTL options:
useModel = 'linear' # modelANOVA or modelLINEAR or modelLINEAR_CROSS
# useModel = modelLINEAR_CROSS # Should work better for Tx specific eQTLs

# The p-value threshold determines which gene-SNP associations are saved.
threshold = 1e-05
# threshold = 0.1
cisThreshold = 0.001
# cisThreshold = 0.90

# cis distance definition:
cis = 1e+06

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

# Next run:
####################

