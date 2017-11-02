#############################
# To be run after 02 normalisation of array data, pheno file processing
# Antonio J Berlanga-Taylor
# 01 March 2016
# Differential expression - GAinS II cohort for pipeline sanity checking
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/public_data_tests.dir/GAinS_2.dir/analysis.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_diff_expression",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is project specific. Check ways of making count comparisons.

# Load results from 02_microarrayxxx file, saved as RData object:
# Re-load a previous R session, data and objects:
#load('R_session_saved_image_probe_filtering.RData', verbose=T)
load('R_session_saved_image_pheno_file_check.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_diff_expression', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:
#install.packages('ellipse')
#install.packages("statmod")


library(limma)
library(ggplot2)
library(ellipse)
library(Hmisc)
library(splines)
library(plyr)
library(statmod)
library(illuminaHumanv4.db)
#############################


#############################

# Read in the file with the sample membership information (group membership, replicate number, 
# treatment status, etc.) to be able to create a design and contrast (done above for PCA plots).

# Sanity check:
# TO DO: Raise error and stop if false:
identical(row.names(membership_file_cleaned), colnames(normalised_filtered))
length(which(row.names(membership_file_cleaned) %in% colnames(normalised_filtered)))

head(membership_file_cleaned)
tail(membership_file_cleaned)
str(membership_file_cleaned)

dim(normalised_filtered)
normalised_filtered[1:5, 1:5]

#############################

#############################

# 6) Differential gene expression analysis
# Linear modelling using weights

# Analyze using linear models in limma. A model is fit for every gene and an empirical Bayes method
# moderates the standard errors of the estimated log-fold changes

# Requires one or two matrices to be specified. A design matrix which has each array/sample per row, 
# and coefficients that describe the RNA sources in each column (eg group, treatment, batch/plate, etc.).
# A contrast matrix can then be used to group the coeffients for comparisons.
# Input data is the ExpressionSet (eset) or the EList class.
# Main functions: lm(), eBayes(), topTable(), etc.

######################### 
# b) Two group comparison of treated vs untreated: 

# Check NAs:
summary(membership_file_cleaned$Characteristics.sepsis.response.signature.group.)
to_remove <- which(is.na(membership_file_cleaned$Characteristics.sepsis.response.signature.group.))
to_remove
normalised_filtered_groups <- normalised_filtered[, -to_remove]
membership_file_cleaned_groups <- membership_file_cleaned[-to_remove, ]
identical(row.names(membership_file_cleaned), colnames(normalised_filtered))

#Compare samples based on groups:
group <- factor(membership_file_cleaned_groups$Characteristics.sepsis.response.signature.group., levels = c('group 1', 'group 2'))
count(group)


#Define design:
design_by_group <- model.matrix(~group)
head(design_by_group)
tail(design_by_group)
count(design_by_group)
count(membership_file_cleaned_groups$Characteristics.sepsis.response.signature.group.)
design_by_group

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_by_group <- lmFit(normalised_filtered_groups, design_by_group)
fit_by_group
names(fit_by_group)
colnames(fit_by_group)
head(row.names(fit_by_group))

fit_by_group_2 <- eBayes(fit_by_group)
fit_by_group_2
head(fit_by_group_2$coefficients)

topTable_groups <- topTable(fit_by_group_2, adjust = 'BH', number = Inf)
count(topTable_groups$adj.P.Val < 10e-3)
count(topTable_groups$adj.P.Val < 0.05 & 2^(topTable_groups$logFC) > 1.49)
count(topTable_groups$adj.P.Val < 0.05 & 2^(topTable_groups$logFC) < 0.68)
count(topTable_groups$adj.P.Val < 0.05 & topTable_groups$logFC > abs(0.57))

# Interpretation: Many significant differences between groups
###################

###################
## Add annotations to topTable:
head(topTable_groups)
illumina_ids  <-  as.character(rownames(topTable_groups))
head(illumina_ids)
length(illumina_ids)

probes_by_ENTREZID_RE  <- unlist(mget(illumina_ids, illuminaHumanv4ENTREZREANNOTATED, ifnotfound = NA))
head(probes_by_ENTREZID_RE)
summary(probes_by_ENTREZID_RE)
str(probes_by_ENTREZID_RE)
count(!is.na(probes_by_ENTREZID_RE))
class(probes_by_ENTREZID_RE)
probes_by_ENTREZID_RE_df <- as.data.frame(probes_by_ENTREZID_RE)
head(probes_by_ENTREZID_RE_df)
head(t(probes_by_ENTREZID_RE))

probes_by_symbol <- unlist(mget(illumina_ids, illuminaHumanv4SYMBOLREANNOTATED, ifnotfound = NA))
summary(probes_by_symbol)
count(!is.na(probes_by_symbol))
probes_by_symbol_df <- as.data.frame(probes_by_symbol)
head(probes_by_symbol_df)

topTable_groups_annot <-  merge(topTable_groups, probes_by_symbol_df)
t

###################

###################
# Write files to disk:
write.table(x=topTable_groups, sep='\t', quote = FALSE,
            col.names = NA, row.names = TRUE,
            file='full_topTable_groups.txt')


#############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for higher level analyses.
#############################
