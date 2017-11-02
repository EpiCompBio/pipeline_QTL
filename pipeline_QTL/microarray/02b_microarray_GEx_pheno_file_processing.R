#############################
# Antonio J Berlanga-Taylor
# 22 Feb 2016
# BEST-D project phenotype/metadata file processing: checks order between array and pheno file and adds 
# columns to pheno file for diff exp analysis for BEST-D.
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_pheno_file_check",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is project specific.

# Load results from 02_microarrayxxx file, saved as RData object:
# Re-load a previous R session, data and objects:
load('R_session_saved_image_probe_filtering.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_pheno_file_check', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:

library(plyr)
#############################

#############################

# 4) Experimental design matrix specification for differential expression analysis in limma.

# Read in the file with the sample membership information (group membership, replicate number, 
# treatment status, etc.) to be able to create a design and contrast (done above for PCA plots).

# Pheno file should have rows for samples and columns for labels/phenotype data. Array file has
# rows for probe expression data and each sample in a column. Pheno file rows must match exactly array file columns.
# If there is a mismatch in numbers it will error (due to NAs for example), but if in wrong order it 
# will not error or give warnings, so take care!

#TO DO: Raise error if numbers don't match to avoid further processing:

dim(membership_file_cleaned)
str(membership_file_cleaned)
head(membership_file_cleaned)
head(row.names(membership_file_cleaned))
head(colnames(normalised_filtered))
head(colnames(normalised_filtered_annotated))
tail(membership_file_cleaned)
dim(normalised_filtered)
str(normalised_filtered)
dim(normalised_filtered_annotated)

# TO DO: Move this section to another file (eg metadata processing) and keep row/column names 
# checks here:
# Add labels for group membership (ie baseline_2000 only):
membership_file_cleaned['group_membership'] <- ifelse(membership_file_cleaned$arm == 0 & membership_file_cleaned$visit_type == 'Randomisation', 'baseline_4000',
                                                      ifelse(membership_file_cleaned$arm == 0 & membership_file_cleaned$visit_type == 'FinalVisit', 'final_4000',
                                                             ifelse(membership_file_cleaned$arm == 1 & membership_file_cleaned$visit_type == 'Randomisation', 'baseline_2000', 
                                                                    ifelse(membership_file_cleaned$arm == 1 & membership_file_cleaned$visit_type == 'FinalVisit', 'final_2000',
                                                                           ifelse(membership_file_cleaned$arm == 2 & membership_file_cleaned$visit_type == 'Randomisation', 'baseline_placebo', 
                                                                                  ifelse(membership_file_cleaned$arm == 2 & membership_file_cleaned$visit_type == 'FinalVisit', 'final_placebo',
                                                                                         NA))))))
count(membership_file_cleaned$group_membership)
head(membership_file_cleaned)
tail(membership_file_cleaned)

#membership_file_cleaned[which(membership_file_cleaned$pt_id == 103485), ]

# Add treatment label to placebo final:
membership_file_cleaned$treatment <- ifelse(membership_file_cleaned$treatment == 'untreated' & membership_file_cleaned$group_membership == 'final_placebo', 'treated_placebo', membership_file_cleaned$treatment)
head(membership_file_cleaned)
tail(membership_file_cleaned)
count(membership_file_cleaned$treatment)

# Add labels for all treated (2000 + 4000) and placebo only:
membership_file_cleaned['all_treated'] <- ifelse(membership_file_cleaned$treatment == 'treated_2000' | membership_file_cleaned$treatment == 'treated_4000',
                                                 'treated_2000+4000', membership_file_cleaned$treatment)
head(membership_file_cleaned)
tail(membership_file_cleaned)
count(membership_file_cleaned$all_treated)

membership_file_cleaned$two_group_Tx <- ifelse(membership_file_cleaned$all_treated == 'treated_placebo', 'untreated', membership_file_cleaned$all_treated)
count(membership_file_cleaned$two_group_Tx)

# Check order of pheno file and array file are the same:
identical(row.names(membership_file_cleaned), colnames(normalised_filtered))
length(which(row.names(membership_file_cleaned) %in% colnames(normalised_filtered)))

# Re-order so that phenotype matches the array:
head(normalised_filtered)[1:5, 1:5]
normalised_filtered_ordered <- normalised_filtered[, order(colnames(normalised_filtered))]
membership_file_cleaned <- membership_file_cleaned[order(row.names(membership_file_cleaned)), ]
head(normalised_filtered_ordered)[1:5, 1:5]
head(membership_file_cleaned)[1:5, 1:5]

#identical(normalised_filtered[, '120000009'], normalised_filtered_ordered[, 2])
identical(row.names(membership_file_cleaned), colnames(normalised_filtered_ordered))
length(which(row.names(membership_file_cleaned) %in% colnames(normalised_filtered_ordered)))

normalised_filtered <- normalised_filtered_ordered
identical(row.names(membership_file_cleaned), colnames(normalised_filtered))

dim(membership_file_cleaned)
dim(normalised_filtered)

############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for differential expression analysis.
#############################
