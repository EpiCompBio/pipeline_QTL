#############################
# To be run after 02 normalisation of array data, pheno file processing
# Antonio J Berlanga-Taylor
# 25 Feb 2016
# BEST-D project differential expression 2, analysis with covariate adjustment
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_diff_expression_2",".txt", sep=""), open='a')
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
#load('R_session_saved_image_diff_expression.RData', verbose=T) # Contains subsetted array data
load('R_session_saved_image_pheno_file_check.RData', verbose=T)


# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_diff_expression_2', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:
#install.packages('ellipse')

library(limma)
library(ggplot2)
library(ellipse)
library(Hmisc)
library(splines)
library(plyr)
library(statmod)
#############################


#############################

# Read in the file with the sample membership information (group membership, replicate number, 
# treatment status, etc.) to be able to create a design and contrast (done above for PCA plots).

#Check dimensions between annotation file with meta-data (must have the same number of rows, otherwise
#errors downstream):
#TO DO: Raise error if not the same.

dim(membership_file_cleaned)
str(membership_file_cleaned)
head(membership_file_cleaned)
tail(membership_file_cleaned)
dim(normalised_filtered)
str(normalised_filtered)
dim(normalised_filtered_annotated)

# Sanity check:
# TO DO: Raise error and stop if false:
identical(row.names(membership_file_cleaned), colnames(normalised_filtered))
length(which(row.names(membership_file_cleaned) %in% colnames(normalised_filtered)))

# Load full phenotype data for covariates adjustment:
phenotype_data <- read.csv('BEST-D_phenotype_file_final.tsv', sep = '\t', 
                           header = TRUE, na.string = c('-99999', "", " ", "NA"))
dim(phenotype_data)
length(which(complete.cases(phenotype_data)))
#View(phenotype_data)
head(phenotype_data)
tail(phenotype_data)
summary(phenotype_data)
str(phenotype_data)
class(phenotype_data)
names(phenotype_data)
#############################

########################

########################
# TO DO: Adjusted for in assoc analysis in plink: 
#vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2
# male, 'vitd6', 'vitd12'

## Multi-level experiment: comparing within and across groups
# This analysis uses paired samples and account for placebo effect.
# Check limma user manul section 9.6.2 Many time points p.49 spline library if variables don't have a linear correlation.

# TO DO: continue from here Feb 22: re-run this as order of GEx file is different to array and indexing will be wrong. 
# TO DO: GEx pheno file also needs subsetting as kit IDs vs pt IDs are arranged differently (ie there are either NAs or 
# duplicated values after the merging of pheno files at the beginning of script).
# TO DO: clean up as this GEx file is no longer used. Check and compare with 03 and 03b.

# Set up labels and groups:
group <- factor(GEx_phenotype_data$group_membership)
levels(group)
count(group)
length(group)

# Pass covariates for adjustment here, consider spline package for non-linear if needed. Cut down to these variables
# looking at correlations (association script for BEST-D). Some variables still correlate (age and gender particularly).
# TO DO: I haven't assessed non-linear relationships.
# Variables adjusted for in association analysis (plink checks co-linearity an excludes some so not final list):
# vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2
# male, 'vitd6', 'vitd12'
incident_fracture <- factor(GEx_phenotype_data$incident_fracture)
levels(incident_fracture)
count(incident_fracture)

incident_resp_infection <- factor(GEx_phenotype_data$incident_resp_infection)
diabetes <- factor(GEx_phenotype_data$diabetes)
heart_disease <- factor(GEx_phenotype_data$heart_disease)
copd_asthma <- factor(GEx_phenotype_data$copd_asthma)
basemed_vitamind <- factor(GEx_phenotype_data$basemed_vitamind)
currsmoker <- factor(GEx_phenotype_data$currsmoker)
bmi0 <- GEx_phenotype_data$bmi0
calendar_age_ra <- GEx_phenotype_data$calendar_age_ra
season_randomisation_2 <- factor(GEx_phenotype_data$season_randomisation_2)
male <- factor(GEx_phenotype_data$male)

covariates_list <- cbind(incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, 
                         currsmoker, bmi0, calendar_age_ra, season_randomisation_2, male)

summary(covariates_list)
count(covariates_list)
covariates_list_df <- as.data.frame(covariates_list)
head(covariates_list_df)
which(!complete.cases(covariates_list_df))
covariates_list_df[244, ]
group[244]

# TO DO: haven't adjusted for cell count. Should be at least partially accounted for in paired analysis and randomisation. Check if
# I can pass top PCs for example.

# R command for limma for covariates would be: 
design_multilevel_covar <- model.matrix(~0 + group +
                                          incident_fracture +
                                          incident_resp_infection +
                                          diabetes +
                                          heart_disease +
                                          copd_asthma +
                                          basemed_vitamind +
                                          currsmoker +
                                          bmi0 +
                                          calendar_age_ra +
                                          season_randomisation_2 +
                                          male
)

# design_multilevel <- model.matrix(~0 + group) 
design_multilevel <- design_multilevel_covar
dim(design_multilevel)

which(!complete.cases(design_multilevel))


count(design_multilevel)
#View(design_multilevel)
dim(design_multilevel)
str(design_multilevel)
colnames(design_multilevel)
head(design_multilevel)
tail(design_multilevel)
summary(design_multilevel)

# Estimate the correlation between measurements made on the same subject. This takes several minutes:
# Without covariates:
# corfit <- duplicateCorrelation(normalised_filtered, design_multilevel, block = membership_file_cleaned$pt_id)
# corfit$consensus

# With covariates:
corfit_covar <- duplicateCorrelation(normalised_filtered, design_multilevel, block = GEx_phenotype_data$pt_id)
corfit_covar$consensus
corfit_covar$cor
corfit_covar$consensus.correlation
head(corfit_covar$atanh.correlations)
#?duplicateCorrelation

dim(design_multilevel)
dim(GEx_phenotype_data)
dim(normalised_filtered)

#Then this inter-subject correlation is input into the linear model fit:
fit_multilevel <- lmFit(normalised_filtered, design_multilevel, block = membership_file_cleaned$pt_id, 
                        correlation = corfit_covar$consensus)
fit_multilevel
fit2_test <- eBayes(fit_multilevel)
fit2_test
str(fit2_test)
fit2_test
topTable(fit2_test, adjust = 'BH', coef = c())

#Specify the comparisons between the experimental conditions:
cont.matrix_groups_before_v_after_multilevel <- makeContrasts(f4000vsb4000=groupfinal_4000-groupbaseline_4000, 
                                                              f2000vsb2000=groupfinal_2000-groupbaseline_2000,
                                                              fplacebovsbplacebo=groupfinal_placebo-groupbaseline_placebo,
                                                              levels=design_multilevel)
cont.matrix_groups_before_v_after_multilevel

cont.matrix_groups_before_v_after_time_multilevel <- makeContrasts(UI4000minusplacebo=(groupfinal_4000-groupbaseline_4000) - (groupfinal_placebo-groupbaseline_placebo), 
                                                                   UI2000minusplacebo=(groupfinal_2000-groupbaseline_2000) - (groupfinal_placebo-groupbaseline_placebo),
                                                                   UI4000minus2000=(groupfinal_4000-groupbaseline_4000) - (groupfinal_2000-groupbaseline_2000),
                                                                   UI4000minus2000minusplacebo=(groupfinal_4000-groupbaseline_4000) - (groupfinal_2000-groupbaseline_2000) - (groupfinal_placebo-groupbaseline_placebo),
                                                                   UI2000minus4000minusplacebo=(groupfinal_2000-groupbaseline_2000) - (groupfinal_4000-groupbaseline_4000) - (groupfinal_placebo-groupbaseline_placebo),
                                                                   levels=design_multilevel)
cont.matrix_groups_before_v_after_time_multilevel

cont.matrix_group_comparisons <- makeContrasts(f4000vsf2000=groupfinal_4000-groupfinal_2000, 
                                               f4000vsfplacebo=groupfinal_4000-groupfinal_placebo,
                                               f2000vsfplacebo=groupfinal_2000-groupfinal_placebo,
                                               levels=design_multilevel)
cont.matrix_group_comparisons

#Compute the contrasts and moderated t-tests:
fit2_multilevel <- contrasts.fit(fit_multilevel, cont.matrix_groups_before_v_after_multilevel)
fit2_multilevel <- eBayes(fit2_multilevel)

fit2_multilevel_time <- contrasts.fit(fit_multilevel, cont.matrix_groups_before_v_after_time_multilevel)
fit2_multilevel_time <- eBayes(fit2_multilevel_time)

fit2_multilevel_group_comparison <- contrasts.fit(fit_multilevel, cont.matrix_group_comparisons)
fit2_multilevel_group_comparison <- eBayes(fit2_multilevel_group_comparison)

# Extract differences for comparison of interest:
topTable(fit2_multilevel, adjust = 'BH')
topTable(fit2_multilevel, adjust = 'BH', coef = 'f4000vsb4000')
topTable(fit2_multilevel, adjust = 'BH', coef = 'f2000vsb2000')
topTable(fit2_multilevel, adjust = 'BH', coef = 'fplacebovsbplacebo')

topTable(fit2_multilevel_time, adjust = 'BH')
topTable(fit2_multilevel_time, adjust = 'BH', coef = 'UI4000minusplacebo')
topTable(fit2_multilevel_time, adjust = 'BH', coef = 'UI2000minusplacebo')
topTable(fit2_multilevel_time, adjust = 'BH', coef = 'UI4000minus2000')
topTable(fit2_multilevel_time, adjust = 'BH', coef = 'UI4000minus2000minusplacebo')
topTable(fit2_multilevel_time, adjust = 'BH', coef = 'UI2000minus4000minusplacebo')

topTable(fit2_multilevel_group_comparison, adjust = 'BH')
topTable(fit2_multilevel_group_comparison, adjust = 'BH', coef = 'f4000vsf2000')
topTable(fit2_multilevel_group_comparison, adjust = 'BH', coef = 'f4000vsfplacebo')
topTable(fit2_multilevel_group_comparison, adjust = 'BH', coef = 'f2000vsfplacebo')

order(fit2_multilevel_group_comparison$F.p.value)[1:10]

# Interpretation: There are no significant differences once covariates are taken into account, neither in final vs baseline
# comparisons or in contrasts when accounting for placebo. TO DO: check I'm not erasing the effect of Tx as raw comparisons of
# before and after are significant.
#############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for higher level analyses.
#############################