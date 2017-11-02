#############################
# To be run after 02 normalisation of array data, pheno file processing
# Antonio J Berlanga-Taylor
# 12 Feb 2016
# BEST-D project differential expression
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")
# setwd('/Users/antoniob/Desktop/BEST_D_12_FEB.DIR/results_1/')

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

# Get functions from other scripts (eg to add annotations to topTable results):
source('/Users/antoniob/Desktop/BEST_D_12_FEB.DIR/scripts/gene_expression_functions.R')
#source('/ifs/devel/antoniob/projects/BEST-D/gene_expression_functions.R')
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
# TO DO: Adjusted for in assoc analysis in plink: 
#vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2
# male, 'vitd6', 'vitd12'

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
## Compare only baseline samples:
# Experimental design matrix specification
#Define design and factors to constrast from annotation file:
head(membership_file_cleaned)
group <- factor(membership_file_cleaned$group_membership)
str(group)
count(group)
head(group)
head(membership_file_cleaned$group_membership)
head(rownames(membership_file_cleaned))
head(colnames(normalised_filtered))

design <- model.matrix(~0+group)
head(design)
head(membership_file_cleaned)
tail(design)
tail(membership_file_cleaned)
design

#Run linear model and set contrasts:
fit <- lmFit(normalised_filtered, design)
head(row.names(fit))
colnames(fit)
str(fit)
names(fit)

cont.matrix <- makeContrasts(b4000vsb2000=groupbaseline_4000-groupbaseline_2000, 
                             b4000vsbplacebo=groupbaseline_4000-groupbaseline_placebo,
                             b2000vsbplacebo=groupbaseline_2000-groupbaseline_placebo,
                             levels=design)
cont.matrix

#Obtain differentially expressed genes based on contrasted factors:
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
fit2$contrasts
fit2
names(fit2)
head(fit2$contrasts)

#Get results and plot:
topTable(fit2, adjust="BH")
topTable(fit2, adjust="BH", coef = 'b4000vsbplacebo')
topTable(fit2, adjust="BH", coef = 'b2000vsbplacebo')
topTable(fit2, adjust="BH", coef = 'b4000vsb2000')
#volcanoplot(fit2)
results <- decideTests(fit2) 
results
vennDiagram(results)

# Interpretation: baseline samples do not differ from each other, as expected.

###############
## Compare only final samples:

#Define design and factors to constrast from annotation file. Use the same definitions as for baseline:
# group <- factor(membership_file_cleaned$group_membership)
# design <- model.matrix(~0+group)

#Run linear model and set contrasts:
fit_final <- lmFit(normalised_filtered, design)
cont.matrix_final <- makeContrasts(f4000vsf2000=groupfinal_4000-groupfinal_2000, 
                             f4000vsfplacebo=groupfinal_4000-groupfinal_placebo,
                             f2000vsfplacebo=groupfinal_2000-groupfinal_placebo,
                             levels=design)
cont.matrix_final

#Obtain differentially expressed genes based on contrasted factors:
fit2_final <- contrasts.fit(fit_final, cont.matrix_final)
fit2_final <- eBayes(fit2_final)
fit2_final$contrasts
fit2_final

#Get results and plot:
topTable(fit2_final, adjust="BH") # This provides ANOVA F tests on the groups
topTable(fit2_final, adjust="BH", coef = 'f4000vsfplacebo') # This provides p-values on the two comparison only
topTable(fit2_final, adjust="BH", coef = 'f4000vsfplacebo', sort.by = 'logFC')
topTable(fit2_final, adjust="BH", coef = 'f4000vsfplacebo', sort.by = 'p')
topTable(fit2_final, adjust="BH", coef = 'f2000vsfplacebo')
topTable(fit2_final, adjust="BH", coef = 'f4000vsf2000')

#volcanoplot(fit2)
results_final <- decideTests(fit2_final) 
results_final
vennDiagram(results_final)

# Interpretation: final samples (after treatment, eg final_2000 vs final_placebo) do not differ from each other,
# there are no significant differences, not as expected.

###############
## Compare before vs after as groups:
# Define the questions of interest:
# Before and after comparisons for three groups:
# Which genes respond to treatment after stimulation at high dose: 4000IU Final Visit vs 4000IU randomisation (baseline)
# Which genes respond to treatment after stimulation at low dose: 2000IU Final Visit vs 2000IU randomisation (baseline)
# Which genes are likely due to noise:  Placebo Final Visit vs Placebo randomisation (baseline)

#Define design and factors to constrast from annotation file. Use the same definitions as for baseline:
# group <- factor(membership_file_cleaned$group_membership)
# design <- model.matrix(~0+group)

#Run linear model and set contrasts:
head(design)
count(design)
head(group)
head(membership_file_cleaned$group_membership)

tail(group)
count(group)
length(group)

fit_groups_before_v_after <- lmFit(normalised_filtered, design)
cont.matrix_groups_before_v_after <- makeContrasts(f4000vsb4000=groupfinal_4000-groupbaseline_4000, 
                                   f2000vsb2000=groupfinal_2000-groupbaseline_2000,
                                   fplacebovsbplacebo=groupfinal_placebo-groupbaseline_placebo,
                                   levels=design)
cont.matrix_groups_before_v_after

#Obtain differentially expressed genes based on contrasted factors:
fit2_groups_before_v_after <- contrasts.fit(fit_groups_before_v_after, cont.matrix_groups_before_v_after)
fit2_groups_before_v_after <- eBayes(fit2_groups_before_v_after)
fit2_groups_before_v_after$contrasts
fit2_groups_before_v_after

#Get results and plot:
topTable_fit2_groups_before_v_after <- topTable(fit2_groups_before_v_after, adjust="BH", number = Inf) # This provides ANOVA F tests on the groups
topTable_fit2_f4000vsb4000 <- topTable(fit2_groups_before_v_after, adjust="BH", coef = 'f4000vsb4000', number = Inf) # This provides p-values on the two comparison only
topTable(fit2_groups_before_v_after, adjust="BH", coef = 'f4000vsb4000', sort.by = 'logFC')
topTable(fit2_groups_before_v_after, adjust="BH", coef = 'f4000vsb4000', sort.by = 'p')
topTable_fit2_f2000vsb2000 <- topTable(fit2_groups_before_v_after, adjust="BH", coef = 'f2000vsb2000', number = Inf)
topTable_fit2_fplacebovsbplacebo <- topTable(fit2_groups_before_v_after, adjust="BH", coef = 'fplacebovsbplacebo', number = Inf)

# Extract results for some comparisons:
head(topTable_fit2_f4000vsb4000)
head(topTable_fit2_fplacebovsbplacebo)
head(topTable_fit2_f2000vsb2000)
topTable_fit2_f2000vsb2000$FC <- 2^topTable_fit2_f2000vsb2000$logFC
str(topTable_fit2_f2000vsb2000)

range(topTable_fit2_f2000vsb2000$AveExpr)
range(topTable_fit2_f2000vsb2000$P.Value)
range(topTable_fit2_f2000vsb2000$adj.P.Val)
summary(topTable_fit2_f2000vsb2000$adj.P.Val < 0.05)
range(topTable_fit2_f2000vsb2000$B)
range(topTable_fit2_f2000vsb2000$FC)

topTable_fit2_f4000vsb4000$FC <- 2^topTable_fit2_f4000vsb4000$logFC
range(topTable_fit2_f4000vsb4000$AveExpr)
range(topTable_fit2_f4000vsb4000$P.Value)
range(topTable_fit2_f4000vsb4000$adj.P.Val)
summary(topTable_fit2_f4000vsb4000$adj.P.Val < 0.05)
range(topTable_fit2_f4000vsb4000$B)
range(topTable_fit2_f4000vsb4000$FC)

topTable_fit2_fplacebovsbplacebo$FC <- 2^topTable_fit2_fplacebovsbplacebo$logFC
range(topTable_fit2_fplacebovsbplacebo$AveExpr)
range(topTable_fit2_fplacebovsbplacebo$P.Value)
range(topTable_fit2_fplacebovsbplacebo$adj.P.Val)
summary(topTable_fit2_fplacebovsbplacebo$adj.P.Val < 0.05)
range(topTable_fit2_fplacebovsbplacebo$B)
range(topTable_fit2_fplacebovsbplacebo$FC)


# Basic volcano plot
# From: http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html
png('volcano_plots_before_vs_after.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
plot(topTable_fit2_fplacebovsbplacebo$FC, -log10(topTable_fit2_fplacebovsbplacebo$adj.P.Val), pch=20, main="Volcano plot: placebo", abline(h=1.30103, v=c(0.9, 1.1)))
plot(topTable_fit2_f2000vsb2000$FC, -log10(topTable_fit2_f2000vsb2000$adj.P.Val), pch=20, main="Volcano plot: low dose", abline(h=1.30103, v=c(0.9, 1.1)))
plot(topTable_fit2_f4000vsb4000$FC, -log10(topTable_fit2_f4000vsb4000$adj.P.Val), pch=20, main="Volcano plot: high dose", abline(h=1.30103, v=c(0.9, 1.1)))
par(mfrow=c(1,1))
dev.off()


png('volcano_plots_fit2_groups_before_v_after.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))

volcanoplot(fit2_groups_before_v_after, coef = 'f2000vsb2000', highlight=0, #names=fit2_groups_before_v_after_time$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='High_dose_before_and_after')

volcanoplot(fit2_groups_before_v_after, coef = 'f2000vsb2000', highlight=0, #names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Low_dose_before_and_after')

volcanoplot(fit2_groups_before_v_after, coef='fplacebovsbplacebo', highlight=0, #names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Placebo_before_and_after')

par(mfrow=c(1,1))
dev.off()

results_groups_before_v_after <- decideTests(fit2_groups_before_v_after) 
results_groups_before_v_after

png('vennDiagram_results_groups_before_v_after_time.png', width = 3, height = 3, units = 'in', res = 300)
vennDiagram(results_groups_before_v_after)
dev.off()

# Estimate correlation from gene expression profiles:
ncol(fit2_groups_before_v_after$coefficients)
colnames(fit2_groups_before_v_after$coefficients)
high_v_placebo_group_corr <- genas(fit2_groups_before_v_after, coef=c(1, 3))
high_v_placebo_group_corr
high_v_low_group_corr <- genas(fit2_groups_before_v_after, coef=c(2, 3))
high_v_low_group_corr

# Interpretation: before vs after comparison by groups (not individuals) do not show differences.

############### 
## Use the same design as above and adjust for group comparisons with before/after but for time (ie placebo) (factorial design):
# I've run a similar analysis below but using pairing (which should be more powerful than this).
# Compare the difference of differences (ie interaction terms): 
# Differences due to high dose, ie 4000IU minus Placebo: (0.FinalVisit vs 0.Randomisation) - (2.FinalVisit vs 2.Randomisation)
# Differences due to low dose, ie 2000IU minus Placebo: (1.FinalVisit vs 1.Randomisation) - (2.FinalVisit vs 2.Randomisation)
# Differences due to high dose only, ie 4000IU minus 2000IU: (0.FinalVisit vs 0.Randomisation) - (1.FinalVisit vs 1.Randomisation) 
# Differences truly due to high dose only, ie (4000IU minus 2000IU minus Placebo:
#(0.FinalVisit vs 0.Randomisation) - (1.FinalVisit vs 1.Randomisation) - (2.FinalVisit vs 2.Randomisation)

# TO DO/CHECK: not sure limma is actually adjusting for placebo group:

# Define contrasts, use design as above:
cont.matrix_groups_before_v_after_time <- makeContrasts(UI4000minusplacebo=(groupfinal_4000-groupbaseline_4000) - (groupfinal_placebo-groupbaseline_placebo), 
                                                   UI2000minusplacebo=(groupfinal_2000-groupbaseline_2000) - (groupfinal_placebo-groupbaseline_placebo),
                                                   UI4000minus2000=(groupfinal_4000-groupbaseline_4000) - (groupfinal_2000-groupbaseline_2000),
                                                   UI4000minus2000minusplacebo=(groupfinal_4000-groupbaseline_4000) - (groupfinal_2000-groupbaseline_2000) - (groupfinal_placebo-groupbaseline_placebo),
                                                   UI2000minus4000minusplacebo=(groupfinal_2000-groupbaseline_2000) - (groupfinal_4000-groupbaseline_4000) - (groupfinal_placebo-groupbaseline_placebo),
                                                   levels=design)
cont.matrix_groups_before_v_after_time

#Obtain differentially expressed genes based on contrasted factors:
fit2_groups_before_v_after_time <- contrasts.fit(fit_groups_before_v_after, cont.matrix_groups_before_v_after_time)
fit2_groups_before_v_after_time <- eBayes(fit2_groups_before_v_after_time)
fit2_groups_before_v_after_time$contrasts


#Get results and plot:
topTable_fit2_groups_before_v_after_time <- topTable(fit2_groups_before_v_after_time, adjust="BH", number = Inf) # This provides ANOVA F tests on the groups
topTable_fit2_UI4000minusplacebo <- topTable(fit2_groups_before_v_after_time, adjust="BH", coef = 'UI4000minusplacebo', number = Inf) # This provides p-values on the two comparison only
topTable(fit2_groups_before_v_after_time, adjust="BH", coef = 'UI4000minusplacebo', sort.by = 'logFC', number = Inf)
topTable(fit2_groups_before_v_after_time, adjust="BH", coef = 'UI4000minusplacebo', sort.by = 'p', number = Inf)
topTable_fit2_UI2000minusplacebo <- topTable(fit2_groups_before_v_after_time, adjust="BH", coef = 'UI2000minusplacebo', number = Inf)
topTable_fit2_UI4000minus2000 <- topTable(fit2_groups_before_v_after_time, adjust="BH", coef = 'UI4000minus2000', number = Inf)
topTable_fit2_UI4000minus2000minusplacebo <- topTable(fit2_groups_before_v_after_time, adjust="BH", coef = 'UI4000minus2000minusplacebo', number = Inf)
topTable_fit2_UI2000minus4000minusplacebo <- topTable(fit2_groups_before_v_after_time, adjust="BH", coef = 'UI2000minus4000minusplacebo', number = Inf)

head(topTable_fit2_UI4000minusplacebo)
head(topTable_fit2_UI2000minusplacebo)
head(topTable_fit2_UI4000minus2000)
head(topTable_fit2_UI4000minus2000minusplacebo)
head(topTable_fit2_UI2000minus4000minusplacebo)

# Get fold changes for each:
topTable_fit2_UI4000minusplacebo$FC <- 2^topTable_fit2_UI4000minusplacebo$logFC
topTable_fit2_UI2000minusplacebo$FC <- 2^topTable_fit2_UI2000minusplacebo$logFC
topTable_fit2_UI4000minus2000$FC <- 2^topTable_fit2_UI4000minus2000$logFC
topTable_fit2_UI4000minus2000minusplacebo$FC <- 2^topTable_fit2_UI4000minus2000minusplacebo$logFC
topTable_fit2_UI2000minus4000minusplacebo$FC <- 2^topTable_fit2_UI2000minus4000minusplacebo$logFC

head(topTable_fit2_UI4000minusplacebo)
range(topTable_fit2_UI4000minusplacebo$AveExpr)
range(topTable_fit2_UI4000minusplacebo$P.Value)
range(topTable_fit2_UI4000minusplacebo$adj.P.Val)
summary(topTable_fit2_UI4000minusplacebo$adj.P.Val < 0.05)
range(topTable_fit2_UI4000minusplacebo$B)
range(topTable_fit2_UI4000minusplacebo$FC)


# Plot the linear model fit for diagnostics, this plots gene-wise residual standard deviations against average log-expression
# to identify mean-variance trends (sigma vs A plot). ?plotSA:
png('sigma_vs_A_plot_lm_factorial_design.png', width = 12, height = 12, units = 'in', res = 300)
plotSA(fit2_groups_before_v_after_time)
dev.off()

# Basic volcano plot
png('volcano_plots_factorial_design.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,2))
plot(topTable_fit2_UI4000minusplacebo$FC, -log10(topTable_fit2_UI4000minusplacebo$adj.P.Val), pch=20, main="UI4000minusplacebo", abline(h=1.30103, v=c(0.9, 1.1)))
plot(topTable_fit2_UI2000minusplacebo$FC, -log10(topTable_fit2_UI2000minusplacebo$adj.P.Val), pch=20, main="UI2000minusplacebo", abline(h=1.30103, v=c(0.9, 1.1)))
plot(topTable_fit2_UI4000minus2000$FC, -log10(topTable_fit2_UI4000minus2000$adj.P.Val), pch=20, main="UI4000minus2000", abline(h=1.30103, v=c(0.9, 1.1)))
plot(topTable_fit2_UI4000minus2000minusplacebo$FC, -log10(topTable_fit2_UI4000minus2000minusplacebo$adj.P.Val), pch=20, main="UI4000minus2000minusplacebo", abline(h=1.30103, v=c(0.9, 1.1)))
plot(topTable_fit2_UI2000minus4000minusplacebo$FC, -log10(topTable_fit2_UI2000minus4000minusplacebo$adj.P.Val), pch=20, main="UI2000minus4000minusplacebo", abline(h=1.30103, v=c(0.9, 1.1)))
par(mfrow=c(1,1))
dev.off()


# Plot all comparisons:
png('volcano_plots_fit2_groups_before_v_after_time.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,2))

volcanoplot(fit2_groups_before_v_after_time, coef='UI4000minusplacebo', highlight=0, #names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='High_dose_minus_Placebo')

volcanoplot(fit2_groups_before_v_after_time, coef='UI2000minusplacebo', highlight=0, #names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Low_dose_minus_Placebo')

volcanoplot(fit2_groups_before_v_after_time, coef='UI4000minus2000', highlight=0, #names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='High_dose_minus_low_dose')

volcanoplot(fit2_groups_before_v_after_time, coef='UI4000minus2000minusplacebo', highlight=0, #names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Truly_high_dose_only')

volcanoplot(fit2_groups_before_v_after_time, coef='UI2000minus4000minusplacebo', highlight=0, #names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Truly_low_dose_only')

par(mfrow=c(1,1))
dev.off()


# Classify results:
results_groups_before_v_after_time <- decideTests(fit2_groups_before_v_after_time) 
results_groups_before_v_after_time

png('vennDiagram_results_groups_before_v_after_time.png', width = 3, height = 3, units = 'in', res = 300)
vennDiagram(results_groups_before_v_after_time)
dev.off()

# Estimate correlation from gene expression profiles:
ncol(fit2_groups_before_v_after_time$coefficients)
colnames(fit2_groups_before_v_after_time$coefficients)
genas(fit2_groups_before_v_after_time, coef=c(1, 2))
genas(fit2_groups_before_v_after_time, coef=c(4, 5))


# Interpretation: Accounting for time (ie placebo before/after) shows no significant probes in treated samples. Genas function
# for correlation between comparisons does not make sense though. It does for the group comparisons above (final treated vs baseline).

########################

##############
# b) Two group comparison of treated vs untreated: 

#Compare samples based on treatment:
treatment <- factor(membership_file_cleaned$two_group_Tx, levels = c('untreated', 'treated_2000+4000'))
count(treatment)

#Define design and set contrasts:
design_by_treatment <- model.matrix(~treatment)
head(design_by_treatment)
tail(design_by_treatment)
count(design_by_treatment)
count(membership_file_cleaned$two_group_Tx)
design_by_treatment

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_by_treatment <- lmFit(normalised_filtered, design_by_treatment)
fit_by_treatment
names(fit_by_treatment)
colnames(fit_by_treatment)
head(row.names(fit_by_treatment))

fit_by_treatment_2 <- eBayes(fit_by_treatment)
fit_by_treatment_2
head(fit_by_treatment_2$coefficients)

topTable(fit_by_treatment_2, adjust = 'BH')

# Interpretation: Broad treated vs untreated shows no differences, so straight up stimulated vs basal is negative.
##############


#####################################
##Compare paired samples (eg based on 'SibShip', here using pt_id) for before vs after according to treatment
# Define factors to constrast from annotation file:
#Pairing:
head(membership_file_cleaned)
tail(membership_file_cleaned)
# kit_id appears once with corresponding pt_id twice (baseline + final visit):
which(membership_file_cleaned == 104018)
membership_file_cleaned[c(275, 560), ]

pairing <- factor(membership_file_cleaned$pt_id)
head(pairing)
length(pairing)

treatment <- factor(membership_file_cleaned$treatment, levels=c('untreated', 'treated_placebo', 'treated_2000', 'treated_4000'))
head(treatment)
summary(treatment)

#Define design and set contrasts:
design_all_pairs <- model.matrix(~pairing+treatment)
head(design_all_pairs)[1:5,1:5]
tail(design_all_pairs)[, -1]

dim(design_all_pairs)
#View(design_all_pairs)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_all_pairs <- lmFit(normalised_filtered, design_all_pairs)
fit_all_pairs_2 <- eBayes(fit_all_pairs)
dim(fit_all_pairs_2)
str(fit_all_pairs_2)
summary(fit_all_pairs_2)
names(fit_all_pairs_2)
head(fit_all_pairs_2$coefficients)#[1:5, 1:5]
head(fit_all_pairs_2$lods)#[1:5, 1:5]
head(fit_all_pairs_2$cov.coefficients)#[1:5, 1:5]
fit_all_pairs_2$coefficients

topTable(fit_all_pairs_2, adjust='BH')
topTable(fit_all_pairs_2, coef="treatmenttreated_4000", adjust='BH')
topTable(fit_all_pairs_2, coef="treatmenttreated_2000", adjust='BH')
topTable(fit_all_pairs_2, coef="treatmenttreated_placebo", adjust='BH')

topTable_pairing_high_dose <- topTable(fit_all_pairs_2, coef="treatmenttreated_4000", adjust='BH', n=Inf)
topTable_pairing_low_dose <- topTable(fit_all_pairs_2, coef="treatmenttreated_2000", adjust='BH', n=Inf)
topTable_pairing_placebo <- topTable(fit_all_pairs_2, coef="treatmenttreated_placebo", adjust='BH', n=Inf)


head(topTable_pairing_high_dose)
head(topTable_pairing_low_dose)
head(topTable_pairing_placebo)

count(topTable_pairing_high_dose$adj.P.Val < 0.05)
count(topTable_pairing_low_dose$adj.P.Val < 0.05)
count(topTable_pairing_placebo$adj.P.Val < 0.05)

#Get FC from logFC:
# TO DO: check this is correct:
topTable_pairing_high_dose$FC <- 2^topTable_pairing_high_dose$logFC
topTable_pairing_low_dose$FC <- 2^topTable_pairing_low_dose$logFC
topTable_pairing_placebo$FC <- 2^topTable_pairing_placebo$logFC

# TO DO: complete subsetting for overlaps between before/after groups. Get gene symbols.
range(topTable_pairing_high_dose$AveExpr)
range(topTable_pairing_high_dose$P.Value)
range(topTable_pairing_high_dose$adj.P.Val)
summary(topTable_pairing_high_dose$adj.P.Val < 0.05)
range(topTable_pairing_high_dose$B)
range(topTable_pairing_high_dose$FC)
count(topTable_pairing_high_dose$adj.P.Val < 0.05 & topTable_pairing_high_dose$FC > 1.1)
count(topTable_pairing_high_dose$adj.P.Val < 0.05 & topTable_pairing_high_dose$FC < 0.9)
high_downReg <- which(topTable_pairing_high_dose$adj.P.Val < 0.05 & topTable_pairing_high_dose$FC < 0.9)
#high_downReg_symbols <- topTable_pairing_high_dose[high_downReg, ]
#high_downReg_symbols$SYMBOL

range(topTable_pairing_low_dose$AveExpr)
range(topTable_pairing_low_dose$P.Value)
range(topTable_pairing_low_dose$adj.P.Val)
summary(topTable_pairing_low_dose$adj.P.Val < 0.05)
range(topTable_pairing_low_dose$B)
range(topTable_pairing_low_dose$FC)
count(topTable_pairing_low_dose$adj.P.Val < 0.05 & topTable_pairing_low_dose$FC > 1.1)
count(topTable_pairing_low_dose$adj.P.Val < 0.05 & topTable_pairing_low_dose$logFC > 0.1365)
count(topTable_pairing_low_dose$adj.P.Val < 0.05 & topTable_pairing_low_dose$FC < 0.9)
count(topTable_pairing_low_dose$adj.P.Val < 0.05 & topTable_pairing_low_dose$logFC < -0.15)
low_downReg <- which(topTable_pairing_low_dose$adj.P.Val < 0.05 & topTable_pairing_low_dose$FC < 0.9)
#low_downReg_symbols <- topTable_pairing_low_dose[low_downReg,]
#low_downReg_symbols$SYMBOL

range(topTable_pairing_placebo$AveExpr)
range(topTable_pairing_placebo$P.Value)
range(topTable_pairing_placebo$adj.P.Val)
summary(topTable_pairing_placebo$adj.P.Val < 0.05)
range(topTable_pairing_placebo$B)
range(topTable_pairing_placebo$FC)
count(topTable_pairing_placebo$adj.P.Val < 0.05 & topTable_pairing_placebo$FC > 1.1)
count(topTable_pairing_placebo$adj.P.Val < 0.05 & topTable_pairing_placebo$FC < 0.9)
placebo_downReg <- which(topTable_pairing_placebo$adj.P.Val < 0.05 & topTable_pairing_placebo$FC < 0.9)
#placebo_downReg_symbols <- topTable_pairing_placebo[low_downReg,]
#placebo_downReg_symbols$SYMBOL


# Plot the linear model fit for diagnostics, this plots gene-wise residual standard deviations against average log-expression
# to identify mean-variance trends (sigma vs A plot). ?plotSA:
png('sigma_vs_A_plot_lm_pairing_treatment.png', width = 12, height = 12, units = 'in', res = 300)
plotSA(fit_all_pairs_2)
dev.off()

#Plot:
png('volcano_plots_pairing_treatment.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
plot(topTable_pairing_placebo$FC, -log10(topTable_pairing_placebo$adj.P.Val), 
     pch=20, main="Placebo paired tests", abline(h=1.30103, v=c(0.9, 1.1), col='red'))
plot(topTable_pairing_low_dose$FC, -log10(topTable_pairing_low_dose$adj.P.Val), 
     pch=20, main="Low dose paired tests", abline(h=1.30103, v=c(0.9, 1.1), col='red'))
plot(topTable_pairing_high_dose$FC, -log10(topTable_pairing_high_dose$adj.P.Val), 
     pch=20, main="High dose paired tests", abline(h=1.30103, v=c(0.9, 1.1), col='red'))
par(mfrow=c(1,1))
dev.off()


# TO DO: Check genas function again, results are positive and seem overly sensitive. Name plot (it's saved as Rplots.pdf)
# Estimate correlation from gene expression profiles:
ncol(fit_all_pairs_2$coefficients)
fit_all_pairs_2$coefficients[1:5, 299:301]
genas(fit_all_pairs_2, coef = c(299, 300), subset = 'logFC', plot = TRUE)
genas(fit_all_pairs_2, coef = c(299, 301), subset = 'logFC', plot = TRUE)
genas(fit_all_pairs_2, coef = c(300, 301), subset = 'logFC', plot = TRUE)


# Clasify results as up-, down- or not significant. See ?decideTests
#propTrueNull(fit_all_pairs_2)
fit_all_pairs_2$coef[1:3, 299:ncol(fit_all_pairs_2$coef)]
fit_all_pairs_2_subset_coef <- fit_all_pairs_2[, 299:301]#fit_all_pairs_2$coef[, 299:301]
results_by_pairing <- decideTests(fit_all_pairs_2_subset_coef, method = 'hierarchical')
head(results_by_pairing)
str(results_by_pairing)

heatDiagram(results_by_pairing, fit_all_pairs_2_subset_coef)#, primary = 'treatmenttreated_4000')
vennDiagram(results_by_pairing)


# Write results to disk:
# write.table(x=head(topTable_pairing_high_dose, n=20), sep='\t', quote = FALSE,
#             col.names = NA, row.names = TRUE,
#             file='top20_topTable_pairing_high_dose.txt')
# 
# write.table(x=head(topTable_pairing_low_dose, n=20), sep='\t', quote = FALSE,
#             col.names = NA, row.names = TRUE,
#             file='top20_topTable_pairing_low_dose.txt')
# 
# write.table(x=head(topTable_pairing_placebo, n=20), sep='\t', quote = FALSE,
#             col.names = NA, row.names = TRUE,
#             file='top20_topTable_pairing_placebo.txt')

write.table(x=topTable_pairing_high_dose, sep='\t', quote = FALSE,
            col.names = NA, row.names = TRUE,
            file='full_topTable_pairing_high_dose.txt')

write.table(x=topTable_pairing_low_dose, sep='\t', quote = FALSE,
            col.names = NA, row.names = TRUE,
            file='full_topTable_pairing_low_dose.txt')

write.table(x=topTable_pairing_placebo, sep='\t', quote = FALSE,
            col.names = NA, row.names = TRUE,
            file='full_topTable_pairing_placebo.txt')

# Interpretation: There are very few significant differences for paired tests (before vs after) for each condition (placebo, 2000, 4000).

########################
# TO DO: extract the FC values for each probe/sample for eQTL analysis.
# Compare joint 2000+4000 vs baseline and placebo with pairing.
# Define factors to constrast from annotation file:
#Pairing:
head(membership_file_cleaned)
tail(membership_file_cleaned)
# kit_id appears once with corresponding pt_id twice (baseline + final visit):
which(membership_file_cleaned == 104018)
membership_file_cleaned[c(275, 560), ]

pairing_joint <- factor(membership_file_cleaned$pt_id)
head(pairing_joint)
length(pairing_joint)
pairing_joint

treatment_joint <- factor(membership_file_cleaned$all_treated, levels = c('untreated', 'treated_2000+4000', 'treated_placebo'))
head(treatment_joint)
str(treatment_joint)
summary(treatment_joint)
treatment_joint

#Define design and set contrasts:
design_all_treated_pairs <- model.matrix(~pairing_joint+treatment_joint)
head(design_all_treated_pairs)[1:5, 1:5]
tail(design_all_treated_pairs)[, -1]
dim(design_all_treated_pairs)
colnames(design_all_treated_pairs)

#########
# TO DO:
# # Estimate the correlation between measurements made on the same subject. This takes several minutes:
# # Without covariates:
# corfit_all_treated <- duplicateCorrelation(normalised_filtered, design_all_treated_pairs, block = membership_file_cleaned$pt_id)
# corfit_all_treated$consensus
# corfit_all_treated$cor
# corfit_all_treated$consensus.correlation
# head(corfit_all_treated$atanh.correlations)
# #?duplicateCorrelation
# # Generated 50 warnings: In sqrt(dfitted.values) : NaNs produced
# Correlation was 1. ??
# #Then this inter-subject correlation is input into the linear model fit:
# fit__all_treated_corfit <- lmFit(normalised_filtered, design_all_treated_pairs, block = membership_file_cleaned$pt_id, 
#                         correlation = corfit_all_treated$consensus)
#########


#Run linear model and obtain differentially expressed genes based on all pairs:
fit_all_treated <- lmFit(normalised_filtered, design_all_treated_pairs)
fit_all_treated_2 <- eBayes(fit_all_treated)
dim(fit_all_treated_2)
str(fit_all_treated_2)
summary(fit_all_treated_2)
names(fit_all_treated_2)
head(fit_all_treated_2$coefficients)#[1:5, 1:5]
head(fit_all_treated_2$lods)#[1:5, 1:5]
head(fit_all_treated_2$cov.coefficients)#[1:5, 1:5]
fit_all_treated_2$coefficients
colnames(fit_all_treated_2)

# Get results:
topTable(fit_all_treated_2, adjust='BH')
topTable(fit_all_treated_2, coef="treatment_jointtreated_2000+4000", adjust='BH')
topTable(fit_all_treated_2, coef="treatment_jointtreated_placebo", adjust='BH')

topTable_pairing_joint_treated <- topTable(fit_all_treated_2, coef="treatment_jointtreated_2000+4000", adjust='BH', n=Inf)
topTable_pairing_joint_placebo <- topTable(fit_all_treated_2, coef="treatment_jointtreated_placebo", adjust='BH', n=Inf)
head(topTable_pairing_joint_treated)
dim(topTable_pairing_joint_treated)
head(topTable_pairing_joint_placebo)
dim(topTable_pairing_joint_placebo)

# Basic counts
count(topTable_pairing_joint_treated$adj.P.Val < 0.05)
count(topTable_pairing_joint_placebo$adj.P.Val < 0.05)

count(topTable_pairing_joint_treated$adj.P.Val < 0.10)
count(topTable_pairing_joint_placebo$adj.P.Val < 0.10)

#count(topTable_pairing_joint_treated$adj.P.Val < 0.05 & topTable_pairing_joint_treated$logFC > 0.1365) # Around 1.1 FC
count(topTable_pairing_joint_treated$adj.P.Val < 0.05 & 2^(topTable_pairing_joint_treated$logFC) > 1.1)
#count(topTable_pairing_joint_treated$adj.P.Val < 0.05 & topTable_pairing_joint_treated$logFC < -0.15) # Around 0.9 FC
count(topTable_pairing_joint_treated$adj.P.Val < 0.05 & 2^(topTable_pairing_joint_treated$logFC) < 0.9)
count(topTable_pairing_joint_placebo$adj.P.Val < 0.05 & topTable_pairing_joint_placebo$logFC > 0.1365) # Around 1.1 FC
count(topTable_pairing_joint_placebo$adj.P.Val < 0.05 & topTable_pairing_joint_placebo$logFC < -0.15) # Around 0.9 FC

# Plot overlaps:
fit_all_treated_2$coef[1:3, 298:ncol(fit_all_treated_2$coef)]
fit_all_treated_2_subset_coef <- fit_all_treated_2[, 299:300]
results_by_pairing_all_treated <- decideTests(fit_all_treated_2_subset_coef)
head(results_by_pairing_all_treated)
str(results_by_pairing_all_treated)
heatDiagram(results_by_pairing_all_treated, fit_all_treated_2_subset_coef)#, primary = 'treatmenttreated_4000')
vennDiagram(results_by_pairing_all_treated)

# Get correlation:
genas(fit_all_treated_2, coef = c(299, 300), subset = 'logFC', plot = TRUE)

# Get annotations for topTable:
topTable_pairing_joint_placebo <- get_gene_symbols(topTable_pairing_joint_placebo)
topTable_pairing_joint_treated <- get_gene_symbols(topTable_pairing_joint_treated)
head(topTable_pairing_joint_treated)
dim(topTable_pairing_joint_treated)
head(topTable_pairing_joint_placebo)
dim(topTable_pairing_joint_placebo)
# View(topTable_pairing_joint_treated)

# Write results to disk:
write.table(x=topTable_pairing_joint_treated, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE,
            file='full_topTable_pairing_all_treated.txt')

write.table(x=topTable_pairing_joint_placebo, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE,
            file='full_topTable_pairing_all_placebo.txt')

# Interpretation: There are significant differences for paired tests (before vs after) for each condition when joining
# treatment groups (placebo and 2000 + 4000).
# Correlation between changes is similar across both groups though (see genas result).

# TO DO: A ?barcodeplot (limma) can also be used to visualise DE genes.
# To provide statistical significance to correlations use roast()

#########################


#############################
# Write files to disk of subsets, needed for eQTL analysis:
# TO DO: write out FC values for each probe for comparisons above (so as to pass for eQTL analysis. This should be 
# in eBayes object, check what to extract (coefficients or log odds?)).


# Placebo at baseline only:
baseline_placebo <- subset(x=membership_file_cleaned, 
                           subset=membership_file_cleaned$visit_type =='Randomisation'
                           & membership_file_cleaned$arm == 2)
dim(baseline_placebo)
head(baseline_placebo)

array_baseline_placebo <- normalised_filtered[, rownames(baseline_placebo)]
array_baseline_placebo[1:5, 1:5]
dim(array_baseline_placebo)

write.table(array_baseline_placebo, 'GEx_baseline_placebo.tsv', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)

# Placebo at final visit only:
finalVisit_placebo <- subset(x=membership_file_cleaned, 
                             subset=membership_file_cleaned$visit_type =='FinalVisit'
                             & membership_file_cleaned$arm == 2)
dim(finalVisit_placebo)

array_finalVisit_placebo <- normalised_filtered[, rownames(finalVisit_placebo)]
array_finalVisit_placebo[1:5, 1:5]
dim(array_finalVisit_placebo)

write.table(array_finalVisit_placebo, 'GEx_finalVisit_placebo.tsv', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)

# 2000 dose at baseline only:
baseline_2000 <- subset(x=membership_file_cleaned, 
                        subset=membership_file_cleaned$visit_type =='Randomisation'
                        & membership_file_cleaned$arm == 1)
dim(baseline_2000)

array_baseline_2000 <- normalised_filtered[, rownames(baseline_2000)]
array_baseline_2000[1:5, 1:5]
dim(array_baseline_2000)

write.table(array_baseline_2000, 'GEx_baseline_2000.tsv', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)

# 2000 dose at final visit only:
finalVisit_2000 <- subset(x=membership_file_cleaned, 
                          subset=membership_file_cleaned$visit_type =='FinalVisit'
                          & membership_file_cleaned$arm == 1)
dim(finalVisit_2000)

array_finalVisit_2000 <- normalised_filtered[, rownames(finalVisit_2000)]
array_finalVisit_2000[1:5, 1:5]
dim(array_finalVisit_2000)

write.table(array_finalVisit_2000, 'GEx_finalVisit_2000.tsv', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)

# 4000 dose at baseline only:
baseline_4000 <- subset(x=membership_file_cleaned, 
                        subset=membership_file_cleaned$visit_type =='Randomisation'
                        & membership_file_cleaned$arm == 0)
dim(baseline_4000)

array_baseline_4000 <- normalised_filtered[, rownames(baseline_4000)]
array_baseline_4000[1:5, 1:5]
dim(array_baseline_4000)

write.table(array_baseline_4000, 'GEx_baseline_4000.tsv', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)

# 4000 dose at final visit only:
finalVisit_4000 <- subset(x=membership_file_cleaned, 
                          subset=membership_file_cleaned$visit_type =='FinalVisit'
                          & membership_file_cleaned$arm == 0)
dim(finalVisit_4000)

array_finalVisit_4000 <- normalised_filtered[, rownames(finalVisit_4000)]
array_finalVisit_4000[1:5, 1:5]
dim(array_finalVisit_4000)

write.table(array_finalVisit_4000, 'GEx_finalVisit_4000.tsv', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)

array_list <- list(array_baseline_placebo, array_finalVisit_placebo,
                    array_baseline_2000, array_finalVisit_2000, 
                    array_baseline_4000, array_finalVisit_4000)

lapply(X=array_list, dim)
# head(array_list)
class(array_list)
str(array_list)

png('boxplots_subsets_intensities.png', width = 12, height = 12, units = 'in', res = 300)
boxplot(array_baseline_placebo, array_finalVisit_placebo,
        array_baseline_2000, array_finalVisit_2000, 
        array_baseline_4000, array_finalVisit_4000)
dev.off()

# All treated (minus final visit placebo):
count(membership_file_cleaned$all_treated)
count(membership_file_cleaned$two_group_Tx)

treated_4000_and_2000 <- subset(x=membership_file_cleaned, 
                          subset=membership_file_cleaned$all_treated == 'treated_2000+4000')
dim(treated_4000_and_2000)
head(treated_4000_and_2000)

array_treated_4000_and_2000 <- normalised_filtered[, rownames(treated_4000_and_2000)]
array_treated_4000_and_2000[1:5, 1:5]
dim(array_treated_4000_and_2000)

write.table(array_treated_4000_and_2000, 'GEx_treated_4000_and_2000.tsv', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)

# Baseline (minus placebo baseline):
count(membership_file_cleaned$all_treated)
count(membership_file_cleaned$two_group_Tx)

baseline_4000_and_2000 <- subset(x=membership_file_cleaned, 
                                subset=membership_file_cleaned$all_treated == 'untreated')
dim(baseline_4000_and_2000)
head(baseline_4000_and_2000)

array_baseline_4000_and_2000 <- normalised_filtered[, rownames(baseline_4000_and_2000)]
array_baseline_4000_and_2000[1:5, 1:5]
dim(array_baseline_4000_and_2000)

write.table(array_baseline_4000_and_2000, 'GEx_baseline_4000_and_2000.tsv', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)
#############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for higher level analyses.
#############################
