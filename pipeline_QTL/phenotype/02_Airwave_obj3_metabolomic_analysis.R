###############
# Airwave phenotype data
# 07 March 2016
# Airwave chronic inflammation proposal
# Objective 3 
# Metabolomic basis of chronic inflammation and NCD risk factors - differential expression analysis
###############

###############
# To DO's:
#   Checklist for project management
#   Checklist for data security
#   Data analysis protocol per objective
#   Checklist exploratory analysis
#   Get metabolome data
#   Get genetics data
#   Query outcome data with Paul
###############


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/Users/antoniob/Desktop/Airwave_inflammation/results_1.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_Airwave_obj3_metabolomics",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is very project specific. Check ways of making count comparisons.

# Re-load a previous R session, data and objects:
# load('R_session_saved_image_Airwave_metabolomics.RData', verbose=T) 

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_Airwave_metabolomics.RData')

#############################


#############################
## Update packages if necessary and load them:
# source("https://bioconductor.org/biocLite.R")
# biocLite('pcaMethods')
# install.packages('pcaMethods')

library(data.table)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(splines)
library(plyr)
library(metabolomics)

# Get script with functions needed:
# source('functions_for_metabolomics.R')
source('moveme.R')
#############################

#############################################
# Run with command line arguments:
options(echo=TRUE) # to see commands in output file. TO DO: check how it works with sink() above.
args <- commandArgs(trailingOnly = TRUE)

#############################################


#############################################
#  Set variables:
# TO DO: pass to configuration file
phenotype_file <- as.character(args[1])
# phenotype_file <- 'screen_data_processed.tsv'
metabolomics_file <- as.chracter(args[2])
# metabolomics_file <- 'Airwave_CPMG_Plasma.txt'
principal_components_file <- as.character(args[3]) 
# principal_components_file <- 'principal_components_normalised_filtered_PC20.tsv'
print(args)
#############################################


#############################
## Load data sets:
## Medical history, biochemical, behaviours, etc:
screen_data <- read.csv(phenotype_file, sep = '\t', header = TRUE, stringsAsFactors=FALSE,
                        na.string = c('NULL', '.', '_1_', '_2_', '_3_', '_4_', '_5_', '_6_',
                                      '', ' ', 'NA', 'NaN', '<0', '<NA>',
                                      '-1', '-2', '-3', '-4' , '-5', '-6'))

# screen_data <- fread(phenotype_file, sep = '\t', header = TRUE, stringsAsFactors=FALSE, 
#                      colClasses = c('DIAG_DIABETES_TYPE_1', 'DIAG_DIABETES_TYPE_2'),
#                         na.string = c('NULL', '.', '_1_', '_2_', '_3_', '_4_', '_5_', '_6_',
#                                       '', ' ', 'NA', 'NaN', '<0', '<NA>',
#                                       '-1', '-2', '-3', '-4' , '-5', '-6'), verbose = TRUE)
# setkey(screen_data, BARCODE)
# ?data.table

# Convert negative values to NAs as these aren't read properly with read.csv:
screen_data <- as.data.frame(lapply(screen_data, function(x){replace(x, x <0, NA)}))
names(screen_data)
dim(screen_data)
class(screen_data)
screen_data[1725:1730, 25:30] # data.table gives warnings for variables types, at col 25 and 28 row 1727
class(screen_data)
dim(screen_data)
names(screen_data)
screen_data[1:5, 1:5]
#View(screen_data)

# Set individual IDs as row names:
row.names(screen_data) <- screen_data$BARCODE

# Read in metabolomics data:
metabolomics_data <- fread(metabolomics_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE) #check if first cell is not empty
class(metabolomics_data)
# setkey(metabolomics_data, 'Row') # data.table doesn't take column numbers
tables()
head(names(metabolomics_data))
head(metabolomics_data[, Row])
metabolomics_data[1:5, 1:5, with = F]
dim(metabolomics_data)

# Set individual IDs as row names:
row.names(metabolomics_data) <- metabolomics_data$Row
# metabolomics_data <- metabolomics_data[, 'Row' := NULL] # deletes this column
metabolomics_data[1:5, 1:5, with = F]
head(row.names(metabolomics_data))
row.names(metabolomics_data)
metabolomics_samples <- row.names(metabolomics_data)
class(metabolomics_samples)
head(metabolomics_samples)
write.table(x = metabolomics_samples, file = sprintf('%s_samples_names.txt', metabolomics_file), sep = '', 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# principal_components_data <- read.csv(principal_components_file, sep = '\t', header = TRUE)
##################


##################
# Check what the IDs correspond to:
# Rows are metabolites and columns are sample IDs as eg 'Plasma_CPMG_NMR_8_998017', 
# last part corresponds to BARCODE:
metabolomics_data[1:5, 1:5, with = F]
class(metabolomics_data)
head(colnames(metabolomics_data))
head(rownames(metabolomics_data))
pheno_to_keep <- which(screen_data$BARCODE %in% metabolomics_data[, Row])
length(pheno_to_keep)
dim(metabolomics_data)

# Subset phenotype data to those with metabolomics data;
pheno_with_metabolomics <- screen_data[pheno_to_keep, ]
head(pheno_with_metabolomics)
dim(pheno_with_metabolomics)
count(pheno_with_metabolomics$BMI_binned)

# Other molecules of interest:
count(pheno_with_metabolomics$BILIRUBIN)
summary(pheno_with_metabolomics$BILIRUBIN)

length(which(!is.na(pheno_with_metabolomics$C_REACTIVE_PROTEIN)))
summary(pheno_with_metabolomics$C_REACTIVE_PROTEIN)

length(which(!is.na(pheno_with_metabolomics$NEUTROPHILS_COUNT)))
summary(pheno_with_metabolomics$NEUTROPHILS_COUNT)

length(which(!is.na(pheno_with_metabolomics$FIBRINOGEN)))
summary(pheno_with_metabolomics$FIBRINOGEN)

length(which(!is.na(pheno_with_metabolomics$GLUCOSE)))
summary(pheno_with_metabolomics$GLUCOSE)

length(which(!is.na(pheno_with_metabolomics$HBA1C_PERCENT)))
summary(pheno_with_metabolomics$HBA1C_PERCENT)

length(which(!is.na(pheno_with_metabolomics$TOTAL_CHOLESTEROL)))
summary(pheno_with_metabolomics$TOTAL_CHOLESTEROL)


# Order files so that rows and columns match, both files have rows as sampe IDs and columns as variables:
pheno_with_metabolomics <- pheno_with_metabolomics[order(pheno_with_metabolomics$BARCODE), ]
metabolomics_data <- metabolomics_data[order(metabolomics_data[, Row]), ]
metabolomics_data[1:5, 1:5, with = F]
# row.names(pheno_with_metabolomics) <- pheno_with_metabolomics$BARCODE

# Sanity check:
# TO DO: Raise error and stop if false:
identical(pheno_with_metabolomics$BARCODE, metabolomics_data[, Row])
length(which(pheno_with_metabolomics$BARCODE %in% metabolomics_data[, Row]))

# Get transposed file as needed downstream:
metabolomics_data_t <- transpose(metabolomics_data)
metabolomics_data_t[1:5, 1:5, with = F]
colnames(metabolomics_data_t) <- as.character(metabolomics_data_t[1, ,])
colnames(metabolomics_data_t)
metabolomics_data_t <- metabolomics_data_t[-1, ,]
metabolomics_data_t[1:5, 1:5, with = F]

metabolomics_data[1:5, 1:5, with = F]
metabolomics_data[nrow(metabolomics_data), 1:5, with = F]
##################

##################
# TO DO: check this as row names get moved around...
# Get object that matches input for metabolomics package:
# metabolomics_data_metab_pack <- as.data.frame(metabolomics_data[1:1000, 1:1000, with = F])
# metabolomics_data_metab_pack[1:5, 1:5]
# class(metabolomics_data_metab_pack)
# head(names(metabolomics_data_metab_pack))
# head(row.names(metabolomics_data_metab_pack))
# row.names(metabolomics_data_metab_pack) <- row.names(pheno_with_metabolomics[1:1000, ])
# bmi_only <- pheno_with_metabolomics[1:1000, c('BMI_binned', 'BARCODE')]
# metabolomics_data_metab_pack <- merge(metabolomics_data_metab_pack, bmi_only)
# metabolomics_data_metab_pack <- metabolomics_data_metab_pack[, moveme(names(metabolomics_data_metab_pack), 
#                                                                       'BMI_binned first')]
# dim(metabolomics_data_metab_pack)
# class(metabolomics_data_metab_pack)
# head(names(metabolomics_data_metab_pack))
# tail(names(metabolomics_data_metab_pack))
# head(row.names(metabolomics_data_metab_pack))
# tail(row.names(metabolomics_data_metab_pack))


##################

##################
# TO DO: continue from here, long computations...

metabolomics_data[1:5, 2:5, with = F]
ncol(metabolomics_data)
nrow(metabolomics_data)

# Get PCs for metabolomics data:
#Plot probes in a multi-dimensional scaling plot:
plot_MDS <- ("MDS_metabolomics.png")
png(plot_MDS, width = 4, height = 4, units = 'in', res = 300)
# plotMDS(metabolomics_data[1:5, 2:100, with = F], pch=1)
plotMDS(metabolomics_data[1:5, 2:ncol(metabolomics_data), with = F], pch=1)
dev.off()

plot_MDS_by_targets <- ("MDS_metabolomics_by_targets.png")
png(plot_MDS_by_targets, width = 4, height = 4, units = 'in', res = 300)
# plotMDS(metabolomics_data[2:100, 2:100, with = F], pch=1, labels = pheno_with_metabolomics$BMI_binned)
plotMDS(metabolomics_data[1:5, 2:ncol(metabolomics_data), with = F], pch=1, labels = pheno_with_metabolomics$BMI_binned)
dev.off()

# Plot PCA of normalised samples:
# Compute the PCs, use transposed expression values:
metabolomics_data_t[, 1:5, with = F]
# pca_test <- prcomp(metabolomics_data_t[, 1:100, with = F], center=TRUE, scale=TRUE)
# pca_test$x
pca_values <- prcomp(metabolomics_data_t, center=TRUE, scale=TRUE)
head(pca_values$x)
tail(pca_values$x)

# TO DO: match PC's and sample names
# Obtain values for all PCs output:
pc <- data.frame(round(pca_values$x, 2))
dim(pc)
head(pc)
pc[1:5, 1:5]
pc[1:5, ncol(pc)-5:ncol(pc)]
pc$sample_id <- row.names(pc)
head(pc$sample_id)

pc <- pc[, moveme(names(pc), 'sample_id first')]
names(pc)[1:10]
class(pc)
write.table(pc, 'principal_components_metabolomics.tsv', quote = FALSE, sep = '\t', row.names = FALSE)

# Explore dimensions and plot first 10 or so components:
dim(pc)
dim(metabolomics_data)
str(pca_values)
head(pc)
# pc[1:5, 1:5]
summary(pca_values)

# Plot PCA results:
plot_PCA <- ('plot_PCA.png')
png(plot_PCA, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
# Histogram of first x PCs:
plot(pca_values, main='Normalised metabolomics values')
# Scatterplot of PC1 and PC2:
biplot(pca_values, main='Normalised metabolomics values')
par(mfrow=c(1,1))
dev.off()

# Check how much of the variance in gene expression values is explained by the first x PCs:
# sum(pca_normalised_filtered$sdev[1:10]^2)/length(normalised_filtered[1,])

# Run PCA analysis by groups of interest: TO DO: Cross files and IDs first
pc_data_top_20 <- data.frame(pc[, 1:20])
str(pc_data_top_20)
head(pc_data_top_20)
row.names(pc_data_top_20)

pca_by_groups <- data.frame(merge(pheno_with_metabolomics, pc_data_top_20, by='row.names'))
head(pca_by_groups)
dim(pca_by_groups)
dim(pc_data_top_20)

head(arrange(pc_data_top_20, PC1), 10)
head(arrange(pca_by_groups, PC1), 10)

plot_PCA_by_groups_1 <- ('plot_PCA_by_groups_1.png')
png(plot_PCA_by_groups_1, width = 13, height = 13, units = 'in', res = 300)
p1 <- ggplot(data=pc_data_top_20, aes(x=PC1, y=PC2, colour=factor(pheno_with_metabolomics$BMI_binned))) + theme(legend.position="bottom")
p1
p2 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(pheno_with_metabolomics$TAKES_INTENSE_EXERCISE)) + theme(legend.position="bottom")
p3 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(pheno_with_metabolomics$IS_SMOKER)) + theme(legend.position="bottom")
# p4 <- qplot(x=PC2, y=PC3, data=pca_by_groups, colour=factor(pheno_with_metabolomics$)) + theme(legend.position="bottom")
# grid.arrange(p1, p2, p3, p4, ncol=2)
grid.arrange(p1, p2, p3, ncol=2)
dev.off()


##################

##################
# Exploratory analyses of metabolomics data:

# Some stats per sample:
class(metabolomics_data) # data.table object so indexing is DT[i, j, by] with cols requiring names not numbers
range(metabolomics_data[3, ])
range(metabolomics_data[100, ])
range(metabolomics_data[1000, ])

# Some stats per metabolite:
summary(metabolomics_data[, list(Plasma_CPMG_NMR_0_500066)])
summary(metabolomics_data[, 100, with = FALSE])
summary(metabolomics_data[, 1000, with = FALSE])

# TO DO: errors, need to get correct input, see above, merging pbs.
# Some plots per metabolite:
# png('Plasma_CPMG_NMR_0_500066_boxplot.png', width = 6, height = 6, units = 'in', res = 300)
# MetBoxPlots(inputdata = as.data.frame(metabolomics_data), metname = 'Plasma_CPMG_NMR_0_500066', main = 'Plasma_CPMG_NMR_0_500066')
# dev.off()

png('densities_metabolomics.png', width = 6, height = 6, units = 'in', res = 300)
limma::plotDensities(metabolomics_data_t, legend = F)
dev.off()

png('MA_plot_metabolomics.png', width = 6, height = 6, units = 'in', res = 300)
limma::plotMA(as.matrix(metabolomics_data[1:5, 2:100, with = F]))
dev.off()


#############################

#############################

# 6) Differential gene expression analysis
# Linear modelling using weights
# b) Two group comparison of treated vs untreated: 

# Check NAs:

# Compare samples based on groups:
# Limma will use the first group defined as the reference group:
pheno_with_metabolomics$BMI_binned <- factor(pheno_with_metabolomics$BMI_binned, 
                                             levels=c('BMI_20_25', 'BMI_25_30', 'BMI_over_30', 'BMI_less_20'))

group <- pheno_with_metabolomics$BMI_binned
count(group)

#Define design:
design_by_group <- model.matrix(~group)
head(design_by_group)
tail(design_by_group)
count(design_by_group)

#Run linear model and obtain differentially expressed metabolites (?) based on all pairs:
head(colnames(metabolomics_data_t))
head(row.names(metabolomics_data_t))

fit_by_group <- lmFit(metabolomics_data_t, design_by_group)
fit_by_group
names(fit_by_group)
colnames(fit_by_group)
head(row.names(fit_by_group))

fit_by_group_2 <- eBayes(fit_by_group)
fit_by_group_2
head(fit_by_group_2$coefficients)

topTable(fit_by_group_2, adjust = 'BH')
topTable(fit_by_group_2, adjust = 'BH', coef = 'groupBMI_over_30')
topTable(fit_by_group_2, adjust = 'BH', coef = 'groupBMI_25_30')
topTable(fit_by_group_2, adjust = 'BH', coef = 'groupBMI_less_20')

topTable_groups_bmi <- topTable(fit_by_group_2, adjust = 'BH')
topTable_groups_v_less_20 <- topTable(fit_by_group_2, adjust = 'BH', coef = 'groupBMI_less_20', number = Inf)
topTable_groups_v_25_30 <- topTable(fit_by_group_2, adjust = 'BH', coef = 'groupBMI_25_30', number = Inf)
topTable_groups_v_over_30 <- topTable(fit_by_group_2, adjust = 'BH', coef = 'groupBMI_over_30', number = Inf)

head(topTable_groups_v_over_30)

head(x = topTable_groups_v_less_20)
count(topTable_groups_v_less_20$adj.P.Val < 10e-3)
count(topTable_groups_v_25_30$adj.P.Val < 10e-3)
count(topTable_groups_v_over_30$adj.P.Val < 10e-3)

# count(topTable_groups$adj.P.Val < 0.05 & 2^(topTable_groups$logFC) > 1.49)
# count(topTable_groups$adj.P.Val < 0.05 & 2^(topTable_groups$logFC) < 0.68)
# count(topTable_groups$adj.P.Val < 0.05 & topTable_groups$logFC > abs(0.57))

colnames(fit_by_group_2)
results_by_groups <- decideTests(fit_by_group_2)
vennDiagram(results_by_groups)
vennDiagram(decideTests(fit_by_group_2[, 3]))

png('volcano_plots_bmi_nmr.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
plot(2^(topTable_groups_v_less_20$logFC), -log10(topTable_groups_v_less_20$adj.P.Val), pch=20, main="Volcano plot: Lean vs BMI under 20", abline(h=2, v=c(0.8, 1.2)))
plot(2^(topTable_groups_v_25_30$logFC), -log10(topTable_groups_v_25_30$adj.P.Val), pch=20, main="Volcano plot: Lean vs BMI 25 - 30", abline(h=2, v=c(0.8, 1.2)))
plot(2^(topTable_groups_v_over_30$logFC), -log10(topTable_groups_v_over_30$adj.P.Val), pch=20, main="Volcano plot: Lean vs BMI over 30", abline(h=2, v=c(0.8, 1.2)))
# volcanoplot(fit_by_group_2, coef = 4)
# volcanoplot(fit_by_group_2, coef = 2)
# volcanoplot(fit_by_group_2, coef = 3)
par(mfrow=c(1,1))
dev.off()

# Interpretation: Many significant differences between groups for BMI. Lots to check first.
# TO DO: Get random subsets of individuals and run diff. abundance analysis, no differences expected.
#############################

#############################
# Write files to disk:

write.table(x=topTable_groups, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE, 
            file='full_topTable_groups.txt')

# Write to disk pheno file for those with NMR data:
write.table(x=pheno_with_metabolomics, sep='\t', quote = FALSE, col.names = TRUE, row.names = FALSE, 
            file='pheno_with_metabolomics_Airwave_CPMG_Plasma.txt')

#############################

#############################
# Question
# Background
# Test
# Result
# Plot
# Interpretation
# Next step
#############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for xxx.
#############################