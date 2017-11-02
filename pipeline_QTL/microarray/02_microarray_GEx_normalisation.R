#############################
# To be run after 01 read and QC array data
# Antonio J Berlanga-Taylor

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd("/ifs/projects/proj043/analysis.dir/gene_expression.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_normalisation",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters:

# Read results from 01_microarrayxxx file (saved as RData object):
# Load a previous R session, data and objects:
load('R_session_saved_image_read_and_QC.RData', verbose=T)
#load('R_session_saved_image_normalisation_full.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_normalisation', '.RData', sep='')
R_session_saved_image_full <- paste('R_session_saved_image_normalisation_full', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:

# TO DO: get rid of unused packages, change library for require

library(limma)
library(beadarray)
library(lumi)
library(illuminaHumanv4.db)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)
library(quantro)
library(flashClust)
library(lumiHumanAll.db)
library(convert)
library(vsn)
library(reshape2)
library(grid)
library(gridExtra)
library(plyr)
library(dendextend)
library(gplots)
library(doParallel)

#############################


#############################

# 3) Background correction, normalisation, transformation and filtering

# Reading the control probe profiles is optional but recommended. If the control probe profiles are available, 
# then the Illumina data can be favorably background corrected and normalized using the neqc or nec functions.

## Test whether quantile normalisation is appropriate.
# The package quantro tests whether quantile normalisation is appropriate for the given dataset:
# http://www.bioconductor.org/packages/release/bioc/vignettes/quantro/inst/doc/quantro-vignette.pdf
# browseVignettes("quantro")

# View the distributions of the samples of interest, matdensity() computes the density for each sample (columns):

png('density_plots_raw_cleaned_quantro.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
matdensity(read_files_cleaned_QC$E, groupFactor = membership_file_cleaned$treatment, col = c(1,2,3), xlab = " ", ylab = "density", 
           main = "Beta Values, by treatment")
#legend('top', c("NeuN_neg", "NeuN_pos"), col = c(2,3), lty= 1, lwd = 3)

matdensity(read_files_cleaned_QC$E, groupFactor = membership_file_cleaned$arm, col = c(1,2,3), xlab = " ", ylab = "density", 
           main = "Beta Values, by arm")

matdensity(read_files_cleaned_QC$E, groupFactor = membership_file_cleaned$visit_type, col = c(1,2), xlab = " ", ylab = "density", 
           main = "Beta Values, by visit type")

par(mfrow=c(1,1))
dev.off()

# matboxplot() orders and colors the samples by a group level variable:
#png('boxplots_raw_cleaned_subset_quantro.png', width = 12, height = 12, units = 'in', res = 300)
#par(mfrow=c(3,1))
#matboxplot(read_files_cleaned_QC$E[1:100,1:200], groupFactor = membership_file_cleaned$treatment, col = c(1,2,3), xaxt = "n", 
#           main = "Beta Values, by treatment")
#matboxplot(read_files_cleaned_QC$E[1:100,1:200], groupFactor = membership_file_cleaned$arm, col = c(1,2,3), xaxt = "n",
#           main = "Beta Values, by arm")
#matboxplot(read_files_cleaned_QC$E[1:100,1:200], groupFactor = membership_file_cleaned$visit_type, col = c(1,2), xaxt = "n",
#           main = "Beta Values, by visit type")
#par(mfrow=c(1,1))
#dev.off()

# Run quantro tests and permutations in parallel with doParallel library:
# TO DO extract number of cores as argument; run complete files for permutations
registerDoParallel(cores=4)
quantro_test_treatment <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$treatment, B=1000)
quantro_test_by_arm <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$arm, B=1000)
quantro_test_visit_type <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$visit_type, B=1000)

quantro_list <- c(quantro_test_treatment, quantro_test_by_arm, quantro_test_visit_type)

# Assess statistical significance:
print(quantro_list)
lapply(quantro_list, summary)
lapply(quantro_list, anova)
lapply(quantro_list, quantroStat)


# Plot the test statistics of the permuted samples. The plot is a histogram of the null test statistics quantroStatPerm 
# from quantro() and the red line is the observed test statistic quantroStat from quantro():

png('quantro_histogram_permutations.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
#lapply(quantro_list, quantroPlot)
quantroPlot(quantro_test_treatment)
quantroPlot(quantro_test_by_arm) 
quantroPlot(quantro_test_visit_type)
par(mfrow=c(1,1))
dev.off()

# a) Background correction
# Non-background corrected, non-normalised, sample and control probe profiles and targets file.
# Check backgroundCorrect and plotFB

class(read_files_cleaned_QC)
summary(read_files_cleaned_QC)
summary(read_files_cleaned_QC$genes$Status)

# Explore foreground vs background intensities and decide which method to use: TO DO: Errors due to no x$Eb present (background intensities)
# plotFB(x=read_files_cleaned_QC, array=1, pch='.')

# Background correct. Single-channel arrays methods incude none, nec or normexp. 
# ?nec, use robust=TRUE? offset set at 16 from suggestion in Ritchie et al PLoS Comp Biol 2011 
# http://www.ploscompbiol.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pcbi.1002276&representation=PDF

range(read_files_cleaned_QC$E)

background_corrected_nec <- nec(x=read_files_cleaned_QC)
background_corrected_nec_robust <- nec(x=read_files_cleaned_QC, robust=TRUE)
background_corrected_nec_16 <- nec(x=read_files_cleaned_QC, offset=16)
background_corrected_nec_16_robust <- nec(x=read_files_cleaned_QC, offset=16, robust=TRUE)

range(background_corrected_nec$E)
range(background_corrected_nec_robust$E)
range(background_corrected_nec_16$E)
range(background_corrected_nec_16_robust$E)

background_corrected_auto <- backgroundCorrect(read_files_cleaned_QC, method='auto', verbose=TRUE)
range(background_corrected_auto$E)

background_corrected_normexp <- backgroundCorrect(read_files_cleaned_QC, method='normexp', verbose=TRUE)
range(background_corrected_normexp$E)


# background_corrected_subtract <- backgroundCorrect(read_files_cleaned_QC, method='subtract', verbose=TRUE)
# range(background_corrected_subtract$E)

# background_corrected_movingmin <- backgroundCorrect(read_files_cleaned_QC, method='movingmin', verbose=TRUE)
# range(background_corrected_movingmin$E)

# background_corrected_edwards <- backgroundCorrect(read_files_cleaned_QC, method='edwards', verbose=TRUE)
# range(background_corrected_edwards$E)


# Visualise and explore the fit of the dataset:
meanSdPlot_background_corrected_nec <- ('mean_SD_plot_background_corrected_nec.png')
png(meanSdPlot_background_corrected_nec, width = 6, height = 12, units = 'in', res = 300)
par(mfrow=c(4,1))
meanSdPlot(background_corrected_nec$E, main='Background corrected expression values - nec')
meanSdPlot(background_corrected_nec_robust$E, main='Background corrected expression values - nec, robust')
meanSdPlot(background_corrected_nec_16$E, main='Background corrected expression values - nec, offset 16')
meanSdPlot(background_corrected_nec_16_robust$E, main='Background corrected expression values - nec, robust, offset 16')
par(mfrow=c(1,1))
dev.off()


meanSdPlot_background_corrected_normexp <- ('mean_SD_plot_background_corrected_normexp_and_auto.png')
png(meanSdPlot_background_corrected_normexp, width = 6, height = 10, units = 'in', res = 300)
par(mfrow=c(2,1))
meanSdPlot(background_corrected_normexp$E, main='Background corrected expression values - normexp')
# meanSdPlot_background_corrected_normexp <- ('mean_SD_plot_background_corrected_normexp.png')
meanSdPlot(background_corrected_auto$E, main='Background corrected expression values - auto')
par(mfrow=c(1,1))
dev.off()


# Mean-difference plots, using limma's plotMA:
# Also check mdplot

plotMA_background_corrected_1 <- ('mean_difference_plots_background_corrected_1.png')
png(plotMA_background_corrected_1, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
limma::plotMA(read_files_cleaned_QC$E, main='Raw expression values')
limma::plotMA(background_corrected_nec$E, main='Background corrected expression values - nec')
limma::plotMA(background_corrected_nec_robust$E, main='Background corrected expression values - nec, robust')
limma::plotMA(background_corrected_nec_16$E, main='Background corrected expression values - nec, offset 16')
par(mfrow=c(1,1))
dev.off()


plotMA_background_corrected_2 <- ('mean_difference_plots_background_corrected_2.png')
png(plotMA_background_corrected_2, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
limma::plotMA(background_corrected_nec_16_robust$E, main='Background corrected expression values - nec, robust, offset 16')
limma::plotMA(background_corrected_normexp$E, main='Background corrected expression values - normexp')
limma::plotMA(background_corrected_auto$E, main='Background corrected expression values - auto')
par(mfrow=c(1,1))
dev.off()


# Plot density distributions:
# ?plotDensity

plotDensity_background_corrected <- ('density_distribution_background_corrected.png')
png(plotDensity_background_corrected, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,2))
plotDensity(read_files_cleaned_QC$E, logMode=F, addLegend=F, main='Density plot of raw expression values')
plotDensity(background_corrected_nec$E, logMode=F, addLegend=F, 
            main='Density distribution of background corrected expression values - nec')
plotDensity(background_corrected_nec_16$E, logMode=F, addLegend=F, 
            main='Density distribution of background corrected expression values - nec, 16')
plotDensity(background_corrected_nec_robust$E, logMode=F, addLegend=F, 
            main='Density distribution of background corrected expression values - nec, robust')
plotDensity(background_corrected_nec_16_robust$E, logMode=F, addLegend=F, 
            main='Density distribution of background corrected expression values - nec, robust, 16')
plotDensity(background_corrected_normexp$E, logMode=F, addLegend=F, main='Background corrected normexp')

par(mfrow=c(1,1))
dev.off()


# TO DO: extract as parameter:
background_corrected = background_corrected_normexp

class(background_corrected)
#head(background_corrected)
str(background_corrected)
summary(background_corrected)
range(background_corrected$E)


#############################
# TO DO: is this needed: clean up:
## Call arrayQualityMetrics
# Transform EListRaw class to ExpressionSet class (as more general and accepted by more downstream methods):

# This extracts expression values (with probes as rownames) and sample names (column headers) and returns as a matrix:
background_corrected_as_matrix <- as.matrix(background_corrected, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
class(background_corrected_as_matrix)
dim(background_corrected_as_matrix)
head(background_corrected_as_matrix)[1:5,1:5]
head(colnames(background_corrected_as_matrix))
head(rownames(background_corrected_as_matrix))


# Remove duplicated probes from expression matrix (these seem to be the negative control probes):
duplicated_probes <- which(duplicated(rownames(background_corrected_as_matrix)))
duplicated_probes

# If probes have duplicated names but for some reason need to be kept: rownames(as_matrix) = make.names(rownames(as_matrix), unique=TRUE)

# Output clean expression matrix:
background_corrected_as_matrix_dedup <- background_corrected_as_matrix[-duplicated_probes,]
head(background_corrected_as_matrix_dedup)
dim(background_corrected_as_matrix_dedup)

# Convert expression matrix to expression set:
background_corrected_minimal_eset <- ExpressionSet(assayData=background_corrected_as_matrix_dedup)
class(background_corrected_minimal_eset)
dim(background_corrected_minimal_eset)

# Call arrayQualityMetrics after processing data and background correction:
#arrayQualityMetrics_background_corrected <- arrayQualityMetrics(expressionset=background_corrected_minimal_eset, 
#                                                                outdir=paste('arrayQualityMetrics_background_corrected.dir'), 
#                                                                force=TRUE, do.logtransform=TRUE)


#############################

# b) Transformation
# log2 transformation is run by neqc, see other approaches.

# c) Normalisation
# Quantile approach is run already by limma s neqc function.
# Also check VST and other approaches.


## The neqc function performs normexp background correction using negative controls, then quantile normalizes 
# and finally log2 transforms [see ]. It also automatically removes the control probes, leaving only the regular 
# probes:
#TO DO: When applying quantile normalization, it is assumed that the distribution in signal should be 
# the same from each array. 

run_neqc <- neqc(read_files_cleaned_QC)
dim(run_neqc)
summary(run_neqc)
class(run_neqc)
range(run_neqc$E)

## Other functions to normalise and  transform: 
#?normalizeBetweenArrays # This is run for single channel arrays. For two-colour arrays normaliseWithinArrays must be run first.
# For single-channel data check scale, quantile and cyclicloess 

# TO DO: extract method as argument

between_arrays <- normalizeBetweenArrays(object=read_files_cleaned_QC, method='scale')
class(between_arrays)
#head(between_arrays)
dim(between_arrays)
range(between_arrays$E)

between_arrays_loess <- normalizeBetweenArrays(object=read_files_cleaned_QC, method='cyclicloess')
class(between_arrays_loess)
#head(between_arrays_loess)
dim(between_arrays_loess)
range(between_arrays_loess$E)

between_arrays_loess_affy <- normalizeBetweenArrays(object=read_files_cleaned_QC, method='cyclicloess', cyclic.method='affy')
class(between_arrays_loess_affy)
dim(between_arrays_loess_affy)
range(between_arrays_loess_affy$E)

# Variance-stabilising normalisation from the 'vsn' package. This background corrects then normalises. 
# Limma's normalizeVSN is an interface to vsn's vsnMatrix function (which run's vsn's vsn2 function, 
# see ?normalizeVSN and ?vsn2 or ?vsnMatrix, vignette('Introduction to vsn'), ?normalizeVSN.

normalize_VSN <- normalizeVSN(read_files_cleaned_QC)
class(normalize_VSN)
dim(normalize_VSN)
range(normalize_VSN$E)
head(normalize_VSN$E)

# Other methods for normalisation (require object as ExpressionSet so transform data object):
# ?normaliseIllumina from beadarray

#transform_vst <- normaliseIllumina(read_files_cleaned_QC, transform='vst')
#run_median <- normaliseIllumina(minimal_eset, method='median')
#run_qspline <- normaliseIllumina(minimal_eset, method='qspline')


#Visualise normalisation methods:
# run_neqc
# between_arrays
# between_arrays_loess
# normalize_VSN

meanSdPlot_normalised <- ('meanSdPlot_normalised.png')
png(meanSdPlot_normalised, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
meanSdPlot(run_neqc$E, main='Normalised expression values - neqc')
meanSdPlot(between_arrays$E, main='Normalised expression values - scale')
meanSdPlot(between_arrays_loess$E, main='Normalised expression values - cyclic loess')
meanSdPlot(normalize_VSN$E, main='Non-background corrected, normalised expression values - VSN')
par(mfrow=c(1,1))
dev.off()


# Mean-difference plots, using limma's plotMA:
# Also check mdplot

plotMA_normalised <- ('mean_difference_plots_normalised.png')
png(plotMA_normalised, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
limma::plotMA(run_neqc$E, main='Normalised expression values - neqc')
limma::plotMA(between_arrays$E, main='Normalised expression values - scale')
limma::plotMA(between_arrays_loess$E, main='Normalised expression values - cyclic loess')
limma::plotMA(normalize_VSN$E, main='Non-background corrected, normalised expression values - VSN')
par(mfrow=c(1,1))
dev.off()


# Plot density distributions:
# ?plotDensity

plotDensity_normalised <- ('density_distribution_normalised.png')
png(plotDensity_normalised, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
plotDensity(run_neqc$E, logMode=F, addLegend=F, main='Density plot of normalised expression values - neqc')
plotDensity(between_arrays$E, logMode=F, addLegend=F, 
            main='Density distribution of normalised expression values - scale')
plotDensity(between_arrays_loess$E, logMode=F, addLegend=F, 
            main='Density distribution of normalised expression values - cyclic loess')
plotDensity(normalize_VSN$E, logMode=F, addLegend=F, 
            main='Density distribution of non-background corrected, normalised expression values - VSN')
par(mfrow=c(1,1))
dev.off()


# TO DO: Turn into parameter set at the beginning of the script .
# Change object names for downstream analysis:
normalised = normalize_VSN

# Explore contents of the ExpressionSet created (as this is the output from normaliseIllumina):
normalised
dim(normalised)
class(normalised)
summary(normalised)
range(normalised$E)


# If looking at an ExpressionSet:
#head(featureNames(object=normalised))
#head(sampleNames(object=normalised))

# Obtain the actual expression values from the ExpressionSet (this returns a matrix):
#expression_values_normalised <- exprs(normalised)

#class(expression_values_normalised)
#head(expression_values_normalised)
#range(expression_values_normalised[,1])


## TO DO: Plot and explore normalised data
# Extract expression values and convert to matrix (needed for some plots):
# normalised_as_matrix <- as.matrix(normalised$E, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
# class(normalised_as_matrix)
# dim(normalised_as_matrix)
# head(normalised_as_matrix)[1:5,1:5]
# head(colnames(normalised_as_matrix))
# head(rownames(normalised_as_matrix))


## Visualise and explore the fit of the dataset:

# Cumulative density plot

# Pairwise sample correlations

# Heatmap



## Summary plots of raw values and chosen background correction and normalisation methods:
# Density distributions:
summary_plots_density <- ('summary_plots_density_distributions.png')
png(summary_plots_density, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
plotDensity(read_files_cleaned_QC$E, logMode=F, addLegend=F, main='Density distribution of raw expression values')
plotDensity(background_corrected$E, logMode=F, addLegend=F, main='Density distribution of background corrected expression values')
plotDensity(normalised$E, logMode=F, addLegend=F, main='Density distribution of normalised expression values')
par(mfrow=c(1,1))
dev.off()

# Plot row standard deviations versus row means:
summary_meanSD_plots <- ('summary_plots_mean_SD.png')
png(summary_meanSD_plots, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
meanSdPlot(read_files_cleaned_QC$E, logMode=F, addLegend=F, main='Density distribution of raw expression values')
meanSdPlot(background_corrected$E, logMode=F, addLegend=F, main='Density distribution of background corrected expression values')
meanSdPlot(normalised$E, logMode=T, addLegend=F, main='Density distribution of normalised expression values')
par(mfrow=c(1,1))
dev.off()


# Plot cumulative distributions, ?plotCDF:
plotCDF_summary <- ('summary_cumulative_distribution_plots.png')
png(plotCDF_summary, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
plotCDF(read_files_cleaned_QC$E, reverse=FALSE, logMode=TRUE, addLegend=FALSE, main='Cumulative distribution of raw expression values')
plotCDF(background_corrected$E, reverse=FALSE, logMode=TRUE, addLegend=FALSE, main='Cumulative distribution of background corrected expression values')
plotCDF(normalised$E, reverse=FALSE, logMode=TRUE, addLegend=FALSE, main='Cumulative distribution of normalised expression values')
par(mfrow=c(1,1))
dev.off()

#############################
# TO DO : skip this as not really needed, acts like a sanity check.
## Call arrayQualityMetrics
# Transform EListRaw class to ExpressionSet class (as more general and accepted by more downstream methods):

# This extracts expression values (with probes as rownames) and sample names (column headers) and returns as a matrix:
normalised_as_matrix <- as.matrix(normalised, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
class(normalised_as_matrix)
dim(normalised_as_matrix)
head(normalised_as_matrix)[1:5,1:5]
head(colnames(normalised_as_matrix))
head(rownames(normalised_as_matrix))


# Remove duplicated probes from expression matrix (these seem to be the negative control probes):
duplicated_probes <- which(duplicated(rownames(normalised_as_matrix)))
duplicated_probes
length(duplicated_probes)

# If probes have duplicated names but for some reason need to be kept: rownames(as_matrix) = make.names(rownames(as_matrix), unique=TRUE)

# Output clean expression matrix:
normalised_as_matrix_dedup <- normalised_as_matrix[-duplicated_probes,]
head(normalised_as_matrix_dedup)
dim(normalised_as_matrix_dedup)

# Convert expression matrix to expression set:
normalised_minimal_eset <- ExpressionSet(assayData=normalised_as_matrix_dedup)
class(normalised_minimal_eset)
dim(normalised_minimal_eset)

# TO DO: Check this section as it's not carried forward and arrayQualityMetrics is not called either...
# Call arrayQualityMetrics after processing data and normalisation:
#arrayQualityMetrics_normalised <- arrayQualityMetrics(expressionset=normalised_minimal_eset, 
#                                                      outdir=paste('arrayQualityMetrics_normalised.dir'), 
#                                                      force=TRUE)


#############################
# Save files required for other packages/programmes:
head(membership_file_cleaned)
write.table(membership_file_cleaned, 'membership_file_cleaned_all.tab', sep='\t', quote=FALSE, na='NA', 
            col.names=NA, row.names=TRUE)
write.csv(membership_file_cleaned, 'membership_file_cleaned_all.csv', quote=FALSE, na='NA')
write.csv(FAILED_QC_unique, 'FAILED_QC_unique.csv', quote=FALSE, col.names = FALSE, row.names = FALSE)


# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image_full, compress='gzip')

objects_to_save <- (c('normalised', 'membership_file_cleaned', 'FAILED_QC_unique'))
save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: run the script for differential gene expression.
#############################
