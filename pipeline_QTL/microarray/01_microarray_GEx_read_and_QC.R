#############################
#---
# Title: "01_microarray)GEx_read_and_QC"
# Author: "Antonio Berlanga-Taylor"
# Date: "8 April 2015"
# Updated: "14 April 2015"
# Output: html_document
#---

# Current input: Illumina HumanHT-12 v4 Expression BeadChip summary-level data.
# Outputs: various plots and tables from QC
# Requires specifying 
  # phenotype file: membership_input_file variable with groupings, treatments, etc.
  # list of sample and probe files named "targets_file.txt" with the control and probe file names for limma
  # (will require) array type (affymetrix or illumina)
  # arguments: phenotype file, samples to exclude, ?, check TO DOs

#############################


#############################

# Notes and references:

# Check Limma s users guide, p. 20, Section 9 Single-Channel Experimental Designs on p.40, 
# and case study in Section 17.3 on p.108:
# http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

# Also check:
# - Blog with simple instructions:
# http://gettinggeneticsdone.blogspot.co.uk/2014/12/importing-illumina-beadarray-data-into-r.html?m=1

# - BeadArray Expression Analysis Using Bioconductor paper in PLoS Comp Bio:
# www.ploscompbiol.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pcbi.1002276&representation=PDF

# - BeadArray vignette (Dunning et al 2014, same author group as ploscompbio paper above):
# http://www.bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf

# - Lumi package vignette:
# http://www.bioconductor.org/packages/release/bioc/vignettes/lumi/inst/doc/lumi.pdf

# NOTE! Different functions output different classes (e.g. EList [epxression list], eset [expression set], etc.)
# These may not work downstream! e.g. neqc() (from limma) requires an EListRaw class as input (ie what is output from read.ilmn.targets) and
# produces an EList but normaliseIllumina() (from beadarray) requires an eset as input. To avoid this use the package 'convert':
# http://www.bioconductor.org/packages/release/bioc/manuals/convert/man/convert.pdf
# use the matrices (ie xxx$E data containing probes and samples only) or convert data to other structures (see below).

# - For limma, see p. 13 for data objects created (eg EListRaw, EList, MArrayLM and TestResults) 
# and p.14 for functions used (eg summary, dim, length, ncol, nrow, dimnames, rownames and colnames).
#############################


#############################

# Get text file data from the array facility including the following:
# - Control probe file for each sample with:
# ProbeID, AVG_Signal, BEAD_STDERR, Avg_NBEADS, Detection Pval
# - Sample probe file for each sample with:
# ProbeID, Symbol, AVG_Signal, BEAD_STDERR, Avg_NBEADS, Detection Pval
# - Annotation columns including:
# SEARCH_KEY, ILMN_GENE, CHROMOSOME, DEFINITION, SYNONYMS

# Illumina gene expression microarray analysis steps (Ritchie et al. 2014 PLoS Comp Bio; 
# Ritchie et al. 2015 Nucleic Acids Res.):

# 1) Data acquisition
# 2) Preprocessing and quality assessment

# 3) Background correction, normalisation, transformation and filtering

# 4) Experimental design matrix specification
# 5) Descriptive statistics
# 6) Differential gene expression analysis

# 7) Higher level analysis: Gene ontology, co-expression, gene set analyses
# 8) Integration with other data
#############################


#############################
# Preliminaries

# Script is intended to run from src location (eg "/ifs/devel/antoniob/projects/BEST-D")
# with files containing data and output deposited in separate location 
# (eg /ifs/projects/proj043/analysis.dir/gene_expression.dir).

# Things to do before running:
# - Set working directories.
# - Provide file names (if running one file (see read.ilmn function), or multiple files (create a summary file 
# first, see read.ilmn.targets function)).
# - Requires Illumina processed files returned as summary level data, no QC for bead level data at the moment. 
# - After QC, arrays will likely need to be removed, this needs intervention 
# (after assessing numbers and plots) and then passing array (Sample IDs) as a list for limma to remove. 

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd("/ifs/projects/proj043/analysis.dir/gene_expression.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_read_and_QC",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

# Load previous session if needed:
#load(file='R_session_saved_image_normalisation_full.RData', verbose=TRUE)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_read_and_QC', '.RData', sep='')
R_session_saved_image_full <- paste('R_session_saved_image_read_and_QC_full', '.RData', sep='')


#############################


#############################
## Update packages if necessary and load them:

# vignette("lumi")

#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
#install.packages('dendextend')

library(limma)
library(beadarray)
library(lumi)
library(illuminaHumanv4.db)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)
library(quantro)
library(lumiHumanAll.db)
library(plyr)

#############################

# 1) Data acquisition:
#   Raw data comprises one observation per pixel, per array.
# Bead-level data comprises one observation per bead, per array.
# Summary-level data comprises one observation per probe type, per sample
# (Ritchie et al. PLos Comp Bio, 2014)

# a) Bead-level data: intensity and location information for each bead on each BeadArray is 
# available (eg txt files with this information plus other optional/recommended files such as .bab 
# [compressed data from the txt files], raw TIFF image data, sdf files, targets file, metrics file, etc., 
# for analysis with beadarray package). See pdf vignette (Dunning et al 2014) above, analysis not included here.

# b) Summary level data (typical): files after processing in the BeadStudio software, for analysis with limma, 
# lumi and beadarray packages. Files can be at different levels of processing (intensities with/without 
# background correction and/or normalised). Files include sample probe profile (text file as data-matrix with 
# 48,000 rowswith non-normalised summary values as output by BeadStudio, required); control probe profile (text 
# file with summarised data for each of the controls of each array used for diagnostics and calibration, recommended);
# targets file (user created text file with information on which sample was hybridised to each array, 
# required if reading in multiple probe profiles). Files with normalised intensities typically end in avg and 
# files with intensites per gene may also be available. Avoid these as well as Illumina background correction and 
# normalisation.

# Get probe summary profiles (containing the intensity data). Best to obtain intensities from 
# GenomeStudio or BeadStudio without background correction or normalization. Probe summary files are 
# tab-delimited and usually arrays processed at one time are written to a single file.
# Also obtain profiles for the control probes from BeadStudio or GenomeStudio processed data.

# After processing with BeadStudio, for each array probe profile files contain: 
#   summarized expression level (AVG_Signal)
# standard error of the bead replicates (BEAD_STDERR)
# number of beads (Avg_NBEADS) 
# detection p-value (DetectionPval), estimation of the probability of a gene being detected 
# above the background level


## Load and read the Illumina probe files (typically one per experiment containing multiple arrays with 
#a separate control probe profile file output). An EListRaw class is created in limma and should hold ~50,127 
#rows and 6 columns. Each row specificies if it is a control or sample probe. 'xxx'$targets data frame specifies 
#which samples were hybridised to each array. 

#Columns will include source, E (with the expression value for each probe), genes, other, Detection Pval, targets.

#If multiple probe files need to be read use the read.ilmn.targets function, 
#(requires a summary file with probe profiles (one, or two columns if control probe files are available).
#Can also just read one file with corresponding control.


##File with names of files with data
#Create a tab-delimited file with two columns ("files" header= sample probe profiles and 
#"ctrlfiles" header= control profiles):
targets_file_input <- "targets_file.txt"
targets_file_input

#File with meta-data (group membership, treatment status, replicate number, etc.)
# TO DO: pass as parameter or leave out as can be done at differential expression step.
#membership_input_file <- 'GEx_BEST-D_targets_file.txt'
#membership_input_file <- 'BEST-D_lab_kit_IDs_for_array_QC.tsv' # This contains all the IDs and comes from the R script merging pheno and lab kit IDs.

membership_input_file <- 'metadata_file.txt'
membership_input_file

# TO DO: pass as file only if exists, otherwise warn and continue:
# List of samples that failed preliminary QC at core facility:
#failed_QC_input_file <- c('120005280', '120005133', '120000312', '120005211', '120000131', '120005098', 
#                          '120000272', '120005145', 'misload')
failed_QC_input_file <- ''
failed_QC_input_file

## Use limma to read and process array data:
targets_file <- readTargets(targets_file_input)
head(targets_file)


#Then read all files:
# TO DO: check this is true: If 'ctrlfiles=NULL' is specified it will not create the $genes$Status column which
# contains the seven types of control probes (important for QC and normalisation). 

read_files <- read.ilmn.targets(targets=targets_file, probeid="PROBE_ID", 
                                annotation=c("SYMBOL", "SEARCH_KEY", "ILMN_GENE", "CHROMOSOME", "DEFINITION", 
                                             "SYNONYMS"), expr="AVG_Signal", other.columns="Detection Pval", 
                                sep="\t", verbose=TRUE)

# read_files <- read.ilmn(files = 'discovery_Davenport_sepsis_Jan2016_raw_270.txt', probeid = 'Probe_ID', expr = 'AVG_Signal.',
#                         other.columns="Detection.", sep = '\t')


# Control probes are labelled: negative, biotin, labeling, cy3_hyb, housekeeping, high_stringency_hyb or low_stringency_hyb.

#############################


#############################

## Read phenotype file and metadata

# Read in the file with the sample membership information (group membership, replicate number, 
# treatment status, etc.) to be able to create a design and contrast (done above for PCA plots).

# For BEST-D: Cross kit_ids from 03_lab_kits file with column names (ie kit id labels) from expression data, 
# this to extract only samples measured and have same number of samples.


#Save sample names to cross with phenotype IDs:
membership_file <- read.csv(membership_input_file, header=TRUE, row.names=1, sep='\t')

dim(membership_file)
head(membership_file)
tail(membership_file)
summary(membership_file)

# TO DO: change this to eg grouping1, or remove, keep separate as not needed here.
# Add columns needed for PCA and other plots and groupings:
#membership_file$treatment <- ifelse(membership_file$arm == 0 & membership_file$visit_type == 'FinalVisit',  'treated_4000', 
#                                    ifelse(membership_file$arm == 1 & membership_file$visit_type == 'FinalVisit',  'treated_2000',
#                                           'untreated'))

#Get sample names (stored as columns in the limma read files):
#sample_IDs_clean <- colnames(read_files_cleaned_QC)

#Cross IDs:
#matched_IDs <- which(rownames(membership_file) %in% sample_IDs_clean)

#Get file with no NAs and only matched samples to those in the Expression Set object:
#membership_file_cleaned <- membership_file[matched_IDs,]

#head(membership_file_cleaned)
#tail(membership_file_cleaned)

#############################


#############################

# 2) Preprocessing and quality assessment

# If available/possible, run per array signal-to-noise value plot (95th percentile of signal divided by 
# the 5th percentile); spatial plots of the intensities across array surface to detect array artefacts; 
# between sample comparison with bloxplots of intensities; and MDS plot to determine between sample differences 
# (biological, batches, etc.) (see Ritchie et al. 2014 Fig 2 for examples)

# Positive controls can be used to identify suspect arrays.
# Negative control probes, which measure background signal on each array, can be used to assess the proportion 
# of expressed probes that are present in a given sample (Beadarray vignette).


#a) Descriptives

#Look at summary information of the probe profiles:

#Array dimensions:
#Number of arrays match number of samples expected? ~570
#Number of rows is expected number of probes? eg ~48,000
dim(read_files)
#All files were included and read?
read_files$targets
read_files$E[1:5, ]

#Check objects class and other attributes:
class(read_files)
class(read_files$E)
str(read_files)
str(read_files$E)

#Number of negative and regular probes. Illumina BeadChip arrays contain 750~1600 negative control probes:
table(read_files$genes$Status)

#View expression values for first 5 columns, first 10 samples:
options(digits=5)
head(read_files$E[1:5,1:10])
#All samples:
head(read_files$other$'Detection Pval')
length(read_files$other$Detection)
length(read_files$other$'Detection Pval')

#Explore contents of file:
dim(read_files)
class(read_files)
summary(read_files)
head(read_files$source)
head(read_files$E)[1:5,1:5]
range(read_files$E)
median(read_files$E)
head(read_files$genes)
summary(read_files$other)
head(read_files$other$'Detection Pval')[1:5,1:5]
range(read_files$other$'Detection Pval')
read_files$targets


#See p-values for detection, these test whether each probe is more intense than the negative control probes.
#Small values indicate that the probe corresponds to true expression:

range(read_files$other$Detection)
mean(read_files$other$Detection)
median(read_files$other$Detection)


#Boxplots of intensities to assess dynamic range from each sample and identify outliers from signal distributions. 
#The intensities vary from about 5 to 14 on the log2 scale:

#Boxplots for x number of samples (run loop for random sets? Plot all separately?):
intensities_plot <- ("boxplots_of_intensities.png")
png(intensities_plot, width = 4, height = 4, units = 'in', res = 300)
boxplot(log2(read_files$E[, 1:20]),range=0,ylab="log2 intensity", xlab="Array")
dev.off()
#All samples (too many to visualise for BEST-D):
#boxplot(log2(read_files$E),range=0,ylab="log2 intensity of probe intensities", xlab="Array")


#Separate boxplots of regular probes and control probes to highlight unusual samples:
#Regular probes:
png("boxplot_intensities_regular_probes.png", width = 4, height = 4, units = 'in', res = 300)
boxplot(log2(read_files$E[read_files$genes$Status == "regular", ]), 
        range= 0, las = 2, xlab = "", ylab =expression(log[2](intensity)), main = "Regular probes")
dev.off()

#Negative probes:
length(which(read_files$genes$Status == 'negative'))
png("boxplot_intensities_negative_probes.png", width = 4, height = 4, units = 'in', res = 300)
boxplot(log2(read_files$E[read_files$genes$Status == "negative", ]), 
        range= 0, las = 2, xlab = "", ylab = expression(log[2](intensity)), main = "Negative control probes")
dev.off()


#MDS plots:
png("MDS_of_intensities.png", width = 4, height = 4, units = 'in', res = 300)
#'Multidimensional scaling plot of probe intensities')
plotMDS(read_files$E, pch=1)
dev.off()

#The 'propexpr' function estimates the proportion of expressed probes in each array by comparing the empirical 
#intensity distribution of the negative control probes with that of the regular probes. A mixture model is fitted 
#to the data from each array to infer the intensity distribution of expressed probes and estimate the expressed 
#proportion.

#Get proportion of expressed probes and descriptives:
proportion <- propexpr(read_files)
head(proportion)
length(proportion)
range(proportion)
mean(proportion)
median(proportion)
quantile(proportion)

png("boxplot_of_propexpr.png", width = 4, height = 4, units = 'in', res = 300)
boxplot(proportion, ylab='Proportion of expressed probes', xlab='All samples')
dev.off()

#Arrays with low expression (which to exclude? what value is min?):
which_low <- which(proportion <= 0.25)
length(which_low)
proportion[which_low]

which_med <- which(proportion > 0.25 & proportion <= 0.75)
length(which_med)
proportion[which_med]

which_high <- which(proportion > 0.75)
length(which_high)
proportion[which_high]

#Which samples to compare if any? Random subsets? Null hypothesis is that 
#t.test(proportion[-which_low], proportion[which_low])
#############################


#############################
## Call arrayQualityMetrics
#The arrayQualityMetrics package collates quality assessment plots for summarized data created by beadarray to 
#identify outlier arrays:

# Transform EListRaw class to ExpressionSet class (as more general and accepted by more downstream methods):
# Different methods, check: 
# ?ExpressionSet
# object<-new(Class=, exprs=as.matrix(normalised_expressed))
# ?new, ?coerce, ?as 
# as(as_matrix, 'ExpressionSet')
# read_files_cleaned_QC_eset <- coerce(read_files_cleaned_QC, to='eset', strict=TRUE)
# convert(read_files_cleaned_QC, 'ExpressionSet')


# This extracts expression values (with probes as rownames) and sample names (column headers) and returns as a matrix:
raw_as_matrix <- as.matrix(read_files, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
class(raw_as_matrix)
dim(raw_as_matrix)
head(raw_as_matrix)[1:5,1:5]
head(colnames(raw_as_matrix))
head(rownames(raw_as_matrix))


# Remove duplicated probes from expression matrix (these seem to be the negative control probes):
# TO DO: if duplicated is 0 then the next call will not error but result in all rows being deleted:
duplicated_probes <- which(duplicated(rownames(raw_as_matrix)))
duplicated_probes

# If probes have duplicated names but for some reason need to be kept: rownames(as_matrix) = make.names(rownames(as_matrix), unique=TRUE)

# Output clean expression matrix:
#raw_as_matrix_dedup <- raw_as_matrix 
raw_as_matrix_dedup <- raw_as_matrix[-duplicated_probes,]
head(raw_as_matrix_dedup)
dim(raw_as_matrix_dedup)

# Convert expression matrix to expression set:
raw_minimal_eset <- ExpressionSet(assayData=raw_as_matrix_dedup)
class(raw_minimal_eset)
dim(raw_minimal_eset)

# Call arrayQualityMetrics before processing data
# This outputs to a folder in the working directory:
arrayQualityMetrics_raw_data <- arrayQualityMetrics(expressionset=raw_minimal_eset, 
                    outdir=paste('arrayQualityMetrics_raw_data.dir'), force=TRUE, do.logtransform=TRUE)

#############################


#############################
# Determine which samples to exclude and output list of sample names for next script (ie from b onwards), save image for 
# loading and list to exclude

# ArrayQualityMetrics (aqm) outputs many plots and tables and allows programmatic access and manipulation of the reporting:
# http://www.bioconductor.org/packages/release/bioc/vignettes/arrayQualityMetrics/inst/doc/aqm.pdf
# There is no agreed way of what to exclude, aqm ouputs a list of arrays that it considers outliers:


#View(arrayQualityMetrics_raw_data$arrayTable)

# Identify where arrayQualityMetrics stores its outlier counts and identification:
head(arrayQualityMetrics_raw_data$arrayTable)

# Transform to data frame:
aqm_outliers <- data.frame(arrayQualityMetrics_raw_data$arrayTable)

# Transform x's to numbers:
aqm_outliers$recode1[aqm_outliers$X.a.href...hm...1..a. =="x"] <- as.numeric('1')
aqm_outliers$recode1[aqm_outliers$X.a.href...hm...1..a. ==" "] <- as.numeric('0')
aqm_outliers$recode2[aqm_outliers$X.a.href...box...2..a. =="x"] <- as.numeric('1')
aqm_outliers$recode2[aqm_outliers$X.a.href...box...2..a. ==" "] <- as.numeric('0')
aqm_outliers$recode3[aqm_outliers$X.a.href...ma...3..a. =="x"] <- as.numeric('1')
aqm_outliers$recode3[ aqm_outliers$X.a.href...ma...3..a. ==" "] <- as.numeric('0')

# R didn't like the variable types so transformed data types again:
aqm_outliers_numeric <- transform(aqm_outliers, recode1= as.numeric(recode1), recode2= as.numeric(recode2), recode3= as.numeric(recode3))

# Seems to have worked:
str(aqm_outliers_numeric)

# Add the columns counting the number of times the sample was identified as an outlier:
aqm_outliers_numeric$count <- rowSums(x=aqm_outliers_numeric[,6:8], na.rm=TRUE)
head(aqm_outliers_numeric)

# Samples identified as outliers zero, one, two or three times:
length(which(aqm_outliers_numeric$count == 0))
length(which(aqm_outliers_numeric$count == 1))
length(which(aqm_outliers_numeric$count == 2))
length(which(aqm_outliers_numeric$count == 3))
length(aqm_outliers_numeric$count)

# Samples with two or more counts as outliers:
# TO DO: pass as argument: This isn't actually used, see below 'failed_aqm'
#aqm_to_exclude <- as.numeric(2)
failed_aqm_1 <- which(aqm_outliers_numeric$count == 1)
failed_aqm_2 <- which(aqm_outliers_numeric$count >= 2)
failed_aqm_3 <- which(aqm_outliers_numeric$count == 3)

length(failed_aqm_1)
length(failed_aqm_2)
length(failed_aqm_3)
length(aqm_outliers_numeric[failed_aqm_2,][,1])
length(aqm_outliers_numeric[failed_aqm_3,][,1])

aqm_outliers_numeric[failed_aqm_1,]
aqm_outliers_numeric[failed_aqm_2,]
aqm_outliers_numeric[failed_aqm_3,]

# Sample array numbers 
failed_aqm_1_sample_IDs <- aqm_outliers_numeric[failed_aqm_1,][,2]
failed_aqm_2_sample_IDs <- aqm_outliers_numeric[failed_aqm_2,][,2]
failed_aqm_3_sample_IDs <- aqm_outliers_numeric[failed_aqm_3,][,2]

# Check which samples also failed the QC from the hybridisation stage:
failed_QC_input_file %in% failed_aqm_1_sample_IDs
failed_QC_input_file %in% failed_aqm_2_sample_IDs
failed_QC_input_file %in% failed_aqm_3_sample_IDs

# TO DO: plot summaries comparing outliers vs others:
#plotMSD
#plotMeanSd
#plotDensity
#boxplots of intensities
#heatmap of distances
#plotMA

# TO DO: pass according to argument of what to exclude (ie >=2  or 3 (of three) outlier counts from ArrayQualityMetrics (aqm)):
#failed_aqm_2_sample_IDs # Samples with two or three outlier counts
#failed_aqm_3_sample_IDs # If only samples with 3 outlier counts

#failed_aqm <- failed_aqm_2_sample_IDs
failed_aqm <- failed_aqm_3_sample_IDs

# To access aqm and extract data use:
#extract_aqm_data <-prepdata(expressionset=raw_minimal_eset, intgroup=c(), do.logtransform=TRUE)
#str(extract_aqm_data)

# aqm_heatmap_outliers <- arrayQualityMetrics_raw_data$modules$heatmap@outliers@which
# aqm_pca_outliers <- arrayQualityMetrics_raw_data$modules$pca@outliers@which
# aqm_boxplot_outliers <- arrayQualityMetrics_raw_data$modules$boxplot@outliers@which
# aqm_density_outliers <- arrayQualityMetrics_raw_data$modules$density@outliers@which
# aqm_meansd_outliers <- arrayQualityMetrics_raw_data$modules$meansd@outliers@which
# aqm_maplot_outliers <- arrayQualityMetrics_raw_data$modules$maplot@outliers@which

#############################


#############################
## b) Remove failed samples from EListRaw object before continuing with analysis. These will mainly be from the hybridisation QC 
# and samples failing the arrayQualityMetrics package metrics:

#Pass samples that failed QC at the hybridisation stage (usually given by microarray facility):
failed_hybridisation_QC <- failed_QC_input_file
failed_hybridisation_QC

# Pass sample labels from arrayQualityMetrics
FAILED_QC <- c(failed_hybridisation_QC, failed_aqm)

length(FAILED_QC)
length(failed_hybridisation_QC)
length(failed_aqm)
# TO DO: check matching options for vectors of differing length:
count(failed_hybridisation_QC %in% failed_aqm)

FAILED_QC_unique <- unique(FAILED_QC)
FAILED_QC_unique
length(FAILED_QC_unique)
#fix(FAILED_QC_unique)

#Get samples IDs and index numbers from EListRaw object:
array_sample_IDs <- colnames(read_files)
to_extract <- match(FAILED_QC_unique, array_sample_IDs)

#Check indexes match ID:
read_files[0,to_extract]
read_files$E[0,to_extract]

#Get clean EListRaw object:
read_files_cleaned_QC <- read_files[,-to_extract]
dim(read_files)
dim(read_files_cleaned_QC)

#Explore contents of file:
head(read_files_cleaned_QC)
class(read_files_cleaned_QC)
summary(read_files_cleaned_QC)
head(read_files_cleaned_QC$source)
head(read_files_cleaned_QC$E)[1:5,1:5]
range(read_files_cleaned_QC$E)
median(read_files_cleaned_QC$E)
head(read_files_cleaned_QC$genes)
summary(read_files_cleaned_QC$other)
head(read_files_cleaned_QC$other$'Detection Pval')[1:5,1:5]
range(read_files_cleaned_QC$other$'Detection Pval')
read_files_cleaned_QC$targets

## Clean membership file. Get sample names (stored as columns in the limma read files):
sample_IDs_clean <- colnames(read_files_cleaned_QC)

#Cross IDs:
matched_IDs <- which(rownames(membership_file) %in% sample_IDs_clean)

#Get file with no NAs and only matched samples to those in the Expression Set object:
membership_file_cleaned <- membership_file[matched_IDs,]

head(membership_file_cleaned)
tail(membership_file_cleaned)
dim(membership_file_cleaned)
dim(read_files_cleaned_QC)

#Separate boxplots of regular probes and control probes to highlight unusual samples in the cleaned file:
#Regular probes:
png("boxplot_intensities_regular_probes_cleaned_QC.png", width = 4, height = 4, units = 'in', res = 300)
boxplot(log2(read_files_cleaned_QC$E[read_files_cleaned_QC$genes$Status == "regular", ]), 
        range= 0, las = 2, xlab = "", ylab =expression(log[2](intensity)), main = "Regular probes")
dev.off()

#Negative probes:
length(which(read_files_cleaned_QC$genes$Status == 'negative'))
png("boxplot_intensities_negative_probes_clean_QC.png", width = 4, height = 4, units = 'in', res = 300)
boxplot(log2(read_files_cleaned_QC$E[read_files_cleaned_QC$genes$Status == "negative", ]), 
        range= 0, las = 2, xlab = "", ylab = expression(log[2](intensity)), main = "Negative control probes")
dev.off()


#############################
## Call arrayQualityMetrics
# Transform EListRaw class to ExpressionSet class (as more general and accepted by more downstream methods):

# This extracts expression values (with probes as rownames) and sample names (column headers) and returns as a matrix:
raw_cleaned_as_matrix <- as.matrix(read_files_cleaned_QC, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
class(raw_cleaned_as_matrix)
dim(raw_cleaned_as_matrix)
head(raw_cleaned_as_matrix)[1:5,1:5]
head(colnames(raw_cleaned_as_matrix))
head(rownames(raw_cleaned_as_matrix))


# Remove duplicated probes from expression matrix (these seem to be the negative control probes):
duplicated_probes <- which(duplicated(rownames(raw_cleaned_as_matrix)))
duplicated_probes

# If probes have duplicated names but for some reason need to be kept: rownames(as_matrix) = make.names(rownames(as_matrix), unique=TRUE)

# Output clean expression matrix:
#raw_cleaned_as_matrix_dedup <- raw_cleaned_as_matrix
raw_cleaned_as_matrix_dedup <- raw_cleaned_as_matrix[-duplicated_probes,]
head(raw_cleaned_as_matrix_dedup)
dim(raw_cleaned_as_matrix_dedup)

# Convert expression matrix to expression set:
raw_cleaned_minimal_eset <- ExpressionSet(assayData=raw_cleaned_as_matrix_dedup)
class(raw_cleaned_minimal_eset)
dim(raw_cleaned_minimal_eset)


# Call arrayQualityMetrics before processing data and after cleaning raw files
# This outputs to a folder in the working directory:
arrayQualityMetrics_raw_cleaned <- arrayQualityMetrics(expressionset=raw_cleaned_minimal_eset, 
                    outdir=paste('arrayQualityMetrics_raw_cleaned.dir'), force=TRUE, do.logtransform=TRUE)

#############################


#############################
# This next run for aqm is a sanity check as samples to exclude have already been identified.
# TO DO: delete this (it's a copy/paste of pre-QC aqm commands above), and only run summary plots for before/after comparison.

#View(arrayQualityMetrics_raw_cleaned$arrayTable)

# Identify where arrayQualityMetrics stores its outlier counts and identification:
head(arrayQualityMetrics_raw_cleaned$arrayTable)

# Transform to data frame:
aqm_outliers_cleaned <- data.frame(arrayQualityMetrics_raw_cleaned$arrayTable)

# Transform x's to numbers:
aqm_outliers_cleaned$recode1[aqm_outliers_cleaned$X.a.href...hm...1..a. =="x"] <- as.numeric('1')
aqm_outliers_cleaned$recode1[aqm_outliers_cleaned$X.a.href...hm...1..a. ==" "] <- as.numeric('0')
aqm_outliers_cleaned$recode2[aqm_outliers_cleaned$X.a.href...box...2..a. =="x"] <- as.numeric('1')
aqm_outliers_cleaned$recode2[aqm_outliers_cleaned$X.a.href...box...2..a. ==" "] <- as.numeric('0')
aqm_outliers_cleaned$recode3[aqm_outliers_cleaned$X.a.href...ma...3..a. =="x"] <- as.numeric('1')
aqm_outliers_cleaned$recode3[ aqm_outliers_cleaned$X.a.href...ma...3..a. ==" "] <- as.numeric('0')

# R didn't like the variable types so transformed data types again:
aqm_outliers_cleaned_numeric <- transform(aqm_outliers_cleaned, recode1= as.numeric(recode1), 
                                          recode2= as.numeric(recode2), recode3= as.numeric(recode3))

# Seems to have worked:
str(aqm_outliers_cleaned_numeric)

# Add the columns counting the number of times the sample was identified as an outlier:
aqm_outliers_cleaned_numeric$count <- rowSums(x=aqm_outliers_cleaned_numeric[,6:8], na.rm=TRUE)
head(aqm_outliers_cleaned_numeric)

# Samples identified as outliers zero, one, two or three times:
length(which(aqm_outliers_cleaned_numeric$count == 0))
length(which(aqm_outliers_cleaned_numeric$count == 1))
length(which(aqm_outliers_cleaned_numeric$count == 2))
length(which(aqm_outliers_cleaned_numeric$count == 3))
length(aqm_outliers_cleaned_numeric$count)

# Samples with two or more counts as outliers:
# TO DO: pass as argument:
aqm_to_exclude_cleaned <- as.numeric(2)
failed_aqm_cleaned_1 <- which(aqm_outliers_cleaned_numeric$count == 1)
failed_aqm_cleaned_2 <- which(aqm_outliers_cleaned_numeric$count >= 2)
failed_aqm_cleaned_3 <- which(aqm_outliers_cleaned_numeric$count == 3)

length(failed_aqm_cleaned_1)
length(failed_aqm_cleaned_2)
length(failed_aqm_cleaned_3)
length(aqm_outliers_cleaned_numeric[failed_aqm_cleaned_2,][,1])
length(aqm_outliers_cleaned_numeric[failed_aqm_cleaned_3,][,1])

aqm_outliers_cleaned_numeric[failed_aqm_cleaned_1,]
aqm_outliers_cleaned_numeric[failed_aqm_cleaned_2,]
aqm_outliers_cleaned_numeric[failed_aqm_cleaned_3,]

# Sample array numbers 
failed_aqm_cleaned_1_sample_IDs <- aqm_outliers_cleaned_numeric[failed_aqm_1,][,2]
failed_aqm_cleaned_2_sample_IDs <- aqm_outliers_cleaned_numeric[failed_aqm_2,][,2]
failed_aqm_cleaned_3_sample_IDs <- aqm_outliers_cleaned_numeric[failed_aqm_3,][,2]

# Check which ones failed hybridisation QC:
failed_QC_input_file %in% failed_aqm_cleaned_1_sample_IDs
failed_QC_input_file %in% failed_aqm_cleaned_2_sample_IDs
failed_QC_input_file %in% failed_aqm_cleaned_3_sample_IDs
#############################


#############################
#The end:
# Remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes)))

# rm(targets_file, raw_as_matrix, duplicated_probes, array_sample_IDs, to_extract, read_files, membership_file, proportion,
#            sample_IDs_clean, matched_IDs, raw_cleaned_as_matrix, raw_as_matrix_dedup, which_high, which_med, which_low,
#            raw_cleaned_as_matrix_dedup, raw_minimal_eset, working_dir,
#    failed_QC_input_file, targets_file_input, membership_input_file) #arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed

# To save R workspace with all objects to use at a later time:
# Also consider dump(), dput() and dget(), see http://thomasleeper.com/Rcourse/Tutorials/savingdata.html

# TO DO: clean up and save only necessary objects, very slow otherwise:
save.image(R_session_saved_image_full, compress='gzip')

# Or save specific objects:
objects_to_save <- c('read_files_cleaned_QC', 'membership_file_cleaned', 'FAILED_QC_unique', 'failed_aqm')
save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: run the script for background correction, normalisation, etc.

#############################
