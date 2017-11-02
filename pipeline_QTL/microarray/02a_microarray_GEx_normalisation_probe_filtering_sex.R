#############################
# To be run after 02 normalisation and filtering script
# Antonio J Berlanga-Taylor

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_probe_filtering",".txt", sep=""), open='a')
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
#load('R_session_saved_image_read_and_QC.RData', verbose=T)
#load('R_session_saved_image_normalisation_full.RData', verbose=T)
load('R_session_saved_image_normalisation.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_probe_filtering', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:

library(limma)
library(illuminaHumanv4.db)
library(plyr)
#############################


#############################

# Once background correction, normalisation and transformation have been performed run probe filtering (ie to exclude multi-mapping probes, probes over
# SNPs, etc.) carry out the following:

## d) Filtering and creating an eset object

# Filter out probes that are not expressed. Keep probes that are expressed in at least three arrays 
# according to a detection p-value of 5% (Limma vignette case study, p. 108):

normalised
dim(normalised)
class(normalised)
summary(normalised)
range(normalised$E)


normalised$other$'Detection Pval'[1:5,1:5]
#normalised$other$Detection[1:5,1:5] # I think these are the inverse of p-values? ie for GAinS data from ArrayExpress
length(normalised$other$'Detection Pval' < 0.05)
length(which(normalised$other$'Detection Pval' < 0.05) == TRUE)
#length(which(normalised$other$Detection < 0.05) == TRUE)
range(normalised$other$'Detection Pval')
dim(normalised$other$'Detection Pval')

#Set filters:
expressed <- which((rowSums(normalised$other$'Detection Pval' < 0.05) >= 3) == TRUE)
#expressed <- which((rowSums(normalised$other$Detection < 0.05) >= 10) == TRUE)
head(expressed, n=10)
length(expressed)

expressed_2 <- rowSums(normalised$other$'Detection Pval' < 0.05) >=3
#expressed_2 <- rowSums(normalised$other$Detection < 0.05) >=10
head(expressed_2, n=10)
head(which(expressed_2 == TRUE), n=10)
length(which(expressed_2 == TRUE))


#Extract data and create an ExpressionSet (eset, or EList object) at the same time (necessary for linear modelling steps):
normalised_expressed <- normalised[expressed,]
#normalised_expressed <- normalised[-expressed,] # I think Array Express data 'Detection' column is the inverse of the p-value

#Explore contents of file:
head(normalised_expressed$E)
dim(normalised_expressed)
class(normalised_expressed)
summary(normalised_expressed)
head(normalised_expressed$source)
head(normalised_expressed$E)[1:5,1:5]
range(normalised_expressed$E)
head(normalised_expressed$genes)
summary(normalised_expressed$other)
head(normalised_expressed$other$'Detection Pval')[1:5,1:5]
range(normalised_expressed$other$'Detection Pval')
normalised_expressed$targets


###################
## Filter based on annotation
# Check http://www.bioconductor.org/packages/release/data/annotation/manuals/illuminaHumanv4.db/man/illuminaHumanv4.db.pdf
# Could also use Bioconductor's Annotation workflows:
# See http://www.bioconductor.org/help/workflows/annotation/annotation/
# Also download Illumina's annotation file:
# https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanht-12/humanht-12_v4_0_r2_15002873_b.txt.zip
# Merge with file above and exclude all unplaced and non-autosomal probes for example.
# TO DO: The R package illuminaHumanv4.db is known to have issues though:


# Remove probes that aren't annotated with a gene symbol:
annotated_probes <- !is.na(normalised_expressed$genes$SYMBOL)
count(annotated_probes)
head(annotated_probes)
tail(annotated_probes)

un_annotated_probes <- which(annotated_probes == FALSE)
length(un_annotated_probes)
head(un_annotated_probes)
tail(un_annotated_probes)

# Remove un-annotated probes:
normalised_expressed_annotated <- normalised_expressed[annotated_probes,]
#normalised_expressed_annotated <- normalised_expressed # for array data without annotation at this point (eg GAinS)

#Explore contents of file and sanity checks:
head(normalised_expressed_annotated$genes$SYMBOL)
dim(normalised_expressed_annotated)
dim(normalised_expressed)
length(annotated_probes)
count(annotated_probes)

#View(normalised_expressed$genes)
###################


###################
# Remove probes based on quality annotations:
# See: http://www.bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf
# for an example from the BeadArray use cases, p. 27.
illumina_ids  <-  as.character(rownames(normalised_expressed_annotated))
head(illumina_ids)
length(illumina_ids)

probes_by_qual  <- unlist(mget(illumina_ids, illuminaHumanv4PROBEQUALITY, ifnotfound = NA))
head(probes_by_qual)
summary(probes_by_qual)
count(probes_by_qual)

# Plot expression signal by probe quality:
aveSignal_probes_by_qual <- rowMeans(normalised_expressed_annotated$E)
png('Average_signal_probes_by_quality.png', width = 4, height = 4, units = 'in', res = 300)
boxplot(aveSignal_probes_by_qual ~ probes_by_qual)
dev.off()

# TO DO: Confirm low quality probes match to multiple locations and should be removed:
#?illuminaHumanv4REPEATMASK

all_probes  <- unlist(mget(illumina_ids, illuminaHumanv4REPEATMASK, ifnotfound = NA))
#View(all_probes)
summary(all_probes)
all_probes_2 <- !is.na(all_probes)
#View(all_probes_2)
summary(all_probes_2)

query_bad_IDs <- names(which(probes_by_qual == "Bad" & aveSignal_probes_by_qual  > 12))
head(query_bad_IDs)
length(query_bad_IDs)

repeat_mask_bad <- unlist(mget(query_bad_IDs, illuminaHumanv4REPEATMASK, ifnotfound = NA))
#View(repeat_mask_bad)
head(repeat_mask_bad)
length(repeat_mask_bad)

repeat_mask_bad_2 <- !is.na(repeat_mask_bad)
#View(repeat_mask_bad_2)
summary(repeat_mask_bad_2)
head(repeat_mask_bad_2)
length(repeat_mask_bad_2)

multiple_matches_bad <- unlist(mget(query_bad_IDs, illuminaHumanv4SECONDMATCHES, ifnotfound = NA))
head(multiple_matches_bad)
length(multiple_matches_bad)

multiple_matches_bad_2 <- !is.na(multiple_matches_bad)
#View(repeat_mask_bad_2)
summary(multiple_matches_bad_2)
head(multiple_matches_bad_2)
length(multiple_matches_bad_2)

# mget("ILMN_1692145", illuminaHumanv4PROBESEQUENCE)

query_perfect_IDs <- names(which(probes_by_qual == "Perfect" & aveSignal_probes_by_qual  > 12))
head(query_perfect_IDs)
length(query_perfect_IDs)

repeat_mask_perfect <- unlist(mget(query_perfect_IDs, illuminaHumanv4REPEATMASK))
head(repeat_mask_perfect)
length(repeat_mask_perfect)

multiple_matches_perfect <- unlist(mget(query_bad_IDs, illuminaHumanv4SECONDMATCHES))
head(multiple_matches_perfect)
length(multiple_matches_perfect)


# Remove probes of low quality:
remove_bad_probes  <- probes_by_qual == "No match" | probes_by_qual == "Bad"
count(remove_bad_probes)
normalised_expressed_annotated_qual <- normalised_expressed_annotated[!remove_bad_probes, ]
dim(normalised_expressed_annotated)
dim(normalised_expressed_annotated_qual)

# Check counts match:
count(!is.na(normalised_expressed_annotated_qual$genes$SYMBOL))
length(rownames(normalised_expressed_annotated_qual))
###################

###################
# Remove probes that overlap SNPs:
probes_by_SNPs  <- !is.na(unlist(mget(as.character(rownames(normalised_expressed_annotated_qual)), 
                                      illuminaHumanv4OVERLAPPINGSNP, ifnotfound = NA)))
head(probes_by_SNPs)
summary(probes_by_SNPs)
#View(probes_by_SNPs)

normalised_expressed_annotated_qual_noSNPs <- normalised_expressed_annotated_qual[!probes_by_SNPs, ]
dim(normalised_expressed_annotated_qual)
dim(normalised_expressed_annotated_qual_noSNPs)
###################


###################
# Keep only probes that have Entrez IDs:
## Get ENTREZ ID mappings:
dim(normalised_expressed_annotated_qual_noSNPs)
probes_by_ENTREZID  <- unlist(mget(as.character(rownames(normalised_expressed_annotated_qual_noSNPs)), 
                                   illuminaHumanv4ENTREZID, ifnotfound = NA))

head(probes_by_ENTREZID)
summary(probes_by_ENTREZID)
str(probes_by_ENTREZID)
probes_without_ENTREZID  <- is.na(probes_by_ENTREZID)
count(probes_without_ENTREZID)
length(which(probes_without_ENTREZID))
#View(probes_by_ENTREZID)

# Remove probes without Entrez IDs:
normalised_expressed_annotated_qual_noSNPs_noID <- normalised_expressed_annotated_qual_noSNPs[!probes_without_ENTREZID, ]
dim(normalised_expressed_annotated_qual_noSNPs)
dim(normalised_expressed_annotated_qual_noSNPs_noID)

# TO DO: Check p. 12 illuminaHumanv4listNewMappings:
# https://www.bioconductor.org/packages/3.3/data/annotation/manuals/illuminaHumanv4.db/man/illuminaHumanv4.db.pdf
# Consider using illuminaHumanv4ENTREZREANNOTATED instead of entrezid re-annotated instead of entrez id.
# Other new mappings are in use but double check.
# The below is currently not used:
probes_by_ENTREZID_RE  <- unlist(mget(illumina_ids, illuminaHumanv4ENTREZREANNOTATED, ifnotfound = NA))
head(probes_by_ENTREZID_RE)
summary(probes_by_ENTREZID_RE)
str(probes_by_ENTREZID_RE)
count(!is.na(probes_by_ENTREZID_RE))
#View(probes_by_ENTREZID_RE)

###################


###################
# Remove probes from X and Y chromosomes and unmapped or with missing values:
## Get chromosome mappings:
# illuminaHumanv4CHR can show duplicate mappings, see: 
# http://rgm.ogalab.net/RGM/R_rdfile?f=illuminaHumanv4.db/man/illuminaHumanv4CHR.Rd&d=R_BC
dim(normalised_expressed_annotated_qual_noSNPs_noID)
probes_by_CHR  <- na.omit(unlist(mget(as.character(rownames(normalised_expressed_annotated_qual_noSNPs_noID)),
                                   illuminaHumanv4CHR, ifnotfound = NA)))

#View(probes_by_CHR)
head(probes_by_CHR)
tail(probes_by_CHR)
length(probes_by_CHR)
summary(probes_by_CHR)
str(probes_by_CHR)
count(probes_by_CHR)

# Lengths differ between chromosome annotation and normalised filtered file (multi-mapping probes when illuminaHumanv4 was written).
# Probe IDs aren't duplicated though and it seems like each probe has one chromosome (for each row).

# Convert files to data frames and merge:
probes_by_CHR <- as.data.frame(probes_by_CHR)
probes_by_CHR['Probe_Id'] <- NA
probes_by_CHR['Probe_Id'] <- row.names(probes_by_CHR)
head(probes_by_CHR)
length(which(duplicated(probes_by_CHR$Probe_Id)))
length(which(!duplicated(probes_by_CHR$Probe_Id)))

# Convert elist object to dataframe (this loses some information but retains what is
# needed (probe expression levels, probe ID, illumina ID, annotations, etc.):
# Informartion lost is source, detection P-values (already used, kept in $other), and $targets (
# which stores the file names read originally for read.illm function):
str(normalised$source)
str(normalised$genes)
str(normalised$other)
str(normalised$targets)

df_normalised_expressed_annotated_qual_noSNPs_noID <- as.data.frame(normalised_expressed_annotated_qual_noSNPs_noID)
#View(df_normalised_expressed_annotated_qual_noSNPs_noID)
head(df_normalised_expressed_annotated_qual_noSNPs_noID)
df_normalised_expressed_annotated_qual_noSNPs_noID['Probe_Id'] <- row.names(df_normalised_expressed_annotated_qual_noSNPs_noID)
head(df_normalised_expressed_annotated_qual_noSNPs_noID$Probe_Id)
length(df_normalised_expressed_annotated_qual_noSNPs_noID$Probe_Id)

# Merge both data frames:
df_normalised_expressed_annotated_qual_noSNPs_noID_chr <- merge(df_normalised_expressed_annotated_qual_noSNPs_noID, probes_by_CHR, by = 'Probe_Id')
df_normalised_expressed_annotated_qual_noSNPs_noID_chr[1:5, 1:5]
head(df_normalised_expressed_annotated_qual_noSNPs_noID_chr)
dim(df_normalised_expressed_annotated_qual_noSNPs_noID_chr)
summary(df_normalised_expressed_annotated_qual_noSNPs_noID_chr$probes_by_CHR)

# Check if annotations between Illumina and illuminaHumanv4.db match:
# Illumina annotations here aret hose from the array itself (ie not the downloaded file from Illumina's webpage, that's another
# option to check).
# This is just sanity as issues have been raised in the past regarding different types 
# of annotations (probe quality, gene symbols, etc.).
illumina_chr <- as.list(df_normalised_expressed_annotated_qual_noSNPs_noID_chr$CHROMOSOME)
head(illumina_chr)
length(illumina_chr)

illuminaHumanv4_chr <- as.list(df_normalised_expressed_annotated_qual_noSNPs_noID_chr$probes_by_CHR)
head(illuminaHumanv4_chr)
length(illuminaHumanv4_chr)

identical(illumina_chr, illuminaHumanv4_chr)
# This test gives false (only takes one mismatch though).

# Check what's happening:
chr_comparison <- df_normalised_expressed_annotated_qual_noSNPs_noID_chr[, c('Probe_Id', 'CHROMOSOME', 'probes_by_CHR')]
head(chr_comparison)
tail(chr_comparison)

# Illumina's annotation has many missing values:
length(which(is.na(chr_comparison$CHROMOSOME)))

# Convert blanks to NAs:
chr_comparison$CHROMOSOME[chr_comparison$CHROMOSOME == ''] <- NA
chr_comparison$probes_by_CHR[chr_comparison$probes_by_CHR == ''] <- NA
#View(chr_comparison)
length(which(is.na(chr_comparison$CHROMOSOME)))
length(which(is.na(chr_comparison$probes_by_CHR)))
count(chr_comparison$CHROMOSOME)
count(chr_comparison$probes_by_CHR)

# TO DO: Check these differences (chromosome counts seem close but don't match...).
# Keeping illuminaHumanv4 package annotations for now (no clear answers on what to use and 
# haven't tested):

# Remove probes mapping to sex chromosomes and unmapped:
# Rename gene expression file to make it manageable:
df_normalised_filtered <- df_normalised_expressed_annotated_qual_noSNPs_noID_chr
remove_sex_probes  <- which(df_normalised_filtered$probes_by_CHR == 'X' | df_normalised_filtered$probes_by_CHR == 'Y' | 
                              df_normalised_filtered$probes_by_CHR == 'Un')

head(which(df_normalised_filtered$probes_by_CHR == 'X'))
head(which(df_normalised_filtered$probes_by_CHR == 'Y'))
head(which(df_normalised_filtered$probes_by_CHR == 'Un'))
# TO DO: these indexes will cause errors when using other files. Convert or comment out:
#df_normalised_filtered[c(39, 930, 16825), c('Probe_Id', 'CHROMOSOME', 'probes_by_CHR')]
#df_normalised_filtered[c(17484, 16809, 16825), c('Probe_Id', 'CHROMOSOME', 'probes_by_CHR')]

#View(remove_sex_probes)
head(remove_sex_probes)
length(remove_sex_probes)
dim(df_normalised_filtered)

df_normalised_filtered <- df_normalised_filtered[-remove_sex_probes, ]
dim(normalised)
dim(normalised_expressed)
dim(normalised_expressed_annotated)
dim(normalised_expressed_annotated_qual)
dim(normalised_expressed_annotated_qual_noSNPs)
dim(normalised_expressed_annotated_qual_noSNPs_noID)
dim(df_normalised_filtered)

which(df_normalised_filtered$probes_by_CHR == 'X')
which(df_normalised_filtered$probes_by_CHR == 'Y')
which(df_normalised_filtered$probes_by_CHR == 'Un')
#df_normalised_filtered[c(39, 930, 16825), c('Probe_Id', 'CHROMOSOME', 'probes_by_CHR')]
###################


##########
# TO DO / CHECK / CLEAN:
#Read Illumina's annotation file downloaded from webpage (alternative to illuminaHumanv4.db package).
# The arrays themselves have columns with annotations. I haven't checked if these are the same and/or which is better
# (illuminaHumanv4 package). All seem to have reported issues.
# Another option is the following:
# http://biorxiv.org/content/biorxiv/early/2015/05/21/019596.full.pdf

# Read Illumina's annotation file:
#illumina_annotation_file <- as.data.frame(read.csv('/ifs/projects/proj043/analysis.dir/HumanHT-12_V4_0_R2_15002873_B.txt', header = TRUE,
#                                     sep = '\t', skip = 8, na.strings=c(""," ","NA")))

#View(illumina_annotation_file)
#class(illumina_annotation_file)
#head(illumina_annotation_file)
#tail(illumina_annotation_file)
#dim(illumina_annotation_file)

#summary(illumina_annotation_file$Chromosome)
#length(which(is.na(illumina_annotation_file$Chromosome)))

#illumina_annotation_file_chr <- illumina_annotation_file[which(!is.na(illumina_annotation_file$Chromosome)), ]
#dim(illumina_annotation_file)
#dim(illumina_annotation_file_chr)
#View(illumina_annotation_file_chr)

#chrs_to_keep <- as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
#                               11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
#                               21, 22))
#chrs_to_keep
#illumina_annotation_file_chr_autosome <- subset(illumina_annotation_file_chr, 
#                                                illumina_annotation_file_chr$Chromosome %in% chrs_to_keep)

#count(illumina_annotation_file_chr$Chromosome)
#length(which(illumina_annotation_file_chr$Chromosome == 'X'))
#length(which(illumina_annotation_file_chr$Chromosome == 'Y'))
#dim(illumina_annotation_file_chr_autosome)

#row.names(illumina_annotation_file_chr_autosome) <- illumina_annotation_file_chr_autosome$Probe_Id
#View(illumina_annotation_file_chr_autosome)
#head(illumina_annotation_file_chr_autosome$Probe_Id)

#normalised_expressed_annotated_qual_noSNPs_noID_nosex <- merge(normalised_expressed_annotated_qual_noSNPs_noID, illumina_annotation_file_chr_autosome, all.x = TRUE)
#Probe_Id <- which(colnames(illumina_annotation_file_chr_autosome) == 'Probe_Id')
#Chromosome <- which(colnames(illumina_annotation_file_chr_autosome) == 'Chromosome')

#probes_to_keep_chr <- illumina_annotation_file_chr_autosome[, c(Probe_Id, Chromosome)]
#head(probes_to_keep_chr)
#dim(probes_to_keep_chr)
############



##############
## Rename object for downstream analysis:
# Rename row IDs:
row.names(df_normalised_filtered) <- df_normalised_filtered$Probe_Id
head(df_normalised_filtered)
#df_normalised_filtered[1:5, 1:5]

# Annotations and expression data:
normalised_filtered_annotated <- df_normalised_filtered

# Keep only expression data for plotting. First 8 columns are annotation, last column is 'probes_by_CHR:
df_normalised_filtered[1:5, 8:12]
normalised_filtered <- df_normalised_filtered[, c(9:(ncol(df_normalised_filtered)-1))]
head(normalised_filtered)

dim(normalised_filtered_annotated)
class(normalised_filtered_annotated)
head(normalised_filtered_annotated$SYMBOL)
#View(normalised_filtered_annotated)

dim(normalised_filtered)
class(normalised_filtered)
str(normalised_filtered)
head(normalised_filtered)
#View(normalised_filtered)

# Sanity checks:
count(!is.na(normalised_filtered_annotated$SYMBOL))
count(!is.na(normalised_filtered_annotated$probes_by_CHR))
count(normalised_filtered_annotated$probes_by_CHR)
count(!is.na(normalised_filtered_annotated$Probe_Id))
count(normalised_filtered_annotated$Status)

#############################



#############################

## Visual inspection of normalised and filtered dataset: 
# TO DO: I changed the normalised_filtered object from an Elist to a dataframe, downstream will now error.
# TO DO: move PCA analysis to a different file separate to probe filtering?
#Plot log2 intensity of regular probes normalised:
boxplot_normalised_filtered_intensities <- ("regular_probes_normalised.png")
png(boxplot_normalised_filtered_intensities, width = 4, height = 4, units = 'in', res = 300)
boxplot(normalised_filtered,range= 0, ylab =expression(log[2](intensity)), 
        las = 2, xlab = "", main = "Regular probes, normalised")
dev.off()

#Plot expressed probes in a multi-dimensional scaling plot:
plot_normalised_filtered <- ("MDS_normalised_filtered.png")
png(plot_normalised_filtered, width = 4, height = 4, units = 'in', res = 300)
plotMDS(normalised_filtered, pch=1)
dev.off()

plot_normalised_filtered_by_targets <- ("MDS_normalised_filtered_by_targets.png")
png(plot_normalised_filtered_by_targets, width = 4, height = 4, units = 'in', res = 300)
plotMDS(normalised_filtered, pch=1, labels=normalised_filtered$targets)
dev.off()

# png("MDS_normalised_filtered_by_sample.png", width = 4, height = 4, units = 'in', res = 300)
# plotMDS(normalised_filtered, pch=1, labels=normalised_filtered$targets) 
# dev.off()

# Plot PCA of normalised samples:
# Compute the PCs, first transpose the expression values and ignore the first column, then run the PCs:
pca_normalised_filtered <- prcomp(t(normalised_filtered), center=TRUE, scale=TRUE)

# Obtain values for all PCs output:
pc <- data.frame(round(pca_normalised_filtered$x, 2))
pc$sample_id <- rownames(pc) 
source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
pc <- pc[, moveme(names(pc), 'sample_id first')]
names(pc)[1:10]
class(pc)
write.table(pc, 'principal_components_normalised_filtered.tsv', quote = FALSE,
            sep = '\t', row.names = FALSE)

# Explore dimensions and plot first 10 or so components:
dim(pc)
dim(normalised_filtered)
str(pca_normalised_filtered)
head(pc)
# pc[1:5, 1:5]
summary(pca_normalised_filtered)


# Plot PCA results:
plot_PCA_normalised_filtered <- ('plot_PCA_normalised_filtered.png')
png(plot_PCA_normalised_filtered, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
# Histogram of first x PCs:
plot(pca_normalised_filtered, main='Normalised expression values')
# Scatterplot of PC1 and PC2:
biplot(pca_normalised_filtered, main='Normalised expression values')

par(mfrow=c(1,1))
dev.off()

# Check how much of the variance in gene expression values is explained by the first x PCs:
# sum(pca_normalised_filtered$sdev[1:10]^2)/length(normalised_filtered[1,])


# Run PCA analysis by groups of interest: TO DO: Cross files and IDs first
head(membership_file_cleaned)
tail(membership_file_cleaned)
pc_data <- data.frame(pca_normalised_filtered$x[,1:13])
str(pc_data)
head(pc_data)

pca_by_groups <- data.frame(merge(membership_file_cleaned, pc_data, by='row.names'))
head(pca_by_groups)
dim(pca_by_groups)
dim(pc_data)
# View(pc_data)
# View(pca_by_groups)
#pc_data['120000222',]
head(arrange(pc_data, PC1), 10)
head(arrange(pca_by_groups, PC1), 10)

plot_PCA_normalised_filtered_by_groups_1 <- ('plot_PCA_normalised_filtered_by_groups_1.png')
png(plot_PCA_normalised_filtered_by_groups_1, width = 13, height = 13, units = 'in', res = 300)
p1 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$arm)) + theme(legend.position="bottom")
p2 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$visit_type)) + theme(legend.position="bottom")
p3 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p4 <- qplot(x=PC2, y=PC3, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

plot_PCA_normalised_filtered_by_groups_2 <- ('plot_PCA_normalised_filtered_by_groups_2.png')
png(plot_PCA_normalised_filtered_by_groups_2, width = 13, height = 13, units = 'in', res = 300)
p5 <- qplot(x=PC3, y=PC4, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p6 <- qplot(x=PC4, y=PC5, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p7 <- qplot(x=PC5, y=PC6, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p8 <- qplot(x=PC6, y=PC7, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
grid.arrange(p5, p6, p7, p8, ncol=2)
dev.off()

plot_PCA_normalised_filtered_by_groups_3 <- ('plot_PCA_normalised_filtered_by_groups_3.png')
png(plot_PCA_normalised_filtered_by_groups_3, width = 13, height = 13, units = 'in', res = 300)
p9 <- qplot(x=PC7, y=PC8, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p10 <- qplot(x=PC8, y=PC9, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p11 <- qplot(x=PC9, y=PC10, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p12 <- qplot(x=PC10, y=PC11, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
grid.arrange(p9, p10, p11, p12, ncol=2)
dev.off()


# Plot dendrogram
# Prepare hierarchical cluster:
correlation <- cor(normalised_filtered, method='pearson')
head(correlation)
# View(correlation)

hc <- flashClust(dist(correlation))
str(hc)
head(hc)

# Convert to a dendogram object:
dendro <- as.dendrogram(hc)

# Add colours by groups
colourCodes <- c(treated_4000="red", treated_2000="green", untreated="blue")

# Assign labels of dendrogram object with new colors:
labels_colors(dendro) <- colourCodes[membership_file_cleaned$treatment][order.dendrogram(dendro)]
head(colourCodes[membership_file_cleaned$treatment][order.dendrogram(dendro)])
head(order.dendrogram(dendro))
head(membership_file_cleaned$treatment)

# Plot simple dendrogram, labels at the same level:
plot_dendrogram_normalised_filtered <- ('plot_dendrogram_normalised_filtered.png')
png(plot_dendrogram_normalised_filtered, width = 13, height = 13, units = 'in', res = 300)
par(mfrow=c(1,2), cex = 1)
plot(dendro, hang = -1)
# Zoom in to a sub tree:
plot(dendro[[1]], horiz = TRUE)
par(mfrow=c(1,1))
dev.off()

# Plot correlation between samples in a heatmap:
cor_normalised_filtered <- melt(correlation)
head(cor_normalised_filtered)
summary(cor_normalised_filtered)
plot_heatmap_correlations <- ('plot_heatmap_correlations_normalised_filtered.png')
png(plot_heatmap_correlations, width = 12, height = 12, units = 'in', res = 300)
p1 <- ggplot(cor_normalised_filtered, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient2(low = 'blue', high = 'yellow')
grid.arrange(p1, ncol=1)
dev.off()


# TO DO: 
# Plot expression values in a heatmap:
#plot_heatmap_correlations <- ('plot_heatmap_correlations_normalised_filtered.png')
#png(plot_heatmap_correlations, width = 12, height = 12, units = 'in', res = 300)
#par(mfrow=c(1,2))
#heatmap(normalised_filtered[,1:10])
#heatmap(normalised_filtered[1:100,1:100], ColSideColors=colourCodes)
# p1 <- ggplot(cor_normalised_filtered, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient2(low = 'blue', high = 'yellow')
# grid.arrange(p1, ncol=1)
#par(mfrow=c(1,1))
#dev.off()

# # Try again with heatmap.2:
# # http://sebastianraschka.com/Articles/heatmaps_in_r.html
# 
# # Create colour palette from red to green:
# my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
# 
# # (optional) defines the color breaks manually for a "skewed" color transition:
# col_breaks = c(seq(-1,0,length=100),  # for red
#                seq(0,0.8,length=100),              # for yellow
#                seq(0.8,1,length=100))              # for green
# 
# # creates a 5 x 5 inch image
# png('heatmap_normalised_filtered_heatmap2.png',    # create PNG for the heat map        
#     width = 5*300,        # 5 x 300 pixels
#     height = 5*300,
#     res = 300,            # 300 pixels per inch
#     pointsize = 8)        # smaller font size
# 
# heatmap.2(normalised_filtered[1:100, 1:100], 
#           cellnote = normalised_filtered[1:100, 1:100],  # same data set for cell labels
#           main = "Correlation", # heat map title
#           notecol="black",      # change font color of cell labels to black
#           density.info="none",  # turns off density plot inside color legend
#           #          trace="none",         # turns off trace lines inside the heat map
#           margins =c(12,9),     # widens margins around plot
#           col=my_palette,       # use on color palette defined earlier 
#           #          breaks=col_breaks,    # enable color transition at specified limits
#           dendrogram="column",     # only draw a column dendrogram
#           #          Colv="NA")            # turn off column clustering
# )
# dev.off()               # close the PNG device


# hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# t <- ExpressionSet(e, AnnotatedDataFrame(tab))
# rv <- rowVars(exprs(t))
# idx <- order(-rv)[1:40]
# heatmap(exprs(t)[idx, ], col = hmcol)


#############################


#############################
# The end:
# Remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes))[1:10])

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))


# Save files required for other packages/programmes:
head(normalised_filtered) # expression values only
head(normalised_filtered_annotated) # expression values plus annotations
#head(membership_file_cleaned)

write.table(normalised_filtered_annotated, 'normalised_filtered_annotated.tab', sep='\t', 
            quote=FALSE, na='NA', col.names=NA, row.names=TRUE)
write.table(normalised_filtered, 'normalised_filtered_expression_values.tab', sep='\t', 
            quote=FALSE, na='NA', col.names=NA, row.names=TRUE)
#col.names=NA, row.names=TRUE is required otherwise it skips the first position and 
#generates a tab file with the first column on rownames (screws up the header).
#write.table(membership_file_cleaned, 'membership_file_cleaned_all.tab', sep='\t', quote=FALSE, na='NA', 
#            col.names=NA, row.names=TRUE)

write.csv(normalised_filtered, 'normalised_filtered_annotated.csv', quote=FALSE, na='NA')
write.csv(normalised_filtered, 'normalised_filtered_expression_values.csv', quote=FALSE, na='NA')

# To save R workspace with all objects to use at a later time:
#save.image(file=R_session_saved_image_full, compress='gzip')

objects_to_save <- (c('normalised_filtered_annotated', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# Print session information:
sessionInfo()
     
q()

# Next: run the script for differential gene expression.
#############################
