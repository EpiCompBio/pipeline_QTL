biomart script

get from probes script

for eg convert illumina probes to ensembl IDs

rsIDs to chr positions

##############
# illumina location file is incomplete, ~12800 QC'd BESTD probes match out of ~16700
# so probably many untested eQTLs
# Use biomaRt See: https://support.bioconductor.org/p/75923/
# Illumina download page is:
# http://support.illumina.com/array/array_kits/humanht-12_v4_expression_beadchip_kit/downloads.html
# Missing about one third of genomic locations from differentially expressed genes in BESTD
# Also see another tool:
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139516
# Reasons? Best to leave excluded if low quality but not filtered by illuminaHumanv4.db? Probe annotations are good or
# perfect according to illuminaHumanv4.db though.
# Can also get from illuminaHumanv4.db but erroring
##############

##############
# Load libraries:
library(biomaRt)
library(illuminaHumanv4.db)
library(plyr)
##############

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Illumina probes with chr locations:
illumina_file <- as.character(args[1])
# illumina_file <- 'illumina_probes_genomic_locations_annot.txt'

hits_file <- as.character(args[2])
# hits_file <- 'full_topTable_pairing_all_treated.txt'

illumina_original_file <- 'HumanHT-12_V4_0_R2_15002873_B.txt'
####################


####################
# Read data:
illumina_data <- read.csv(illumina_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
class(illumina_data)
str(illumina_data)
head(illumina_data)
dim(illumina_data)

hits_data <- read.csv(hits_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
class(hits_data)
str(hits_data)
head(hits_data)
dim(hits_data)

illumina_original <- read.csv(illumina_original_file, sep = '\t', skip = 8, header = TRUE, 
                              stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c('', ' ', 'NA'))
class(illumina_original)
str(illumina_original)
head(illumina_original)
tail(illumina_original)
dim(illumina_original)

illumina_data[1:10, 1]
hits_data[1:10, 2]

# Rename headers:
colnames(hits_data)[2] <- 'probe_IDs'
colnames(hits_data)
####################

####################
# Explore Illumina original file:
colnames(illumina_original)
# View(illumina_original)
length(which(is.na(illumina_original$Chromosome)))
length(which(is.na(illumina_original$ILMN_Gene)))
length(which(is.na(illumina_original$Probe_Coordinates)))
# TO DO: check what processing I did to the original file. For now look at other annotation databases.
####################


####################
# Basic checks:
length(unique(hits_data$probe_IDs))
which(as.character(hits_data$probe_IDs) %in% as.character(illumina_data$probe_IDs))
length(which(as.character(hits_data$probe_IDs) %in% as.character(illumina_data$probe_IDs)))
dim(illumina_data)
dim(hits_data)
####################

####################
# Sanity check on QC filtering script, test whether QC'd and normalised probes are of good quality:
# Get probe IDs:
rownames(hits_data) <- hits_data$probe_IDs
head(hits_data)
illumina_ids  <-  as.character(rownames(hits_data))
class(illumina_ids)
illumina_ids

# Get quality annotations:
probes_by_quality <- unlist(mget(illumina_ids, illuminaHumanv4PROBEQUALITY, ifnotfound = NA))
probes_by_quality <- as.data.frame(probes_by_quality)
head(probes_by_quality)
dim(probes_by_quality)
summary(probes_by_quality)
# All probes are of good or perfect quality.
# TO DO: contact Illumina, try to get re-annotations/locations from other pipeline.
####################

####################
# Run biomaRt to get missing locations:
# See: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf

# Check what's available:
listMarts()
# Set the source (database) and dataset to use:
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Check what information is available
filters <- listFilters(ensembl)
head(filters)
# View(listAttributes(ensembl))
# Search for keyword:
grep('illumina', listAttributes(ensembl), ignore.case = TRUE)

# Obtain values from the database and dataset selected:
# illumina_ids (obtained from QC'd probes in hits_file) and mart were already set:
biomart_attributes <- getBM(attributes=c('illumina_humanht_12_v4', 'hgnc_symbol', 'ensembl_gene_id', 'entrezgene',
                                         'chromosome_name', 'start_position', 'end_position'), 
                            filters = 'illumina_humanht_12_v4', values = illumina_ids, mart = ensembl)

head(biomart_attributes)
tail(biomart_attributes)
dim(biomart_attributes)
length(illumina_ids)
dim(hits_data)
dim(illumina_data)
####################

####################
# Write results to disk:
write.table(biomart_attributes, 'biomart_QCd_probes_genomic_locations_annot.txt', row.names = FALSE, quote = FALSE, sep = '\t', col.names = TRUE)

# Generate file for MatrixeQTL:
system('cat biomart_QCd_probes_genomic_locations_annot.txt | cut -f1,5,6,7 | awk -v OFS='\t' '$2="chr"$2' - | sed '1d' - > biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt')
####################

