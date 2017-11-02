###############
# Airwave data extraction
# Oct 31, 2016
# Airwave phenotype data - missing gender variable
###############


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_airwave_extracted_variables",".txt", sep=""), open='a')
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
# load('R_session_saved_image_Airwave_exploratory.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_output_airwave_extracted_variables.RData')

#############################


#############################
## Update packages if necessary and load them:
# source("https://bioconductor.org/biocLite.R")
# biocLite()

# library(ggplot2)
library(data.table)
#############################

#############################
# Read data:

## Auto-enrollment data:
# TO DO: Change questionnaire codes to labels...
pheno_file = 'part_025_290216.csv'
pheno_data <- read.csv(pheno_file, sep = ',', header = TRUE, 
                      na.string = c(-Inf, 'NULL', NULL, '.', '_3_', '_4_', '', ' ', 'NA', 'NaN'))

class(pheno_data)
dim(pheno_data)
names(pheno_data)
# View(pheno_data)

# part data only has substudy ID, add BARCODE
pheno_file_2 = 'screen_025_290216.csv'
screen_data <- read.csv(pheno_file_2, sep = ',', header = TRUE, stringsAsFactors=FALSE,
                        na.string = c('NULL', '.', '_1_', '_2_', '_3_', '_4_', '_5_', '_6_',
                                      '', ' ', 'NA', 'NaN', '<0',
                                      '-1', '-2', '-3', '-4' , '-5', '-6'))
# Convert negative values to NAs as these aren't read properly with read.csv:
screen_data <- as.data.frame(lapply(screen_data, function(x){replace(x, x <0, NA)}))
names(screen_data)[1:10]

pheno_data <- merge(pheno_data, screen_data[, c("substudy_part_id", "BARCODE")])
dim(pheno_data)
head(pheno_data)

# Keep only gender and BARCODE:
pheno_data <- pheno_data[, c("BARCODE", "sex")]
dim(pheno_data)
head(pheno_data)

# Metabolomics data:
# Read in metabolomics data:
metabolomics_file <- '../metabolomics_data.dir/Airwave_CPMG_Plasma.txt'
metabolomics_data <- fread(metabolomics_file, sep = '\t', header = TRUE, 
                           stringsAsFactors = FALSE, select = c("Row")) #check if first cell is not empty
metabolomics_data <- as.data.frame(metabolomics_data)
class(metabolomics_data)
dim(metabolomics_data)
metabolomics_data[1:5, ]

#############
# Get subset of those only with metabolomics data:
# 'BARCODE' is the ID, var 'Row' in metabolomics data, rename:
names(metabolomics_data)[1] <- "BARCODE"
metabolomics_data[1:5, ]

samples_w_metab <- which(pheno_data[, "BARCODE"] %in% metabolomics_data[, "BARCODE"])
length(samples_w_metab)
head(samples_w_metab)
tail(samples_w_metab)

pheno_w_metabolomics <- pheno_data[samples_w_metab, ]
head(pheno_w_metabolomics)
tail(pheno_w_metabolomics)
grepl('47083', metabolomics_data)
grepl('55518', metabolomics_data)
grepl('55517', metabolomics_data)

dim(metabolomics_data)
dim(pheno_data)
dim(pheno_w_metabolomics)
#View(pheno_w_metabolomics)
#############################

#############################
# Write to file:
file_name <- sprintf('subset_%s', pheno_file)
file_name
setwd('~/Desktop/')
write.table(pheno_w_metabolomics, file_name, row.names = FALSE, quote = FALSE, sep = ',',
            na = 'NA', col.names = TRUE)
#############################

#############################
#The end:

# To save R workspace with all objects to use at a later time:
# save.image(file=R_session_saved_image, compress='gzip')
sessionInfo()

q()
#############################