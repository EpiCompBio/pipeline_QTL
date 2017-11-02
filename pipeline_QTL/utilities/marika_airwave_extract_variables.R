###############
# Airwave data extraction
# Oct 11, 2016
# Airwave phenotype data
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

## Medical history, biochemical, behaviours, etc:
pheno_file = 'screen_025_290216.csv'
screen_data <- read.csv(pheno_file, sep = ',', header = TRUE, stringsAsFactors=FALSE,
                        na.string = c('NULL', '.', '_1_', '_2_', '_3_', '_4_', '_5_', '_6_',
                                      '', ' ', 'NA', 'NaN', '<0',
                                      '-1', '-2', '-3', '-4' , '-5', '-6'))

# Convert negative values to NAs as these aren't read properly with read.csv:
screen_data <- as.data.frame(lapply(screen_data, function(x){replace(x, x <0, NA)}))

class(screen_data)
dim(screen_data)
names(screen_data)
#View(screen_data)
var_classes <- sapply(screen_data, class)
as.data.frame(var_classes)
head(screen_data)
summary(screen_data$DIAG_STROKE_AGE)
summary(screen_data$C_REACTIVE_PROTEIN)
length(which(!is.na(screen_data$DIAG_STROKE_AGE)))
which(!is.na(screen_data$DIAG_DIABETES_COMMENTS))
screen_data[1:5, 1:5]

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
# Variables needed for study:
vars_file <- '~/Desktop/Downloads_to_delete/screen_selected.txt'
vars_selected <- read.csv(vars_file, header = F, stringsAsFactors = F)
str(vars_selected)

# Convert to upper case:
# vars_selected_df <- as.data.frame(sapply(vars_selected, toupper))
# vars_selected_df <- as.data.frame(as.character(vars_selected_df[, 1]))
#View(vars_selected_df)
vars_selected
dim(vars_selected)

# Check how many:
# length(vars_selected)
length(vars_selected[, 1])
length(which(names(screen_data) %in% vars_selected[, 1]))

#############
# Variables to extract:
to_extract <- which(names(screen_data) %in% vars_selected[, 1])
vars_extracted <- screen_data[, to_extract]

# Check:
#View(vars_extracted)
dim(vars_extracted)
names(vars_extracted)
length(which(as.character(vars_selected[, 1]) %in% as.character(names(vars_extracted))))

# Missing:
which(!as.character(vars_selected[, 1]) %in% as.character(names(vars_extracted)))
missing_vars <- which(!as.character(vars_selected[, 1]) %in% as.character(names(vars_extracted)))
missing_vars <- sapply(as.character(vars_selected[as.integer(missing_vars), ]), tolower)
missing_vars

# which(sapply(missing_vars, grepl, names(screen_data), ignore.case = TRUE))
# Not in screen file.

# Get subset of those only with metabolomics data:
# 'BARCODE' is the ID, var 'Row' in metabolomics data, rename:
names(metabolomics_data)[1] <- "BARCODE"
metabolomics_data[1:5, ]

samples_w_metab <- which(vars_extracted[, "BARCODE"] %in% metabolomics_data[, "BARCODE"])
length(samples_w_metab)

pheno_w_metabolomics <- vars_extracted[samples_w_metab, ]

dim(metabolomics_data)
dim(screen_data)
dim(vars_selected)
dim(vars_extracted)
dim(pheno_w_metabolomics)
# View(pheno_w_metabolomics)

# TO DO/CHECK:
# Sanity:
sanity_numbers <- pheno_w_metabolomics[1:5, "BARCODE"]
sanity_numbers
screen_data[which((sanity_numbers %in% screen_data[, "BARCODE"]) == TRUE), "BARCODE"]
pheno_w_metabolomics[which((sanity_numbers %in% pheno_w_metabolomics[, "BARCODE"]) == TRUE), "BARCODE"]
length(which(pheno_w_metabolomics[, "BARCODE"] %in% screen_data[, "BARCODE"]))


#############################

#############################
# Write to file:
file_name <- sprintf('subset_%s', pheno_file)
file_name
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