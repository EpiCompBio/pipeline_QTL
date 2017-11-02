###############
# Airwave phenotype data
# 29 Feb 2016
# Airwave chronic inflammation proposal
# Objective 1
# What is chronic inflammation?
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
#   Check effect of season, see:
# http://www.nature.com/ncomms/2015/150512/ncomms8000/pdf/ncomms8000.pdf
###############


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/Users/antoniob/Desktop/Airwave_inflammation/results_1.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_Airwave_obj1_exploratory",".txt", sep=""), open='a')
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
R_session_saved_image <- paste('R_session_saved_image_Airwave_exploratory.RData')

#############################


#############################
## Update packages if necessary and load them:
# source("https://bioconductor.org/biocLite.R")
# biocLite()
#install.packages('foreign')

library(ggplot2)
library(Hmisc)
library(splines)
library(plyr)
source('moveme.R')
#############################

#############################
## Load data sets:
## Medications:
trt_025_data <- read.csv('trt_025_290216.csv', sep = ',', header = TRUE, 
                         na.string = c(-Inf, 'NULL', NULL, '.', '_3_', '_4_', '', ' ', 'NA', 'NaN'))


class(trt_025_data)
dim(trt_025_data)
head(trt_025_data)
tail(trt_025_data)
names(trt_025_data)
str(trt_025_data)

## Demographics and region data:
# TO DO: check north v. south for InfBiom
part_data <- read.csv('part_025_290216.csv', sep = ',', header = TRUE, 
                      na.string = c(-Inf, 'NULL', NULL, '.', '_3_', '_4_', '', ' ', 'NA', 'NaN'))

class(part_data)
dim(part_data)
names(part_data)

## Medical history, biochemical, behaviours, etc:
screen_data <- read.csv('screen_025_290216.csv', sep = ',', header = TRUE, stringsAsFactors=FALSE,
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

## Auto-enrollment data:
# TO DO: Change questionnaire codes to labels...
auto_data <- read.csv('Auto_025_290216.csv', sep = ',', header = TRUE, stringsAsFactors=FALSE, fileEncoding="latin1",
                      na.string = c('NULL', '.', '_1_', '_2_', '_3_', '_4_', '_5_', '_6_',
                                      '', ' ', 'NA', 'NaN', '<0',
                                      '-1', '-2', '-3', '-4' , '-5', '-6'))
class(auto_data)
dim(auto_data)
names(auto_data)
#View(auto_data)


## ECG data:
# TO DO: get interpretation codes
ecg_data <- read.csv('Ecg_025_290216.csv', sep = ',', header = TRUE, stringsAsFactors=FALSE,
                      na.string = c('NULL', '.', '_1_', '_2_', '_3_', '_4_', '_5_', '_6_',
                                    '', ' ', 'NA', 'NaN', '<0',
                                    '-1', '-2', '-3', '-4' , '-5', '-6'))
class(ecg_data)
dim(ecg_data)
names(ecg_data)
#View(ecg_data)
summary(ecg_data$interpretation_system)
count(ecg_data$interpretation_system)

## ECG detail data:
# TO DO: difference with ecg data?
ecg_detail_data <- read.csv('ecg_detail_025_290216.csv', sep = ',', header = TRUE, stringsAsFactors=FALSE,
                     na.string = c('NULL', '.', '_1_', '_2_', '_3_', '_4_', '_5_', '_6_',
                                   '', ' ', 'NA', 'NaN', '<0',
                                   '-1', '-2', '-3', '-4' , '-5', '-6'))
class(ecg_detail_data)
dim(ecg_detail_data)
names(ecg_detail_data)
#View(ecg_detail_data)
summary(ecg_detail_data$interpretation_system)
count(ecg_detail_data$interpretation_system)

## Cognitive data?:
# TO DO: Get question labels...
precogn_data <- read.csv('pre_cogn_025_290216.csv', sep = ',', header = TRUE, stringsAsFactors=FALSE,
                            na.string = c('NULL', '.', '_1_', '_2_', '_3_', '_4_', '_5_', '_6_',
                                          '', ' ', 'NA', 'NaN', '<0',
                                          '-1', '-2', '-3', '-4' , '-5', '-6'))
class(precogn_data)
dim(precogn_data)
names(precogn_data)
#View(precogn_data)
summary(ecg_detail_data$interpretation_system)
##################


##################
# Exploratory analyses of screen data:
summary(screen_data$FIBRINOGEN)
screen_data$FIBRINOGEN
summary(screen_data$C_REACTIVE_PROTEIN)
summary(screen_data$BILIRUBIN)
dim(screen_data)

qplot(screen_data$FIBRINOGEN); ggsave('FIBRINOGEN.pdf')
qplot(screen_data$C_REACTIVE_PROTEIN); ggsave('C_REACTIVE_PROTEIN.pdf')
qplot(screen_data$age); ggsave('age.pdf')
qplot(screen_data$HBA1C_PERCENT); ggsave('HBA1C_PERCENT.pdf')
qplot(screen_data$IS_PREGNANT)
qplot(screen_data$IS_SMOKER); ggsave('IS_SMOKER.pdf')
qplot(screen_data$BODY_MASS_INDEX); ggsave('histogram_bmi.pdf')
qplot(screen_data$WAIST_HIP)
qplot(screen_data$TAKES_INTENSE_EXERCISE); ggsave('TAKES_INTENSE_EXERCISE.pdf')
qplot(screen_data$BILIRUBIN); ggsave('BILIRUBIN.pdf')
qplot(screen_data$GLUCOSE); ggsave('GLUCOSE.pdf')
qplot(screen_data$NEUTROPHILS_COUNT); ggsave('NEUTROPHILS_COUNT.pdf')
qplot(screen_data$WHITE_BLOOD_CELL_COUNT); ggsave('WHITE_BLOOD_CELL_COUNT.pdf')
qplot(screen_data$FAT_PERCENTAGE); ggsave('FAT_PERCENTAGE.pdf')
summary(screen_data$BODY_TYPE)

ggplot(screen_data, aes(FIBRINOGEN, age)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_FIBRINOGEN_age.pdf')

ggplot(screen_data, aes(FIBRINOGEN, WHITE_BLOOD_CELL_COUNT)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_FIBRINOGEN_WHITE_BLOOD_CELL_COUNT.pdf')

ggplot(screen_data, aes(FIBRINOGEN, BODY_MASS_INDEX)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_FIBRINOGEN_BODY_MASS_INDEX.pdf')

ggplot(screen_data, aes(y=FIBRINOGEN, as.factor(IS_SMOKER))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggsave('boxplot_FIBRINOGEN_IS_SMOKER.pdf')

ggplot(screen_data, aes(y=C_REACTIVE_PROTEIN, as.factor(IS_SMOKER))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggsave('boxplot_C_REACTIVE_PROTEIN_IS_SMOKER.pdf')

ggplot(screen_data, aes(y=NEUTROPHILS_COUNT, as.factor(IS_SMOKER))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggsave('boxplot_NEUTROPHILS_COUNT_IS_SMOKER.pdf')

ggplot(screen_data, aes(y=LYMPHOCYTES_COUNT, as.factor(IS_SMOKER))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

ggplot(screen_data, aes(y=BILIRUBIN, as.factor(IS_SMOKER))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggsave('boxplot_BILIRUBIN_IS_SMOKER.pdf')

ggplot(screen_data, aes(BILIRUBIN, age)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_BILIRUBIN_age.pdf')

ggplot(screen_data, aes(BILIRUBIN, FIBRINOGEN)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_BILIRUBIN_FIBRINOGEN.pdf')

ggplot(screen_data, aes(BILIRUBIN, C_REACTIVE_PROTEIN)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_BILIRUBIN_CRP.pdf')

# Hardly any match between CRP and bilirubin measurements:
summary(screen_data$C_REACTIVE_PROTEIN)
length(which(!is.na(screen_data$C_REACTIVE_PROTEIN) & !is.na(screen_data$BILIRUBIN)))
length(which(is.na(screen_data$C_REACTIVE_PROTEIN) & is.na(screen_data$BILIRUBIN)))
screen_data[, c('BILIRUBIN', 'C_REACTIVE_PROTEIN')]

ggplot(screen_data, aes(BILIRUBIN, C_REACTIVE_PROTEIN)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_BILIRUBIN_BMI.pdf')

ggplot(screen_data, aes(BILIRUBIN, BODY_MASS_INDEX)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_BILIRUBIN_BMI.pdf')

ggplot(screen_data, aes(y=BILIRUBIN, as.factor(TAKES_INTENSE_EXERCISE))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggsave('scatter_BILIRUBIN_EXERCISE.pdf')

ggplot(screen_data, aes(y=BILIRUBIN, as.factor(BODY_TYPE))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggsave('scatter_BILIRUBIN_BODY_TYPE.pdf')

summary(screen_data$GLUCOSE)
ggplot(screen_data, aes(BILIRUBIN, GLUCOSE)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_BILIRUBIN_GLUCOSE.pdf')

summary(screen_data$TOTAL_CHOLESTEROL)
ggplot(screen_data, aes(BILIRUBIN, TOTAL_CHOLESTEROL)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_BILIRUBIN_TOTAL_CHOLESTEROL.pdf')

summary(screen_data$C_PEPTIDE)
ggplot(screen_data, aes(BILIRUBIN, C_PEPTIDE)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_BILIRUBIN_C_PEPTIDE.pdf')

summary(screen_data$HBA1C_PERCENT)
ggplot(screen_data, aes(BILIRUBIN, HBA1C_PERCENT)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggsave('scatter_BILIRUBIN_HBA1C_PERCENT.pdf')


ggplot(screen_data, aes(C_REACTIVE_PROTEIN, age)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(screen_data, aes(FIBRINOGEN, C_REACTIVE_PROTEIN)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(screen_data, aes(y=FIBRINOGEN, x=1)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

ggplot(screen_data, aes(y=GLUCOSE, x=1)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

ggplot(screen_data, aes(y=WHITE_BLOOD_CELL_COUNT, x=1)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

ggplot(screen_data, aes(y=WHITE_BLOOD_CELL_COUNT, x=1)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

ggplot(screen_data, aes(y=NEUTROPHILS_COUNT, x=1)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

ggplot(screen_data, aes(y=LYMPHOCYTES_COUNT, x=1)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

ggplot(screen_data, aes(y=MONOCYTES_COUNT, x=1)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

ggplot(screen_data, aes(y=BASOPHILS_COUNT, x=1)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

ggplot(screen_data, aes(y=EOSINOPHILS_COUNT, x=1)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

summary(screen_data$large_unstain_cells_count)
ggplot(screen_data, aes(y=large_unstain_cells_count, x=1)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

#############################

#############################
## Question
# Subset data for measures of interest:
#   Biochemical measurements: fibrinogen, CRP, neutrophil counts, cell count ratios
#   Nested case v. control definitions on NCD risk factors: obese v. lean, athletic v. inactive, smokers v. non-smokers,
#     unhealthy diet v. healthy diet, alcohol v. no alcohol intake

# TO DO: process so that all v. all can be done for correlations:
rcorr(as.matrix(screen_data))
length(which(complete.cases(screen_data)))

# TO DO: Check other biochem traits and relationship
# TO DO: Subset for current disease, medications, etc
names(screen_data)

# Physical activity:
summary(screen_data$TAKES_INTENSE_EXERCISE)
summary(screen_data$BODY_TYPE)

chi_exercise <- chisq.test(screen_data$TAKES_INTENSE_EXERCISE, screen_data$BODY_TYPE)
str(chi_exercise)
chi_exercise$observed
chi_exercise$p.value

phys_active <- subset(screen_data, subset = (screen_data$TAKES_INTENSE_EXERCISE == 'YES'))
phys_others <- subset(screen_data, subset = (screen_data$TAKES_INTENSE_EXERCISE == 'NO'))

summary(phys_active$FIBRINOGEN)
summary(phys_others$FIBRINOGEN)
boxplot(phys_active$FIBRINOGEN, phys_others$FIBRINOGEN)

qqnorm(phys_active$FIBRINOGEN)
qqnorm(phys_others$FIBRINOGEN)
t.test(phys_active$FIBRINOGEN, phys_others$FIBRINOGEN)
t.test(phys_active$age, phys_others$age)
wilcox.test(phys_active$FIBRINOGEN, phys_others$FIBRINOGEN)
wilcox.test(phys_active$age, phys_others$age)

boxplot(phys_active$C_REACTIVE_PROTEIN, phys_others$C_REACTIVE_PROTEIN)
qqnorm(phys_active$C_REACTIVE_PROTEIN)
qqnorm(phys_others$C_REACTIVE_PROTEIN)
summary(phys_others$C_REACTIVE_PROTEIN)
summary(phys_active$C_REACTIVE_PROTEIN)
wilcox.test(phys_active$C_REACTIVE_PROTEIN, phys_others$C_REACTIVE_PROTEIN)

boxplot(phys_active$NEUTROPHILS_COUNT, phys_others$NEUTROPHILS_COUNT)
qqnorm(phys_active$NEUTROPHILS_COUNT)
qqnorm(phys_others$NEUTROPHILS_COUNT)
summary(phys_others$NEUTROPHILS_COUNT)
summary(phys_active$NEUTROPHILS_COUNT)
wilcox.test(phys_active$NEUTROPHILS_COUNT, phys_others$NEUTROPHILS_COUNT)


boxplot(phys_active$BILIRUBIN, phys_others$BILIRUBIN)
qqnorm(phys_active$BILIRUBIN)
qqnorm(phys_others$BILIRUBIN)
summary(phys_others$BILIRUBIN)
summary(phys_active$BILIRUBIN)
wilcox.test(phys_active$BILIRUBIN, phys_others$BILIRUBIN)
t.test(phys_active$BILIRUBIN, phys_others$BILIRUBIN)


## Alcohol
summary(screen_data$BEER_DRUNK)
summary(screen_data$BEER_TYPE)
summary(screen_data$OTHER_ALCOHOL_UNITS)

qqnorm(screen_data$OTHER_ALCOHOL_UNITS)
qqnorm(screen_data$BEER_DRUNK)
rcorr(screen_data$BEER_DRUNK, screen_data$OTHER_ALCOHOL_UNITS, type = 'spearman')

ggplot(screen_data, aes(BEER_DRUNK, OTHER_ALCOHOL_UNITS)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(screen_data, aes(BEER_DRUNK, CIGARETTES_PER_DAY)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(screen_data, aes(OTHER_ALCOHOL_UNITS, CIGARETTES_PER_DAY)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(screen_data, aes(FIBRINOGEN, BEER_DRUNK)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(screen_data, aes(FIBRINOGEN, OTHER_ALCOHOL_UNITS)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(screen_data, aes(C_REACTIVE_PROTEIN, BEER_DRUNK)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(screen_data, aes(C_REACTIVE_PROTEIN, OTHER_ALCOHOL_UNITS)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(screen_data, aes(NEUTROPHILS_COUNT, BEER_DRUNK)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(screen_data, aes(NEUTROPHILS_COUNT, OTHER_ALCOHOL_UNITS)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

## Smoking:
summary(screen_data$CIGARETTES_PER_DAY)
summary(screen_data$IS_SMOKER)
qqnorm(screen_data$CIGARETTES_PER_DAY)

smokes_yes <- subset(screen_data, subset = (screen_data$IS_SMOKER == 'YES'))
smokes_no <- subset(screen_data, subset = (screen_data$IS_SMOKER == 'NO'))
dim(smokes_yes)
dim(smokes_no)
summary(smokes_yes$IS_SMOKER)
summary(smokes_no$IS_SMOKER)

ggplot(screen_data, aes(CIGARETTES_PER_DAY, FIBRINOGEN)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)
ggplot(screen_data, aes(IS_SMOKER, y=FIBRINOGEN)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
wilcox.test(smokes_yes$FIBRINOGEN, smokes_no$FIBRINOGEN)
t.test(smokes_yes$FIBRINOGEN, smokes_no$FIBRINOGEN)

ggplot(screen_data, aes(CIGARETTES_PER_DAY, C_REACTIVE_PROTEIN)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)
ggplot(screen_data, aes(IS_SMOKER, y=C_REACTIVE_PROTEIN)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
summary(smokes_no$C_REACTIVE_PROTEIN)
summary(smokes_yes$C_REACTIVE_PROTEIN)
summary(screen_data$CIGARETTES_PER_DAY)
count(screen_data$CIGARETTES_PER_DAY)
smoke_and_CRP <- screen_data[, c('C_REACTIVE_PROTEIN', 'CIGARETTES_PER_DAY')]
head(smoke_and_CRP)
wilcox.test(smokes_yes$C_REACTIVE_PROTEIN, smokes_no$C_REACTIVE_PROTEIN)
t.test(smokes_yes$C_REACTIVE_PROTEIN, smokes_no$C_REACTIVE_PROTEIN)

ggplot(screen_data, aes(CIGARETTES_PER_DAY, NEUTROPHILS_COUNT)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)
ggplot(screen_data, aes(IS_SMOKER, y=NEUTROPHILS_COUNT)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
wilcox.test(smokes_yes$NEUTROPHILS_COUNT, smokes_no$NEUTROPHILS_COUNT)
t.test(smokes_yes$NEUTROPHILS_COUNT, smokes_no$NEUTROPHILS_COUNT)

## Age and main biochemical markers of inflammation:
qqnorm(screen_data$age)

ggplot(screen_data, aes(age, FIBRINOGEN)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)
rcorr(screen_data$age, screen_data$FIBRINOGEN, type = 'spearman')
rcorr(screen_data$age, screen_data$FIBRINOGEN)

ggplot(screen_data, aes(age, C_REACTIVE_PROTEIN)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)
rcorr(screen_data$age, screen_data$C_REACTIVE_PROTEIN, type = 'spearman')
rcorr(screen_data$age, screen_data$C_REACTIVE_PROTEIN)

ggplot(screen_data, aes(age, NEUTROPHILS_COUNT)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)
rcorr(screen_data$age, screen_data$NEUTROPHILS_COUNT, type = 'spearman')
rcorr(screen_data$age, screen_data$NEUTROPHILS_COUNT)

## BMI
# TO DO: check fat percentage, impedance, WC, body water, etc.

qqnorm(screen_data$BODY_MASS_INDEX)
ggplot(screen_data, aes(BODY_MASS_INDEX, FIBRINOGEN)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)
rcorr(screen_data$BODY_MASS_INDEX, screen_data$FIBRINOGEN, type = 'spearman')
rcorr(screen_data$BODY_MASS_INDEX, screen_data$FIBRINOGEN)

ggplot(screen_data, aes(BODY_MASS_INDEX, C_REACTIVE_PROTEIN)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)
rcorr(screen_data$BODY_MASS_INDEX, screen_data$C_REACTIVE_PROTEIN, type = 'spearman')
rcorr(screen_data$BODY_MASS_INDEX, screen_data$C_REACTIVE_PROTEIN)

ggplot(screen_data, aes(BODY_MASS_INDEX, NEUTROPHILS_COUNT)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)
rcorr(screen_data$BODY_MASS_INDEX, screen_data$NEUTROPHILS_COUNT, type = 'spearman')
rcorr(screen_data$BODY_MASS_INDEX, screen_data$NEUTROPHILS_COUNT)

# Subset BMI:
screen_data$BMI_binned <- ifelse(screen_data$BODY_MASS_INDEX < 20, 'BMI_less_20',
                                 ifelse(screen_data$BODY_MASS_INDEX >= 20 & screen_data$BODY_MASS_INDEX < 25, 'BMI_20_25',
                                        ifelse(screen_data$BODY_MASS_INDEX >= 25 & screen_data$BODY_MASS_INDEX < 30, 'BMI_25_30',
                                               ifelse(screen_data$BODY_MASS_INDEX >= 30, 'BMI_over_30',
                                                      screen_data$BODY_MASS_INDEX
                                               ))))
screen_data$BMI_binned <- factor(screen_data$BMI_binned, levels=c('BMI_less_20', 'BMI_20_25', 'BMI_25_30', 'BMI_over_30'))
count(screen_data$BMI_binned)
screen_data$BMI_binned

ggplot(screen_data, aes(BMI_binned, y=BODY_MASS_INDEX)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
dev_to_pdf()
ggsave('BMI_binned_BODY_MASS_INDEX.pdf')

ggplot(screen_data, aes(BMI_binned, y=FIBRINOGEN)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggsave('BMI_binned_FIBRINOGEN.pdf')

ggplot(screen_data, aes(BMI_binned, y=NEUTROPHILS_COUNT)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggsave('BMI_binned_NEUTROPHILS_COUNT.pdf')

ggplot(screen_data, aes(BMI_binned, y=C_REACTIVE_PROTEIN)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggsave('BMI_binned_C_REACTIVE_PROTEIN.pdf')

ggplot(screen_data, aes(BMI_binned, y=age)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggsave('BMI_binned_age.pdf')



## Diet
names(screen_data)

# TO DO/CHECK: check very high outliers for several variables of interest;
# age doesn't correlate strongly with InfBioms. Beer is protective : )
# Outcome data is needed...

# Background
# Test
# Result
# Plot
# Interpretation
# Next step
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
# Write files to disk:
write.table(screen_data, 'screen_data_processed.tsv', row.names = FALSE, quote = FALSE, sep = '\t', 
            na = 'NA', col.names = TRUE)

# Write file for plink:
# Plink ID:
screen_data$BARCODE_plink <- screen_data$BARCODE
screen_data[1:5, 1:5]
screen_data[1:5, c('BARCODE_plink')]

# c1.rle <- rle(screen_data$BARCODE)
# c1.rle
# screen_data$BARCODE_plink <- paste0(rep(c1.rle$values, times = c1.rle$lengths), "_", 
#                                     unlist(lapply(c1.rle$values, seq_len)))

write.table(screen_data, 'screen_data_processed_plink.tsv', row.names = FALSE, quote = FALSE, sep = '\t', 
            na = '-9', col.names = TRUE)

# Write out minimal file for plink with age, gender, arbitrary case v. control
class(screen_data)
minimal_covariates_airwave <- screen_data[, c('BARCODE_plink', 'substudy_part_id', 'BARCODE',
                                              'BODY_MASS_INDEX', 'TAKES_INTENSE_EXERCISE', 
                                              'WHITE_BLOOD_CELL_COUNT', 'C_REACTIVE_PROTEIN',
                                              'age')]
head(minimal_covariates_airwave)
# Add gender:
dim(part_data)
dim(screen_data)
names(part_data)

minimal_covariates_airwave <- merge(minimal_covariates_airwave, part_data[, c('substudy_part_id', 'sex')])
head(minimal_covariates_airwave)

# Add cases (=2) and controls (=1) for plink tests, these are random:
# https://www.cog-genomics.org/plink2/formats#fam
# Sex code ('1' = male, '2' = female, '0' = unknown)
# Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
case_labels <- sample(c(1,2), size = nrow(screen_data), replace = TRUE)
minimal_covariates_airwave <- cbind(minimal_covariates_airwave, 'phenotype' = case_labels)

# Recode gender:
minimal_covariates_airwave$sex_cat <- ifelse(minimal_covariates_airwave$sex == 'MALE', 1, 
                                             ifelse(minimal_covariates_airwave$sex == 'FEMALE', 2, 
                                                    0))
# Rename and move columns  FID and IID for plink covariates file:
# https://www.cog-genomics.org/plink2/input#pheno
head(minimal_covariates_airwave)
names(minimal_covariates_airwave)[2] <- 'FID'
names(minimal_covariates_airwave)[3] <- 'IID'
names(minimal_covariates_airwave)

minimal_covariates_airwave <- minimal_covariates_airwave[moveme(names(minimal_covariates_airwave), "FID first")]
minimal_covariates_airwave <- minimal_covariates_airwave[moveme(names(minimal_covariates_airwave), "IID after FID")]
head(minimal_covariates_airwave)

tail(minimal_covariates_airwave)
dim(minimal_covariates_airwave)

summary(minimal_covariates_airwave$sex)
summary(minimal_covariates_airwave$phenotype)

# Write to file:
write.table(minimal_covariates_airwave, 'minimal_covariates_airwave.tsv', row.names = FALSE, quote = FALSE, sep = '\t', 
            na = '-9', col.names = TRUE)
#############################

#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for higher level analyses.
#############################