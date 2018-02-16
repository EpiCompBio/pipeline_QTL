##########
library(qvalue)
library(data.table)
##########

##########
"The program takes a list of p-values and computes their estimated π0 - the proportion of 
features that are truly null - based on their distribution (the assumption used is that p-values of 
truly alternative cases tend to be close to zero, while p-values of null features will be uniformly 
distributed among [0,1])."
##########

##########
browseVignettes("qvalue")
setwd('~/Documents/quickstart_projects/BEST_D_molecular.p_q/results/reviewers')
# Read first file:
input_name_1 <- '2000+4000-baseline-1.eQTL_cis'
input_data_1 <- fread(input_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
input_data_1
# p-value column:
# input_data_1_pval <- as.data.frame(input_data_1[, c('SNP', 'gene', 'p-value')])
# input_data_1_pval

# Read second file:

# p-value column:

##########

##########
# Obtain q-values:
qobj_1 <- qvalue(p = input_data_1$`p-value`)
str(qobj_1)
summary(qobj_1)
# FDR:
qobj_1$qvalues[1:10]
input_data_1$FDR[1:10]
# Estimate of the proportion of null p-values:
qobj_1$pi0
##########

##########
"The quantity π1 = 1−π0 estimates the lower bound of the proportion of truly 
alternative features, i.e. the proportion of true positives (TP). 
Replication and sharing between two samples is reported as the proportion of TP (π1) estimated 
from the p-value distribution of independent eQTLs discovered in sample 1 in the second sample 
(exact SNP-probe combinations are tested)."

##########