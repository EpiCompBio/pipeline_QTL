#############################################
# code blurb, this is andy's from BEST-D additional tests
# using here for some sva tests

library(sva)
library(SmartSVA)
sink( 'run_sva.Rout' )
load( 'R_session_saved_image_pheno_file_check.RData' )

time	<- as.numeric( membership_file_cleaned[,'visit_type'] == 'FinalVisit' )
vit_D	<- 2-as.numeric( membership_file_cleaned[,'arm'] ) 
pt_id	<- membership_file_cleaned[,'pt_id']

unique( time  )
unique( vit_D )

expr		<- t(normalised_filtered)
dim( expr )
expr[1:5, 1:10]
ngenes	<- nrow( expr )


expr_sv	<- scale(expr)
t(expr_sv[1:5, 1:5])
t(expr_sv[nrow(t(expr_sv)):-5, ncol(t(expr_sv)):-5])

X			<- cbind( time, time*vit_D )
mod0  <- model.matrix(~1, data=data.frame(expr_sv) )
mod   <- model.matrix(~1+X, data=data.frame(expr_sv) )
svs1  <- sva( dat=t(expr_sv), n.sv=10, mod=mod0, method="two-step" )$sv
svs2  <- smartsva.cpp( dat=t(expr_sv), n.sv=10, mod=mod	, mod0=mod0 )$sv

head(svs2)

#top SVs are correlated between sup/unsup implementations
sapply( 1:10, function(i) cor( svs1[,i], svs2[,i] ) )


pvals0	<- sapply( 1:ngenes, function(i)
  summary( lm( as.numeric( expr[,i] ) ~ time + time:vit_D ) )$coef[3,4]
)
pvals1	<- sapply( 1:ngenes, function(i)
  summary( lm( as.numeric( expr[,i] ) ~ time + time:vit_D + svs1 ) )$coef[3,4]
)
pvals2	<- sapply( 1:ngenes, function(i)
  summary( lm( as.numeric( expr[,i] ) ~ time + time:vit_D + svs2 ) )$coef[3,4]
)

pdf( 'sva_pvals.pdf', width=15, height=7 )
par( mfrow=c(1,3) )

hist( pvals0, xlim=0:1, breaks=30, freq=F, main='Uncorrected' )
hist( pvals1, xlim=0:1, breaks=30, freq=F, main='UnSup SVA Corrected (K=10)' )
hist( pvals2, xlim=0:1, breaks=30, freq=F, main='Sup SVA Corrected (K=10)' )

dev.off()
#############################################


#############################################
# code blubs for sva tests from docs from elsewhere

## Methylation M values (CpG by Sample, 27K by 1,000)
require(sva)
require(SmartSVA)
Y <- matrix(rnorm(100*270), 270, 100)
head(Y)
dim(Y)
df <- data.frame(pred = gl(2, 50))
head(df)
dim(df)
summary(df)

## Determine the number of SVs
Y.r <- t(resid(lm(t(Y) ~ pred, data = df)))
n.sv <- 50
mod_test <- model.matrix( ~ pred, df)
head(mod_test)
tail(mod_test)
sv.obj1 <- smartsva.cpp(Y, mod = mod_test, mod0 = NULL, B = 5, alpha = 1, VERBOSE = TRUE, n.sv = n.sv)

str(sv.obj1)
head(sv.obj1$sv)
#############################################



#############################################
# Imputation blurb for sva tests

library(mice)
library(miceadds)


#############################################
##########
# Impute missing values
# TO DO: move to separate script
# See:
# https://www.r-bloggers.com/imputing-missing-data-with-r-mice-package/
# https://www.jstatsoft.org/article/view/v045i03
# How to formally check data is missing at random or not?

# Check proportion of missing data, usually <5%:
prop_NA <- function(x) {sum(is.na(x)) / length(x) * 100}
# Individuals with more than 5% of missing variables:
apply(all_data, 1, prop_NA) # by rows
length(all_data[which(apply(all_data, 1, prop_NA) > 10), 'pt_id'])
# By columns:
apply(all_data[, covars_list], 2, prop_NA)
apply(all_data[, cytokines_0], 2, prop_NA)
apply(all_data[, cytokines_12], 2, prop_NA)

# See pattern using VIM and mice libraries
# Plot missing values:
# View(md.pattern(all_data))
png('missing_data_plots_vars_interest.png', width = 7.3, height = 5, units = 'in', res = 600)
aggr_plot <- aggr(all_data[, c(covars_list, cytokines_0, cytokines_12)],
                  only.miss = T, # Plot only missing variables
                  col = c('lightgrey', 'red'), # 1 colour for missing data, 2 observed, 3 imputed
                  numbers = T, sortVars = T,
                  labels = names(all_data[, c(covars_list, cytokines_0, cytokines_12)]),
                  cex.axis = 0.4,
                  gap = 2,
                  ylab = c('Proportion of missing data', 'Pattern'))
dev.off()
##########

###########
# Impute missing data
# Easy tutorial using library mice:
# http://web.maths.unsw.edu.au/~dwarton/missingDataLab.html

names(all_data)
all_data[1:5, 1:5]
all_data_interest <- all_data[, c(covars_list, cytokines_0, cytokines_12)]
# Keep samples IDs but exclude them from imputation (non missing and would be used to estimate
# imputation if left as column):
rownames(all_data_interest) <- all_data[, 'pt_id']
identical(all_data[, 'pt_id'], rownames(all_data_interest))
head(all_data_interest)

# Run imputation
# Roughly one imputation per percent of incomplete data (White et al.,2011),
# but the more the better, 100 can easily be run on small datasets on a laptop
# Roughly 20-30 iterations should be enough, use plot() to check convergence:
# http://stats.stackexchange.com/questions/219013/how-do-the-number-of-imputations-the-maximum-iterations-affect-accuracy-in-mul
imp_all_data <- mice(all_data_interest,
                     m = 50, # Number of imputed datasets, 5 is default
                     maxit = 50, 
                     # meth = 'pmm', # predictive mean matching, leave empty for 
                     # auto selection depending on variable type
                     diagnostics = T,
                     seed = 500)
summary(imp_all_data)
# Plot convergence of imputed data, only plots the last 3 variables:
png('convergence_plot_imputations.png', width = 7.3, height = 5, units = 'in', res = 600)
plot(imp_all_data)
dev.off()

# Check the imputed data, e.g.:
imp_all_data$method
imp_all_data$imp$vitd12 # Each column is the number of datasets ran
png('missing_data_scatterplots_VD12_cyto12.png', width = 7.3, height = 5, units = 'in', res = 600)
xyplot(imp_all_data, 
       vitd12 ~ Ln_IFNgamma12 + 
         Ln_IL10_12 +
         Ln_IL6_12 +
         Ln_IL8_12 +    
         Ln_TNFalpha12,
       # pch = 1, cex = 1, strip = T, 
       ylab = '25OHD levels at 12 months',
       xlab = '',
       strip = strip.custom(factor.levels = labels),
       type = c('p')) # Magenta are imputed, blue observed
dev.off()
##########

###########
# Further exploratory plots:
png('missing_data_densityplots.png', width = 7.3, height = 5, units = 'in', res = 600)
densityplot(imp_all_data) # Plots all numerical variables with 2 or more missing values 
dev.off()

png('missing_data_bwplots.png', width = 7.3, height = 5, units = 'in', res = 600)
bwplot(imp_all_data)
dev.off()

# Stripplots look better:
png('missing_data_stripplots.png', width = 7.3, height = 5, units = 'in', res = 600)
stripplot(imp_all_data,
          subset = (.imp == 1 | .imp == 2 | .imp == 3 | .imp == 4 | .imp == 5 |
                      .imp == 6 | .imp == 7 | .imp == 8 | .imp == 9 | .imp == 10),
          # col = mdc(1:2), #col = mdc(1:2), pch=20, cex=1.5,
          pch = 1, cex = 0.7,
          strip = strip.custom(par.strip.text = list(cex = 0.7)))
# Magenta are imputed
dev.off()
##########

###########
# Sanity check observed and imputed data
# Variables used as predictors for imputation of each incomplete variable:
imp_all_data$pred
# Implausible results for specific variables:
which(imp_all_data$imp$vitd12 <= 1 | imp_all_data$imp$vitd12 >= 250)
which(imp_all_data$imp$calendar_age_ra <= 60 | imp_all_data$imp$calendar_age_ra >= 100)

# Create dataset with both observed and imputed data:
# imp_all_data_completed <- complete(imp_all_data, action = 'repeat', include = T)
# 'repeated' includes original data (colnames "xxx.0") and all imputations ("xxx.1, etc"). 
imp_all_data_completed <- complete(imp_all_data, include = T)
head(imp_all_data_completed)
dim(imp_all_data_completed)
length(complete.cases(imp_all_data_completed) == TRUE)
# View(imp_all_data_completed)

complete.cases(all_data_interest)
complete.cases(imp_all_data$data)
complete.cases(imp_all_data_completed)
# If mice::complete(action = 'repeat') then this should be TRUE, else FALSE like here:
# action returns the first completed imputed dataset
identical(complete.cases(all_data_interest), complete.cases(imp_all_data_completed))

sapply(all_data_interest, function(x) sum(is.na(x)))
sapply(imp_all_data$data, function(x) sum(is.na(x)))
sapply(imp_all_data_completed, function(x) sum(is.na(x)))

# These should all be TRUE:
identical(all_data$vitd12, all_data_interest$vitd12)
identical(all_data_interest$vitd12, imp_all_data$data$vitd12)
# identical(imp_all_data$data$vitd12, imp_all_data_completed$vitd12.0) # 'xxx.0' is the original data

# Next step: use mice and miceadds libraries to run regression models (bottom of script)
# with imputed dataset and pooled values from multiple imputed datasets and observed data.
##########


###########
# TO DO:
# Consider library(Hmisc) with aregImpute() using additive regression, 
# bootstrapping, and predictive mean matching (pmm)
# It also adapts the method based on variable type automatically
# PMM (with mice or others) for numerical variables
# For categorical in mice use
# polyreg(Bayesian polytomous regression) for factor variables with >= 2 levels
# proportional odds model for ordered variables with >= 2 levels
##########
#############################################
