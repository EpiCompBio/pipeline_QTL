# Update sample IDs for Airwave when formatting for MatrixEQTL

setwd('/Users/antoniob/Desktop/Airwave_inflammation/results_3.dir')

# Create file with IDs:
system("cat chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome.fam | cut -d ' ' -f1 > sample_IDs.txt")

# Read IDs:
sample_IDs <- read.csv('sample_IDs.txt', header = FALSE, stringsAsFactors = FALSE, sep = ' ')
head(sample_IDs)
str(sample_IDs)

# Split plink double IDs:
sample_IDs_split <- strsplit(sample_IDs$V1, split = '[_]')
# Create data frame:
sample_IDs_split_df <- data.frame(matrix(unlist(sample_IDs_split), 
                                         nrow = length(sample_IDs_split), 
                                         byrow = T), 
                                  stringsAsFactors=FALSE)

dim(sample_IDs_split_df)
str(sample_IDs_split_df)

# Create plink file:
# See: 
# https://www.cog-genomics.org/plink2/data#update_indiv
new_sample_IDs <- cbind(sample_IDs, sample_IDs, sample_IDs_split_df)
dim(new_sample_IDs)
head(new_sample_IDs)
tail(new_sample_IDs)

# Write to file:
write.table(new_sample_IDs, 'new_sample_IDs.txt', row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = '\t')

# Check:
system('head new_sample_IDs.txt')
system('wc -l new_sample_IDs.txt')

# Update plink:
# plink2 --bfile chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome --update-ids new_sample_IDs.txt --make-bed --out chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome
# bash 00_eQTL_genotype_process.sh chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome.matrixQTL chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome.A-transpose chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome.A-transpose.matrixQTL.geno
q()
