# File renaming code blurb and convention
# saving for re-working later on
# this code doesn't do anything, saving as notes

# Naming convention for input files
# =================================
#   
#   Output files get named based on the input files. The script assumes "cohort" is the same for input files (but only takes it from the genotype file).
# 
# Please rename your files in the following way:
#   
#   Infile: cohort-platform-other_descriptor.suffix
# 
# Outfile: cohort-platform_infile1-descriptor1-platform_infile2-descriptor2.new_suffix
# 
# For example:
#   
#   genotype file: airwave-illumina_exome-all_chrs.geno
# 
#   phenotype file: airwave-NMR-blood.pheno
# 
# and depending on the input and arguments you might get:
#   
#   airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL
# airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.cis
# airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.trans
# airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.qqplot.svg
# airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.degrees_condition.txt
# airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.log
# 
# File names can get long so use abbreviations or short versions.
# 
# You can also override this and simply choose your outfile prefix.

##########

infile_1 <- 'cohort-platform1-descriptor1.suffix1'
infile_2 <- 'cohort-platform2-descriptor2.suffix2'
SNP_file <- infile_1
expression_file <- infile_2



# Set output file names:
if (is.null(args[['-O']])) {
  print('Output file name prefix not given. Using:')
  # Separate names componenets based on convention for input in this script:
  infile_1 <- strsplit(SNP_file, '-')
  cohort <- infile_1[[1]][1]
  platform1 <- infile_1[[1]][2]
  descriptor_1 <- strsplit(infile_1[[1]][3], '\\.')[[1]][1]
  # Do the same for the phenotype file, take cohort only from genotype file:
  infile_2 <- strsplit(expression_file, '-')
  platform2 <- infile_2[[1]][2]
  descriptor_2 <- strsplit(infile_2[[1]][3], '\\.')[[1]][1]
  output_file_name <- sprintf('%s-%s-%s-%s-%s.MxEQTL', cohort, platform1, descriptor_1, platform2, descriptor_2)
  print(output_file_name)
} else {
  output_file_name <- as.character(args[['-O']])
  # output_file_name <- 'testing'
  output_file_name <- sprintf('%s.MxEQTL', output_file_name)
  print(sprintf('Output file names will contain %s', output_file_name))
}
##########