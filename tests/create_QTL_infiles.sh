#!/usr/bin/env bash
  
# Set bash script options:
# https://kvz.io/blog/2013/11/21/bash-best-practices/
set -o errexit
set -o pipefail
set -o nounset


# Create softlinks and rename files to use as input tests for pipeline
# Run as:
#bash create_QTL_infiles.sh ~/Documents/github.dir/EpiCompBio/pipeline_QTL/tests/ 

dirname=$1
sim_dir=simulated_data
# Simulated files with pipeline starting at plink and pheno with rows as vars
# and cols as samples:

#####
echo 'Simulating pheno dataset: '
mkdir -p ${sim_dir}
cd ${sim_dir}

python ${dirname}/simulate_cont_var.py --createDF \
                                       --sample-size=1000 \
                                       --var-size=1500 \
                                       --sd=1.0 \
                                       --mean=0.0 \
                                       --lower-bound=-5.0 \
                                       --upper-bound=5.0
echo 'Done'
#####

#####
echo 'Simulating genotype dataset and running QC: '
bash ${dirname}/simulate_plink.sh dummy_binary
echo 'Done'
#####

#####
# Create links in results directory (one above):
echo 'Creating softlinks and changing names for pipeline input...'
cd ..

# Change filenames:
ln -fs ${sim_dir}/continuous_var_simulation.tsv simulated-NMR-blood.pheno
ln -fs ${sim_dir}/dummy_binary.QC.bed simulated-dummy_binary-QC.bed
ln -fs ${sim_dir}/dummy_binary.QC.bim simulated-dummy_binary-QC.bim
ln -fs ${sim_dir}/dummy_binary.QC.fam simulated-dummy_binary-QC.fam
ln -fs ${dirname}/SNP_exclusion_regions.txt .

# Files from MatrixQTL which don't need plink processing and ared ready for
# PCA:
ln -fs ${dirname}/Covariates_ext.txt mxqtl_ext-illumina_exome-all_chrs-NMR-blood.cov
ln -fs ${dirname}/SNP_ext.txt mxqtl_ext-illumina_exome-all_chrs.geno
ln -fs ${dirname}/GE_ext.txt mxqtl_ext-NMR-blood.pheno

ln -fs ${dirname}/Covariates.txt airwave-illumina_exome-all_chrs-NMR-blood.cov
ln -fs ${dirname}/SNP.txt airwave-illumina_exome-all_chrs.geno
ln -fs ${dirname}/GE.txt airwave-NMR-blood.pheno

ln -fs ${dirname}/Covariates_ext.txt other-illumina_exome-all_chrs-NMR-urine.cov
ln -fs ${dirname}/SNP_ext.txt other-illumina_exome-all_chrs.geno
ln -fs ${dirname}/GE_ext.txt other-NMR-urine.pheno
echo 'Done'
#####

#####
#echo 'Removing intermediate files...'
#rm -rf ${dirname}.bed \
#       ${dirname}.bim \
#       ${dirname}.fam

#echo 'Done'
#####
