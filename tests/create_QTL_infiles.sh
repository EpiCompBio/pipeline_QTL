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
ln -fs ${sim_dir}/continuous_var_simulation.tsv simulated1-NMR-blood.pheno
ln -fs ${sim_dir}/dummy_binary.QC.bed simulated1-dummy_binary-all_chrs.QC.bed
ln -fs ${sim_dir}/dummy_binary.QC.bim simulated1-dummy_binary-all_chrs.QC.bim
ln -fs ${sim_dir}/dummy_binary.QC.fam simulated1-dummy_binary-all_chrs.QC.fam

# Create an artificial second set just with different symlink names:
ln -fs ${sim_dir}/continuous_var_simulation.tsv simulated2-NMR-blood.pheno
ln -fs ${sim_dir}/dummy_binary.QC.bed simulated2-dummy_binary-all_chrs.QC.bed
ln -fs ${sim_dir}/dummy_binary.QC.bim simulated2-dummy_binary-all_chrs.QC.bim
ln -fs ${sim_dir}/dummy_binary.QC.fam simulated2-dummy_binary-all_chrs.QC.fam

cp -f ${dirname}/SNP_exclusion_regions.txt .

# Files from MatrixQTL which don't need plink processing and ared ready for
# PCA:
# Rename the covariates so PCA isn't run and this file gets picked up instead:
# TO DO: Fix this so an additional non molecular pheno non-geno can be passed (age, gender, etc):
# Also: Ruffus updates tasks if input files are newer than output files and re-writes even if symlinked?
#ln - sf ${dirname}/Covariates.txt airwave-illumina_exome-all_chrs.geno.mx_qtl.pcs.tsv
#cp -f ${dirname}/Covariates.txt airwave-illumina_exome-all_chrs.geno.mx_qtl.pcs.tsv
#cp -f ${dirname}/SNP.txt airwave-illumina_exome-all_chrs.geno
#cp -f ${dirname}/GE.txt airwave-NMR-blood.pheno


# Create a second set of files to test multiple sets in pipeline:
#cp -f ${dirname}/Covariates.txt other-illumina_exome-all_chrs.geno.mx_qtl.pcs.tsv
#cp -f ${dirname}/SNP.txt other-illumina_exome-all_chrs.geno
#cp -f ${dirname}/GE.txt other-NMR-urine.pheno
echo 'Done'
#####

#####
#echo 'Removing intermediate files...'
#rm -rf ${dirname}.bed \
#       ${dirname}.bim \
#       ${dirname}.fam

#echo 'Done'
#####
