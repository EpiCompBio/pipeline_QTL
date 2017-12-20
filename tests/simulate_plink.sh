#!/usr/bin/env bash
  
# Set bash script options:
# https://kvz.io/blog/2013/11/21/bash-best-practices/
set -o errexit
set -o pipefail
set -o nounset


# Intended to generate a dummy file with test genotypes. Easier for uploading
# to GitHub

# See also:
# https://raw.githubusercontent.com/chrchang/plink-ng/master/1.9/tests/test_setup.sh

# Runs on plink 1.9 as plink

#plink --dummy [sample ct] [SNP ct] {missing geno freq} {missing pheno freq} <acgt | 1234 | 12> <scalar-pheno>

# where (From plink --help):
# This generates a fake input dataset with the specified number of samples
#    and SNPs.  By default, the missing genotype and phenotype frequencies are
#    zero, and genotypes are As and Bs (change the latter with
#    'acgt'/'1234'/'12').  The 'scalar-pheno' modifier causes a normally
#    distributed scalar phenotype to be generated instead of a binary one.


# Variables:
plink_file=$1
#dummy_binary # dummy_scalar

# Generate only one set of unrelated cases:
plink --dummy 1000 10000 0 0 12 --out ${plink_file}
#plink --dummy 1000 10000 0.02 0.03 12 scalar-pheno --out ${plink_file}

# Run a basic QC:
plink --bfile ${plink_file} \
      --maf 0.05 --geno 0.01 --mind 0.01 --hwe 5e-6 \
      --filter-founders --autosome \
      --make-bed --out ${plink_file}.QC

