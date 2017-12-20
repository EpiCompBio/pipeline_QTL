#!/bin/bash

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


# Generate only one set of unrelated cases:
plink --dummy 1000 10000 0.02 0.03 12 --out dummy_binary_pheno
plink --dummy 1000 10000 0.02 0.03 12 scalar-pheno --out dummy_scalar_pheno

