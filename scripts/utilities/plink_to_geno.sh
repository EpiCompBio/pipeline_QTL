#!/usr/bin/env bash

##########
# Run as eg:
# bash plink_to_geno.sh plink_bfile= out_oxford_recode= traw_file= out_geno_file=
# TO DO: move to ruffus although will need to run on Imperial as well.

# eQTL file processing script for genotype data
# Antonio J Berlanga-Taylor
# 05 Feb 2016

# Input: quality controlled genotyping data after processing with plink

# Outputs: genotype data in a format ready for MatrixEQTL

# Inputs are:
# QCd plink genotype file

# e.g.
# ln -s /ifs/projects/proj043/analysis.dir/genotypes.dir/P140343-Results_FinalReport_clean-individuals_and_SNPs* .
# Run the 'bash ln eQTL xxx' script first in the working directory (contains the commands for ln links for files and scripts).
# After this script run the 01_eQTL_xxx MatrixeQTL script for processing and matching covariates, gene expression and genotype files.
##########


##########
# Set bash script options:
# https://kvz.io/blog/2013/11/21/bash-best-practices/
set -o errexit
set -o pipefail
set -o nounset


# Set variables:
plink_bfile=$1
#P140343-Results_FinalReport_clean_SNPs_autosome
out_oxford_recode=$2
#P140343-Results_FinalReport_clean_SNPs_autosome.matrixQTL
traw_file=$3
#P140343-Results_FinalReport_clean_SNPs_autosome.A-transpose (which plink will convert to xxx.traw)
out_geno_file=$4
#P140343-Results_FinalReport_clean_SNPs_autosome.A-transpose.matrixQTL.geno
##########

##########
# Run commands:
# Convert plink binary (bed) to text (ped) with:
plink --bfile $plink_bfile --recode tab --out $plink_bfile


# Genotype files need to be processed for use in Matrix eQTL. First convert plink formats to oxford format:
# https://www.cog-genomics.org/plink2/data#recode
# https://www.cog-genomics.org/plink2/formats#gen
plink --bfile $plink_bfile --recode oxford --out $out_oxford_recode

# Then transpose and convert SNP coding using plink1.9 with --recode A-transpose:
# https://www.cog-genomics.org/plink2/data#recode
# https://www.cog-genomics.org/plink2/formats#raw
plink --bfile $plink_bfile --recode A-transpose --out $traw_file

# Finally cut (keep) columns required.
# Columns in the xxx.traw file are: CHR	SNP	(C)M	POS	COUNTED	ALT	sample_ID1
# So keep columns 2, 7 and onwards (which will give SNP and all samples):
cat $traw_file.traw | cut -f2,7- > $out_geno_file

# Remove intermediate files:
rm -f *ped *hh *nosex *traw *map *sample *gen

echo 'Done converting plink binary to MatrixEQTL format'
##########
