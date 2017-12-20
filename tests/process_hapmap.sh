#!/usr/bin/env bash
  
# Set bash script options:
# https://kvz.io/blog/2013/11/21/bash-best-practices/
set -o errexit
set -o pipefail
set -o nounset


# See https://github.com/gabraham/flashpca/tree/master/HapMap3 :
##wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
##wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2

# There are ~1194 individuals, reduce to the first 100:
##head -n 100 hapmap3_r2_b36_fwd.consensus.qc.poly.ped


# Reduce the number of SNPs
plink_file=hapmap3_r2_b36_fwd.MEX.qc.poly

# Get a smaller set, pick the smallest file from the various 1000G pop:
wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.qc.poly.tar.bz2
bunzip2 hapmap3_r2_b36_fwd.qc.poly.tar.bz2
tar xvfz hapmap3_r2_b36_fwd.qc.poly.tar
mv hapmap3_pop/${plink_file}.* .
head ${plink_file}.ped | cut -f1-10

# Create a map file:
cat ${plink_file}.ped | cut -f1-6 > ${plink_file}.map

# Convert to binary format:
plink --file ${plink_file} --make-bed --out ${plink_file}

# Run a basic QC:
plink --bfile ${plink_file} \
      --maf 0.05 --geno 0.01 --mind 0.01 --hwe 5e-6 \
      --filter-founders --autosome \
      --make-bed --out ${plink_file}.QC

# Delete files no longer needed (~6 Gb):
rm -rf hapmap3_pop/
rm -rf hapmap3_r2_b36_fwd.qc.poly.tar

rm -rf ${plink_file}.ped \
       ${plink_file}.map \
       ${plink_file}.b* \
       ${plink_file}.fam \
       ${plink_file}.log \
       ${plink_file}.hh
