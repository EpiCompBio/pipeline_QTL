#!/bin/bash

Rscript="/Library/Frameworks/R.framework/Resources/bin/Rscript"
R_script="eQTL_plotting_airwave.R"

set -e


###########
# Plots at baseline:
eQTL_file="me_25_covariates.eQTL"
geno_file="cut_chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome.geno_matched.tsv"
gex_file="cut_Airwave_CPMG_Plasma.txt_matched.tsv"
PC_file='PCs_to_adjust_for_cut_22_air_clean_SNPs.geno_matched.tsv.txt'
PCs_to_correct='25'

SNP="22:18906839"
PROBE="Plasma_CPMG_NMR_4_127146"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="22:18911333"
PROBE="Plasma_CPMG_NMR_4_127146"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="22:18889966"
PROBE="Plasma_CPMG_NMR_4_127146"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="22:18889969"
PROBE="Plasma_CPMG_NMR_4_127146"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="22:18910479"
PROBE="Plasma_CPMG_NMR_4_127146"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="22:18890037"
PROBE="Plasma_CPMG_NMR_4_127146"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="22:18905964"
PROBE="Plasma_CPMG_NMR_4_127146"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="22:18887537"
PROBE="Plasma_CPMG_NMR_4_127483"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="22:18906839"
PROBE="Plasma_CPMG_NMR_4_127483"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="22:18887537"
PROBE="Plasma_CPMG_NMR_4_127146"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

