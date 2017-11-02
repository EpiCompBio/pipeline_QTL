#!/bin/bash

set -e

qsub -v GENO="cut_22_air_clean_SNPs.geno_matched.tsv" \
	-v EXPR="cut_air_probes.txt_matched.tsv" \
	-v threshold="1e-5" \
	-v cisThreshold="0.001" \
	qsub_mePipe_air.sh



