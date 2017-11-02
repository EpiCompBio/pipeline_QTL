#!/bin/bash

#Use current working directory:
#$ -cwd

#Run on n processors:
#$ -pe dedicated 2

# select all.q queue:
#$ -q all.q

#Memory option, each thread specified above will use the amount of memory specified here:
#$ -l mem_free=10G

#Save standard error and out to files:
#$ -e stderr_eQTl_mePipe.file
#$ -o stdout_eQTL_mePipe.file

#Run the job - program must have full path:
/ifs/apps/apps/R-3.1.3/bin/Rscript run_mePipe_air.R \
					$GENO \
					$EXPR \
					$threshold \
					$cisThreshold
				#       $PC_seq_to_test \
