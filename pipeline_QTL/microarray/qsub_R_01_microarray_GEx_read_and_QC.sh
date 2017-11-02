#!/bin/bash

#Use current working directory:
#$ -cwd

#Run on n processors:
#$ -pe dedicated 3 

# select all.q queue:
#$ -q all.q

#Memory option, each thread specified above will use the amount of memory specified here:
#$ -l mem_free=25G

#Save standard error and out to files:
#$ -e stderr_QC.file
#$ -o stdout_QC.file

#Run the job - program must have full path:
/ifs/apps/apps/R-3.1.3/bin/Rscript 01_microarray_GEx_read_and_QC.R 
