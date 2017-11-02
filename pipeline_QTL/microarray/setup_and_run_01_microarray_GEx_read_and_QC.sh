!#/bin/bash

# Script to set up files in order to run the gene expression analysis scripts

# Copy files with expression data (these already include a target file here for read.ilmn function):
ln -s /ifs/projects/sftp/backup/proj043/BESTD_all/galahad_data_backup_project_BESTD_25_Nov_2014.dir/BESTD/expression/raw/* .

# Copy file with membership information for gene expression samples:
#ln -s /ifs/projects/proj043/analysis.dir/GEx_BEST-D_targets_file.txt .
ln -s /ifs/projects/proj043/analysis.dir/BEST-D_lab_kit_IDs_for_array_QC.tsv .

# Leave only sample and control probe profiles files:
rm -f *docx md5sum* Plate_* *Table*

# Copy scripts needed to run analysis:
ln -s /ifs/devel/antoniob/projects/BEST-D/*microarray* .

# Rename files as the original ones have (two) spaces:
rename " " _ *
rename " " _ *
rename " " _ *

# Run analysis first step of QC:
#qsub qsub_R_01_microarray_GEx_read_and_QC.sh

