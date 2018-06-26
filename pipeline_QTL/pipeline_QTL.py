'''
pipeline_QTL
===================

:Author: Antonio Berlanga-Taylor
:Release: |version|
:Date: |today|


Overview
========

This is a pipeline that wraps tools such as matrixEQTL.

|long_description|


Purpose
=======


Ruffus pipeline which runs quantitative trait analysis.

Originally thought for human QTls on gene expression but depending on the tool
can run on other quantitative molecular phenotypes.


Usage and options
=================

For command line help type:

    pipeline_QTL --help


Configuration
=============

This pipeline is built using a Ruffus/CGAT approach. You need to have Python,
Ruffus, CGAT core tools and any other specific dependencies needed for this
script.

A configuration file needs to be created running the config option.

Use this to extract any arbitrary parameters that could be changed in future
re-runs of the pipeline.


Input files
===========

Plink formatted binary files (bed, bim, fam) and a molecular phenotype file
(continuous variables) with samples (individuals) as columns and variables
(features, phenotypes) as rows.

Both files have to be named as explained below.

Phenotype and genotype data must have been quality controlled already.

Optional files:

- covariates
- error covariance matrix
- SNP position file
- probe position file

Annotation files are needed for cis vs trans analysis (SNPs positions and
probe positions and any associated annotations).

Files need to be in the same format as MatrixEQTL requires:

SNPs/genes/covariates in rows, individuals in columns with (dummy) headers and
row names (the first column and first row are skipped when read).

Missing values must be set as "NA".

SNP and probe position files, e.g.

   - snp146Common_MatrixEQTL_snp_pos.txt
   - biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt

must be tab separated.

See:

http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

A bash script to create test files can be found in:

    https://github.com/EpiCompBio/pipeline_QTL/tests

Run as eg:

    bash /path_to/create_QTL_infiles.sh /path_to_further_scripts/pipeline_QTL/tests

Naming convention for input files
=================================

Please rename your files in the following way (use soft links to rename for
example):

Infile: cohort-platform-other_descriptor.QC.suffix

Outfile: cohort-platform_infile1-descriptor1-platform_infile2-descriptor2.new_suffix

For example:

genotype file: airwave-illumina_exome-all_chrs.QC.bed/bim/fam

phenotype file: airwave-NMR-blood.pheno

covariates file: airwave-illumina_exome-all_chrs-NMR-blood.cov

Please note principal components are run for all files and don't need to be
provided though.

Starting files must be plink formatted and must be named as above with the
suffix '.QC.bed/bim/fam'.

If your genotype file is processed and ready make sure it is named:

genotype file: airwave-illumina_exome-all_chrs.geno

with the corresponding pheno file as above.

Output files get named based on the input files. The script assumes "cohort" is
the same for input files and matches pairs of files accordingly.

Depending on the input and arguments you might get as output:

airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.cis
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.trans
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.qqplot.svg
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.degrees_condition.txt
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.log

File names can get long so use abbreviations or short versions.

In the pipeline.yml file you can try to override output file names and
pass these as options for some of the scripts called.

If you do not use an outfile name and your files do not follow the naming above
you will get errors in the pipeline or something like:

"SNP.txt-NA-NA-NA-NA.MxEQTL"


Pipeline output
===============

Various plots (qqplots, boxplots, pca plots, etc.) and the results from the QTL
association analyses.

All files are saved in the working directory.


Requirements
============

See requirements files and Dockerfile for full information on the requirements.


Documentation
=============

    For more information see:

        |url|

'''
################
# Get modules needed:
import sys
import os
import re
import subprocess
import pandas as pd

# Pipeline:
from ruffus import *
# Ruffus combinatorics doesn't get imported:
from ruffus.combinatorics import *

# Database:
import sqlite3

# CGAT tools:
import CGATCore.Pipeline as P
import CGATCore.Experiment as E
import CGATCore.IOTools as IOTools


# TO DO: this doesn't get used at the moment and with entry point pipeline_QTL it is not
# recognised
# Import this project's module, uncomment if building something more elaborate:
#import pipeline_QTL.PipelineQTL as QTL
import PipelineQTL as QTL

# Import additional packages:
# Set path if necessary:
# os.system('''export PATH="~/xxxx/xxxx:$PATH"''')
################

################
# Load options from the config file
# Pipeline configuration

#PARAMS = P.Parameters.get_params()
#PARAMS = P.Parameters.get_parameters(getParamsFiles()) # works

P.get_parameters(
        ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
            "../pipeline.yml",
            "pipeline.yml"],
        )

PARAMS = P.PARAMS

# Set the python executable in case matplotlib is used in a Mac (where pythonw
# might be needed):
def get_py_exec():
    '''
    Look for the python executable. This is only in case of running on a Mac
    which needs pythonw for matplotlib for instance.
    '''

    try:
        if str('python') in PARAMS["general"]["py_exec"]:
            py_exec = '{}'.format(PARAMS["general"]["py_exec"])
    except NameError:
        E.warn('''
               You need to specify the python executable, just "python" or
               "pythonw" is needed in pipeline.yml.
               ''')
    return(py_exec)

def getINIpaths():
    '''
    Get the path to scripts for this project, e.g.
    project_xxxx/code/project_xxxx/:
    e.g. my_cmd = "%(scripts_dir)s/bam2bam.py" % P.Parameters.get_params()
    This is a legacy function now. Scripts can be run from the CLI if set in
    setup.py and P.PARAMS from cgat-core will find and load paths from yml
    files.
    pipeline_scriptsdir from cgat-core should also contain the location
    '''
    try:
        project_scripts_dir = '{}/'.format(PARAMS['general']['project_scripts_dir'])
        E.info('''
                   Location set for the projects scripts is:
                   {}
                   '''.format(project_scripts_dir)
                   )
    except KeyError:
        # Use the ini location if variable is set manually:
        project_scripts_dir = QTL.getDir()
        E.info('''
                   Location set for the projects scripts is:
                   {}
                   '''.format(project_scripts_dir)
                   )
    return(project_scripts_dir)
################


################
# TO DO: this is a workaround for now as pipeline_QTL CLI doesn't work and
# setting eg tools = PARAMS XXXX directly errors.
def populate_params():
    '''
    Access the ini/yml file and populate values needed.
    '''
    # Get command line tools to run:
    #tools = PARAMS['pipeline']['tools']
    tools = PARAMS['pipeline_tools']

    # Get the location of the pipeline specific scripts:
    project_scripts_dir = str(getINIpaths())

    # Set the name of this pipeline (for report softlinks):
    project_name = PARAMS['metadata']['project_name']

    # Set if running many input files:
    many_infiles = PARAMS['pipeline']['many_infiles']

    return(tools,
           project_scripts_dir,
           project_name,
           many_infiles)

################


################
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"]["name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations"]["database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh
################


################
##########
# Run PCA on genotype file with plink binary files as input
@jobs_limit(1) # plink outputs plink.prune.in and plink.prune.out which may
               # cause overwriting
@transform('*.QC.bim',
           regex('(.+).bim'),
           [r'\1.plink.prune.in',
            r'\1.plink.prune.out'],
           'SNP_exclusion_regions.txt')
def prune_SNPs(infile, outfiles, exclude):
    '''
    Prune/thin genotype data by LD using plink.
    Requires a bim plink file as input and a bed file with regions to exclude
    specifically named 'SNP_exclusion_regions.txt'
    '''
    # Add any options passed to the ini file for flashpca:
    tool_options = PARAMS['plink']['prune_options']
    if tool_options == None:
        tool_options = ''
    else:
        pass

    # Split at the last suffix separated by '.' to keep only prefix that plink
    # needs:
    infile = infile.rsplit('.', 1)[0]
    outfile1 = outfiles[0]
    outfile2 = outfiles[1]

    statement = '''
                plink \
                --bfile %(infile)s \
                --exclude range %(exclude)s \
                %(tool_options)s ;
                mv plink.prune.in %(outfile1)s ;
                mv plink.prune.out %(outfile2)s ;
                mv plink.log %(infile)s.plink.prune_SNPs.log
                '''
    P.run(statement)

@follows(prune_SNPs)
@transform('*.plink.prune.in',
           regex('(.+).plink.prune.in'),
           r'\1.plink.extracted.touch',
           r'\1',
           r'\1.plink.extracted')
def extract_SNPs(infile, outfile, plink_infile_prefix, plink_outfile_prefix):
    '''
    Extract SNPs after pruning using plink. Output can then be passed to
    Flashpca2.
    '''
    tool_options = PARAMS['plink']['extract_options']
    if tool_options == None:
        tool_options = ''
    else:
        pass

    statement = '''
                plink \
                --bfile %(plink_infile_prefix)s \
                --extract %(infile)s \
                --make-bed \
                --out %(plink_outfile_prefix)s \
                %(tool_options)s ;
                touch %(outfile)s
                '''
                # plink option with --out generate a log with that prefix, no
                # need to rename it
    P.run(statement)

    # TO DO: Needs a post-task touch otherwise tasks are re-run
    # Remove intermediate files:
    #E.info('''
    #       Deleting intermediate files:
    #       *.plink.prune.out
    #       *.plink.prune.in
    #       ''')
    #statement = '''
    #            rm -f *.plink.prune.out ;
    #            rm -f *.plink.prune.in
    #            '''
    #P.run(statement)


@follows(extract_SNPs)
@transform('*.QC.plink.extracted.bim',
           regex('(.+).bim'),
           r'\1.flashpca.touch',
           r'\1')
def PC_geno(infile, outfile, plink_infile_prefix):
    '''
    Run Flashpca2 on genotype data.
    Infile must be a plink bim file after QC and pruning.
    It must have the suffix '.QC.plink.pruned.bim'
    '''
    # Add any options passed to the ini file for flashpca:
    tool_options = PARAMS['flashpca']['options']
    if tool_options == None:
        tool_options = ''
    else:
        pass

    prefix1 = '.flashpca'
    prefix2 = '.tsv'

    statement = '''
                flashpca \
                --bfile %(plink_infile_prefix)s \
                --verbose \
                --outload %(plink_infile_prefix)s%(prefix1)s.loadings%(prefix2)s \
                --outval %(plink_infile_prefix)s%(prefix1)s.eigenvalues%(prefix2)s \
                --outvec %(plink_infile_prefix)s%(prefix1)s.eigenvectors%(prefix2)s \
                --outpc %(plink_infile_prefix)s%(prefix1)s.pcs%(prefix2)s \
                --outpve %(plink_infile_prefix)s%(prefix1)s.pve%(prefix2)s \
                %(tool_options)s ;
                touch %(outfile)s
                '''
    P.run(statement)

    # TO DO: Needs a post-task touch otherwise tasks are re-run
    # Remove intermediate files:
    #E.info('''
    #       Deleting intermediate files:
    #       *{}.bim/bed/fam
    #       '''.format(%(plink_infile_prefix)s))
    #statement = '''
    #            rm -f *%(plink_infile_prefix)s.bim *%(plink_infile_prefix)s.bed *%(plink_infile_prefix)s.fam
    #            '''
    #P.run(statement)

@follows(PC_geno)
@transform('*.QC.plink.extracted.flashpca.pcs.tsv',
           regex('(.+).pcs.tsv'),
           add_inputs(r'\1.pve.tsv'),
           r'\1.svg.touch')
def plot_PC_geno(infiles, outfile):
    '''
    Plot the results from flashpca2 on genotype data.
    '''

    # Plot flashpca results:
    tool_options = PARAMS['flashpca_plot']['options']
    if tool_options == None:
        tool_options = ''
    else:
        pass

    pcs = infiles[0]
    pve = infiles[1]

    statement = '''
                plot_flashpca.R \
                --pcs %(pcs)s \
                --pve %(pve)s \
                %(tool_options)s ;
                touch %(outfile)s
                '''
    P.run(statement)
##########


##########
# Get PCs for molecular pheno data:
@follows(plot_PC_geno)
@transform('*.pheno',
           suffix('.pheno'),
           '.pheno.pca.tsv')
def PC_and_plot_pheno(infile, outfile):
    '''
    Run PCA on molecular phenotype data using run_PCA.R from stats_utils.
    Inputs must be tab separated, rows must be features (phenotypes, variables)
    and columns must be samples (individuals).
    Outputs a tsv file and plot.
    See stats_utils and run_PCA.R -h for more info.
    '''
    # Add any options passed to the ini file for :
    tool_options = PARAMS['run_PCA']['options']
    if tool_options == None:
        tool_options = ''
    else:
        pass

    statement = '''
                run_PCA.R \
                -I %(infile)s \
                -O %(outfile)s \
                %(tool_options)s ;
                '''
    P.run(statement)
##########


##########
# Prepare files for MxEQTL
# Transform plink to MxEQTL format:
# TO DO: add:
#@active_if('matrixeqtl' in tools)
# which seems to cause errors with PARAMS can't accessing on --help for example
@follows(PC_and_plot_pheno)
@transform('*.QC.bim',
           suffix('.QC.bim'),
           '.geno')
# Files to process are the QC'd plink files.
def plink_to_geno(infile, outfile):
    '''
    Process the plink genotype files to use as input for MatrixEQTL using the
    plink_to_geno.sh script and convert plink double IDs to single with
    plink_double2singleID.R
    '''

    # Split at the last suffix separated by '.':
    infile = infile.rsplit('.', 1)[0]
    # If setup.py scripts option doesn't work use eg:
    #project_scripts_dir = str(getINIpaths() + '/matrixQTL/')
    #Rscript %(project_scripts_dir)s/run_matrixEQTL.R \
    #%(project_scripts_dir)s/

    statement = '''
                bash plink_to_geno.sh \
                        %(infile)s \
                        %(infile)s.matrixQTL \
                        %(infile)s.A-transpose \
                        %(outfile)s ;
                rm -rf *.sample *.gen *.traw *.ped *.hh *.nosex *.map ;
                '''
    P.run(statement)

    statement = '''
                plink_double2singleID.R -I %(outfile)s ;
                mv IID_%(outfile)s %(outfile)s ;
                '''
    P.run(statement)
##########

##########
# Order and match samples between geno, pheno and covariates of each:
@follows(plink_to_geno)
@transform('*geno',
           formatter('(?P<path>.+)/(?P<cohort>.+)-(?P<platform>.+)-(?P<descriptor>.+).geno'),
           add_inputs('{cohort[0]}*.pheno'),
           '{cohort[0]}.matched_geno_pheno.touch')
def orderAndMatch1(infiles, outfile):
    '''
    Order and match genotype and phenotype files using the script
    order_and_match_QTL.R
    See order_and_match_QTL.R -h
    '''
    geno = infiles[0]
    pheno = infiles[1]
    print(geno, pheno)

    statement = '''
                order_and_match_QTL.R \
                        --file1 %(geno)s \
                        --file2 %(pheno)s ;
                touch %(outfile)s
                '''
    P.run(statement)

# TO DO: transpose both files first, e.g.
# Rscript /Users/antoniob/anaconda/envs/r_test/lib/python3.5/site-packages/pipeline_QTL-0.1.0-py3.5.egg/pipeline_QTL/utilities/transpose_metabolomics.R -I matched_airwave-illumina_exome-all_chrs.flashpca.pcs.tsv
@follows(orderAndMatch1)
@transform('*.pcs.tsv',
           suffix('.pcs.tsv'),
           add_inputs('*.pheno.pca.tsv'),
           '.matched_covs.touch')
def orderAndMatch2(infiles, outfile):
    '''
    Order and match the covariates files to each other using the script
    order_and_match_QTL.R
    '''

    cov_geno = infiles[0]
    cov_pheno = infiles[1]

    statement = '''
                order_and_match_QTL.R \
                        --file1 %(cov_geno)s \
                        --file2 %(cov_pheno)s ;
                touch %(outfile)s
                '''
    P.run(statement)


# TO DO: continue here
@follows(orderAndMatch2)
@transform('matched*.pcs.tsv',
           formatter('(?P<path>.+)/matched_(?P<cohort>.+)-(?P<platform>.+)-(?P<descriptor>.+).pcs.tsv'),
           add_inputs('matched*.geno'),
           '{cohort[0]}-{platform[0]}.merged_covs')
def mergeCovs(infiles, outfile, PCs_keep_geno, PCs_keep_pheno):
    '''
    Merge covariate files from geno and pheno principal component data
    using the script merge_dataframes.R 
    '''

    cov_geno = infiles[0]
    cov_pheno = infiles[1]
    # TO DO: replace with loop
    PCs_keep_geno = 10
    PCs_keep_pheno = 35
    statement = '''
                merge_dataframes.R \
                        --file1 %(cov_geno)s \
                        --file2 %(cov_pheno)s \
                        --file1-PCs %(PCs_keep_geno)s \
                        --file1-PCs %(PCs_keep_pheno)s \
                        -O %(outfile)s
                '''
    P.run(statement)

#@follows(mergeCovs)
#@transform('*merged_covs',
#           formatter('(?P<path>.+)/matched_(?P<cohort>.+)-(?P<platform>.+)-(?P<descriptor>.+).pcs.tsv'),
#           add_inputs('matched*.pheno.pca.tsv'),
#           '{cohort[0]}-{platform[0]}.merged_covs')

##########


##########
#####
# TO DO: continue here
# Confounder loop with components
# Make optional, so that cov is passed directly or runs with loop.

# first geno, then geno plus pheno cov
# Subset geno and pheno for loop? CCs on full set but loop on subset?
# Run full set each time unless it takes days
# Run a pragmatic approach:
# Submit in batches, not loop, with step of 10 for 0, 10, 20, 30, 40, 50, 60
# Then compare results, pick highest, loop with e.g. step of 2 within
#####

#####
# Call mxeqtl geno PCs, pick most significant

#####

#####
# Join geno PCs to pheno PCs

# Run mxeqtl in batches of step e.g. 10


# pick most significant

# Run within for step of 2, pick most significant


# Delete intermediary files
#####

#####
# Run order and match of geno, pheno, covs

#####
##########


##########
# Run matrixeqtl
@follows(populate_params)
#@active_if('matrixeqtl' in tools)
@follows(mergeCovs)
@transform('*.geno',
           formatter('(?P<path>.+)/(?P<cohort>.+)-(?P<platform>.+)-(?P<descriptor>.+).geno'),
           add_inputs(['*.pheno',
                       '*.cov',
                       ]),
           '{cohort[0]}-{platform[0]}.MxEQTL.touch')
def run_MxEQTL(infiles, outfile):
    '''
    Run MatrixEQTL wrapper script for a single set of files.
    '''
    # Set up infiles:
    geno_file = infiles[0]
    pheno_file = infiles[1][0]
    cov_file = infiles[1][1]

    # Check there is an actual covariates file present, otherwise run without:
    if cov_file:
        pass
    else:
        cov_file = None

    tool_options = PARAMS['matrixeqtl']['options']
    #project_scripts_dir = str(getINIpaths() + '/matrixQTL/')
    #Rscript %(project_scripts_dir)s/run_matrixEQTL.R \
    #%(project_scripts_dir)s/

    statement = '''
                run_matrixEQTL.R \
                --gex %(pheno_file)s \
                --geno %(geno_file)s \
                --cov %(cov_file)s \
                %(tool_options)s ;
                touch %(outfile)s
                '''
    P.run(statement)


@follows(run_MxEQTL)
@transform('*.tsv', suffix('.tsv'), '.tsv.load')
def load_MxEQTL(infile, outfile):
    '''
    Load the results of run_MxEQTL() into an SQL database.
    '''
    P.load(infile, outfile)
##########


##########
# Run some other tool:


##########
################


################
# Copy to log enviroment from conda:
# TO DO: add to pipeline.py?
def conda_info():
    '''
    Print to screen conda information and packages installed.
    '''

    statement = '''conda info -a ;
                   conda list -e > conda_packages.txt ;
                   conda list --show-channel-urls ;
                   conda env export > environment.yml
                '''
    P.run(statement)

################

################
# Specify function to create reports pre-configured with sphinx-quickstart:

# Convert any svg files to PDF if needed:
@transform('*.svg', suffix('.svg'), '.pdf')
def svgToPDF(infile, outfile):
    '''
    Simple conversion of svg to pdf files with inkscape
    '''
    statement = '''
                inkscape --without-gui \
                         --export-area-drawing \
                         --export-margin=2 \
                         --file=%(infile)s \
                         --export-pdf=%(outfile)s
                '''
    P.run(statement)


# Build the report:
@follows(conda_info, svgToPDF)
def make_report():
    ''' Generates html and pdf versions of restructuredText files
        using sphinx-quickstart pre-configured files (conf.py and Makefile).
        Pre-configured files need to be in a pre-existing report directory.
        Existing reports are overwritten.
    '''
    report_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               'pipeline_report'
                                               ))
    print('Copying report templates from: {}'.format(report_path))

    if (os.path.exists(report_dir) and
            os.path.isdir(report_dir) and not
            os.listdir(report_dir)):
        statement = '''cp %(report_path)s/* pipeline_report ;
                       cd pipeline_report ;
                       make html ;
                       ln -sf _build/html/report_pipeline_pq_example.html . ;
                       make latexpdf ;
                       ln -sf _build/latex/pq_example.pdf .
                    '''
        E.info("Building pdf and html versions of your rst files.")
        P.run(statement)

    elif (os.path.exists(report_dir) and
            os.path.isdir(report_dir) and
            os.listdir(report_dir)):
        sys.exit(''' {} exists, not overwriting. You can manually run:
                       make html ;
                       ln -sf _build/html/report_pipeline_pq_example.html . ;
                       make latexpdf ;
                       ln -sf _build/latex/pq_example.pdf .
                       Or delete the folder and re-run make_report
                 '''.format(report_dir))

    else:
        sys.exit(''' The directory "pipeline_report" does not exist.
                     Are the paths correct?
                     Template files were tried to be copied from:
                     {}
                     You can also manually copy files and run "make html" or
                     "make latexpdf".
                 '''.format(report_path))
    return
################

################
# Create the "full" pipeline target to run all functions specified
@follows(load_MxEQTL, conda_info)
def full():
    pass
################

################
# TO DO:
# Check if docopt and argparse can play to show my_pipeline options and P.py
# options

def main():
    sys.exit(P.main(sys.argv))

#def main(argv=None):
#    if argv is None:
#        argv = sys.argv
#    P.main(argv)
################

################
if __name__ == "__main__":
    main()
    #sys.exit(P.main(sys.argv))
################
