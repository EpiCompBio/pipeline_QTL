'''
|project_name|
===================

:Author: |author_name|
:Release: |version|
:Date: |today|


Overview
========

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

A configuration file was created at the same time as this script.

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


Naming convention for input files
=================================

Please rename your files in the following way (use soft links to rename for
example):

Infile: cohort-platform-other_descriptor.suffix

Outfile: cohort-platform_infile1-descriptor1-platform_infile2-descriptor2.new_suffix

For example:

genotype file: airwave-illumina_exome-all_chrs.geno

phenotype file: airwave-NMR-blood.pheno

covariates file: airwave-illumina_exome-all_chrs-NMR-blood.cov


Output files get named based on the input files. The script assumes "cohort" is
the same for input files (but only takes it from the genotype file).

Depending on the input and arguments you might get as output:

airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.cis
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.trans
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.qqplot.svg
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.degrees_condition.txt
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.log

File names can get long so use abbreviations or short versions.

In the pipeline_XXXX.ini file you can try to override output file names and
pass these as options for some of the scripts called.

If you do not use an outfile name and your files do not follow the naming above
you might errors in the pipeline or something like:

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

# Try getting CGAT:
try:
    # import CGAT.IOTools as IOTools
    import CGATPipelines.Pipeline as P
    import CGAT.Experiment as E

except ImportError:
    print('\n', '''Warning: Couldn't import CGAT modules,
                   these are required. Exiting...''')
    raise

# required to make iteritems python2 and python3 compatible
# from builtins import dict

# Import this project's module, uncomment if building something more elaborate:
try:
    import pipeline_QTL as QTL
except ImportError:
    print("Could not import this project's module, exiting")
    raise

# Import additional packages:
# Set path if necessary:
# os.system('''export PATH="~/xxxx/xxxx:$PATH"''')
################

################
# Load options from the config file
# Pipeline configuration
ini_paths = [os.path.abspath(os.path.dirname(sys.argv[0])),
             "../",
             os.getcwd(),
             ]


def getParamsFiles(paths = ini_paths):
    '''
    Search for python ini files in given paths, append files with full
    paths for P.getParameters() to read.
    Current paths given are:
    where this code is executing, one up, current directory
    '''
    p_params_files = []
    for path in ini_paths:
        for f in os.listdir(os.path.abspath(path)):
            ini_file = re.search(r'pipelin(.*).ini', f)
            if ini_file:
                ini_file = os.path.join(os.path.abspath(path),
                                        ini_file.group())
                p_params_files.append(ini_file)
    return(p_params_files)


P.getParameters(getParamsFiles())

PARAMS = P.PARAMS
# Print the options loaded from ini files and possibly a .cgat file:
# pprint.pprint(PARAMS)
# From the command line:
# python ../code/pq_example/pipeline_pq_example/pipeline_pq_example.py printconfig


# Set global parameters here, obtained from the ini file
# e.g. get the cmd tools to run if specified:
# cmd_tools = P.asList(PARAMS["cmd_tools_to_run"])


def get_py_exec():
    '''
    Look for the python executable. This is only in case of running on a Mac
    which needs pythonw for matplotlib for instance.
    '''
    try:
        PARAMS["py_exec"]
        py_exec = '{}'.format(PARAMS['py_exec'])
    except NameError:
        E.warn('''
               You need to specify the python executable, just "python" or
               "pythonw" is needed. Trying to guess now...
               ''')
    else:
        test_cmd = subprocess.check_output(['which', 'pythonw'])
        sys_return = re.search(r'(.*)pythonw', str(test_cmd))
        if sys_return:
            py_exec = 'pythonw'
        else:
            py_exec = 'python'
    return(py_exec)


def getINIpaths():
    '''
    Get the path to scripts for this project, e.g.
    project_xxxx/code/project_xxxx/:
    e.g. my_cmd = "%(scripts_dir)s/bam2bam.py" % P.getParams()
    '''
    try:
        project_scripts_dir = '{}/'.format(PARAMS['project_scripts_dir'])
        if project_scripts_dir == str('/'):
            # dir not set in ini file so use installation directory:
            project_scripts_dir = QTL.getDir()
            E.info('''
                   Location set for the projects scripts is:
                   {}
                   '''.format(project_scripts_dir)
                   )
        else:
            # Use the ini location if variable is set manually:
            project_scripts_dir = '{}/'.format(PARAMS['project_scripts_dir'])
            E.info('''
                   Location set for the projects scripts is:
                   {}
                   '''.format(project_scripts_dir)
                   )
    except KeyError:
        E.warn('''
               Could not set project scripts location, this needs to be
               specified in the project ini file.
               ''')
        raise

    return(project_scripts_dir)
################


################
# Get command line tools to run:
tools = P.asList(PARAMS["pipeline_tools"])

# Get the location of the pipeline specific scripts:
project_scripts_dir = str(getINIpaths())

# Set the name of this pipeline (for report softlinks):
project_name = PARAMS['metadata_project_name']
# 'pipeline_QTL'

# Set if running many input files:
many_infiles = P.asList(PARAMS["pipeline_many_infiles"])
################


################
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh
################


################
##########
# Naming:
# Infile: cohort-platform-other_descriptor.suffix
# Outfile: cohort-platform1-descriptor1-platform2-descriptor2.new_suffix

# TO DO: some R scripts and bash scripts have the paths hardcoded, correct
# these.

# Process phenotype and genotype files
# Run PCA on each type of file

#@posttask(touch_file("prune_SNPs.touch"))
@transform('*.bim',
           suffix('.bim'),
           '.pruned', 'SNP_exclusion_regions.txt')
def prune_SNPs(infile, outfile, exclude):
    '''
    Prune genotype data using plink.
    Requires a bim plink file as input and a bed file with regions to exclude
    named 'exclusion<any other characters>.txt'
    '''
    # Add any options passed to the ini file for flashpca:
    tool_options = P.substituteParameters(**locals())["plink_exclude_options"]

    # Split at the last suffix separated by '.':
    infile = infile.rsplit('.', 1)[0]
    statement = '''
                plink \
                --bfile %(infile)s \
                --exclude range %(exclude)s \
                %(tool_options)s ;
                checkpoint
                '''
    P.run()

    tool_options = P.substituteParameters(**locals())["plink_prune_options"]
    statement = '''
                plink \
                --bfile %(infile)s \
                --extract plink.prune.in \
                --make-bed \
                --out %(outfile)s \
                %(tool_options)s ;
                '''
    P.run()

#@posttask(touch_file("PC_geno.touch"))
@follows(prune_SNPs)
@transform('*.bim',
           suffix('.bim'),
           '.PC.touch')
def PC_geno(infile, outfile):
    '''
    Run Flashpca2 on genotype data.
    Infile must be a plink bim file after QC and pruning.
    The file must have the suffix '.pruned.bim' as generated by
    prune_SNPs() function in this script.
    '''
    # Add any options passed to the ini file for flashpca:
    tool_options = P.substituteParameters(**locals())["flashpca_options"]

    infile = infile.rsplit('.', 1)[0]
    prefix1 = '.flashpca'
    prefix2 = '.tsv'
    statement = '''
                flashpca \
                --bfile %(infile)s \
                --verbose \
                --outload %(infile)s%(prefix1)s.loadings%(prefix2)s \
                --outval %(infile)s%(prefix1)s.eigenvalues%(prefix2)s \
                --outvec %(infile)s%(prefix1)s.eigenvectors%(prefix2)s \
                --outpc %(infile)s%(prefix1)s.pcs%(prefix2)s \
                --outpve %(infile)s%(prefix1)s.pve%(prefix2)s \
                %(tool_options)s ;
                checkpoint ;
                touch %(outfile)s
                '''
    P.run()

    # Remove intermediate files:
    #statement = '''
    #            rm -f %(infile)s.pruned.* ;
    #            rm -rf plink*
    #            '''
    #P.run()
# TO DO: touch and posttask don't work on PC_geno and prune_SNPs....

@follows(PC_geno)
@transform('*.pcs.tsv',
           regex('(.+).(.+).(.+).tsv'),
           add_inputs('*.pve.tsv'),
           r'\1.svg.touch')
def plot_PC_geno(infile, outfile):
    '''
    Plot the results from flashpca2 on genotype data.
    '''

    # Plot flashpca results:
    tool_options = P.substituteParameters(**locals())["flashpca_plot_options"]
    project_scripts_dir = str(getINIpaths() + '/scripts/utilities/')

    pcs = infile[0]
    pve = infile[1]

    statement = '''
                Rscript %(project_scripts_dir)s/plot_flashpca.R \
                --pcs %(pcs)s \
                --pve %(pve)s \
                %(tool_options)s ;
                checkpoint
                '''
    P.run()


@follows(plot_PC_geno)
@transform('*.bim',
           suffix('.bim'),
           '.A-transpose.matrixQTL.geno')
def plink_to_geno(infile, outfile):
    '''
    Process the plink genotype files to use as input for MatrixEQTL.
    '''
    # Add any options passed to the ini file for :
    #tool_options = P.substituteParameters(**locals())["_options"]

    project_scripts_dir = str(getINIpaths() + '/scripts/utilities/')

    # Split at the last suffix separated by '.':
    infile = infile.rsplit('.', 1)[0]

    statement = '''
                bash %(project_scripts_dir)s/plink_to_geno.sh \
                        %(infile)s \
                        %(infile)s.matrixQTL \
                        %(infile)s.A-transpose \
                        %(outfile)s ;
                checkpoint
                '''
    P.run()

    # TO DO: delete intermediary files

# TO DO, add:
# bash /Users/antoniob/Documents/github.dir/EpiCompBio/pipeline_QTL/scripts/utilities/plink_to_geno.sh airwave-illumina_exome-all_chrs airwave-illumina_exome-all_chrs.matrixQTL airwave-illumina_exome-all_chrs.A-transpose airwave-illumina_exome-all_chrs.A-transpose.matrixQTL.geno
# Rscript /Users/antoniob/Documents/github.dir/EpiCompBio/pipeline_QTL/scripts/utilities/plink_double2singleID.R -I airwave-illumina_exome-all_chrs.A-transpose.matrixQTL.geno
# mv IID_airwave-illumina_exome-all_chrs.A-transpose.matrixQTL.geno airwave-illumina_exome-all_chrs.A-transpose.matrixQTL.geno 
# Rscript /Users/antoniob/Documents/github.dir/EpiCompBio/pipeline_QTL/scripts/utilities/order_and_match_QTL.R --file1 airwave-illumina_exome-all_chrs.A-transpose.matrixQTL.geno --file2 AIRWAVE-CPMG_BatchCorrected_log_Var_Data_Sample-plasma.transposed.tsv

@transform('*.pheno',
           suffix('.pheno'),
           '.pheno.pca.touch')
def PC_pheno(infile, outfile):
    '''
    Run PCA on molecular phenotype data.
    '''
    # Add any options passed to the ini file for :
    #tool_options = P.substituteParameters(**locals())["_options"]

    project_scripts_dir = str(getINIpaths() + '/scripts/utilities/')
    tool_options = P.substituteParameters(**locals())["run_PCA_options"]
    statement = '''
                Rscript %(project_scripts_dir)s/run_PCA.R \
                -I %(infile)s \
                -O %(outfile)s \
                %(tool_options)s ;
                checkpoint
                '''
    P.run()
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
@active_if('matrixeqtl' in tools)
@follows(PC_pheno, plink_to_geno)
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

    tool_options = P.substituteParameters(**locals())["matrixeqtl_options"]
    project_scripts_dir = str(getINIpaths() + '/scripts/matrixQTL/')

    statement = '''
                Rscript %(project_scripts_dir)s/run_matrixEQTL.R \
                --gex %(pheno_file)s \
                --geno %(geno_file)s \
                --cov %(cov_file)s \
                %(tool_options)s ;
                checkpoint ;
                touch %(outfile)s
                '''
    P.run()


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
    P.run()

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
    P.run()


# Build the report:
@follows(conda_info, svgToPDF)
def make_report():
    ''' Generates html and pdf versions of restructuredText files
        using sphinx-quickstart pre-configured files (conf.py and Makefile).
        Pre-configured files need to be in a pre-existing report directory.
        Existing reports are overwritten.
    '''
    if os.path.exists('pipeline_report'):
        statement = ''' cd pipeline_report ;
                        checkpoint ;
                        make html ;
                        checkpoint ;
                        make latexpdf
                    '''
        E.info("Building pdf and html versions of your rst files.")

    else:
        E.info(''' The directory "pipeline_report" does not exist. Did you run the config
                   option? This should copy across templates for easier
                   reporting of your pipeline.
                   If you changed the dir names, just go in and run "make html" or
                   "make latexpdf" or follow Sphinx docs.
                ''')
        sys.exit()

    if (os.path.exists('pipeline_report/_build/html/index.hmtl') and
       os.path.exists(os.path.join('pipeline_report/_build/latex/',
                                   project_name, '.pdf'))):
        statement = '''
                    ln -s pipeline_report/_build/html/index.hmtl %(project_name)s.html ;
                    ln -s pipeline_report/_build/latex/%(project_name)s.pdf .
                    '''
        E.info('''Done, linkts to the pdf and html versions of your rst files are in the main
               folder.''')
        P.run()

    else:
        E.info('''
               The html and/or latex/pdf files did not build correctly. See the
               logs and go into pipeline_report to find out. You can also try
               building the report manually with make html and make latexpdf.
               ''')
        sys.exit()

    return
################

################
# Create the "full" pipeline target to run all functions specified
#@follows(load_MxEQTL, conda_info, make_report)
@follows(run_MxEQTL, conda_info)
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
