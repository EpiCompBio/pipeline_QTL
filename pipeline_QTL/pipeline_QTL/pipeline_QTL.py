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

These are based on CGATPipelines_ and Ruffus_, not docopt.

.. _CGATPipelines: https://github.com/CGATOxford/CGATPipelines

.. _Ruffus: http://www.ruffus.org.uk/


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

At least correctly formatted files with genotype and phenotype data.

Depending on the tool run covariate files, SNP position and probe position
files can be added.


Quality controlled (molecular) phenotype (e.g. gene expression) and genotyping data.

Optionally covariates and error covariance matrix.

Annotation files are also needed for cis vs trans analysis (SNPs positions and probe positions and any associated annotations).

Files need to be in the same format as MatrixEQTL requires:

SNPs/genes/covariates in rows, individuals in columns with (dummy) headers and row names (the first column and first row are skipped when read).

Missing values must be set as "NA".

SNP and probe position files, e.g.
   - snp146Common_MatrixEQTL_snp_pos.txt
   - biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt
must be tab separated.

See:

http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/


Naming convention for input files
=================================

Output files get named based on the input files. The script assumes "cohort" is the same for input files (but only takes it from the genotype file).

Please rename your files in the following way (use soft links to rename for example):

Infile: cohort-platform-other_descriptor.suffix

Outfile: cohort-platform_infile1-descriptor1-platform_infile2-descriptor2.new_suffix

For example:

genotype file: airwave-illumina_exome-all_chrs.geno

phenotype file: airwave-NMR-blood.pheno

covariates file: airwave-NMR-blood.cov

and depending on the input and arguments you might get:

airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.cis
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.trans
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.qqplot.svg
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.degrees_condition.txt
airwave-illumina_exome-all_chrs-NMR-blood.MxEQTL.log

File names can get long so use abbreviations or short versions.

You can also override this and simply choose your outfile prefix.

If you do not use an outfile name and your files do not follow the naming above you might get something like:

"SNP.txt-NA-NA-NA-NA.MxEQTL"

If running in a pipeline .geno, .pheno and .cov are required.

Pipeline output
===============

Namely a qqlot and tables of genotype molecular phenotype associations. These are saved in the working directory.


Requirements
============

See requirements files and Dockerfile for full information. At the least you'll need:

* CGATCore
* R >= 3.2
* Python >= 3.5
* r-matrixeqtl
* r-docopt
* r-data.table
* r-ggplot2


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

# Pipeline:
from ruffus import *

# Database:
import sqlite3

# Try getting CGAT: 
try:
    import CGAT.IOTools as IOTools
    import CGATPipelines.Pipeline as P
    import CGAT.Experiment as E

except ImportError:
    print('\n', "Warning: Couldn't import CGAT modules, these are required. Exiting...")
    raise

# required to make iteritems python2 and python3 compatible
from builtins import dict

# Import this project's module, uncomment if building something more elaborate: 
try:
    import pipeline_QTL as QTL
except ImportError:
    print("Could not import this project's module, exiting")
    raise

# Import additional packages: 
# Set path if necessary:
#os.system('''export PATH="~/xxxx/xxxx:$PATH"''')
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
                ini_file = os.path.join(os.path.abspath(path), ini_file.group())
                p_params_files.append(ini_file)
    return(p_params_files)

P.getParameters(getParamsFiles())

PARAMS = P.PARAMS
# Print the options loaded from ini files and possibly a .cgat file:
#pprint.pprint(PARAMS)
# From the command line:
#python ../code/pq_example/pipeline_pq_example/pipeline_pq_example.py printconfig


# Set global parameters here, obtained from the ini file
# e.g. get the cmd tools to run if specified:
#cmd_tools = P.asList(PARAMS["cmd_tools_to_run"])

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
        if project_scripts_dir != str('/'):
            E.info('''
                   Location set for the projects scripts is:
                   {}
                   '''.format(project_scripts_dir)
                   )
        else:
            project_scripts_dir = QTL.getDir()
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
project_scripts_dir = str(getINIpaths() + '/matrixQTL/')
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
#####
# Run matrixeqtl:
@active_if('matrixeqtl' in tools)
@mkdir('MatrixEQTL')
@transform('*.geno',
           regex(r'(.+).geno'),
           add_inputs([r'\1.pheno',
                       r'\1.cov',
                       ]),
           r'\1.touch')
def run_MxQTL(infiles, outfiles):
    '''
    Run MatrixEQTL wrapper script.
    '''
    # Set up infiles:
    geno_file = infiles[0]
    pheno_file = infiles[1][0]
    cov_file = infiles[1][1]
    print(infiles)
    # Check there is an actual covariates file present, otherwise run without:
    if cov_file:
        cov_file = cov_file
    else:
        cov_file = None

    tool_options = P.substituteParameters(**locals())["matrixeqtl_options"]

    statement = '''
                cd MatrixEQTL ;
                Rscript %(project_scripts_dir)s/run_matrixEQTL.R \
                --gex %(pheno_file)s \
                --geno %(geno_file)s \
                --cov %(cov_file)s \
                %(tool_options)s
                '''
                #-O %(outfile)s
                #--model modelANOVA
                #--pvOutputThreshold 0.05
                #--snpspos snpsloc.txt
                #--genepos geneloc.txt
                #--pvOutputThreshold.cis 0.1
                #--session mxqtl
                #--condition cond_test

    P.run()

@follows(run_MxQTL)
@transform('*.tsv', suffix('.tsv'), '.tsv.load')
def load_MxQTL(infile, outfile):
    '''
    Load the results of run_MxQTL() into an SQL database.
    '''
    P.load(infile, outfile)
#####

#####
# Run some other tool:
#####
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
# Create the "full" pipeline target to run all functions specified
@follows(load_MxQTL, conda_info)
def full():
    pass
################

################
# Specify function to create reports pre-configured with sphinx-quickstart:
@follows(full)
def make_report():
    ''' Generates html and pdf versions of restructuredText files
        using sphinx-quickstart pre-configured files (conf.py and Makefile).
        Pre-configured files need to be in a pre-existing report directory.
        Existing reports are overwritten.
    '''
    if os.path.exists('report'):
        statement = ''' cd pipeline_report ;
                        checkpoint ;
                        make html ;
                        ln -s _build/html/index.hmtl . ;
                        checkpoint ;
                        make latexpdf ;
                        ln -s _build/latex/pq_example.pdf .
                    '''
        E.info("Building pdf and html versions of your rst files.")
        P.run()

    else:
        E.stop(''' The directory "report" does not exist. Did you run the config
                   option? This should copy across templates for easier
                   reporting of your pipeline.
                   If you changed the dir names, just go in and run "make html" or
                   "make latexpdf" or follow Sphinx docs.
                ''')

    return
################

def main():
    sys.exit(P.main(sys.argv))

################
if __name__ == "__main__":
    main()
    #sys.exit(P.main(sys.argv))
################
