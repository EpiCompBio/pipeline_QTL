.. include:: substitution_vars.rst

.. copy across your travis "build..." logo so that it appears in your Github page

.. image:: https://travis-ci.org/EpiCompBio/pipeline_QTL.svg?branch=master
    :target: https://travis-ci.org/EpiCompBio/pipeline_QTL

.. do the same for ReadtheDocs image:

.. image:: https://readthedocs.org/projects/pipeline-qtl/badge/?version=latest
    :target: http://pipeline-qtl.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. .. Zenodo gives a number instead, this needs to be put in manually here:
.. .. image:: https://zenodo.org/badge/#######.svg
      :target: https://zenodo.org/badge/latestdoi/#####


-----

**UNDER ACTIVE DEVELOPMENT AND NOT READY FOR USE**

-----


################################################
pipeline_QTL
################################################


.. The following is a modified template from RTD
    http://www.writethedocs.org/guide/writing/beginners-guide-to-docs/#id1

.. For a discussion/approach see 
    http://tom.preston-werner.com/2010/08/23/readme-driven-development.html

Python based pipeline for quantitative trait loci analysis.

For the general framework and workflow see project_quickstart_ docs and cgat-core_.

.. _project_quickstart: https://github.com/AntonioJBT/project_quickstart

.. _cgat-core: https://github.com/cgat-developers/cgat-core

Features
--------

- Currently runs MatrixEQTL on an arbitrary number of inputs and outputs association files, basic plots such as qqplots.


Requirements
------------

See requirements files and Dockerfile for full information. At the least you'll need:

* cgat-core_
* R >= 3.2
* Python >= 3.5
* r-matrixeqtl
* r-docopt
* r-data.table
* r-ggplot2
* stats_utils_

.. _stats_utils: https://github.com/EpiCompBio/stats_utils


Installation
------------

.. code-block:: bash
   
   pip install git+git://github.com/EpiCompBio/pipeline_QTL.git


To use
------

.. code-block:: bash 

    # Create a folder or a whole data science project, e.g. project_quickstart -n QTL_project
    cd QTL_project/results
    mkdir tests ; cd tests
    # Copy or generate test files if needed. See scripts in this repository:
    # https://github.com/EpiCompBio/pipeline_QTL/tree/master/tests
    # You can also download test files from e.g.:
    # wget -nH -np -r --cut-dirs=4 -A .txt http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/
    # https://github.com/gabraham/flashpca/tree/master/HapMap3 :
    # wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
    
    # Once you have QC'd binary files in plink format and a molecular phenotype file with samples in columns and variables in rows, run: 
    pipeline_QTL --help
    pipeline_QTL config
    # Edit pipeline_QTL.ini to adjust the parameters you want, this is essential
    pipeline_QTL show full -v 5
    pipeline_QTL printconfig
    pipeline_QTL plot full -v 5
    pipeline_QTL make full --local -v 5
    pipeline_QTL make make_report --local -v 5
    open pipeline_report/_build/latex/pipeline_QTL.pdf


Bugs and Contributions
-------------------------

- Pull requests welcome!
- Please report bugs using the `Issue Tracker`_

.. _`Issue Tracker`: https://github.com/EpiCompBio/pipeline_QTL/issues
