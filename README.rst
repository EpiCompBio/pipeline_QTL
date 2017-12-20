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

For the general framework and workflow see project_quickstart_ docs and CGAT_.

.. _project_quickstart: https://github.com/AntonioJBT/project_quickstart

.. _CGAT: https://github.com/CGATOxford/CGATCore

Features
--------

- Currently runs MatrixEQTL on an arbitrary number of inputs and outputs association files, basic plots such as qqplots.


Requirements
------------

See requirements files and Dockerfile for full information. At the least you'll need:

* CGATCore
* R >= 3.2
* Python >= 3.5
* r-matrixeqtl
* r-docopt
* r-data.table
* r-ggplot2


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
    # Download test files, e.g.:
    wget -nH -np -r --cut-dirs=4 -A .txt http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/
    wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
    # See https://github.com/gabraham/flashpca/tree/master/HapMap3 :
    wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
    plink --bfile hapmap3_r2_b36_fwd.consensus.qc.poly --maf 0.01 --geno 0.01 --mind 0.01 --hwe 5e-6 --filter-founders --autosome
    pipeline_QTL --help
    pipeline_QTL config
    # Edit pipeline_QTL.ini to adjust the parameters you want, this is essential
    pipeline_QTL show full
    pipeline_QTL printconfig
    pipeline_QTL plot full
    pipeline_QTL make full --local
    pipeline_QTL make make_report --local
    open pipeline_report/_build/latex/pipeline_QTL.pdf


Bugs and Contributions
-------------------------

- Pull requests welcome!
- Please report bugs using the `Issue Tracker`_

.. _`Issue Tracker`: https://github.com/EpiCompBio/pipeline_QTL/issues
