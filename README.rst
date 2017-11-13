.. include:: README_substitution_vars.rst

.. copy across your travis "build..." logo so that it appears in your Github page

.. image:: https://travis-ci.org/EpiCompBio/pipeline_QTL.svg?branch=master
    :target: https://travis-ci.org/EpiCompBio/pipeline_QTL

.. do the same for ReadtheDocs image:

.. image:: https://readthedocs.org/projects/pipeline_QTL/badge/?version=latest
    :target: http://pipeline_QTL.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. .. Zenodo gives a number instead, this needs to be put in manually here:
.. .. image:: https://zenodo.org/badge/#######.svg
      :target: https://zenodo.org/badge/latestdoi/#####

################################################
pipeline_QTL
################################################


.. The following is a modified template from RTD
    http://www.writethedocs.org/guide/writing/beginners-guide-to-docs/#id1

.. For a discussion/approach see 
    http://tom.preston-werner.com/2010/08/23/readme-driven-development.html

Python based pipeline for quantitative trait loci analysis.

Features
--------

- Python based pipeline for QTL analysis
- Currently runs MatrixEQTL on an arbitrary number of inputs and outputs association files, basic plots such as qqplots.


Installation
------------

.. code-block:: bash
    
    pip install git+git://github.com/EpiCompBio/pipeline_QTL.git


To use
------

.. code-block:: bash 

    pip install git+git://github.com/EpiCompBio/pipeline_QTL.git
    # Create a folder or a whole data science project, e.g. project_quickstart -n QTL_project
    cd QTL_project/results
    mkdir tests ; cd tests
    # Download test files, e.g.:
    wget -nH -np -r --cut-dirs=4 -A .txt http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/
    python pipeline_QTL --help
    python pipeline_QTL config
    # Edit pipeline_QTL.ini to adjust the parameters you want, this is essential
    python pipeline_QTL show full
    python pipeline_QTL printconfig
    python pipeline_QTL plot full
    python pipeline_QTL make full --local
    python pipeline_QTL make make_report --local
    open pipeline_report/_build/latex/pipeline_QTL.pdf


Contribute
----------

- Issue Tracker: github.com/EpiCompBio/pipeline_QTL/issues
- Source Code: github.com/EpiCompBio/pipeline_QTL

Support
-------

If you are having issues, please let us know.
