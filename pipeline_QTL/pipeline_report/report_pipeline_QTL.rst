###########################
pipeline_QTL report
###########################

Authored by:

|all_author_names|

Date: |date|

Keywords: |keywords|

version = |version|

Licensed as |license|

Check the project at:

|project_url|

Correspondence to: |author_email|


.. See rst-basics_ for webpages and tutorials.

.. .. _rst-basics: https://github.com/EpiCompBio/welcome/blob/master/rst_basics.rst


Introduction
############


QTL pipeline report

|short_description|

This file was created using Sphinx_. To modify it you'll need to change the
following files:

- Makefile
- conf.py
- index.rst
- report_pipeline_pq_example.rst
- report_substitution_vars.rst

See the `Sphinx tutorial`_ to get a better idea of what's happening.

.. _Sphinx: http://www.sphinx-doc.org
.. _`Sphinx tutorial`: http://www.sphinx-doc.org/en/stable/tutorial.html


Results
#######

Figures
============

.. figure:: ../top_10_PCs_simulated2-NMR-blood.pheno.transposed.pcs.tsv.*
   :scale: 100 %
   :alt: pheno_pca.svg

   **Principal component analysis of molecular phenotypes.** Principal components (PC) for molecular phenotype data. Top-left: Proportion of variance explained for the top 100 PCs. Other panels show scatterplots of the top PCs.


.. figure:: ../top_10_PCs_simulated2-dummy_binary-all_chrs.QC.plink.extracted.flashpca.*
   :scale: 100 %
   :alt: geno_pca.svg

   **Principal component analysis of genotypes.** Principal components (PC) for molecular phenotype data. Top-left: Proportion of variance explained for the top 100 PCs. Other panels show scatterplots of the top PCs.


.. figure:: ../simulated2-dummy_binary-all_chrs-NMR-blood.MxEQTL.qqplot.*
   :scale: 100 %
   :alt: qqplot.svg

   **Quantile-quantile plot of fold change QTLs.** Expected versus observed p-values of cis (red, local) and trans (blue, distant) QTLs. If cis was not specified only trans will appear.


Tables
============


.. csv-table::
   :file: ../simulated2-dummy_binary-all_chrs-NMR-blood.MxEQTL.counts.tsv
   :delim: tab



-----



.. csv-table:: 
   :file: ../XXXXXX.MxEQTL.cis.tsv
   :delim: tab



-----



.. csv-table:: 
   :file: ../XXXX.MxEQTL.trans.tsv
   :delim: tab



-----



.. csv-table::
   :file: ../simulated2-dummy_binary-all_chrs-NMR-blood.MxEQTL.tsv
   :delim: tab



References
##########

Code used is available at `pipeline_QTL`_.

.. _`pipeline_QTL`: https://github.com/EpiCompBio/pipeline_QTL

References, e.g. :cite:`RN2398`, are pulled from a bibtex file which must be
specified at the bottom of the page.

See the sphinxcontrib-bibtex_ extension for details.

There are other ways of including citations, see the `citation directive`_ for a simple approach.

.. _sphinxcontrib-bibtex: https://github.com/mcmtroffaes/sphinxcontrib-bibtex

.. _`citation directive`: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#citations


.. bibliography:: scipy_references.bib

