###########################
pipeline_QTL
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


.. figure:: ../airwave-NMR-blood-NMR-blood.MxEQTL.qqplot.*
   :scale: 100 %
   :alt: qqplot.svg

   **Quantile-quantile plot of fold change eQTLs.** Expected versus observed p-values of cis (red, local) and trans (blue, distant) QTLs.


Tables
============


.. csv-table::
   :file: ../airwave-NMR-blood-NMR-blood.MxEQTL.counts.tsv
   :delim: tab



-----



.. csv-table:: 
   :file: ../airwave-NMR-blood-NMR-blood.MxEQTL.cis.tsv
   :delim: tab



-----



.. csv-table:: 
   :file: ../airwave-NMR-blood-NMR-blood.MxEQTL.trans.tsv
   :delim: tab



-----



.. csv-table::
   :file: ../airwave-NMR-blood-NMR-blood.MxEQTL.tsv
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

