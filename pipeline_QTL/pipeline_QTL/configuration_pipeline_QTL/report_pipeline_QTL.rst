See rst-basics_ for webpages and tutorials.

.. _rst-basics: https://github.com/EpiCompBio/welcome/blob/master/rst_basics.rst


##############
|project_name|
##############

-----

author1, author2, author3 …

affiliation1, affiliation2, affiliation3 …

Correspondence should be addressed to:

Keywords:

Running title:


-----


Abstract
########

Background: 

Methods: 

Findings: 

Interpretation:

Funding: 

Copyright: Open access article under the terms of CC BY.

Introduction
############

Include other rst files::

  .. toctree::
      :maxdepth: 2
      :numbered:
      :titlesonly:
      :glob:
      :hidden:

      intro.rst
      chapter1.rst
      chapter2.rst

See the toctree_ directive for full info.

.. _toctree: http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html#include-other-rst-files-with-the-toctree-directive


Methods
#######
      
Briefly, the primary objectives of 

Briefly, 

All data and results were handled according to the institutional guidelines in secure servers within 

-----

Laboratory procedure1

Laboratory procedure2

Laboratory procedure3

-----

Statistical analysis

A priori power for 

Method 1

.. Glob all methods*.rst files:

.. toctree::
    :glob:
    
    methods*rst

Method 2

Method 3


Statistical power calculations for 

We estimated that 

Taking a two-sided alpha of 0.05 t-test, beta of 10%, standard deviation of 0.7, fold-change of 2 and equal size group estimated a sample size of for each group.

Calculations were performed using the base package in R and pwr (v1.1-3).

-----

Quality control, normalisation and xxx

Quality control was carried out according to 

Quality assessment of 

We removed samples that failed 

We found that xxx samples 

We excluded 

Pre-processing and filtering included  

Linear regression association tests were carried out using 

We corrected for 

To explore 

The primary comparison was a 

We tested for linear or quadratic effects of 

We also performed a linear mixed model analysis with 

We used an additive linear model as implemented in the R package 
R packages were run using R 3.2.4 (R Core Team, 2016). We used ggplot2 (v. 2.1.0), package specific functions (pwr, xxx) as well as scripts from xxx for xxx and data.table (v. 1.9.6) for processing.

-----

Role of the funding source

The study received funding from xxx grant number xxx.

The funders had no role in data collection, analysis, interpretation or writing of the report. 

All authors had access to all the data in the study. 

Registration:


Results
#######

Result 1

.. Glob all legends*.rst files:

.. toctree::
    :glob:
    
    legends*rst
    interpretation*rst    

Include an image::

  .. image:: images/ball1.gif
  
Or::

  .. image:: images/xxx.png
    :height: 100
    :width: 200
    :scale: 50
    :alt: alternate text

See image_ directive full markup.

.. _image: http://docutils.sourceforge.net/docs/ref/rst/directives.html#images

Or import a figure which can have a caption and whatever else you add::

  .. figure:: xxx.jpg
      :width: 200px
      :align: center
      :height: 100px
      :alt: alternate text
      :figclass: align-center
      
      a caption would be written here as plain text. You can add more with eg::
  
    .. code-block:: python

        import image

Include a simple csv table::

  .. csv-table:: a title
     :header: "name", "firstname", "age"
     :widths: 20, 20, 10
     
     "Smith", "John", 40
     "Smith", "John, Junior", 20

See csv-table_ directive for example.

.. _csv-table: http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html#the-csv-table-directive

Result 2

Result 3

-----

We had xx% power to detect xxx or greater change in xxx with a sample size of xxx per arm in xxx

Considering all samples

We next formally tested for the

Our primary pre-defined comparison sought to define

We then considered

We next increased power by

These results would thus require validation in larger studies.

However, neither of these analyses

Given that our primary comparison

We hypothesised that 

Discussion
##########

Overall, we found

To our knowledge, this is the

Although preliminary, our results suggest that

An important question in the field is whether 

Our study has several limitations. We did not carry out 

Indeed, other studies have observed

We cannot address whether

This study shows that 

Future studies will need to 

Our study highlights 


Research in context
###################

Evidence before this study

Added value of this study

Implications of all the available evidence


Funding and acknowledgements
############################
We would like to thank all the study participants, 

XYZ was funded by xxx (Grant xxx) 

We thank the xxx with grant xxx for the generation of data.


Data access
###########
xxx data are available through ArrayExpress (xxx). 

xxx, phenotype and xxx data are available through the European Genome-Phenome Archive (EGA, request through EGASxxxx). 

Code used is available at  https://github.com/xxx .

Figure legends
##############

Figure 1:

Figure 2:

Figure 3:

Supplementary information, figures and tables
#############################################

Appendix 1: Data analysis protocol

Supplementary Figure 1:

Supplementary Figure 2:

Supplementary Figure 3:

Supplementary Table 1:

Supplementary Table 2:

Supplementary Table 3:


References
##########

References, e.g. [CIT2002]_ are defined at the bottom of the page as::

  .. [CIT2002] A citation

and called with::

  [CIT2002]_


