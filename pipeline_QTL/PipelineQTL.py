'''PipelineQTL.py - Utility functions for QTL analysis
==============================================================

QTL intro

different sources of input (:term:`fastq` files, :term:`sra`
files) of different data (single end, paired end) with different
mapping algorithms (bowtie, tophat, stampy).

This module provides utility functions to abstract some of these variations.

The pipeline does not know what kind of data it gets (a :term:`sra`
archive might contain single end or paired end data or both).

A pipeline might get several input data (:term:`fastq` and :term:`sra`
formatted files at the same time).

It implements:
   * .sra: paired-end and single-end
   * .fastq: paired-end and single-end
   * .csfasta: colour-space, single-end

The basic class :class:`SequenceCollectionProcessor`.  The aim of this
class and its derivatives is to build a sequence of command line
statements that can be send to a single node on a cluster to process
the input data.

The basic usage inside a pipeline task is as such::
    @transform()
    def mapReads(infile, outfile):
        # initialize the Tool
        m = PipelineMapping.Hisat(
             executable=P.substituteParameters(**locals())["hisat_executable"],
             strip_sequence=PARAMS["strip_sequence"])
        # build the command line statement
        statement = m.build((infile,), outfile)
        P.run()

When implementing a tool, avoid specifying algorithmic options as
class variables. Instead use an option string that can be set in
:file:`pipeline.ini`. The only arguments to a tool constructor should
pertain to pipeline integration, such as filenames, index locations,
threading and in general processing options that change the tools
input/output, as these need to be tracked by the pipeline.
The benefit of this approach is to provide compleeat control to the
user and is likely to work with different versions of a tool, i.e., if
a command line option to a tool changes, just the configuration file
needs to be changed, but the code remains the same.

Requirements:

* cufflinks >= 2.2.1
* fastq-dump >= 2.1.7
* fastqc >= 0.11.2
* fastq_screen >= 0.4.4
* sailfish >= 0.6.3
* picardtools >= 1.106
* samtools >= 1.1
* tophat >= 2.0.13 (optional)
* bowtie >= 1.0.0 (optional)
* bowtie2 >= 2.2.3 (optional)
* bwa >= 0.7.8
* gsnap >= 2014-01-21 (optional)
* star >= 2.3.0e (optional)
* bismark >= 0.12.5 (optional)
* stampy >= 1.0.23 (optional)
* butter >= 0.3.2 (optional)
* hisat >= 0.1.4 (optional)
* shortstack >3.4 (optional)

Reference
---------

'''

