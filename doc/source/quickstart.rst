
Quick start
===========

A typical pipeline for analysis of capture c data has three main steps:


Build genome index and restriction enzyme fragment list.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since capC-MAP uses bowtie [Langmead2009]_ for sequence
alignment, a bowtie index for the reference genome must be built. This
requires a single fasta file containing the reference genome. For
example

::

     bowtie-build mygenome.fasta mygeneome

The same fasta file must then be used to generate a list of restriction
enzyme fragments for the genome. This is done using the ``genomedigest``
function in capC-MAP.

::

   capC-MAP genomedigest -f mygenome.fasta -r DpnII \
                         -o mygenome_dpnII_fragments.bed

where the options are as follows:

+-----------------------+-----------------------+
| ``-f mygenome.fasta`` | specifies the fasta   |
|                       | for the reference     |
|                       | genome                |
+-----------------------+-----------------------+
| ``-r DpnII``          | specifies the         |
|                       | restriction enzyme    |
|                       | used in the           |
|                       | experiment            |
+-----------------------+-----------------------+
| ``-o mygenome_dpnII_f | specifies the output  |
| ragments.bed``        | bed file              |
+-----------------------+-----------------------+

Note that the ``\`` character means that a single command is broken across lines.

Pre-built bowtie indexes for many genomes are available for download on the bowtie website (http://bowtie-bio.sourceforge.net), though a fasta file is required to build the list of restriction enzyme fragments. It is essential to ensure that the index and fragments list are built from the same reference genome (for this reason we recommend building you own index).


Perform quality control on fastq data files.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


We recommend performing standard quality control on the fastq files, for
example using the `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc>`_ software. Since the
core aim in capture c data it to enrich the library for specific
‘target’ fragments, this may get flagged up the the FastQC report. Also,
since the capture c protocol recommends fragments are sonicated to have
an average length of 200-300bp, depending on the sequencing read length, it
could be that there is a significant proportion of
read-through into the adapters – again this may be flagged
up by FastQC. By default the capC-MAP pipeline includes an adapter
trimming step.



Run capC-MAP analysis pipeline.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The main capC-MAP pipeline is then run using the command

::

   capC-MAP pipeline -c config_file.txt \
                     -o mycaptureCexperiment

where the options are as follows:

+-----------------------------+-----------------------+
| ``-c config_file.txt``      | specifies the         |
|                             | experiment            |
|                             | configuration file    |
+-----------------------------+-----------------------+
| ``-o mycaptureCexperiment`` | specifies a directory |
|                             | where all output will |
|                             | be saved              |
+-----------------------------+-----------------------+

All capC-MAP options are specified in the configuration text file. We
recommend using the example configuration file provided with capC-MAP as
a template. The output directory must not exist. Required inputs
(specified in the configuration file) are:

-  a pair of fastq files;

-  a bed file containing a list of target restriction enzyme fragments;

-  the bowtie index for the reference genome; and

-  a bed file containing a genome wide list of restriction enzyme
   fragments for the reference genome.

The following output files will be generated in the output directory:

-  an *in silico* digested version of the input sequencing data in a
   single fastq file.

-  a SAM file containing mapped reads; the same data is also given in
   the compressed BAM format.

-  a capC-MAP report file giving details of mapped fragments.

-  for each target specified in the targets bed file, a 'validpairs' file
   containing a list of all valid intrachromosomal interactions, and a
   'validinterchrom' files containing a list of valid interchromosomal
   interactions.

-  for each target specified in the targets bed file, and depending on
   options specified in the configuration file, a bedGraph showing
   binned, smoothed, and normalized to reads per million genome wide
   interactions.

Once the analysis is complete, some of these outputs are not needed for
most down-stream analyses. For example the *in silico* digested fastq
file, SAM files and the raw pairs files can usually be deleted to save
disk space.

In order to generate further bedGraphs with different binning or normalization
options without re-running the full analysis, capC-MAP can be run in
``postprocess`` mode (see section :ref:`sectionpostprocess`).

