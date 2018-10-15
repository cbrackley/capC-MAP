
.. only:: latex

   Introduction
   ============

capC-MAP is a software package for the analysis of sequencing data from
Capture-C experiments [Hughes2014]_ [Davies2016]_. It is actually a
suit of programs written in C++ and Python 2.7, along with a Python
wrapper script which allows a full analysis pipeline to be run using a
single command on Unix-based systems. capC-MAP was written so as to be
as easy to use as possible, but allow maximum customisation and separate
use of the component programs by advanced users. While other other tools
which can analyse Capture-C data are available, these are mostly
extensions to software designed with HiC data in mind. To our knowledge
capC-MAP is the first software dedicated to, and optimised for Capture-C
data; it is designed to be used by beginners and ‘C-method’ experts
alike.

The aim of a Capture-C experiment is to obtain an interaction profile for a
set of "target" genomic loci. This is similar to the aim of a 4C experiment,
but the method allows multiple targets to be probes in a single experiment.
This is achieved using oligo capture technology. A frequently cutting
restriction enzyme is used to fragment the DNA so as to obtain interactions at
high-resolution. Oligos are designed against a set of restriction enzyme
fragments of interest. Throughout this manual we use the term "target" to refer
to restriction enzyme fragments for which oligos have been designed, and
"reporter" for fragments which are found ligated to target fragments. The term
"target" is synonymous with the "viewpoint" or "bait" fragment in 4C.
