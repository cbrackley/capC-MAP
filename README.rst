capC-MAP : Analysis software for Capture-C data
***********************************************

.. image:: https://readthedocs.org/projects/capc-map/badge/?version=latest
	   :target: https://capc-map.readthedocs.io/en/latest/?badge=latest
	   :alt: Documentation Status

capC-MAP is a software tool designed to automate the processing of Capture-C data. It allows a user to run the whole analysis pipe line using a single command line. It comprises a set of C++ programs and a wrapper script written in Python. It has been designed so that the component programs can also be run in isolation for advanced usage.

Citation
========

We have a paper presenting capC-MAP currently under review. In the mean time, please take a look at our pre-print, and cite it if you use capC-MAP in your work. Pre-print : https://doi.org/10.1101/456160


Documentation
=============

Full documentation is available at http://capc-map.readthedocs.io/

Also, a small example data set with instructions on how to run an analysis is included in the ``example\`` directory, and a set of tutorials on how to perform downstream analysis and plotting of capC-MAP output are provided in the ``downstream\`` directory.


Requirements
============

The following software must be installed on your system and visible on the path

* Python (>=2.7)
* biopython package for python
* bowtie (>=1.1.1 note that bowtie2 is a distinct piece of software and is not supported)
* samtools (>=1.3.1)
* cutadapt (>=1.11)
* a C++ compiler

Version numbers are those which have been tested, and older versions may also work.


Installation
============

There are several ways to install capC-MAP. One of the easiest is if the bioconda/conda packaging system is available on your system (see https://bioconda.github.io/#install-bioconda for details on how to get it). In this case capC-MAP and all of its requirements can be installed with the single command.
::

  conda install capc-map


If you do not have conda on your system, you will need to install the requirements listed above separately. Then install capC-MAP with the following steps:

1. Download the source code or clone from github.

2. capC-MAP consists of a set of programs written in C++ and a Python package, 
   which both need to be installed on your system. If you have root privileges 
   and pip is available on your system you can compile and install the software 
   by running the following commands in the capC-MAP package directory
   ::

     ./configure
     make
     make install
     pip install .

   If you do not have pip installed, instead run
   ::

     ./configure
     make
     make install
     python setup.py install 

   If you do not have root privileges, you can install a local 
   copy of capC-MAP in your home directory using, for example, the commands
   ::

      ./configure --prefix=${HOME}/.local/
      make
      make install
      pip install --user .

   If you have installed capC-MAP in your home directory, you will need to 
   ensure the binaries directory is present on the system PATH. For example, 
   by adding the following line to your ``~/.bashrc`` file
   ::

     export PATH=$PATH:~/.local/bin


