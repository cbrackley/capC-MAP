capC-MAP : Analysis software for Capture-C data
***********************************************

.. image:: https://readthedocs.org/projects/capc-map/badge/?version=latest
	   :target: https://capc-map.readthedocs.io/en/latest/?badge=latest
	   :alt: Documentation Status

capC-MAP is a software tool designed to automate the processing of Capture-C data. It allows a user to run the whole analysis pipe line using a single command line. It comprises a set of C++ programs and a wrapper script written in Python. It has been designed so that the component programs can also be run in isolation for advanced usage.

Documentation
=============

Full documentation is available at http://capc-map.readthedocs.io/


Requirements
============

The following software must be installed on your system and visible on the path

* Python (>=2.7; Python 3 is not currently supported)
* biopython package for python
* bowtie (>=1.1.1 note that bowtie2 is a distinct piece of software and is not supported)
* samtools (>=1.3.1)
* cutadapt (>=1.11)
* a C++ compiler

Version numbers are those which have been tested, and older versions may also work.


Installation
============

To install capC-MAP perform the following steps:

1. Download the source code or clone from github.

2. capC-MAP consists of a set of programs written in C++ and a Python package, 
   which both need to be installed on your system. If you have root privileges 
   on your system you can compile and install the C++ programs by running the 
   following commands in the capC-MAP package directory
   ::

     ./configure
     make
     make install
     
   and you can install the Python package using pip with the command
   ::

      pip install .

   when in the capC-MAP package directory.

   Alternatively, if you do not have root privileges, you can install a local 
   copy of capC-MAP in your home directory using, for example, the commands
   ::

      ./configure --prefix=${HOME}/.local/
      make
      make install

   and
   ::

      pip install --user .

   If you do not have pip available on your system you can instead install the 
   python package using the command
   ::

      python setup.py install 

   or without root privileges
   ::

      python setup.py install --user

3. If you have installed capC-MAP in your home directory, you will need to 
   ensure the binaries directory is present on the system PATH. For example, 
   by adding the following line to your ``~/.bashrc`` file
   ::

     export PATH=$PATH:~/.local/bin


