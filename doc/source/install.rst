
Installation
============


Requirements
************

capC-MAP requires a C++ compiler and Python (2.7); also the following software should be installed and visible on the system path:

-  cutadapt (>=1.11)

-  bowtie (>=1.1.1 - note that capC-MAP is not compatible with bowtie2)

-  samtools (>=1.3.1)

and the following python packages should be installed:

-  BioPython

Version numbers are those which have been tested, and other versions may also work.


Installation
************

To install capC-MAP perform the following steps:

1. Download the software and unzip into a directory in your home directory.
   Or clone from git using the command:
   ::

     git clone https://git.ecdf.ed.ac.uk/cbrackle/capC-MAP.git
    

2. capC-MAP consists of a set of programs written in C++ and a Python package, 
   which both need to be installed on your system. If you have root privileges 
   on your system you can compile and install the C++ programs by running the 
   following commands in the capC-MAP root directory
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
