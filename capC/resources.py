# capC-MAP - a software package for analysis of Capture-C data
# Copyright (C) 2018, University of Edinburgh
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 
#  Written by Chris Brackley <C.Brackley@ed.ac.uk>
# 


# Reminder - here use bed file convention that bp start counting from 0
#          - SAM files use a different convention, and start from 1

import sys
import os
import re

# List of dependencies
dependencies=["cutadapt","samtools","bowtie"]

# command names for external capC-MAP executables
capCmap_extern = { "main": "capCmain",
                   "digestfastq": "capCdigestfastq",
                   "pair2bg": "capCpair2bg",
                   "pileup2binned": "capCpileup2binned",
                   "location2fragment": "capClocation2fragment"
               }


# Bowtie arguments CONSERVATIVE or RELAXED
bowtie_args_C_nocombine="-m 1 --best --strata --sam --chunkmb 256"
bowtie_args_C_combine="--best --strata --sam --chunkmb 256" # -m option to be added
bowtie_args_R="-m 10 --best --strata --sam --chunkmb 256"

# adapter sequences
adapter_FOR = "AGATCGGAAGAGC"
adapter_REV = "AGATCGGAAGAGC"

# Get the list of restriction enzymes
this_dir, this_filename = os.path.split(__file__)
enzymesfile = os.path.join(this_dir, "resources", "restriction_enzymes.txt")
if not ( os.path.isfile(enzymesfile) ):
    sys.stdout.write("ERROR : Cannot find the enzymes resources file %s\n"%enzymesfile)
    sys.exit()
enzyme_list = {}
with open(enzymesfile) as f:
    for line in f:
        ll = line.split()
        if (ll[0][0] != "#"):
            if not len(ll) == 3:
                # test there are three fields
                sys.stdout.write("Error parsing enzyme file. Line: %s\n"%line)
                sys.exit()                   
            if not re.match(r'^[ATCG]+$',ll[1]):
                # test ll[1] is only ATGC
                sys.stdout.write("Error parsing enzyme file. Line: %s\n"%line)
                sys.exit()
            if not ll[2].isdigit():
                # test ll[2] is int 
                sys.stdout.write("Error parsing enzyme file. Line: %s\n"%line)
                #print ll[2],ll[2].isdigit()
                sys.exit()
            if int(ll[2])>len(ll[1]):
                # test ll[2] is less than len(ll[1])
                sys.stdout.write("Error parsing enzyme file. Line: %s\n"%line)
                sys.exit()
            # add as upper case to dictionart
            enzyme_list[ll[0].upper()] = ( ll[1].upper(), int(ll[2]) )


