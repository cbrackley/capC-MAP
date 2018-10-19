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
import os.path

def get_chrom_sizes(fullpathrestfragfile,outfile=None):

    chromsizes = {}
    
    with open(fullpathrestfragfile) as rf:
        for line in rf:
            word = line.split()
            if (not word[0] in chromsizes) or (word[0] in chromsizes and chromsizes[word[0]]<int(word[2])):
                chromsizes[word[0]]=int(word[2])

    if not outfile is None:
        ouf = open(outfile, 'w')
        for C in chromsizes:
            ouf.write( "%s\t%i\n"%(C,chromsizes[C]) )
        ouf.close()

    return chromsizes



def do_get_chromsizes(args):
    """ Generate chrom.sizes from fragments list """

    fullpathfragfile = os.path.realpath(args.fragmentsfile[0])
    if args.outfile is None:
        outfile = "chrom.sizes"
    else:
        outfile = args.outfile[0]
    fullpathoutfile = os.path.realpath( os.getcwd()+"/"+outfile )

    sys.stdout.write("\n####################################################\n")
    sys.stdout.write("## Generating chromosome sizes file.\n\n")
  
    #####################################################################################  
    # Check files
    try:
        if not os.path.isfile(fullpathfragfile):
            raise RuntimeError("Cannot find file %s.\n"%fullpathfragfile)
        if os.path.isfile(fullpathoutfile):
            raise RuntimeError("File %s already exists.\n"%fullpathoutfile)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code

    #####################################################################################  
    # Generate chromosome sizes
    get_chrom_sizes(fullpathfragfile,outfile=fullpathoutfile)

    sys.stdout.write("Done.\n\n")
