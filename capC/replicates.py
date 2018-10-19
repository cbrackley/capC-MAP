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


import sys
import os.path
from argparse import Namespace

from .pipeline import ParseConfigFile,ParseTargets
from .postprocess import do_postprocess


def report2readcounts(reportfile,target):
    """ Get a list of the read counts for each target """

    Tflag = False
    found = False

    total_reads = {}

    with open(reportfile) as rf:
        for line in rf:
            word = line.split()
            if Tflag and len(word)==5:
                if word[1] in target:
                    total_reads[word[1]]=[int(word[2]),int(word[3]),int(word[4])]
                else:
                    raise RuntimeError("Unknown target %s in report file.\n"%word[1])
            if len(word)>2 and word[1]=="target" and word[2]=="name":
                Tflag = True

    # Check all targets were found
    for T in target:
        if not T in total_reads:
            raise RuntimeError("Did not find totalreads for target %s in report file.\n"%T)

    return total_reads



def do_combine_reps(args):
    """ Combine replicates """

    experimentname = args.outdir[0]
    fulloutputdirpath = os.path.realpath(args.outdir[0])
    paramsfilefullpath = os.path.realpath(args.configfile[0])
    inputdirfullpath = []

    for i in range(len(args.indir)):
        inputdirfullpath.append( os.path.realpath(args.indir[i][0]) )

    #####################################################################################
    # Check directories are there
    try:
        if ( len(inputdirfullpath) <2 ):
            raise RuntimeError("More than one replicate directory required to combine.\n")
        for i in range( len(inputdirfullpath) ):
            if ( not os.path.isdir(inputdirfullpath[i]) ):
                raise RuntimeError("Cannot find directory %s .\n"%inputdirfullpath[i])
        if ( os.path.isdir(fulloutputdirpath) ):
            raise RuntimeError("Output directory %s already exists.\n"%fulloutputdirpath)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################
    # read parameters file
    try:
        sys.stdout.write("Parsing configuration file %s ..."%paramsfilefullpath)
        params = ParseConfigFile( paramsfilefullpath )
        sys.stdout.write(" Done.\n")
    except IOError as e:
        sys.stdout.write("Error : can't open file %s\n"%paramsfilefullpath)
        sys.exit(1) # exit with error code
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################       
    # get targets list
    try:
        if not os.path.isfile(params.targfile):
            raise RuntimeError("Cannot find file %s .\n"%params.targfile)
        fullpathtargfile = os.path.abspath(params.targfile)
        target = ParseTargets(fullpathtargfile)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################
    # Create combine direcotry
    os.makedirs(fulloutputdirpath)
    sys.stdout.write("Creating output directory %s ... Done.\n"%fulloutputdirpath)

    #####################################################################################       
    # Check all files are present
    pairsfile = "captured"
    for reps in range( len(inputdirfullpath) ):
        for T in target:
            repfile = os.path.join(inputdirfullpath[reps],"%s_validpairs_%s.pairs"%(pairsfile,T))
            if not os.path.isfile(repfile):
                raise RuntimeError("Cannot find file %s .\n"%repfile)

            # also do the interchrom files
            if params.save_inter:
                repfile = os.path.join(inputdirfullpath[reps],"%s_validinterchom_%s.pairs"%(pairsfile,T))
                if not os.path.isfile(repfile):
                    raise RuntimeError("Cannot find file %s .\n"%repfile)

    #####################################################################################       
    # Change into output directory and cat all the valid pairs files together
    owd = os.getcwd()
    os.chdir(fulloutputdirpath)

    for T in target:
        outfile = "%s_validpairs_%s.pairs"%(pairsfile,T)
        with open(outfile,'w') as ouf:
            for reps in range( len(inputdirfullpath) ):
                repfile = os.path.join(inputdirfullpath[reps],"%s_validpairs_%s.pairs"%(pairsfile,T))
                with open(repfile) as inf:
                    for line in inf:
                        ouf.write(line)

        # also do the interchrom files
        if params.save_inter:
            outfile = "%s_validinterchom_%s.pairs"%(pairsfile,T)
            with open(outfile,'w') as ouf:
                for reps in range( len(inputdirfullpath) ):
                    repfile = os.path.join(inputdirfullpath[reps],"%s_validinterchom_%s.pairs"%(pairsfile,T))
                    with open(repfile) as inf:
                        for line in inf:
                            ouf.write(line)

    #####################################################################################       
    # Generate a captured report file
    
    combined_total_reads = {}
    for T in target:
        combined_total_reads[T] = [0,0,0]
    for reps in range( len(inputdirfullpath) ):
        reportfile = os.path.join(inputdirfullpath[reps],"%s_report.dat"%pairsfile)
        total_reads = report2readcounts(reportfile,target) 
        for T in target:
            combined_total_reads[T] = [sum(x) for x in zip(combined_total_reads[T],total_reads[T])]

    reportfile = "%s_report.dat"%pairsfile
    with open(reportfile,'w') as ouf:
        ouf.write("##################################################################################\n")
        ouf.write("###\n")
        ouf.write("### Report for combining replicates for Capture C Analysis\n")
        ouf.write("###\n")
        ouf.write("##################################################################################\n")
        ouf.write("###\n")
        ouf.write("### Replicates combined from the following directories\n")
        ouf.write("###\n")
        for reps in range( len(inputdirfullpath) ):
            reportfile = os.path.join(inputdirfullpath[reps],"%s_report.dat"%pairsfile)
            ouf.write("### "+reportfile+"\n")
        ouf.write("###\n")
        ouf.write("##################################################################################\n")
        ouf.write("###\n")
        ouf.write("### Per target information \n")
        ouf.write("###\n")
        ouf.write("###   target name           intra       inter       total\n")
        for T in target:
            ouf.write("###   "+T+"\t\t\t%i\t%i\t%i\n"%(combined_total_reads[T][0],
                                                combined_total_reads[T][1],
                                                combined_total_reads[T][2]) )

    #####################################################################################       
    # Do all of the binning etc. - can call do_postprocess

    os.chdir(owd)

    ppargs = Namespace()
    ppargs.outdir = args.outdir
    ppargs.configfile = args.configfile
    ppargs.pileup = True

    do_postprocess(ppargs)
