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

from Bio import SeqIO,Seq 

from . import resources as rs

# e.g. split site for dpnII is +3, so that fragments start with TC and end
# with GA.

# biopython counts bases from 0
# bed files count bases from 0, and are half open (the end position is not
# included)
# so here need to convert split site to 0-based

def do_genomedigest(args):
    # take a fasta file of the geneome, and generate list of 
    # restriction enzyme fragments

    fastafilefullpath = os.path.realpath(args.inputfasta[0])
    outputbedfullpath = os.path.realpath(args.outputbed[0])
    outputbedfilenameonly = os.path.basename(args.outputbed[0])
    enzyme_name = args.enzymename[0]

    # load enzyme list from resource file
    if not enzyme_name.upper() in rs.enzyme_list:
        sys.stdout.write("Unknown restriction enzyme %s ... \n"%enzyme_name)
        sys.stdout.write("New enzymes can be added to the file\n")
        sys.stdout.write("     %s\n"%rs.enzymesfile)
        sys.stdout.write("... Exiting.\n")
        sys.exit()  
    site,split = rs.enzyme_list[enzyme_name.upper()]
    split = split - 1 # convert to 0 based

    # Check input file is there
    if ( (not os.path.isfile(fastafilefullpath))  ):
        sys.stdout.write("Cannot find file %s ... Exiting.\n"%fastafilefullpath)
        sys.exit()

    # Check output file is not there, and create it
    if ( os.path.isfile(outputbedfullpath) ):
        sys.stdout.write("File %s already exists. ... Exiting.\n"%outputbedfullpath)
        sys.exit()
    ouf = open(outputbedfullpath, 'w')

    # Also create a log file
    logfilename=os.path.dirname(outputbedfullpath)+"/log_%s.log"%outputbedfilenameonly
    if ( os.path.isfile(logfilename) ):
        sys.stdout.write(" Warning : Cannot create log file %s (already exists).\n"%logfilename
                         +"Log file will not be created.\n")
        logflag = 0
    else:
        logflag = 1
        logouf = open(logfilename, 'w')

    sys.stdout.write("Digesting geneome in file %s using %s (%s) and output to file %s.\n\n"%(fastafilefullpath,enzyme_name,site,outputbedfilenameonly))

    for record in SeqIO.parse(fastafilefullpath,"fasta"):
        sites = []
        a = 0
        recordUC = record.upper()  # make sure we are searching upper case
        while ( a != -1 ):
            a = recordUC.seq.find(site,a+len(site) )
            if ( a != -1 ):
                sites.append(a)

        if ( len(sites) > 0 ):
            ouf.write( "%s\t%i\t%i\n"%(record.name,0,sites[0]+split) )
            if ( len(sites) > 1 ):
                for i in range(0,len(sites)-1):
                    ouf.write( "%s\t%i\t%i\n"%(record.name,sites[i]+split,sites[i+1]+split) )
            ouf.write( "%s\t%i\t%i\n"%(record.name,sites[-1]+split,len(record)) )
        else:
            # there are no restriction enzyme sites in this chromosome
            warningline="Warning: chromosome %s (len=%i) does not contain any instance of %s."%(record.name,len(record.seq),site)
            #sys.stdout.write('\r  '+warningline)
            sys.stdout.flush()
            if ( logflag == 1 ):
                logouf.write( warningline+"\n" )
            ouf.write( "%s\t%i\t%i\n"%(record.name,0,len(record.seq)) )

    ouf.close();
    if ( logflag == 1 ):
        logouf.close()

    return 0



