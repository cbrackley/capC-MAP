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
import subprocess
import re
from shutil import copyfile

from . import resources as rs
from .resfrags2chromsize import get_chrom_sizes 

class Parameters:
    def __init__(self):
        # required parameters
        self.fastq1 = ""
        self.fastq2 = ""
        self.targfile = ""
        self.indexpath = ""
        self.restfragfile = ""
        self.enzyme = ""

        # optional parameters -- defaults set
        self.processors = 1
        self.biningcount = 0
        self.binlist = []
        self.bowtie_args = rs.bowtie_args_C_nocombine
        self.exclusion = 500
        self.trim_adapters = True
        self.save_inter = False
        self.normalize = False
        self.combinemode = False
        self.dryrun = False
        self.combinecount = 1
        self.alignmode = "CONSERVATIVE"
        self.aligncustom = ""

        # flags to check if parameter is set
        self.flag_fastq1 = False
        self.flag_fastq2 = False
        self.flag_targfile = False
        self.flag_indexpath = False
        self.flag_restfragfile = False
        self.flag_enzyme = False
        


def ParseTargets(targetsfile):
    """ Read a bed file of targets and return a dictionary """
    targs = {}

    with open(targetsfile) as tf:
        for line in tf:
            word = line.split()
            if len(word) < 4 or not (word[1].isdigit() and word[2].isdigit()):
                raise RuntimeError("cannot read targets from file %s .\n"%targetsfile)
            targs[word[3]]=[word[0],int(word[1]),int(word[2])]

    return targs

        
def ParseConfigFile(configfile):
    """ Read from the config file """
    params = Parameters()

    # read the file
    with open(configfile) as cf:
        for line in cf:
            if line[0] == "#" or line[0] == "\n" or line.isspace():
                continue
            else:
                word = line.split()
                if word[0] == "FASTQ1":
                    if len(word) < 2 or params.flag_fastq1:
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    params.fastq1 = word[1]
                    params.flag_fastq1 = True

                elif word[0] == "FASTQ2":
                    if len(word) < 2 or params.flag_fastq2:
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    params.fastq2 = word[1]
                    params.flag_fastq2 = True

                elif word[0] == "TARGETS":
                    if len(word) < 2 or params.flag_targfile:
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    params.targfile = word[1]
                    params.flag_targfile = True

                elif word[0] == "INDEX":
                    if len(word) < 2 or params.flag_indexpath:
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    params.indexpath = word[1]
                    params.flag_indexpath = True

                elif word[0] == "RESTFRAGS":
                    if len(word) < 2 or params.flag_restfragfile:
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    params.restfragfile = word[1]
                    params.flag_restfragfile = True

                elif word[0] == "ENZYME":
                    if len(word) < 2 or params.flag_enzyme:
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    params.enzyme = word[1].upper()
                    params.flag_enzyme = True

                elif word[0] == "EXCLUDE":
                    if len(word) < 2 or ( not word[1].isdigit() ):
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    params.exclusion = int(word[1])

                elif word[0] == "BIN":
                    if len(word)<3 :
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)   
                    if ( not word[1].isdigit() ) or ( not word[2].isdigit() ):
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    step = int(word[1])
                    wind = int(word[2])
                    if wind<step :
                        raise RuntimeError("Window size must be greater or equal to step size.\n")
                    params.binlist.append( [step,wind] )
                    params.biningcount += 1
                    
                elif word[0] == "INTERCHROM":
                    if len(word) < 2 :
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    if word[1].upper() == "TRUE":
                        params.save_inter = True

                elif word[0] == "NORMALIZE" or word[0] == "NORMALISE":
                    if len(word) < 2 :
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    if word[1].upper() == "TRUE":
                        params.normalize = True

                elif word[0] == "COMBINEMODE":
                    if len(word) < 2 :
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    if word[1].upper() == "TRUE":
                        params.combinemode = True

                elif word[0] == "COMBINECOUNT":
                    if len(word) < 2 or ( not word[1].isdigit() ):
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    params.combinecount = int(word[1])

                elif word[0] == "ALIGNMODE":
                    if len(word) < 2 and (not word[1] in ["CONSERVATIVE","RELAXED","CUSTOM"]):
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    if word[1]=="CUSTOM" and len(word) < 3:
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    params.alignmode = word[1]
                    if word[1]=="CUSTOM":
                        params.aligncustom =  " ".join(word[2:])

                elif word[0] == "PARALLEL":
                    if len(word) < 2 or ( not word[1].isdigit() ):
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    params.processors = int(word[1])

                elif word[0] == "DRYRUN":
                    if len(word) < 2 :
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    if word[1].upper() == "TRUE":
                        params.dryrun = True
                    
                elif word[0] == "TRIMADAPTERS":
                    if len(word) < 2 and (not word[1] in ["TRUE","FALSE","true","false"]):
                        raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
                    if word[1].upper() == "TRUE":
                        params.trim_adapters = True
                    else:
                        params.trim_adapters = False                      

                else:
                    raise RuntimeError("Error reading configuration file line:\n      %s\n"%line)
           
    # test that all required parameters are set
    if not (params.flag_fastq1 and params.flag_fastq2 and params.flag_targfile \
            and params.flag_indexpath and params.flag_restfragfile \
            and params.flag_enzyme):
        raise RuntimeError("Some required parameters are not set in configuration file.\n")
                    
    # check combinedcout is set if in combinemode
    if params.combinemode and params.alignmode=="CONSERVATIVE" and params.combinecount<=1:
        raise RuntimeError("Combine mode has been selected, but COMBINECOUNT is unset or set <=1.\n"%line)

    # set bowtie arguments
    if params.alignmode != "CUSTOM":
        if params.alignmode=="CONSERVATIVE":
            if params.combinemode: 
                params.bowtie_args = rs.bowtie_args_C_combine+" -m %i"%params.combinecount
            else:
                params.bowtie_args = rs.bowtie_args_C_nocombine
        if params.alignmode=="RELAXED":
                params.bowtie_args = rs.bowtie_args_R

    return params

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def check_depends():
    for i in rs.dependencies:
        if ( not which(i) ):
            raise RuntimeError("Error : Cannot find %s on PATH. ... Exiting.\n"%i)
    for i in rs.capCmap_extern:
        if ( not which( rs.capCmap_extern[i] ) ):
            raise RuntimeError("Error : capC-MAP is not properly installed - cannot find %s on PATH. ... Exiting.\n"%i)


def report2totreads(reportfile,target):
    """ Get a list of the total reads for each target """

    Tflag = False
    found = False

    total_reads = {}

    with open(reportfile) as rf:
        for line in rf:
            word = line.split()
            if Tflag and len(word)==5:
                if word[1] in target:
                    total_reads[word[1]]=int(word[4])
                else:
                    raise RuntimeError("Unknown target %s in report file.\n"%word[1])
            if len(word)>2 and word[1]=="target" and word[2]=="name":
                Tflag = True

    # Check all targets were found
    for T in target:
        if not T in total_reads:
            raise RuntimeError("Did not find totalreads for target %s in report file.\n"%T)

    return total_reads

def dryrun_report2totreads(target):
    """ Make a dummy list of the total reads for each target """
    
    total_reads = {}
        
    for T in target:
        total_reads[T] = 0

    return total_reads
        

def getCombinedTargets(target):
    """ Any targets to be combined should have _C* at the end of their name """
    combined_targets = {}

    tnames = target.keys()

    tocombine = [T for T in tnames if "_C" in T]

    for T in tocombine:
        front = T.split("_C",1)[0]+"_combined"
        if front in combined_targets:
            combined_targets[front].append(T)
        else:
            combined_targets[front] = [T]

    return combined_targets


def normalize_pileup(infile, outfile, total_reads):
    """ Read in a rawpileup file and generate a normalizedpileup file """
    
    fullpathinfile = os.path.realpath(infile)
    fullpathoutfile = os.path.realpath(outfile)

    ouf= open(fullpathoutfile,'w') 
    
    with open(fullpathinfile, 'r') as inf:
        for line in inf:
            if line.startswith("track"):
               ouf.write(line)
            else:
                line = line.split()
                normed = 1.0e6*float(line[3])/float(total_reads)
                ouf.write("%s\t%s\t%s\t%f\n"%(line[0],line[1],line[2],normed))
           
    ouf.close()


def run_pipeline(args):
    """ run the pipeline """

    experimentname = args.outdir[0]
    fulloutputdirpath = os.path.realpath(args.outdir[0])
    paramsfilefullpath = os.path.realpath(args.configfile[0])
    
    sys.stdout.write("\n####################################################\n")
    sys.stdout.write("## Running capC-MAP full analysis\n\n")
    

    #####################################################################################
    # Check directory is not there, and create it
    try:
        if ( os.path.isdir(fulloutputdirpath) ):
            raise RuntimeError("Directory %s already exists.\n"%fulloutputdirpath)
        if ( os.path.isfile(fulloutputdirpath) ):
            raise RuntimeError("%s already exists and is a file.\n"%fulloutputdirpath)
        os.makedirs(fulloutputdirpath)
        sys.stdout.write("Creating output directory %s ... Done.\n"%fulloutputdirpath)
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
    # Check required programs are on the PATH
    try:
        for i in rs.dependencies:
            if ( not which(i) ):
                raise RuntimeError("Error : Cannot find %s on PATH. ... Exiting.\n"%i)

        command = ["samtools", "--version-only"]
        version = subprocess.check_output(command)
        if ( not version.decode(sys.stdout.encoding).startswith("1.") ):
            raise RuntimeError("Error : Unsupported version of samtools %s ... Exiting.\n")

        for i in rs.capCmap_extern:
            if ( not which( rs.capCmap_extern[i] ) ):
                raise RuntimeError("Error : capC-MAP is not properly installed - cannot find %s on PATH. ... Exiting.\n"%i)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code
    except subprocess.CalledProcessError as e:
        sys.stdout.write("\nError : Unsupported version of samtools ... Exiting.\n")
        sys.exit(1) # exit with error code


    #####################################################################################
    # Check input files are there
    try:
        if ( not os.path.isfile(params.fastq1) ):
            raise RuntimeError("Cannot find file %s .\n"%params.fastq1)
        if ( not os.path.isfile(params.fastq2) ):
            raise RuntimeError("Cannot find file %s .\n"%params.fastq2)
        if ( params.fastq1 == params.fastq2 ):
            raise RuntimeError("Fastq files must be different.\n")
        if not os.path.isfile(params.targfile):
            raise RuntimeError("Cannot find file %s .\n"%params.targfile)
        if not os.path.isfile(params.restfragfile):
            raise RuntimeError("Cannot find file %s .\n"%params.restfragfile)
        fullpathfastq1 = os.path.abspath(params.fastq1)
        fullpathfastq2 = os.path.abspath(params.fastq2)
        fullpathtargfile = os.path.abspath(params.targfile)
        fullpathrestfragfile = os.path.abspath(params.restfragfile)
        if not os.path.isfile("%s.1.ebwt"%params.indexpath):
            raise RuntimeError("Cannot find bowtie index %s.\n"%params.indexpath)
        fullpathindex = os.path.abspath("%s.1.ebwt"%params.indexpath)[:-7]
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code
        

    #####################################################################################       
    # get targets list
    try:
        target = ParseTargets(fullpathtargfile)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################       
    # Get any targets which are to be combined
    try:
        if params.combinemode:
            combined_targets = getCombinedTargets(target)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################       
    # Get enzyme
    try:
        if ( params.enzyme in rs.enzyme_list ):
            cut_seq = rs.enzyme_list[params.enzyme]
        elif set(params.enzyme) <= set('ACTG') and len(params.enzyme) >=4:
            cut_seq = (params.enzyme,3)
        else:
            raise RuntimeError("Error: '%s' is not a valid restriction enzyme or cutting sequence.\n"%params.enzyme)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################       
    # Change into output directory, open some log files, and run the pipeline
    os.chdir(fulloutputdirpath)

    mainlogfile = open("capC-MAP.commands.log","w")
    mainerrorfile = open("capC-MAP.errors.log","w")

    #####################################################################################       
    # Run cutadapt to remove any adapter contamination.
    try:
        trimmedfile1 = "trimmed.1.fastq"
        trimmedfile2 = "trimmed.2.fastq"
        if params.trim_adapters:
            sys.stdout.write("\nRunning cutadapt to remove adapter sequence. \n"+
                             "\t see cutadapt.stdout.log and cutadapt.stderr.log for details.\n...")
            sys.stdout.flush()
            stdoutfile = open("cutadapt.stdout.log",'w')
            stderrfile = open("cutadapt.stderr.log",'w')
            command = ["cutadapt","-a",rs.adapter_FOR,"-A",rs.adapter_REV,
                       "-o",trimmedfile1,"-p",trimmedfile2,
                       fullpathfastq1,fullpathfastq2]
            mainlogfile.write(subprocess.list2cmdline(command)+"\n")
            if not params.dryrun:
                returncode = subprocess.call(command,
                                             stdout=stdoutfile,stderr=stderrfile)
                if ( not returncode==0 ):
                    raise RuntimeError("Error in cutadapt...Exiting...\n")
            sys.stdout.write(" \t Done.\n")
            stdoutfile.close()
            stderrfile.close()
        else:
            # skipping trimming step, just copy the files
            sys.stdout.write("\nSkipping step to remove adapter sequence. \n")
            copyfile(fullpathfastq1,trimmedfile1)
            copyfile(fullpathfastq2,trimmedfile2)

    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################    
    # Do in silico dpnII digestion of the resulting fastq files
    try:
        sys.stdout.write("\nRunning in silico digest for restriction enzyme cutting at %s\n"%cut_seq[0]+
                         "\t see digestfastq.stdout.log and digestfastq.stderr.log for details.\n...")
        sys.stdout.flush()
        stdoutfile = open("digestfastq.stdout.log",'w')
        stderrfile = open("digestfastq.stderr.log",'w')
        digestedfile = "digested.fastq"
        command = [rs.capCmap_extern["digestfastq"],"-1",trimmedfile1,"-2",trimmedfile2,
                   "-e",cut_seq[0],"-p","%s"%cut_seq[1],"-o",digestedfile]
        mainlogfile.write(subprocess.list2cmdline(command)+"\n")
        if not params.dryrun:
            returncode = subprocess.call(command,
                                         stdout=stdoutfile,stderr=stderrfile)
            if ( not returncode==0 ):
                raise RuntimeError("Error in digestion...Exiting...\n")
        sys.stdout.write(" \t Done.\n")
        stdoutfile.close()
        stderrfile.close()
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################  
    # Align with bowtie as single end reads
    try:
        sys.stdout.write("\nAligning to genome using bowtie\n...")
        sys.stdout.flush()
        stdoutfile = open("bowtie.stdout.log",'w')
        stderrfile = open("bowtie.stderr.log",'w')
        alignedsam = "aligned.sam"
        command = ("bowtie -p %s %s %s %s %s"%(params.processors,params.bowtie_args,fullpathindex,digestedfile,alignedsam)).split()
        mainlogfile.write(subprocess.list2cmdline(command)+"\n")
        if not params.dryrun:
            returncode = subprocess.call(command,
                                         stdout=stdoutfile,stderr=stderrfile)
            if ( not returncode==0 ):
                raise RuntimeError("Error in alignment...Exiting...\n")
        sys.stdout.write(" \t Done.\n")
        stdoutfile.close()
        stderrfile.close()
        # some useful information is put into the stderrfile
        #   -- copy this into stdout file
        with open("bowtie.stderr.log") as inf:
            with open("bowtie.stdout.log", "w") as ouf:
                for line in inf:
                    if not line.startswith("Warning:"):
                        ouf.write(line)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################  
    # Clean up some intermediate files at this point
    if not params.dryrun:
        os.remove(trimmedfile1)   # Remove the trimmed fastq files
        os.remove(trimmedfile2)
        os.remove(digestedfile)   # Remove the in silico digested fastq
        mainlogfile.write("# rm %s \n"%trimmedfile1)
        mainlogfile.write("# rm %s \n"%trimmedfile2)
        mainlogfile.write("# rm %s \n"%digestedfile)

    #####################################################################################  
    # Sort the output SAM file by read name
    try:
        sys.stdout.write("\nSorting SAM file by read name\n...")
        sys.stdout.flush()
        alignedbam = "aligned.bam"
        command = ["samtools", "view", "-S", "-b", "-@", "%i"%params.processors, "-o", alignedbam, alignedsam]
        mainlogfile.write(subprocess.list2cmdline(command)+"\n")
        if not params.dryrun:
            returncode = subprocess.call(command)
            if ( not returncode==0 ):
                raise RuntimeError("Error while converting sam -> bam ...Exiting...\n")
            # can now remove the aligned sam
            os.remove(alignedsam)

        sortedbam = "srt_aligned.bam"
        command = ["samtools", "sort", "-n", "-@", "%i"%params.processors, "-o", sortedbam, "-T", "tempsort", alignedbam]
        mainlogfile.write(subprocess.list2cmdline(command)+"\n")
        if not params.dryrun:
            returncode = subprocess.call(command)
            if ( not returncode==0 ):
                raise RuntimeError("Error while sorting bam ...Exiting...\n")
            # can now remove the aligned bam
            os.remove(alignedbam)

        sortedsam = "srt_aligned.sam"
        command = ["samtools", "view", "-h", "-@", "%i"%params.processors, "-o", sortedsam, sortedbam]
        mainlogfile.write(subprocess.list2cmdline(command)+"\n")
        if not params.dryrun:
            returncode = subprocess.call(command)
            if ( not returncode==0 ):
                raise RuntimeError("Error while converting bam -> sam ...Exiting...\n")

        # we now have a sorted BAM and a sorted SAM file
        sys.stdout.write(" \t Done.\n")
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################  
    # Run the main capture c analysis program
    try:
        sys.stdout.write("\nDoing the main capture-c analysis\n")
        sys.stdout.flush()
        pairsfile = "captured"
        stdoutfile = open("capCmain.stdout.log",'w')
        stderrfile = open("capCmain.stderr.log",'w')
        command = [rs.capCmap_extern["main"],"-r",fullpathrestfragfile,"-t",fullpathtargfile,
                   "-s",sortedsam,"-o",pairsfile,"-e","%s"%params.exclusion]
        # always save interchromosomal interactions
        command.extend(["-i"])
        mainlogfile.write(subprocess.list2cmdline(command)+"\n")
        if not params.dryrun:
            returncode = subprocess.call(command,stderr=stderrfile)
            if ( not returncode==0 ):
                raise RuntimeError("Error in main processing step...Exiting...\n")
        sys.stdout.write("... \t Done.\n")
        stdoutfile.close()
        stderrfile.close()
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################  
    # Pile-up pairs files into bedgraphs
    try:
        sys.stdout.write("\nPiling up valid pairs into a bedgraph file for each target\n")
        stdoutfile = open("pileup.stdout.log",'w')
        stderrfile = open("pileup.stderr.log",'w')
        reportfile = "%s_report.dat"%pairsfile
        if not params.dryrun:
            total_reads = report2totreads(reportfile,target)
        for T in target:
            infile = "%s_validpairs_%s.pairs"%(pairsfile,T)
            outfile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
            loc = "%s:%i-%i"%(target[T][0],target[T][1],target[T][2])
            sys.stdout.write("...target %s\n"%T)
            command = [rs.capCmap_extern["pair2bg"],"-i",infile,"-o",outfile,"-n",T,"-t",loc]
            mainlogfile.write(subprocess.list2cmdline(command)+"\n")
            if not params.dryrun:
                returncode = subprocess.call(command,
                                             stdout=stdoutfile,stderr=stderrfile)
                if ( not returncode==0 ):
                    raise RuntimeError("Error piling up pairs file %s .\n"%infile)
            if params.normalize:
                infile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
                outfile = "%s_normalizedpileup_%s.bdg"%(pairsfile,T)
                if not params.dryrun: 
                    normalize_pileup(infile, outfile, total_reads[T])

        # now do the same for interchromosomal
        if params.save_inter:
            sys.stdout.write("\nPiling up interchromosomal pairs into a bedgraph file for each target\n")
            reportfile = "%s_report.dat"%pairsfile
            if not params.dryrun:
                total_reads = report2totreads(reportfile,target)
            for T in target:
                infile = "%s_validinterchom_%s.pairs"%(pairsfile,T)
                outfile = "%s_rawpileup_interchom_%s.bdg"%(pairsfile,T)
                loc = "%s:%i-%i"%(target[T][0],target[T][1],target[T][2])
                sys.stdout.write("...target %s\n"%T)
                command = [rs.capCmap_extern["pair2bg"],"-i",infile,"-o",outfile,"-n",T,"-t",loc,"--interchrom"]
                mainlogfile.write(subprocess.list2cmdline(command)+"\n")
                if not params.dryrun:
                    returncode = subprocess.call(command,
                                                 stdout=stdoutfile,stderr=stderrfile)
                    if ( not returncode==0 ):
                        raise RuntimeError("Error piling up pairs file %s .\n"%infile)
                if params.normalize:
                    infile = "%s_rawpileup_interchom_%s.bdg"%(pairsfile,T)
                    outfile = "%s_normalizedpileup_interchom_%s.bdg"%(pairsfile,T)
                    if not params.dryrun: 
                        normalize_pileup(infile, outfile, total_reads[T])

        stdoutfile.close()
        stderrfile.close()
        sys.stdout.write("... \t Done.\n")

    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code  


    #####################################################################################  
    # Do windowed binning
    try:

        if params.biningcount > 0 :

            # generate chromosome sizes file from restriction frags list
            if not params.dryrun:
                sys.stdout.write("\nGenerating chrom.sizes file from restriction frags file\n...")
                get_chrom_sizes(fullpathrestfragfile,outfile="chrom.sizes")

            if  params.normalize:
                if not params.dryrun:
                    # binned profiles are to be normalized by reads per million
                    # need to get these values from report file
                    reportfile = "%s_report.dat"%pairsfile
                    total_reads = report2totreads(reportfile,target)
                else:
                    total_reads = dryrun_report2totreads(target)
                    
            stdoutfile = open("windowbin.stdout.log",'w')
            stderrfile = open("windowbin.stderr.log",'w')
            for B in params.binlist:
                step = B[0]
                window = B[1]
                sys.stdout.write("Generating binned interaction profiles with step %ibp and window %ibp ..."%(step,window))
                for T in target:
                    infile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
                    if params.normalize:
                        outfile = "%s_bin_%i_%i_RPM_%s.bdg"%(pairsfile,step,window,T)
                        command = [rs.capCmap_extern["pileup2binned"],"-i",infile,"-o",outfile,"-c","chrom.sizes",
                                   "-t",T,"-b","%i"%step,"%i"%window,
                                   "-n","%i"%total_reads[T]]
                    else:
                        outfile = "%s_bin_%i_%i_%s.bdg"%(pairsfile,step,window,T)
                        command = [rs.capCmap_extern["pileup2binned"],"-i",infile,"-o",outfile,"-c","chrom.sizes",
                                   "-t",T,"-b","%i"%step,"%i"%window]
                    mainlogfile.write(subprocess.list2cmdline(command)+"\n")
                    if not params.dryrun:
                        returncode = subprocess.call(command,stdout=stdoutfile,stderr=stderrfile)
                        if ( not returncode==0 ):
                            raise RuntimeError("Error binning pile-up files %s .\n"%infile)
                sys.stdout.write(" Done.\n")
            stdoutfile.close()
            stderrfile.close()

    except RuntimeError as e:
        sys.stdout.write("\nError : "+str(e))
        sys.exit(1) # exit with error code  
       

    #####################################################################################  
    # If any targets need to be combined, do that here
    try:
        if params.combinemode:
            # some of the targets need to be combined
            combined_targets = getCombinedTargets(target)

            stdoutfile = open("pileup.stdout.log",'a')
            stderrfile = open("pileup.stderr.log",'a')
            for T in combined_targets.keys():
                # first do pile-up
                sys.stdout.write("\nPiling up valid pairs into a bedgraph file for combined targets %s\n"%T)
                infiles = [ "%s_validpairs_%s.pairs"%(pairsfile,i) for i in combined_targets[T] ]
                locs = [ "%s:%i-%i"%(target[i][0],target[i][1],target[i][2]) for i in combined_targets[T] ]
                outfile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
                command = [rs.capCmap_extern["pair2bg"],"-o",outfile,"-n",T]
                for i in range(len(infiles)):
                    command.extend(["-i",infiles[i]])
                    command.extend(["-t",locs[i]])
                mainlogfile.write(subprocess.list2cmdline(command)+"\n")
                if not params.dryrun:
                    returncode = subprocess.call(command,stdout=stdoutfile,stderr=stderrfile)
                    if ( not returncode==0 ):
                        raise RuntimeError("Error piling up combined pairs file %s .\n"%infile)
                if params.normalize:
                    infile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
                    outfile = "%s_normalizedpileup_%s.bdg"%(pairsfile,T)
                    totalreadsT = 0
                    for i in combined_targets[T]:
                        totalreadsT = totalreadsT + total_reads[i]
                    if not params.dryrun: 
                        normalize_pileup(infile, outfile, totalreadsT)
                    
            # and do the same for interchrom
            if params.save_inter:
                sys.stdout.write("\nPiling up interchromosomal pairs into a bedgraph file for combined targets\n")
                reportfile = "%s_report.dat"%pairsfile
                total_reads = report2totreads(reportfile,target)
                for T in combined_targets.keys():
                    infiles = [ "%s_validinterchom_%s.pairs"%(pairsfile,i) for i in combined_targets[T] ]
                    locs = [ "%s:%i-%i"%(target[i][0],target[i][1],target[i][2]) for i in combined_targets[T] ]
                    outfile = "%s_rawpileup_interchom_%s.bdg"%(pairsfile,T)
                    command = [rs.capCmap_extern["pair2bg"],"-o",outfile,"-n",T,"--interchrom"]
                    for i in range(len(infiles)):
                        command.extend(["-i",infiles[i]])
                        command.extend(["-t",locs[i]])
                    mainlogfile.write(subprocess.list2cmdline(command)+"\n")
                    if not params.dryrun:
                        returncode = subprocess.call(command,stdout=stdoutfile,stderr=stderrfile)
                        if ( not returncode==0 ):
                            raise RuntimeError("Error piling up interchromosomal combined pairs file %s .\n"%infile)
                    if params.normalize:
                        infile = "%s_rawpileup_interchom_%s.bdg"%(pairsfile,T)
                        outfile = "%s_normalizedpileup_interchom_%s.bdg"%(pairsfile,T)
                        totalreadsT = 0
                        for i in combined_targets[T]:
                            totalreadsT = totalreadsT + total_reads[i]
                        if not params.dryrun: 
                            normalize_pileup(infile, outfile, totalreadsT)

            stdoutfile.close()
            stderrfile.close()
            sys.stdout.write("... \t Done.\n")

    except RuntimeError as e:
        sys.stdout.write("\nError : "+str(e))
        sys.exit(1) # exit with error code  


    #####################################################################################  
    # Now do the binning for any targets to be combined
    try:
        if params.combinemode:
            # then, if needed, do binning -- chrom sizes file will alread be there, and total_reads already obtained
            if params.biningcount > 0 :
                stdoutfile = open("windowbin.stdout.log",'a')
                stderrfile = open("windowbin.stderr.log",'a')
                for T in combined_targets.keys():
                    if  params.normalize:
                        totalreadsT = 0
                        for i in combined_targets[T]:
                            totalreadsT = totalreadsT + total_reads[i]
                    for B in params.binlist:
                        step = B[0]
                        window = B[1]
                        sys.stdout.write("\nGenerating binned interaction profiles for combined targets %s with step %ibp and window %ibp\n"%(T,step,window))
                        infile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
                        if params.normalize:
                            outfile = "%s_bin_%i_%i_RPM_%s.bdg"%(pairsfile,step,window,T)
                            command = [rs.capCmap_extern["pileup2binned"],"-i",infile,"-o",outfile,"-c","chrom.sizes",
                                       "-t",T,"-b","%i"%step,"%i"%window,
                                       "-n","%i"%totalreadsT]
                        else:
                            outfile = "%s_bin_%i_%i_%s.bdg"%(pairsfile,step,window,T)
                            command = [rs.capCmap_extern["pileup2binned"],"-i",infile,"-o",outfile,"-c","chrom.sizes",
                                       "-t",T,"-b","%i"%step,"%i"%window]
                        mainlogfile.write(subprocess.list2cmdline(command)+"\n")
                        if not params.dryrun:
                            returncode = subprocess.call(command,stdout=stdoutfile,stderr=stderrfile)
                            if ( not returncode==0 ):
                                raise RuntimeError("Error binning combined pairs file %s .\n"%infile)
                stdoutfile.close()
                stderrfile.close()
                sys.stdout.write("... \t Done.\n")

    except RuntimeError as e:
        sys.stdout.write("\nError : "+str(e))
        sys.exit(1) # exit with error code  



    #####################################################################################
    # remove SAM (we already made a BAM)
    try:
        if not params.dryrun:
            os.remove(sortedsam)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        # does not exit



    sys.stdout.write("\nAll done.\n")
    sys.stdout.write("\n####################################################\n")

    
    #####################################################################################       
    # close the log files

    mainlogfile.write("Done.\n")
    
    mainlogfile.close()
    mainerrorfile.close()
    
    return 0
