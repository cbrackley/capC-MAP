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
import subprocess

from . import resources as rs

from .pipeline import Parameters,ParseConfigFile,check_depends,ParseTargets,getCombinedTargets,get_chrom_sizes,report2totreads,normalize_pileup

def do_postprocess(args):
    """" Do some binning / smoothing / normalization """

    experimentname = args.outdir[0]
    fulloutputdirpath = os.path.realpath(args.outdir[0])
    paramsfilefullpath = os.path.realpath(args.configfile[0])
    
    sys.stdout.write("\n####################################################\n")
    sys.stdout.write("## Running capC-MAP post-processing of raw pile-ups.\n\n")
  

    #####################################################################################
    # Check directory is there
    try:
        if ( not os.path.isdir(fulloutputdirpath) ):
            raise RuntimeError("Cannot find directory %s .\n"%fulloutputdirpath)
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

    # postprocess does not work in dryrun mode
    if params.dryrun:
        sys.stdout.write("DRYRUN flag is set in configuration file.\n"
                         +"capC-MAP postprocess will not do anything.\n"
                         +"Set DRYRUN FALSE and try again.\n")
        sys.exit(1) # exit without error

    #####################################################################################
    # Check required programs are on the PATH
    try:
        check_depends()
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code


    #####################################################################################
    # Check input files are there
    try:
        if not os.path.isfile(params.targfile):
            raise RuntimeError("Cannot find file %s .\n"%params.targfile)
        fullpathtargfile = os.path.abspath(params.targfile)
        if not os.path.isfile(params.restfragfile):
            raise RuntimeError("Cannot find file %s .\n"%params.restfragfile)
        fullpathrestfragfile = os.path.abspath(params.restfragfile)
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
    # Change into output directory
    os.chdir(fulloutputdirpath)
    pairsfile = "captured"


    #####################################################################################  
    # Do pile-ups -- always attempt this
    try:
        sys.stdout.write("\nChecking if pile-ups are present. If so they will not be overwritten ...")
        stdoutfile = open("postprocess.pileup.stdout.log",'a')
        stderrfile = open("postprocess.pileup.stderr.log",'a')
        stdoutfile.write("\n# ####### Appended to log during postprocess run:\n")
        stderrfile.write("\n# ####### Appended to log during postprocess run:\n")
        reportfile = "%s_report.dat"%pairsfile
        total_reads = report2totreads(reportfile,target)
        for T in target:
            infile = "%s_validpairs_%s.pairs"%(pairsfile,T)
            outfile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
            loc = "%s:%i-%i"%(target[T][0],target[T][1],target[T][2])
            command = [rs.capCmap_extern["pair2bg"],"-i",infile,"-o",outfile,"-n",T,"-t",loc]
            if not os.path.isfile(outfile):
                sys.stdout.write("\n   Pile-up for target %s not found. Generating file %s\n"%(T,outfile))
                returncode = subprocess.call(command,
                                             stdout=stdoutfile,stderr=stderrfile)
                if ( not returncode==0 ):
                    raise RuntimeError("Error piling up pairs file %s .\n"%infile)
            if params.normalize:
                infile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
                outfile = "%s_normalizedpileup_%s.bdg"%(pairsfile,T)
                if not os.path.isfile(outfile):
                    sys.stdout.write("   Normalized pile-up for target %s not found. Generating file %s\n"%(T,outfile))
                    normalize_pileup(infile, outfile, total_reads[T])

        # now do the same for interchromosomal
        if params.save_inter:
            sys.stdout.write("\nChecking if interchromosomal pile-ups are present. If so they will not be overwritten ...")
            reportfile = "%s_report.dat"%pairsfile
            total_reads = report2totreads(reportfile,target)
            for T in target:
                infile = "%s_validinterchom_%s.pairs"%(pairsfile,T)
                outfile = "%s_rawpileup_interchom_%s.bdg"%(pairsfile,T)
                loc = "%s:%i-%i"%(target[T][0],target[T][1],target[T][2])
                command = [rs.capCmap_extern["pair2bg"],"-i",infile,"-o",outfile,"-n",T,"-t",loc,"--interchrom"]
                if not os.path.isfile(outfile):
                    sys.stdout.write("\n   Interchrom pile-up for target %s not found. Generating file %s\n"%(T,outfile))
                    returncode = subprocess.call(command,
                                                 stdout=stdoutfile,stderr=stderrfile)
                    if ( not returncode==0 ):
                        raise RuntimeError("Error piling up pairs file %s .\n"%infile)
                if params.normalize:
                    infile = "%s_rawpileup_interchom_%s.bdg"%(pairsfile,T)
                    outfile = "%s_normalizedpileup_interchom_%s.bdg"%(pairsfile,T)
                    if not os.path.isfile(outfile):
                        sys.stdout.write("   Normalized interchrom pile-up for target %s not found. Generating file %s\n"%(T,outfile))
                        normalize_pileup(infile, outfile, total_reads[T])

        stdoutfile.close()
        stderrfile.close()
        sys.stdout.write(" Done.\n\n")

    except RuntimeError as e:
        sys.stdout.write("\nError : "+str(e))
        sys.exit(1) # exit with error code  


    #####################################################################################  
    # Do windowed binning
    try:

        if params.biningcount > 0 :

            # generate chromosome sizes file from restriction frags list
            if not os.path.isfile("chrom.sizes"):
                sys.stdout.write("\nGenerating chrom.sizes file from restriction frags file\n...")
                get_chrom_sizes(fullpathrestfragfile,outfile="chrom.sizes")

            if  params.normalize:
                # binned profiles are to be normalized by reads per million
                # need to get these values from report file
                reportfile = "%s_report.dat"%pairsfile
                total_reads = report2totreads(reportfile,target)

            stdoutfile = open("postprocess.windowbin.stdout.log",'a')
            stderrfile = open("postprocess.windowbin.stderr.log",'a')
            stdoutfile.write("\n# ####### Appended to log during postprocess run:\n")
            stderrfile.write("\n# ####### Appended to log during postprocess run:\n")
            for B in params.binlist:
                step = B[0]
                window = B[1]
                sys.stdout.write("Generating binned interaction profiles with step %ibp and window %ibp ...\n"%(step,window))
                sys.stdout.write("... Existing files will not be overwritten ...")
                for T in target:
                    infile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
                    if os.path.isfile(infile):
                        if params.normalize:
                            outfile = "%s_bin_%i_%i_RPM_%s.bdg"%(pairsfile,step,window,T)
                            command = [rs.capCmap_extern["pileup2binned"],"-i",infile,"-o",outfile,"-c","chrom.sizes",
                                       "-t",T,"-b","%i"%step,"%i"%window,
                                       "-n","%i"%total_reads[T]]
                        else:
                            outfile = "%s_bin_%i_%i_%s.bdg"%(pairsfile,step,window,T)
                            command = [rs.capCmap_extern["pileup2binned"],"-i",infile,"-o",outfile,"-c","chrom.sizes",
                                       "-t",T,"-b","%i"%step,"%i"%window]
                        if not os.path.isfile(outfile):
                            returncode = subprocess.call(command,stdout=stdoutfile,stderr=stderrfile)
                            if ( not returncode==0 ):
                                raise RuntimeError("Error binning pile-up files %s .\n"%infile)
                    else:
                        sys.stdout.write("   Warning : cannot find file %s. Skipping. No output will be generated.\n"%infile)
                sys.stdout.write(" Done.\n")
            stdoutfile.close()
            stderrfile.close()

    except RuntimeError as e:
        sys.stdout.write("\nError : "+str(e))
        sys.exit(1) # exit with error code  


    #####################################################################################  
    # If any targets need to be combined, do the pile ups
    try:
        if params.combinemode:
            combined_targets = getCombinedTargets(target)

            stdoutfile = open("postprocess.pileup.stdout.log",'a')
            stderrfile = open("postprocess.pileup.stderr.log",'a')
            for T in combined_targets.keys():
                sys.stdout.write("\nChecking if pile-ups are present for combined targets %s ..."%T)
                infiles = [ "%s_validpairs_%s.pairs"%(pairsfile,i) for i in combined_targets[T] ]
                locs = [ "%s:%i-%i"%(target[i][0],target[i][1],target[i][2]) for i in combined_targets[T] ]
                outfile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
                command = [rs.capCmap_extern["pair2bg"],"-o",outfile,"-n",T]
                for i in range(len(infiles)):
                    command.extend(["-i",infiles[i]])
                    command.extend(["-t",locs[i]])
                if not os.path.isfile(outfile):
                    returncode = subprocess.call(command,stdout=stdoutfile,stderr=stderrfile)
                    if ( not returncode==0 ):
                        raise RuntimeError("Error piling up combined pairs file %s .\n"%infile)
                
                # now normalize
                outfile = "%s_normalizedpileup_%s.bdg"%(pairsfile,T)
                if params.normalize and ( not os.path.isfile(outfile) ):
                    infile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
                    totalreadsT = 0
                    for i in combined_targets[T]:
                        totalreadsT = totalreadsT + total_reads[i]
                    normalize_pileup(infile, outfile, totalreadsT)

            sys.stdout.write(" Done.\n")

            # including interchrom
            if params.save_inter:
                reportfile = "%s_report.dat"%pairsfile
                total_reads = report2totreads(reportfile,target)
                for T in combined_targets.keys():
                    sys.stdout.write("\nChecking if interchromosomal pile-ups are present for combined targets %s ..."%T)
                    infiles = [ "%s_validinterchom_%s.pairs"%(pairsfile,i) for i in combined_targets[T] ]
                    locs = [ "%s:%i-%i"%(target[i][0],target[i][1],target[i][2]) for i in combined_targets[T] ]
                    outfile = "%s_rawpileup_interchom_%s.bdg"%(pairsfile,T)
                    command = [rs.capCmap_extern["pair2bg"],"-o",outfile,"-n",T,"--interchrom"]
                    for i in range(len(infiles)):
                        command.extend(["-i",infiles[i]])
                        command.extend(["-t",locs[i]])
                    if not os.path.isfile(outfile):
                        returncode = subprocess.call(command,stdout=stdoutfile,stderr=stderrfile)
                        if ( not returncode==0 ):
                            raise RuntimeError("Error piling up combined pairs file %s .\n"%infile)

                    # now normalize
                    outfile = "%s_normalizedpileup_interchom_%s.bdg"%(pairsfile,T)
                    if params.normalize and ( not os.path.isfile(outfile) ):
                        infile = "%s_rawpileup_interchom_%s.bdg"%(pairsfile,T)
                        totalreadsT = 0
                        for i in combined_targets[T]:
                            totalreadsT = totalreadsT + total_reads[i]
                        normalize_pileup(infile, outfile, totalreadsT)

            stdoutfile.close()
            stderrfile.close()
            sys.stdout.write(" Done.\n")

    except RuntimeError as e:
        sys.stdout.write("\nError : "+str(e))
        sys.exit(1) # exit with error code  


    #####################################################################################  
    # Does any binning need to be done for combined
    try:
        if params.combinemode and params.biningcount > 0:
            combined_targets = getCombinedTargets(target)
            stdoutfile = open("postprocess.windowbin.stdout.log",'a')
            stderrfile = open("postprocess.windowbin.stderr.log",'a')
            for T in combined_targets.keys():
                if  params.normalize:
                    totalreadsT = 0
                    for i in combined_targets[T]:
                        totalreadsT = totalreadsT + total_reads[i]
                for B in params.binlist:
                    step = B[0]
                    window = B[1]
                    sys.stdout.write("\nGenerating binned interaction profiles for combined targets %s with step %ibp and window %ibp\n"%(T,step,window))
                    sys.stdout.write("... Existing files will not be overwritten ...")
                    infile = "%s_rawpileup_%s.bdg"%(pairsfile,T)
                    if os.path.isfile(infile):
                        if params.normalize:
                            outfile = "%s_bin_%i_%i_RPM_%s.bdg"%(pairsfile,step,window,T)
                            command = [rs.capCmap_extern["pileup2binned"],"-i",infile,"-o",outfile,"-c","chrom.sizes",
                                       "-t",T,"-b","%i"%step,"%i"%window,
                                       "-n","%i"%totalreadsT]
                        else:
                            outfile = "%s_bin_%i_%i_%s.bdg"%(pairsfile,step,window,T)
                            command = [rs.capCmap_extern["pileup2binned"],"-i",infile,"-o",outfile,"-c","chrom.sizes",
                                       "-t",T,"-b","%i"%step,"%i"%window]
                        if not os.path.isfile(outfile):
                            returncode = subprocess.call(command,stdout=stdoutfile,stderr=stderrfile)
                            if ( not returncode==0 ):
                                raise RuntimeError("Error binning combined pairs file %s .\n"%infile)
                    else:
                        sys.stdout.write("   Warning : cannot find file %s. Skipping. No output will be generated.\n"%infile)
                    sys.stdout.write(" Done.\n")
            stdoutfile.close()
            stderrfile.close()
    except RuntimeError as e:
        sys.stdout.write("\nError : "+str(e))
        sys.exit(1) # exit with error code  


    return 0

