
import sys
import os.path
import numpy as np


def read_report(reportfile):
    """ 
    Function to read numbers from a capC-MAP report file.

    Argument : file name for 'captured_report.dat'
    Returns  : a dictionary containing the following keys with integer values :
                  reads             - total number of paired end reads
                  duplicates        - number of reads removed as duplicates
                  reads_after_rmdup - number of remaining reads
                  nonmapped         - number of reads where no fragments mapped
                  notarget          - number of reads without a target
                  noreporter        - number of reads without a reporter
                  multitarget       - number of reads discared because there were multiple targets
                  exclusion         - number of reads discared because they mapped to an exclusion zone
                  validinter        - number of valid interchromosomal interactions
                  validintra        - number of valid intrachromosomal interactions
                  totalvalid        - total number of valid reads
    """


    txt = {}
    txt['reads'] = "### total number of read sets"
    txt['duplicates'] = "###   number of duplicate sets removed"
    txt['nonmapped'] = "###   number of sets where no fragments mapped"
    txt['notarget'] = "###   number of sets with no targets fragments"
    txt['noreporter'] = "###   number of sets with no reporters"
    txt['multitarget'] = "###   number of sets with multiple targets"
    txt['exclusion'] = "###   number of interactions within exclusion region"
    txt['validinter'] = "###   interchromosomal "
    txt['validintra'] = "###   intrachromosomal"

    
    try:
        fullreportpath = os.path.realpath(os.path.expanduser(reportfile))
        if ( not os.path.isfile(fullreportpath) ):
            raise RuntimeError("Cannot find file %s .\n"%fullreportpath)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code
    except:
        sys.stdout.write("Error reading report file.\n")
        sys.exit(1) # exit with error code

    report_data = {}

    try:
        with open(fullreportpath) as rf:
            for line in rf:

                if line.startswith(txt['reads']):
                    word = line.split()
                    report_data['reads'] = int(word[-1])

                if line.startswith(txt['duplicates']):
                    word = line.split()
                    report_data['duplicates'] = int(word[-3])

                if line.startswith(txt['nonmapped']):
                    word = line.split()
                    report_data['nonmapped'] = int(word[-3])

                if line.startswith(txt['notarget']):
                    word = line.split()
                    report_data['notarget'] = int(word[-3])

                if line.startswith(txt['noreporter']):
                    word = line.split()
                    report_data['noreporter'] = int(word[-3])

                if line.startswith(txt['multitarget']):
                    word = line.split()
                    report_data['multitarget'] = int(word[-3])

                if line.startswith(txt['exclusion']):
                    word = line.split()
                    report_data['exclusion'] = int(word[-3])

                if line.startswith(txt['validinter']):
                    word = line.split()
                    report_data['validinter'] = int(word[-3])

                if line.startswith(txt['validintra']):
                    word = line.split()
                    report_data['validintra'] = int(word[-3])

        for key in txt:
            if not (key in report_data):
                raise RuntimeError("Error parsing report. Could not find %s \n"%key)

    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code

    report_data['reads_after_rmdup'] = report_data['reads'] - report_data['duplicates']
    report_data['totalvalid'] = report_data['validinter'] + report_data['validintra']

    return report_data




def read_interactioncounts(ptfile):
    """ 
    Function to read per-target stats from capC-MAP output.

    Argument : file name for 'captured_interactioncounts.dat'
    Returns  : a dictionary of numpy arrays containing the following keys :
                  targets           - list of target names
                  intrachrom        - array with per-target intrachromosomal read count
                  interchrom        - array with per-target interchromosomal read count
                  within5Mb         - per-target intrachromosomal reads within 5Mb of target
                  within1Mb         - per-target intrachromosomal reads within 1Mb of target
                  total             - array with per-target total read count
    """
    
    try:
        fullptfilepath = os.path.realpath(os.path.expanduser(ptfile))
        if ( not os.path.isfile(fullptfilepath) ):
            raise RuntimeError("Cannot find file %s .\n"%fullptfilepath)
    except RuntimeError as e:
        sys.stdout.write("Error : "+str(e))
        sys.exit(1) # exit with error code
    except:
        sys.stdout.write("Error reading file.\n")
        sys.exit(1) # exit with error code

    ## first count targets
    try:
        N = 0
        with open(fullptfilepath) as ptf:
            for line in ptf:
                 if not line.startswith('#'): 
                     N = N + 1
    except:
        sys.stdout.write("Error reading file (2).\n")
        sys.exit(1) # exit with error code

    # get arrays
    per_target = {}

    per_target['targets'] = []
    per_target['intrachrom'] = np.zeros(N,dtype=int)
    per_target['interchrom'] = np.zeros(N,dtype=int)
    per_target['within5Mb'] = np.zeros(N,dtype=int)
    per_target['within1Mb'] = np.zeros(N,dtype=int)
    per_target['total'] =  np.zeros(N,dtype=int)

    # read file
    try:
        i = 0
        with open(fullptfilepath) as ptf:
            for line in ptf:
                 if not line.startswith('#'): 
                     a = line.split()
                     per_target['targets'].append(a[0])

                     per_target['intrachrom'][i] = int(a[1])
                     per_target['total'][i] = int(a[2])
                     per_target['interchrom'][i] = int(a[3])
                     per_target['within5Mb'][i] = int(a[4])
                     per_target['within1Mb'][i] = int(a[5])

                     i = i + 1
    except:
        sys.stdout.write("Error reading file (3).\n")
        sys.exit(1) # exit with error code

    # done

    return per_target

