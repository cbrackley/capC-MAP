/* capC-MAP - a software package for analysis of Capture-C data
 * Copyright (C) 2018, University of Edinburgh
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Written by Chris Brackley <C.Brackley@ed.ac.uk>
 *
 */

///////////////////////////////////////////////////////////////////////////////
//
// Main capture C analysis program
// 
// Loads a group of fragments.
// Removes duplicates.
// Identifies original restriction enzyme fragments.
// Keeps only reads with one target and one reporter.
// Splits intra and interchromosomal interactions.
// Removes target-reporter pairs which are too close to another target.
// Counts different kinds of read, and outputs statistics.
//
///////////////////////////////////////////////////////////////////////////////


// Reminder - here use bed file convention that bp start counting from 0
//          - SAM files use a different convention, and start from 1


#include "main_process.h"
#include "genome.h"
#include "parse_sam.h"
#include "bedfiles.h"
#include "messages.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <stdexcept>

using namespace CAPCMAIN_NS;


parameters::parameters() {
  // constructor for parameters structure
  // set default values here
  exclusion = 500;
}



int main(int argc, char *argv[]) {
  
  genome gnm;
  
  filenames fname;
  parameters params;

  
  // parse command line
  try {
    parse_command_line(argc,argv,fname,params);
  } catch (const std::runtime_error& e) {
    std::cerr<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR in main processing stage : An unknown error has occured."
	     <<std::endl;
    return EXIT_FAILURE;
  }


  // initialize genome -- load restriction enzyme fragments
  try {
    gnm.load_rest_frags(fname.restfrags);
  } catch (const std::runtime_error& e) {
    std::cerr<<"ERROR in main processing stage : "<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR in main processing stage : An unknown error has occured."
	     <<std::endl;
    return EXIT_FAILURE;
  }


  
  // initialize targets and set up counters
  try {
    gnm.load_targets(fname.targets);
    gnm.count.setup();
  } catch (const std::runtime_error& e) {
    std::cerr<<"ERROR in main processing stage : "<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR in main processing stage : An unknown error has occured."
	     <<std::endl;
    return EXIT_FAILURE;
  }



  // parse SAM file
  try {
    parse_sam_file(gnm,fname.samfile,fname.outfile,params);
  } catch (const std::runtime_error& e) {
    std::cerr<<"ERROR in main processing stage : "<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR in main processing stage : An unknown error has occured."
	     <<std::endl;
    return EXIT_FAILURE;
  }

 
  
  // output counts and report
  try {
    gnm.count.output_interchrom(fname.outfile+"_interactioncounts.dat");
    gnm.count.output_report(fname.outfile+"_report.dat");
  } catch (const std::runtime_error& e) {
    std::cerr<<"ERROR outputing report : "<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR in main processing stage : An unknown error has occured."
	     <<std::endl;
    return EXIT_FAILURE;
  }


  // Done!
  return EXIT_SUCCESS;
  
}


void CAPCMAIN_NS::parse_command_line(const int &argc,  char **argv,
				     filenames &fname, parameters &params) {
  // parse the command line

  const std::string usage_message ="\nUsage :\n"
    "   capCmain -r frag_file -t targ_file -s sam_file -o name [-e N] [-i]\n"
    "\n"
    "   Required arguments :\n"
    "       -r  frag_file   is a bed file of restriction enzyme fragments genome wide\n"
    "       -t  targ_file   is a bed file of capture targets\n"
    "       -s  sam_file    is a SAM file containing groups of aligned\n"
    "                       digested fragments, sorted by name\n"
    "       -o  name        is the first part of the output file name\n"
    "\n"
    "   Options :\n"
    "       -e N            exclusion zone; reporter fragments mapping within N bp of\n"
    "                       a target fragment are discarder. Default N=500.\n"
    "       -i              save interchromosomal. If present, interchomosomal\n"
    "                       interactions will be saved as well as counted.\n"
    "\n";
  
  std::string exclusion;
  unsigned short int narg = 4,       // number of required arguments
    resflag = 0,                     // flags for required arguments
    targflag = 0,
    samflag = 0,
    outflag = 0,
    saveIflag = 0;
  unsigned short int   excflag = 0;  // flags for optional arguments
  
  int argi=1;

  while (argi < argc) {

    if ( std::string(argv[argi]) == "-r" ) {
      // restriction fragments file
      if (!(argi+1 < argc) || resflag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      fname.restfrags = std::string(argv[argi+1]);
      resflag++;
      argi += 2;
      
    } else if ( std::string(argv[argi]) == "-t" ) {
      // targets file
      if (!(argi+1 < argc) || targflag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      fname.targets = std::string(argv[argi+1]);
      targflag++;
      argi += 2;
	  
    } else if ( std::string(argv[argi]) == "-s" ) {
      // sam file
      if (!(argi+1 < argc) || samflag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      fname.samfile = std::string(argv[argi+1]);
      samflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "-o" ) {
      // output file
      if (!(argi+1 < argc) || outflag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      fname.outfile = std::string(argv[argi+1]);
      outflag++;
      argi += 2;
      
    } else if ( std::string(argv[argi]) == "-e" ) {
      // exclusion zone
      if (!(argi+1 < argc) || excflag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      exclusion = std::string(argv[argi+1]);
      excflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "-i" ) {
      // exclusion zone
      if (!(argi < argc) || saveIflag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      saveIflag++;
      argi ++;
      
    } else {
      throw std::runtime_error("Unknown option "+std::string(argv[argi])+"\n"
			       +usage_message);
    }
    
  }

  // Check required parameters are there
  if ( resflag + targflag + samflag + outflag != narg ) {
      throw std::runtime_error(usage_message);
  }
    
  // Chech optional parameters
  if ( excflag == 1 ) {

    if ( exclusion.find_first_not_of("0123456789") != std::string::npos ) {
	throw std::runtime_error("Error parsing command line : option -e "
				 "requires positive integer");
    }
    std::istringstream(exclusion) >> params.exclusion;

    if ( params.exclusion < 1 ) {
	throw std::runtime_error("Error parsing command line : option -e "
				 "requires integer >0");
    }
    
  }

  if ( saveIflag == 1 ) {
    params.save_inter = true;
  } else {
    params.save_inter = false;
  }
  
}
