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
// Program which reads a bedGraph of the pileup of interactions
// Normalize to reads per million genome wide
// Put into bins with sliding windows
//
///////////////////////////////////////////////////////////////////////////////

#include "pileup2binned.h"
#include "postprocess.h"
#include "binprofile.h"

#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <string>
#include <map>
#include <sstream>

using namespace POSTPROCESS_NS;

int main(int argc, char *argv[]) {

  parameters params;

  // parse command line
  try {
    parse_pileup2binned_command_line(argc,argv,params);
  } catch (const std::runtime_error& e) {
    std::cerr<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR converting pairs to bedgraph : An unknown error has occured."<<std::endl;
    return EXIT_FAILURE;
  }

  if (params.binflag == 0) {
    std::cerr<<"Sorry : normalization without binning has not yet been implemented"<<std::endl;
    return EXIT_FAILURE;
  }

  // load chrom sizes
  std::map<std::string, long int> chromsizes;
  try {
    chromsizes = load_chrom_sizes(params.chromsizesname);
  } catch (const std::runtime_error& e) {
    std::cerr<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR converting pairs to bedgraph : An unknown error has occured."<<std::endl;
    return EXIT_FAILURE;
  }
  
  
  // do binning
  binprofile binned(params.binwidth,params.windowwidth);
  try {
    do_binning(params,binned,chromsizes);
  } catch (const std::runtime_error& e) {
    std::cerr<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR converting pairs to bedgraph : An unknown error has occured."<<std::endl;
    return EXIT_FAILURE;
  }
  
  // Done!
  return EXIT_SUCCESS;

}



void POSTPROCESS_NS::parse_pileup2binned_command_line(const int &argc,
						      char **argv,
						      parameters &params) {
  // parse the command line

 const std::string usage_message ="\nUsage :\n"
   "   capCpileup2binned -i pileupfile -o outfile -c chromsizes -t target [-b bin wind] [-n totalreads]\n"
   "\n"
   "   Required arguments :\n"
   "       -i  pileupfile  is the input pile-up file name\n"
   "       -o  outfile     is the file name for the output bedGraph\n"
   "       -c  chromsizes  is the file name for the list of chromosome sizes\n"
   "       -t  target      is the name of the target\n"
   "\n"
   "   Options : one or more of the following optional arguments must be present\n"
   "       -b  bin wind    pile-up will be up into sliding window bins with step\n"
   "                       size of 'bin' and window width of 'wind'\n"
   "       -n  totalreads  pile-up will be normalized to reads per million\n"
   "                       genome wide; requires total number of reads (available\n"
   "                       from capC main process report file; includes both inter\n"
   "                       and intra chromosomal).\n"  
   "\n";

  unsigned short int inflag = 0,  // flags for required arguments
    nameflag = 0,
    chromsizesflag = 0,
    oflag= 0;

  unsigned short int binflag = 0, // flags for optional arguments
    normflag =0;

  std::string binwidth,
    windowwidth,
    total_reads;
  
  int argi=1;
  
  while (argi < argc) {
    
    if ( std::string(argv[argi]) == "-i" ) {
      // input files
      if (!(argi+1 < argc) || inflag!=0 ) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.infile.push_back( std::string(argv[argi+1]) );
      inflag++;
      argi += 2;
      
    } else if ( std::string(argv[argi]) == "-o" ) {
      // input files
      if (!(argi+1 < argc) || oflag!=0 ) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.outfile = std::string(argv[argi+1]);
      oflag++;
      argi += 2;
      
    } else if ( std::string(argv[argi]) == "-t" ) {
      // input files
      if (!(argi+1 < argc)) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.targetname = std::string(argv[argi+1]);
      nameflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "-c" ) {
      // input files
      if (!(argi+1 < argc)) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.chromsizesname = std::string(argv[argi+1]);
      chromsizesflag++;
      argi += 2;
      
    } else if ( std::string(argv[argi]) == "-b" ) {
      // input files
      if (!(argi+1 < argc) || binflag!=0 ) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      binwidth = std::string(argv[argi+1]);
      windowwidth = std::string(argv[argi+2]);
      binflag++;
      argi += 3;
      
    } else if ( std::string(argv[argi]) == "-n" ) {
      // input files
      if (!(argi+1 < argc) || normflag!=0 ) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      total_reads = std::string(argv[argi+1]);
      normflag++;
      argi += 2;
      
    } else {
      throw std::runtime_error("Unknown option "+std::string(argv[argi])+"\n"+usage_message);
    }

  }

  // Check required parameters are there
  if ( inflag!=1 || nameflag!=1 || oflag!=1 || chromsizesflag==0 ) {
      throw std::runtime_error(usage_message);
  }

  // Check one of the options is there
  if ( binflag + normflag < 1 ) {
    throw std::runtime_error(usage_message);
  }

  params.binflag = binflag;
  params.normflag = normflag;

  // Convert bin params to ints
  if ( binflag > 0 ) {
    if ( binwidth.find_first_not_of("0123456789") != std::string::npos ) {
	throw std::runtime_error("Error parsing command line : bin step width "
				 "must be positive integer");
    }
    std::istringstream(binwidth) >> params.binwidth;
    if ( windowwidth.find_first_not_of("0123456789") != std::string::npos ) {
	throw std::runtime_error("Error parsing command line : window width "
				 "must be positive integer.");
    }
    std::istringstream(windowwidth) >> params.windowwidth;

    if ( params.binwidth > params.windowwidth ) {
	throw std::runtime_error("Error parsing command line : window width "
				 "must be greater or equal to bin step width.");
    }
    
  }
  if ( normflag > 0 ) {
    if ( total_reads.find_first_not_of("0123456789") != std::string::npos ) {
	throw std::runtime_error("Error parsing command line : totalreads "
				 "must be positive integer");
    }
    std::istringstream(total_reads) >> params.total_reads;
  }

  
}
