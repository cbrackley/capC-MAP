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
// Program which reads a bed file of genomic windows
// and converts these to restriction enzyme fragments
//
///////////////////////////////////////////////////////////////////////////////

#include "location2fragment.h"
#include "genome.h"
#include "bedfiles.h"
#include "messages.h"

#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <set>
#include <vector>

using namespace LOC2FRAG_NS;

int main(int argc, char *argv[]) {

  parameters params;
  std::vector<CAPCMAIN_NS::rest_fragment> locations;
  std::vector<CAPCMAIN_NS::rest_fragment>  restfrags;

  std::ifstream inf;
  std::ofstream ouf;

  CAPCMAIN_NS::genome gnm;

  // parse command line
  try {
    parse_loc2frag_command_line(argc,argv,params);
  } catch (const std::runtime_error& e) {
    std::cerr<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR converting genomic location to restriction fragment : An unknown error has occured."<<std::endl;
    return EXIT_FAILURE;
  }

  // initialize genome -- load restriction enzyme fragments
  try {
    gnm.load_rest_frags(params.fragfile);
  } catch (const std::runtime_error& e) {
    std::cerr<<"ERROR : "<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR : An unknown error has occured."
	     <<std::endl;
    return EXIT_FAILURE;
  }


  // get list of locations
  try {
    if (params.bedfileflag>0) {
      locations = load_locations(params);
    } else {
      locations.push_back( CAPCMAIN_NS::rest_fragment(params.location_chrom,
						   params.location_start,
						   params.location_end) );
    }
  }catch (const std::runtime_error& e) {
    std::cerr<<"ERROR : "<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR : An unknown error has occured."
	     <<std::endl;
    return EXIT_FAILURE;
  }

  // convert locations to fragments
  try {
    restfrags = locs2frags(gnm,locations);
  }catch (const std::runtime_error& e) {
    std::cerr<<"ERROR : "<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR : An unknown error has occured."
	     <<std::endl;
    return EXIT_FAILURE;
  }


  // output
  try {
    if (params.outfileflag) {
      inf.open( params.outfile.c_str() );
      if ( inf.good() ) {
	throw std::runtime_error("file "+params.outfile+" already exists (will not "
				 "overwrite).");
      }
      inf.close();
      ouf.open( params.outfile.c_str() );
      for (int i=0;i<restfrags.size();i++) {
	ouf<<restfrags[i].chrom<<"\t"
		 <<restfrags[i].start<<"\t"
		 <<restfrags[i].end<<"\t"
		 <<restfrags[i].restofline<<std::endl;
      }
      ouf.close();
    } else {
      for (int i=0;i<restfrags.size();i++) {
	std::cout<<restfrags[i].chrom<<"\t"
		 <<restfrags[i].start<<"\t"
		 <<restfrags[i].end<<"\t"
		 <<restfrags[i].restofline<<std::endl;
      }
    }
  }catch (const std::runtime_error& e) {
    std::cerr<<"ERROR in main processing stage : "<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR in main processing stage : An unknown error has occured."
	     <<std::endl;
    return EXIT_FAILURE;
  }

  // Done!
  return EXIT_SUCCESS;

}


void LOC2FRAG_NS::parse_loc2frag_command_line(const int &argc,
					      char **argv,
					      parameters &params) {
  // parse the command line

 const std::string usage_message ="\nUsage :\n"
   "   capClocation2fragment -r restfragfile -o outfile [-i inputfile | -l location]\n"
   "\n"
   "   Required arguments :\n"
   "       -r  restfragfile    filename for bed file containing the list of restriction fragments\n"
   "       -o  outfile         filename for output bed file (if not present output with be sent to stdout)\n"
   "\n"
   "   Options : exactly one the following optional arguments must be present\n"
   "       -i  inputfile       filename for bed file containing genomic locations\n"
   "       -l  location        is a single genomic location in format chr1:1234-5678 \n"

   "\n";

 // flags for required arguments
 params.rfileflag = 0;  

 // flags for optional arguments
 params.bedfileflag = 0;
 params.locflag = 0;
 params.outfileflag = 0;

 std::string bp;

  int argi=1;
  while (argi < argc) {

    if ( std::string(argv[argi]) == "-r" ) {
      // input files
      if (!(argi+1 < argc) || params.rfileflag!=0 ) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.fragfile = std::string(argv[argi+1]);
      params.rfileflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "-i" ) {
      // input files
      if (!(argi+1 < argc) || params.bedfileflag!=0 ) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.bedfile = std::string(argv[argi+1]);
      params.bedfileflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "-o" ) {
      // input files
      if (!(argi+1 < argc) || params.outfileflag!=0 ) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.outfile = std::string(argv[argi+1]);
      params.outfileflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "-l" ) {
      // input files
      if (!(argi+1 < argc) ) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.location_string = std::string(argv[argi+1]);
      params.locflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "--version" ) {
      // version -- overrides all other option
      COMMON_NS::print_version();
      std::exit(EXIT_SUCCESS);
      
    } else {
      throw std::runtime_error("Unknown option "+std::string(argv[argi])+"\n"+usage_message);
    }

  }

  // Check required parameters are there
  if ( params.rfileflag!=1 ) {
      throw std::runtime_error(usage_message);
  }

  // Check one and only of the options is there
  if ( params.locflag + params.bedfileflag != 1 ) {
    throw std::runtime_error("Must enter exactly one of -l and -i options.\n"+usage_message);
  }

  // If -l option is true, parse the location string
  if ( params.locflag >0 ) {
    params.location_chrom = 
      params.location_string.substr(0, 
				    params.location_string.find_first_of(":") );
    
    bp = params.location_string.
      substr(params.location_string.find_first_of(":")+1,
	     params.location_string.find_first_of("-")
	     -params.location_string.find_first_of(":")-1);

    if ( bp.find_first_not_of("0123456789") != std::string::npos ) {
      throw std::runtime_error("Error parsing command line : genomic location is not correctly formed.");
    }
    std::istringstream(bp) >> params.location_start;
    
    bp = params.location_string.
      substr(1+params.location_string.find_first_of("-"));

    if ( bp.find_first_not_of("0123456789") != std::string::npos ) {
      throw std::runtime_error("Error parsing command line : genomic location is not correctly formed.");
    }
    std::istringstream(bp) >> params.location_end;

    if (params.location_end<params.location_start) {
      throw std::runtime_error("Error parsing command line : genomic location is not correctly formed.");

    }

  }


}


std::vector<CAPCMAIN_NS::rest_fragment
	 > LOC2FRAG_NS::load_locations(const parameters &params) {
  // load the list of locations from a bed file

  std::vector<CAPCMAIN_NS::rest_fragment> locs;
  std::ifstream inf;
  std::string line;
  std::stringstream sline;
  std::string bpss,
    bpee;
  CAPCMAIN_NS::rest_fragment newloc;

  inf.open( params.bedfile.c_str() );
  if ( ! inf.good() ) {
    throw std::runtime_error("cannot open file "+params.bedfile);
  }

  while ( getline(inf,line) )  {
    if ( line.compare(0,5,"track")!=0 && 
	 line.compare(0,7,"browser")!=0 ) { // ignore track definition lines
      sline.clear(); sline.str( line );
      sline>>newloc.chrom>>bpss>>bpee;

      if ( bpss.find_first_not_of("0123456789") != std::string::npos ||
	   bpee.find_first_not_of("0123456789") != std::string::npos ) {
	throw std::runtime_error("unreadable entry in bedgraph file");
      }
      newloc.start=atoi( bpss.c_str() );
      newloc.end=atoi( bpee.c_str() );
      getline(sline,newloc.restofline);

      locs.push_back(newloc);
    }
  }

  inf.close();

  return locs;

}

std::vector<CAPCMAIN_NS::rest_fragment
	 > LOC2FRAG_NS::locs2frags(const CAPCMAIN_NS::genome &gnm, 
				   const std::vector<CAPCMAIN_NS::rest_fragment> &loc) {

  std::vector<CAPCMAIN_NS::rest_fragment> frags;
  CAPCMAIN_NS::rest_fragment a_location;
  CAPCMAIN_NS::genome::const_it_rest_map achromosome;
  CAPCMAIN_NS::genome::const_it_rest_set afragment;

  // loop round the locations
  for (int i=0;i<loc.size();i++) {
    
    // set the location to be in the middle of the input location
    a_location = loc[i];
    a_location.start = int(0.5*(a_location.start+a_location.end));
    a_location.end = a_location.start+1;

    achromosome = gnm.restriction_fragments.find(a_location.chrom);
    if ( achromosome == gnm.restriction_fragments.end() ) {
      throw std::runtime_error("cannot find chromosome "+a_location.chrom+
			       " in the restriction fragments file");
    }

    // find the fragment
    afragment = --( achromosome->second.upper_bound(a_location) );
    if ( !(a_location.start >= afragment->start &&
	   a_location.end <= afragment->end ) ) {
      
      throw std::runtime_error("something has gone wrong in deterimning"
			       " the fragment.");
    } else {
      frags.push_back( *afragment );
      frags.back().restofline=loc[i].restofline;
    }
    
  }

  return frags;

}

