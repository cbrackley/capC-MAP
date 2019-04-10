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
// Program which reads a pairs file output from capC main process program
// Generates a pile-up of interactions and outputs a bedgraph file
// Does not normalize
//
///////////////////////////////////////////////////////////////////////////////

#include "pair2bg.h"
#include "postprocess.h"
#include "bedfiles.h"
#include "messages.h"

#include<string>
#include<iostream>
#include<fstream>
#include<stdexcept>
#include<cstdlib>
#include<map>

using namespace POSTPROCESS_NS;

int main(int argc, char *argv[]) {

  parameters params;

  // parse command line
  try {
    parse_pair2bg_command_line(argc,argv,params);
  } catch (const std::runtime_error& e) {
    std::cerr<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR converting pairs to bedgraph : An unknown error has occured."<<std::endl;
    return EXIT_FAILURE;
  }

  // get pile-up and output to bedgraph
  try {
      do_pileup(params);
  } catch (const std::runtime_error& e) {
    std::cerr<<"ERROR converting pairs to bedgraph : "<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR converting pairs to bedgraph : An unknown error has occured."<<std::endl;
    return EXIT_FAILURE;
  }

  
  // Done!
  return EXIT_SUCCESS;
  
}



void POSTPROCESS_NS::parse_pair2bg_command_line(const int &argc, char **argv,
				        parameters &params) {
  // parse the command line

 const std::string usage_message ="\nUsage :\n"
   "   capCpair2bg -i pairsfile -o bgfile -n targetname -t chr:start-end [--interchrom]\n"
   "\n"
   "   Required arguments :\n"
   "       -i  pairfile       is the input file name; can use this option more\n"
   "                          than once to combine multiple targets into one\n"
   "       -o  bgfile         is the file name for the output bedGraph\n"
   "       -n  targetname     is the name of the target\n"
   "       -t  chrom:start-end  is the genomic location of the target site; can\n"
   "                          use this option more than once if multiple pair\n"
   "                          files are specified.\n"
   "  Options  :\n"
   "       --interchrom       flag to specify interchromosomal interactions are present\n"
   "\n";
 
  unsigned short int infcount = 0,  // flags for required arguments
    loccount = 0,
    nameflag = 0,
    oflag= 0;
  
  unsigned short int flagichrom = 0;  // flags for required arguments

  int argi=1;
  
  while (argi < argc) {
    
    if ( std::string(argv[argi]) == "-i" ) {
      // input files
      if (!(argi+1 < argc) ) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.infile.push_back( std::string(argv[argi+1]) );
      infcount++;
      argi += 2;
      
    } else if ( std::string(argv[argi]) == "-o" ) {
      // output file
      if (!(argi+1 < argc) || oflag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.outfile = std::string(argv[argi+1]);
      oflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "-n" ) {
      // target name
      if (!(argi+1 < argc) || nameflag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.targetname = std::string(argv[argi+1]);
      nameflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "-t" ) {
      // target location
      if (!(argi+1 < argc)) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.location.push_back( std::string(argv[argi+1]) );
      loccount++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "--interchrom" ) {
      // are interchromosomal interactions present
      if (!(argi < argc)) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      flagichrom++;
      argi += 1;

    } else if ( std::string(argv[argi]) == "--version" ) {
      // version -- overrides all other option
      COMMON_NS::print_version();
      std::exit(EXIT_SUCCESS);
      
    } else {
      throw std::runtime_error("Unknown option "+std::string(argv[argi])+"\n"+usage_message);
    }
    
  }

  // Check required parameters are there
  if ( infcount==0 || loccount==0 || nameflag!=1 || oflag!=1 ) {
      throw std::runtime_error(usage_message);
  }

  // Set if interchromosomal interactions present
  if ( flagichrom > 0 ) {
    params.interchromflag=1;
  } else {
    params.interchromflag=0;
  }

  // Check locations and targets have same count
  if ( loccount != infcount ) {
    throw std::runtime_error("Number of input files does not equal number of locations");
  }

  // Get list of chroms from locations
  for (int i=0;i<params.location.size();i++) {
    std::size_t colon = params.location[i].find(":"),
      dash = params.location[i].find("-");
    if ( colon == std::string::npos ||
	 params.location[i].find(":",colon+1) != std::string::npos ||
	 dash == std::string::npos ||
	 params.location[i].find(":",dash+1) != std::string::npos ||
	dash<colon ) {
      // there should be one : and one - in that order
      throw std::runtime_error("Incorrect target location format "+
			       params.location[i]);
    }

    params.chrom.insert( params.location[i].substr(0,colon) );
    
  }
  

}


void POSTPROCESS_NS::do_pileup(const parameters &params) {
  
  std::map<std::string,
	   std::map<CAPCMAIN_NS::rest_fragment,unsigned int> > pileup;
  std::map<std::string,
	   std::map<CAPCMAIN_NS::rest_fragment,unsigned int> >::iterator itpileup;
  std::map<CAPCMAIN_NS::rest_fragment,unsigned int>::iterator itfrags;
  
  std::ifstream inf;
  std::ofstream ouf;
  std::string line;
  
  // have a pilup for each chromosome
  for (parameters::itchrom C=params.chrom.begin() ; C!=params.chrom.end() ;
       ++C )  {
    pileup[ *C ] = std::map<CAPCMAIN_NS::rest_fragment,unsigned int>() ;
  }
  
  // test output file
  inf.open( params.outfile.c_str() );
  if ( inf.good() ) {
    throw std::runtime_error("file "+params.outfile+" already exists.");
  }
  inf.close();
  
  // test all input files
  for (int F=0; F<params.infile.size(); F++) {
    inf.open( params.infile[F].c_str() );
    if ( ! inf.good() ) {
      throw std::runtime_error("cannot open file "+params.infile[F]);
    }
    inf.close();
  }
  
  
  // loop round input files
  for (int F=0; F<params.infile.size(); F++) {
    
    // open file
    inf.open( params.infile[F].c_str() );
    
    while ( getline(inf,line) )  {
      // load line by line
      
      itpileup = pileup.find( line.substr(0,line.find("\t")) );
      if ( itpileup != pileup.end() ) {
	itpileup->second.operator[]( CAPCMAIN_NS::bed_feature::line2bed_feature(line) )++;
      } else if ( params.interchromflag ) { 
	// Interchromosomal interaction to new chromosome
	pileup[ line.substr(0,line.find("\t")) ] = std::map<CAPCMAIN_NS::rest_fragment,unsigned int>() ;
	itpileup = pileup.find( line.substr(0,line.find("\t")) );
	itpileup->second.operator[]( CAPCMAIN_NS::bed_feature::line2bed_feature(line) )++;
      } else {
	throw std::runtime_error("interchromosomal interaction in pairs file");
      }
      
    }
    
    inf.close();
    
  }
  
  
  // output bedgraph
  ouf.open( params.outfile.c_str() );
  
  // output a header line
  ouf<<"track type=bedGraph name=\""<<params.targetname<<"\" description=\""
     <<"Target: "<<params.targetname;
  if ( params.location.size() > 1 ) {
    ouf<<" at locations";
  } else {
    ouf<<" at location";
  }
  for (int L=0;L<params.location.size();L++) {
    ouf<<" "<<params.location[L];
  }
  ouf<<"\""<<std::endl;
    
  // now output pilups
  for (itpileup=pileup.begin() ; itpileup!=pileup.end() ; ++itpileup) {
    
    for (itfrags=itpileup->second.begin() ; itfrags!=itpileup->second.end() ;
	 ++itfrags) {
      
      ouf<<itfrags->first.chrom<<"\t"
	 <<itfrags->first.start<<"\t"
	 <<itfrags->first.end<<"\t"
	 <<itfrags->second<<std::endl;
      
    }
    
  }
  
  ouf.close();
  
  
}
