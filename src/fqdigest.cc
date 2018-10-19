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
// Program which reads a pair of fastq files and performs an in silico 
// restriction enzyme digestion
//
///////////////////////////////////////////////////////////////////////////////


#include "fqdigest.h"
#include "fastq.h"
#include "messages.h"

#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

#include <vector>
#include <fstream>

using namespace DIGEST_NS;


parameters::parameters() {
  // constructor for parameters structure
  // set default values here
  long_flag = 0;
}



int main(int argc, char *argv[]) {

  filenames fname;
  parameters params;
  counters count;

  // parse command line
  try {
    parse_command_line(argc,argv,fname,params);
  } catch (const std::runtime_error& e) {
    std::cerr<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR doing in silico digestion : An unknown error has occured."<<std::endl;
    return EXIT_FAILURE;
  }

  // do digestion
  try {
    if (params.long_flag) {
      do_digestion_long(fname,params,count);
    } else {
      do_digestion(fname,params,count);
    }
  } catch (const std::runtime_error& e) {
    std::cerr<<"ERROR doing in silico digestion : "<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR doing in silico digestion : An unknown error has occured."<<std::endl;
    return EXIT_FAILURE;
  }


  // output a log file
  try {
    write_log(fname,params,count);
  } catch (const std::runtime_error& e) {
    std::cerr<<e.what()<<std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr<<"ERROR doing in silico digestion : An unknown error has occured."<<std::endl;
    return EXIT_FAILURE;
  }

  // Done!
  return EXIT_SUCCESS;
  
}


void DIGEST_NS::parse_command_line(const int &argc,  char **argv, filenames &fname, parameters &params) {
  // parse the command line

  const std::string usage_message ="\nUsage :\n"
    "   capCdigestfastq -1 first_fq -2 second_fq -o output_fq -e SEQ -p X [--long]\n"
    "\n"
    "   Required arguments :\n"
    "       -1  first_fq    is the first of the pair of fastq files\n"
    "       -2  second_fq   is the second of the pair of fastq files\n"
    "       -o  output_fq   is the name of the output fastq file\n"
    "       -e  SEQ         is the sequence of the restriction enzyme\n"
    "                       must be characters ACGT only\n"
    "       -p  X           is the bp position within SEQ where the\n"
    "                       cut will occur (first base is 1; Xth base\n"
    "                       will be the start of the right hand\n"
    "                       fragment)\n"
    "\n"
    "   Options :\n"
    "       --long          option switches on 'long' mode, where only the\n"
    "                       longest of the restriction fragments in each of\n"
    "                       the pairs is kept\n"
    "\n";
  
 
  unsigned short int narg = 5,       // number of required arguments
    fq1flag = 0,                     // flags for required arguments
    fq2flag = 0,
    outflag = 0,
    enzymeflag = 0,
    posflag = 0;
  unsigned short int longflag = 0;   // flags for optional arguments

  std::string position;
  
  int argi=1;

  while (argi < argc) {

    if ( std::string(argv[argi]) == "-1" ) {
      // first fastq file
      if (!(argi+1 < argc) || fq1flag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      fname.fq1 = std::string(argv[argi+1]);
      fq1flag++;
      argi += 2;
      
    } else if ( std::string(argv[argi]) == "-2" ) {
      // second fastq file
      if (!(argi+1 < argc) || fq2flag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      fname.fq2 = std::string(argv[argi+1]);
      fq2flag++;
      argi += 2;
      
    } else if ( std::string(argv[argi]) == "-o" ) {
      // output fastq file
      if (!(argi+1 < argc || outflag!=0)) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      fname.outfile = std::string(argv[argi+1]);
      outflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "-e" ) {
      // second fastq file
      if (!(argi+1 < argc) || enzymeflag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      params.enzyme = std::string(argv[argi+1]);
      enzymeflag++;
      argi += 2;

    } else if ( std::string(argv[argi]) == "-p" ) {
      // second fastq file
      if (!(argi+1 < argc) || posflag!=0) {
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      position = std::string(argv[argi+1]);
      posflag++;
      argi += 2;
      
    } else if ( std::string(argv[argi]) == "--long" ) {
      // switch for long mode
      if (!(argi < argc)) { // allow multiple instances of '--long'
	throw std::runtime_error("Error parsing command line.\n"+usage_message);
      }
      longflag++;
      argi += 1;

    } else {
      throw std::runtime_error("Unknown option "+std::string(argv[argi])+"\n"+usage_message);
    }

  }


  // Check required parameters are there
  if ( fq1flag + fq2flag + outflag + enzymeflag + posflag != narg ) {
      throw std::runtime_error(usage_message);
  }

  // check enzyme sequence
  if ( params.enzyme.find_first_not_of("ATCG") != std::string::npos ) {
	throw std::runtime_error("Error parsing command line : enzyme sequence must be characters ACGT only.");
  }
  
  // check position
  if ( position.find_first_not_of("0123456789") != std::string::npos ) {
	throw std::runtime_error("Error parsing command line : enzyme cut point must be a positive integer");
  }
  std::istringstream(position) >> params.cutpoint;
  if ( params.cutpoint<1 || params.cutpoint>params.enzyme.size() ) {
    throw std::runtime_error("Error parsing command line : enzyme cut point must be a positive integer between 1 and the sequence length");
  }

  // Chech optional parameters
  if ( longflag == 1 ) {
    params.long_flag = 1;
  } else {
    params.long_flag = 0;
  }
  
}


void DIGEST_NS::do_digestion(const filenames &fname, const parameters &params, counters &count) {
  // Actually do the digestion

  std::ifstream inf[2];
  std::ofstream ouf;

  fastq_entry fastq[2];

  std::vector<fastq_entry> digested;
  std::vector<unsigned int> cutpoints;
  std::size_t cp;
  
  // open output file, checking it doesn't exist
  inf[0].open( fname.outfile.c_str() );
  if ( !inf[0].fail() ) {
    throw std::runtime_error("File "+fname.outfile+" already exists.");
  }
  inf[0].close();
  ouf.open( fname.outfile.c_str() );

  // Check and open fastq files
  inf[0].open( fname.fq1.c_str() );
  inf[1].open( fname.fq2.c_str() );
  if ( inf[0].fail() ) {
    throw std::runtime_error("Cannot open file "+fname.fq1+".");
  }
  if ( inf[1].fail() ) {
    throw std::runtime_error("Cannot open file "+fname.fq2+".");
  }

  // give some output messages
  COMMON_NS::message("\nDoing an in silico digest of the fastq files "
		     +fname.fq1+" and "+fname.fq2+" ...\n"
		     "... writting to file "+fname.outfile+" ...");
  
  
  // parse FASTQ
  while ( fastq[0].load_new(inf[0]) && fastq[1].load_new(inf[1]) ) {

    count.pairs++;
    
    if ( fastq[0].name != fastq[1].name ) {
      throw std::runtime_error("Names in fastq files do not match.\n"
			       "All entries must have a pair.\n"
			       "Fastq pair files must be in the same order.\n");
    }

    digested.clear();   // clear the list of digested fragments
    
    for (int i=0;i<2;i++) {
      
      cutpoints.clear();  // clear the list of cut points
      cp = fastq[i].seq.find(params.enzyme);

      // always start at the start
      cutpoints.push_back(0);
      
      while ( cp != std::string::npos ) {
	cutpoints.push_back(cp+params.cutpoint-1);
	cp = fastq[i].seq.find(params.enzyme,cp+params.enzyme.size());
      }

      // always end at the end
      cutpoints.push_back( fastq[i].seq.size() );

      // this will keep start and end points which are just bits of the
      // cutting sequence
      
      for (int cc=1;cc<cutpoints.size();cc++) {
	// add to the list of digested fragments
	digested.push_back( fastq[i].extract(cutpoints[cc-1],cutpoints[cc]-1) );
      }
      
    } // end the loop round i

    // Now output all the fragments
    for (int i=0;i<digested.size();i++) {
      ouf<<digested[i].name<<"DIGEST"<<i+1<<std::endl
	 <<digested[i].seq<<std::endl
	 <<digested[i].thirdline<<std::endl
	 <<digested[i].phred<<std::endl;      
    }

    // count the number of cuts made
    count.cuts += digested.size()-2;
    count.outfrags += digested.size();
    
  }
  
  // close the file
  inf[0].close();
  inf[1].close();
  ouf.close();

  // Output message
  COMMON_NS::message("... Done.");
  
}




void DIGEST_NS::do_digestion_long(const filenames &fname, const parameters &params, counters &count) {
  // Actually do the digestion --long mode

  std::ifstream inf[2];
  std::ofstream ouf[2];

  fastq_entry fastq[2];

  std::vector<fastq_entry> digested;
  std::vector<unsigned int> cutpoints;  
  std::size_t cp;

  std::string outfiles[2];

  // get output file names
  if ( fname.outfile.find_last_of(".") != std::string::npos ) {
    outfiles[0] = fname.outfile.substr(0,fname.outfile.find_last_of("."))+
      "_1"+fname.outfile.substr(fname.outfile.find_last_of("."));
    outfiles[1] = fname.outfile.substr(0,fname.outfile.find_last_of("."))+
      "_2"+fname.outfile.substr(fname.outfile.find_last_of("."));
  } else {
   outfiles[0] = fname.outfile+"_1";
   outfiles[1] = fname.outfile+"_2";
  }

  
  // open output file, checking it doesn't exist
  for (int i=0;i<2;i++) {
    inf[i].open( outfiles[i].c_str() );
    if ( !inf[i].fail() ) {
      throw std::runtime_error("File "+outfiles[i]+" already exists.");
    }
    inf[i].close();
    ouf[i].open( outfiles[i].c_str() );
  }

  
  // Check and open fastq files
  inf[0].open( fname.fq1.c_str() );
  inf[1].open( fname.fq2.c_str() );
  if ( inf[0].fail() ) {
    throw std::runtime_error("Cannot open file "+fname.fq1+".");
  }
  if ( inf[1].fail() ) {
    throw std::runtime_error("Cannot open file "+fname.fq2+".");
  }


  // give some output messages
  COMMON_NS::message("\n'--long' option has been used.\n"
		     "Only the longest fragment from each read pair is "
		     "kept.\nDigested fragments will be written to two "
		     "files.\n");
  COMMON_NS::warning_message("This output is not suitable for downstream "
			     "analysis with capCmain which uses the unmapped \n"
			     " short fragments when identifying duplicates.\n");
  
  COMMON_NS::message("\nDoing an in silico digest of the fastq files "
		     +fname.fq1+" and "+fname.fq2+" ...\n"
		     "... writting to output files "+outfiles[0]+" and "
		     +outfiles[1]+" ...");
  
  
  // parse FASTQ
  while ( fastq[0].load_new(inf[0]) && fastq[1].load_new(inf[1]) ) {

    count.pairs++;
    
    if ( fastq[0].name != fastq[1].name ) {
      throw std::runtime_error("Names in fastq files do not match.\n"
			       "All entries must have a pair.\n"
			       "Fastq pair files must be in the same order.\n");
    }

    
    for (int i=0;i<2;i++) {
      
      cutpoints.clear();  // clear cutpoints
      digested.clear();   // clear digested
      cp = fastq[i].seq.find(params.enzyme);

      // always start at the start
      cutpoints.push_back(0);
      
      while ( cp != std::string::npos ) {
	cutpoints.push_back(cp+params.cutpoint-1);
	cp = fastq[i].seq.find(params.enzyme,cp+params.enzyme.size());
      }

      // always end at the end
      cutpoints.push_back( fastq[i].seq.size() );

      // this will keep start and end points which are just bits of the
      // cutting sequence
      
      for (int cc=1;cc<cutpoints.size();cc++) {
	// add to the list of digested fragments
	digested.push_back( fastq[i].extract(cutpoints[cc-1],cutpoints[cc]-1) );
      }

      // find the largest one
      int LARGEST = 0;
      for (int ll=1;ll<digested.size();ll++) {
	if ( digested[ll].seq.size() > digested[LARGEST].seq.size() ) {
	  LARGEST = ll;
	}
      }

      // now output
      ouf[i]<<digested[LARGEST].name<<"/"<<i+1<<std::endl
	    <<digested[LARGEST].seq<<std::endl
	    <<digested[LARGEST].thirdline<<std::endl
	    <<digested[LARGEST].phred<<std::endl; 

      count.cuts += digested.size()-1;
      
    } // end of i loop

    count.outfrags += 2;
    
  };
  

  // close the files
  for (int i=0;i<2;i++) {
    ouf[i].close();
    inf[i].close();
  }

  // Output message
  COMMON_NS::message("... Done.");
  
}


void DIGEST_NS::write_log(const filenames &fname, const parameters &params,
			  const counters &count) {
  // Write a log file with counts
  
  std::string filename;
  std::ifstream inf;
  std::ofstream ouf;

  // test and open file
  filename = "digestlog_"+fname.outfile+".log";
  inf.open( filename.c_str() );
  if ( !inf.fail() ) {
    COMMON_NS::warning_message("Log file "+filename+" already exists. Will not overwrite.");
    inf.close();
    return;
  }
  ouf.open( filename.c_str() );

  // write log
  ouf<<"###################################################"<<std::endl;
  ouf<<"# Performed in silico digest on fastq files."<<std::endl;
  ouf<<"# INPUT FILES  : "<<fname.fq1<<std::endl
     <<"#                "<<fname.fq2<<std::endl;
  if (params.long_flag) {
    ouf<<"# OUTPUT FILES : "
       <<fname.outfile.substr(0,fname.outfile.find_last_of("."))<<"_1"
       <<fname.outfile.substr(fname.outfile.find_last_of("."))
       <<std::endl<<"#                "
       <<fname.outfile.substr(0,fname.outfile.find_last_of("."))<<"_2"
       <<fname.outfile.substr(fname.outfile.find_last_of("."))
       <<std::endl;
  } else {
    ouf<<"# OUTPUT FILE  : "<<fname.outfile<<std::endl;
  }
  ouf<<"#"<<std::endl;
  ouf<<"# Number of input read pairs       :    "<<count.pairs<<std::endl;
  ouf<<"# Number of cut sites              :    "<<count.cuts<<std::endl;
  ouf<<"# Number of output fragments       :    "<<count.outfrags<<std::endl;
  ouf<<"###################################################"<<std::endl;

  // send a message
  COMMON_NS::message("Log file written to "+filename+"\n");
  
  // close the file
  ouf.close();

}
