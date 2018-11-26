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


#include "parse_sam.h"
#include "bedfiles.h"
#include "samfragments.h"
#include "genome.h"
#include "targets.h"
#include "main_process.h"
#include "messages.h"

#include <map>
#include <string>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <algorithm>

using namespace CAPCMAIN_NS;


bool CAPCMAIN_NS::peakheader(std::ifstream& inf) {
  // peak at the next line in the SAM file
  // return TRUE if it is a header
  // return FLASE if it is not a header
  // throw an exception if reached the end of file

  int pp = inf.peek();

  if ( pp == EOF ) {
    throw std::runtime_error("samfile does not contain any entries");
  } else if ( pp == '@') {
    return 1;
  } else {
    return 0;
  }

}


bool CAPCMAIN_NS::peaksamline(std::ifstream& inf, const std::string& setname) {
  // look at the next line in the sam file
  // compare name with current name
  // return TRUE if they match
  // return FALSE if they match

  std::string line,
    nextsetname;
  std::stringstream sline;
  
  // peak at next line - return 0 if it is EOF
  if ( inf.peek() == EOF ) {
    return 0;
  }
    
  // save the file position
  std::streampos fileposition = inf.tellg();
  getline(inf,line);

  sline.clear(); sline.str(line);
  sline>>nextsetname;
  nextsetname = name2setname(nextsetname);

  // return to that file position
  inf.seekg( fileposition );
  
  if ( nextsetname==setname ) {
    // the next line if from the same read set
    return 1;
  } else {
    // the next line is a new read set
    return 0;
  }
  
}


  


void CAPCMAIN_NS::parse_sam_file(genome& gnm, const std::string& samfile,
				 const std::string& fname_out,
				 const parameters &params) {
  // Main function for parsing the sam file
  
  std::ifstream infsam;
  std::map< std::string,std::ofstream* > oufpairs,
    oufinter;
  std::string astring,
    line;
  std::vector<samfrag> current_sams;
  std::set<rest_fragment> current_frags;
  std::vector<rest_fragment> set_of_interchroms;
  target current_target;
  samfrag asam;
  int currentNtargs,
    nonAdjacent,
    n_mapped;

  // set up output files
  for ( genome::it_targs T=gnm.targets.begin(); T != gnm.targets.end(); ++T ) {
    // test file does not exist
    astring = fname_out + "_validpairs_" + T->name + ".pairs";
    infsam.open( astring.c_str() );
    if ( infsam.good() ) {
      throw std::runtime_error("file "+astring+" already exists (will not "
			       "overwrite).");
    }
    infsam.close();

    // open file for writting
    oufpairs[ T->name ] = new  std::ofstream( astring.c_str() );
  }
  
  // if save_inter is true, set up inter out files
  if ( params.save_inter ) {
    for ( genome::it_targs T=gnm.targets.begin(); T != gnm.targets.end(); ++T ) {
      // test file does not exist
      astring = fname_out + "_validinterchom_" + T->name + ".pairs";
      infsam.open( astring.c_str() );
      if ( infsam.good() ) {
	throw std::runtime_error("file "+astring+" already exists (will not "
				 "overwrite).");
      }
      infsam.close();
      
      // open file for writting
      oufinter[ T->name ] = new  std::ofstream( astring.c_str() );
    }
  }
  
  // Open SAM files
  infsam.open( samfile.c_str() );
  if ( ! infsam.good() ) {
    throw std::runtime_error("cannot open file "+samfile+".");
  }

  // throw away headers
  while ( peakheader(infsam) ) {
    getline(infsam,line);
  }

  
  // parse rest of sam file
  do {
    
    current_frags.clear();
    current_sams.clear();
    
    // read a group of fragments
    getline(infsam,line);
    
    current_sams.push_back( samfrag::samline2samfrag(line)  );
    gnm.count.total_read_frags++;
    while ( peaksamline(infsam,current_sams.back().setname) ) {
      // if the next line has the same name, get it
      getline(infsam,line);
      current_sams.push_back( samfrag::samline2samfrag(line)  );
      gnm.count.total_read_frags++;
    }

    gnm.count.total_read_sets++;


    // count number mapped, discard if none mapped
    n_mapped = 0;
    for (int i=0; i<current_sams.size(); i++) {
      if ( current_sams[i].chrom != "*" ) {
	n_mapped++;
      }
    }
    if ( n_mapped==0 ) {
      gnm.count.none_mapped++;
      continue;
    }
    
    // check for duplicates -- do this after discarding none mapped
    if ( gnm.is_duplicate(current_sams) ) {
      gnm.count.duplicates_removed++;
      continue;
    }

    
    // expand to restfrags
    for (int i=0; i<current_sams.size(); i++) {
      if ( current_sams[i].chrom != "*" ) {
	// expand mapped fragment to restriction fragment
	// if the same frag is inserted more than once, we only get one copy
	current_frags.insert( current_sams[i].expand_to_restfrag(gnm) );
      }
    }

    // discard none mapped
    //if ( current_frags.size()==0 ) {
    //  gnm.count.none_mapped++;
    //  continue;
    //}

    
    // count targets
    currentNtargs = 0;
    for (genome::it_rest_set F=current_frags.begin(); F!=current_frags.end();
	 ++F) {
      if ( F->is_target ) {
	currentNtargs++;
      }
    }
   
    // discard notargets
    if ( currentNtargs == 0 ) {
      gnm.count.no_targets++;
      continue;
    }

    // discard multitargets
    if ( currentNtargs > 1 ) {
      gnm.count.multiple_targets++;
      continue;
    }

    
    // get the current target
    for (genome::it_rest_set F =current_frags.begin(); F!=current_frags.end() ;
	 ) {
      if ( F->is_target ) {
	genome::it_targs T;
	T = gnm.targets.lower_bound( target(F->targetname) );
	if ( T->name != F->targetname ) {
	  throw std::runtime_error("something went wrong "
				   "identifying a target.");
	}
	current_target = *(T);
	current_frags.erase(F++);
	F = current_frags.end();
      } else {
	++F;
      }
    }
  
    // discard no reporters
    if ( current_frags.size() == 0 ) {
      gnm.count.no_reporters++;
      continue;
    }


    // discard due to exclusion around targets
    for (genome::it_targs T=gnm.targets.begin()  ; T!=gnm.targets.end()  ;
	 ++T) {
      for (genome::it_rest_set F=current_frags.begin(); F!=current_frags.end();
	   ) {
	if ( F->chrom==T->chrom &&
	     std::abs(0.5*(F->start+F->end) - 0.5*(T->start+T->end))
	     <=params.exclusion ) {
	  current_frags.erase(F++);
	} else {
	  ++F;
	}
      } 
    }
    if ( current_frags.size() == 0 ) {
      gnm.count.exclusion++;
      continue;
    }
    
    // count and then discard if only interchrom
    set_of_interchroms.clear();
    for (genome::it_rest_set F=current_frags.begin(); F!=current_frags.end();
	 ) {
      if ( F->chrom != current_target.chrom ) {
	set_of_interchroms.push_back(*F);
	current_frags.erase(F++);
	// this is the correct way to remove without invalidating the iterator
	// passes a copy of F, then increments the original F
	// before the current one is erased
      } else {
	++F;
      }
    }
    if ( current_frags.size() == 0 ) {
      gnm.count.onlyInter[current_target.name]++;
      gnm.count.total_interchrom++;
      if ( params.save_inter ) {
	// choose one at random by shuffling the vector
	std::random_shuffle ( set_of_interchroms.begin(), set_of_interchroms.end() );
	*(oufinter.find(current_target.name)->second)<<
	  set_of_interchroms[0].chrom<<"\t"<<
	  set_of_interchroms[0].start<<"\t"<<
	  set_of_interchroms[0].end<<//"\t"<<     // uncomment these lines to 
	  //current_sams.back().setname<<         // also output the read name
	  std::endl;
      }
      continue;
    }

    
    // discard multi nonadjacent reporters
    nonAdjacent = 0;
    if ( current_frags.size() > 1 ) {
      genome::it_rest_set F = current_frags.begin();
      int lastend = F->end;
      for ( ++F ; F != current_frags.end() ; ++F ) {
	if ( F->start != lastend ) {
	  nonAdjacent++;
	}
	lastend = F->end;
      }
    }
    if ( nonAdjacent > 0 ) { 
      gnm.count.multiple_reporters++;
      continue;
    }


    // it must be a valid read!!!!
    gnm.count.total_validPairs++;
    gnm.count.validPairs[current_target.name]++;

    // is it within 5Mb of the target? (lets use the start coords of the first)
    if (abs(current_frags.begin()->start-current_target.start)<=5e6) {
      gnm.count.within5Mb[current_target.name]++;
      // is it within 1Mb of the target?
      if (abs(current_frags.begin()->start-current_target.start)<=1e6) {
	gnm.count.within1Mb[current_target.name]++;
      }
    }


    // Choice here to give only the middle of any adjacent set of frags
    {
      genome::it_rest_set F = current_frags.begin();
      std::advance(F, std::floor( 0.5*double(current_frags.size()) ));
      if ( oufpairs.find(current_target.name)==oufpairs.end() ) {
	// check if found one : this should never happen
	throw std::runtime_error("target error");
      }
      *(oufpairs.find(current_target.name)->second)<<F->chrom<<"\t"
						   <<F->start<<"\t"
						   <<F->end//<<"\t" // uncomment these lines to 
	//<<current_sams.back().setname  // also output the read name
						   <<std::endl;
    }

    
  } while ( infsam.peek() != EOF );

  // Now close the files and tidy up
  for (std::map< std::string,std::ofstream* >::iterator it  = oufpairs.begin() ;
       it != oufpairs.end() ; ++it) {
    it->second->close();
    delete it->second;
  }

  if ( params.save_inter ) {
    for (std::map< std::string,std::ofstream* >::iterator it  = oufinter.begin() ;
	 it != oufinter.end() ; ++it) {
      it->second->close();
      delete it->second;
    }
  }


  // Output a message
  std::stringstream mymessage;
  mymessage<<"...Parsed "<<gnm.count.total_read_sets
	   <<" reads from SAM file "<<samfile;
  COMMON_NS::message( mymessage.str() );
  
}

