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


#include "genome.h"
#include "messages.h"
#include "bedfiles.h"
#include "samfragments.h"
#include "targets.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>

using namespace CAPCMAIN_NS;


void genome::load_targets(const std::string &filename) {
  // Load a bed file of targets

  std::ifstream inf;
  std::string line;
  target newtarget;
  int id;
  it_rest_set pntr_rest_frag;
  it_rest_map pntr_chrom;
  std::pair<it_targs,bool> pair_ptr_target;
  std::pair<std::set<std::string>::iterator,bool> pair_ptr_string;

  if ( are_restfrags_loaded==0 ) {
    throw std::runtime_error("attempted to load targets without first loading restriction fragments");
  }
  
  inf.open( filename.c_str() );
  if ( ! inf.good() ) {
    throw std::runtime_error("cannot open file "+filename);
  }

  id = 0;
  while ( getline(inf,line) )  {
    newtarget = bed_feature::line2bed_feature(line);
    newtarget.target_id = id;
    id++;

    // check name is not empty
    if ( newtarget.name == "" ) {
      throw std::runtime_error("an entry in the targets list has no name - all"
			       " target fragments must have a unique name");
    }
      
    pntr_chrom = restriction_fragments.find( newtarget.chrom );
    if ( pntr_chrom==restriction_fragments.end() ) {
      // chromsome is not in fragments list
      throw std::runtime_error("unknown target chromosome "+newtarget.chrom);
    }
   
    pntr_rest_frag = pntr_chrom->second.lower_bound( rest_fragment(newtarget) );
    if (pntr_rest_frag->chrom != newtarget.chrom
	|| pntr_rest_frag->start != newtarget.start
	|| pntr_rest_frag->end != newtarget.end) {
      // this target does not match the restiction fragment
      throw std::runtime_error("target "+newtarget.name+" is not a restriction"
			       "enzyme fragment");
    }

    if ( pntr_rest_frag->is_target==1 ) {
      // a target at this location already existed
      COMMON_NS::warning_message("more than one target defined at same location - "
		      "ignoring the dupliacte entry '"+newtarget.name+"'.");    
    } else {
      // all seems in order - try to add to the targets list
      pair_ptr_target = targets.insert( newtarget );
      if ( pair_ptr_target.second==false ) {
	// this target name was alreasy in use
	throw std::runtime_error("multiple targets with same name - all target"
				 "fragments must have a unique name");
      }
      pntr_rest_frag->is_target = 1;
      pntr_rest_frag->targetname = newtarget.name;
    }

      
  }
    
  inf.close();

  // count the targets
  Ntargets = targets.size();
  std::stringstream mymessage;
  mymessage<<"...Loaded list of "<<Ntargets<<" targets from file "<<filename;
  COMMON_NS::message( mymessage.str() );
  
  are_targets_loaded=1;
  
}


void genome::load_rest_frags(const std::string &filename) {
  // Load a bed file of restriction fragments

  std::ifstream inf;
  std::string line,
    thischrom;
  rest_fragment newfrag;
  
  restriction_fragments.clear();

  inf.open( filename.c_str() );
  if ( ! inf.good() ) {
    throw std::runtime_error("cannot open file "+filename);
  }

  // read fragments from bed file
  while ( getline(inf,line) )  {
    newfrag = bed_feature::line2bed_feature(line);
    restriction_fragments[newfrag.chrom].insert( newfrag );    
  }

  inf.close();

  Nchrom = restriction_fragments.size();

  int laststart,
    lastend,
    nonabuting;

  // test restriction fragments do not overlap
  for (it_rest_map chr=restriction_fragments.begin();
       chr != restriction_fragments.end(); ++chr) {
    it_rest_set fr = chr->second.begin();
    thischrom = fr->chrom;
    laststart=fr->start;
    lastend=fr->end;

    if (laststart != 0) {
      COMMON_NS::warning_message("first restriction fragment in chromosome "
		      +thischrom+" does not start at 0");
    }
    
    nonabuting = 0;
    ++fr; // move to next fragment in this chromosome
    for ( ; fr != chr->second.end(); ++fr) {
      // test that start is not less then lastend
      if ( fr->start < lastend ) {
	throw std::runtime_error("restriction fragments are overlapping in "
				 +thischrom);
      }
      if ( fr->start != lastend ) {
	nonabuting++;
      }
      laststart=fr->start;
      lastend=fr->end;
    }

    if ( nonabuting>0 ) {
      COMMON_NS::warning_message("restriction fragments do not cover a contiguous region"
		      " in "+thischrom);
    }
    
  }
  
  std::stringstream mymessage;
  mymessage<<"...Loaded restriction enzymes from file "<<filename;
  COMMON_NS::message( mymessage.str() );
  are_restfrags_loaded=1;

}


bool genome::is_duplicate(const std::vector<samfrag> &fragset) {
  // check fragset agaist list to see if it is a duplicate

  std::stringstream sline;
  std::pair< std::map<std::string,int>::iterator , bool > is_dup;
     
  for (int i=0;i<fragset.size();i++) {
    if (fragset[i].chrom=="*") {
      // if it didn't map, have to compare sequence
      sline<<fragset[i].sequence<<" "; 
    } else {
      sline<<fragset[i].chrom<<" "<<fragset[i].start<<" "
	   <<fragset[i].length<<" ";
    }
  }

  // insert if new; does not insert if already in list
  is_dup = list_for_duplicates.insert ( std::pair<std::string,int>( sline.str() ,1) );

  if ( is_dup.second == true ) {
    // it did not already exist, and was added to the map, return false
    return !is_dup.second;
  } else {
    // it did already already exist - count it and return true
    is_dup.first->second++;
    return !is_dup.second;
  }

    
}

void genome::counters::output_interchrom(const std::string &filename) const {
  // Ouput interchromosomal counts

  std::ofstream ouf;
  std::ifstream inf;

  inf.open( filename.c_str() );
  if ( inf.good() ) {
    throw std::runtime_error("file "+filename+" already exists (will not "
			     "overwrite).");
  }
  inf.close();

  ouf.open( filename.c_str() );
  ouf<<"# target name, total intrachromosomal interactions, total "
    "interchromosomal interactions, total interactions, "
    "interactions within 5Mb, interactions within 1Mb"<<std::endl;
  for (it_targs T=me.targets.begin(); T != me.targets.end(); ++T) {
    ouf<<T->name<<"\t"
       <<validPairs.find(T->name)->second<<"\t"
       <<onlyInter.find(T->name)->second<<"\t"
       <<validPairs.find(T->name)->second+onlyInter.find(T->name)->second<<"\t"
       <<within5Mb.find(T->name)->second<<"\t"
       <<within1Mb.find(T->name)->second
       <<std::endl;;
  }

  ouf.close();
  
}


void genome::counters::output_report(const std::string &filename) const {
  // Ouput report

  std::ofstream ouf;
  std::ifstream inf;

  double remainingPC = 100./double(total_read_sets-duplicates_removed);

  inf.open( filename.c_str() );
  if ( inf.good() ) {
    throw std::runtime_error("file "+filename+" already exists (will not "
			     "overwrite).");
  }
  inf.close();

  ouf.open( filename.c_str() );

  ouf<<"##################################################################################"<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"### Report for Capture C Analysis "<<filename<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"##################################################################################"<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"### total number of read fragments                    :  "<<std::setw(13)<<std::right<<total_read_frags<<std::endl;
  ouf<<"### total number of read sets                         :  "<<std::setw(13)<<std::right<<total_read_sets<<std::endl;
  ouf<<"###"<<std::endl; 
  ouf<<"### Duplicates :"<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"###   number of duplicate sets removed                :  "<<std::setw(13)<<std::right<<duplicates_removed
     <<std::setw(3)<<"("<<std::fixed<<std::setprecision(2)<<100.*duplicates_removed/total_read_sets<<" %)"<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"### Of the remaining "<<total_read_sets-duplicates_removed<<" read sets :"<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"###   number of sets where no fragments mapped        :  "<<std::setw(13)<<std::right<<none_mapped
     <<std::setw(3)<<"("<<none_mapped*remainingPC<<" %)"<<std::endl;
  ouf<<"###   number of sets with no targets fragments        :  "<<std::setw(13)<<std::right<<no_targets
     <<std::setw(3)<<"("<<no_targets*remainingPC<<" %)"<<std::endl;
  ouf<<"###   number of sets with multiple targets            :  "<<std::setw(13)<<std::right<<multiple_targets
     <<std::setw(3)<<"("<<multiple_targets*remainingPC<<" %)"<<std::endl;
  ouf<<"###   number of sets with no reporters                :  "<<std::setw(13)<<std::right<<no_reporters
     <<std::setw(3)<<"("<<no_reporters*remainingPC<<" %)"<<std::endl;
  ouf<<"###   number of interactions within exclusion region  :  "<<std::setw(13)<<std::right<<exclusion
     <<std::setw(3)<<"("<<exclusion*remainingPC<<" %)"<<std::endl;
  ouf<<"###   number of interactions with multiple non-adjacent reporters  :  "<<std::setw(13)<<std::right<<multiple_reporters
     <<std::setw(3)<<"("<<multiple_reporters*remainingPC<<" %)"<<std::endl;
  ouf<<"###"<<std::endl; 
  ouf<<"### Valid Interactions :"<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"###   interchromosomal                                :  "<<std::setw(13)<<std::right<<total_interchrom
     <<std::setw(3)<<"("<<total_interchrom*remainingPC<<" %)"<<std::endl;
  ouf<<"###   intrachromosomal                                :  "<<std::setw(13)<<std::right<<total_validPairs
     <<std::setw(3)<<"("<<total_validPairs*remainingPC<<" %)"<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"##################################################################################"<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"### Overall yield : "<<std::endl;
  ouf<<"###                      total read sets : "<<std::setw(13)<<total_read_sets<<std::endl;
  ouf<<"###                    valid interaction : "<<std::setw(13)<<total_interchrom+total_validPairs<<std::endl;
  ouf<<"###                            % overall : "<<std::setw(13)<<100.*double(total_interchrom+total_validPairs)/double(total_read_sets)<<" %"<<std::endl;
  ouf<<"###           % after duplicates removed : "<<std::setw(13)<<100.*double(total_interchrom+total_validPairs)/double(total_read_sets-duplicates_removed)<<" %"<<std::endl;
  ouf<<"### % of valid which are intrachomosomal : "<<std::setw(13)<<100.*total_validPairs/double(total_interchrom+total_validPairs)<<" %"<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"##################################################################################"<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"### Per target information "<<std::endl;
  ouf<<"###"<<std::endl;
  ouf<<"###   "<<std::setw(15)<<std::left<<"target name"
     <<"  "<<std::setw(10)<<std::right<<"intra"
     <<"  "<<std::setw(10)<<"inter"
     <<"  "<<std::setw(10)<<"total"
     <<std::endl;
  for (it_targs T=me.targets.begin(); T != me.targets.end(); ++T) {
    ouf<<"###   "<<std::setw(15)<<std::left<<T->name
       <<"  "<<std::setw(10)<<std::right<<validPairs.find(T->name)->second
       <<"  "<<std::setw(10)<<std::right<<onlyInter.find(T->name)->second
       <<"  "<<std::setw(10)<<std::right<<validPairs.find(T->name)->second+onlyInter.find(T->name)->second
       <<std::endl;;
  }
  ouf<<"###"<<std::endl;
  ouf<<"##################################################################################"<<std::endl;
  
  ouf.close();

  // write a message
  std::stringstream mymessage;
  mymessage<<"...Report written to file "<<filename;
  COMMON_NS::message( mymessage.str() );
  
}



void genome::counters::setup() {
  // set all counters to zero

  if ( me.are_restfrags_loaded == 0 ) {
      throw std::runtime_error("attempted to setup counters before targets "
			       "were loaded.");
  }

  // set global counters
  total_read_frags = 0;
  total_read_sets = 0;
  duplicates_removed = 0;
  none_mapped = 0;
  no_targets = 0;
  multiple_targets = 0;
  no_reporters = 0;
  exclusion = 0;
  multiple_reporters = 0;
  total_interchrom = 0;
  total_validPairs = 0;
  
  // set per target counters
  for (it_targs T=me.targets.begin() ; T != me.targets.end() ; ++T ) {
    validPairs[ T->name ] = 0;
    onlyInter[ T->name ] = 0;
    within1Mb[ T->name ] = 0;
    within5Mb[ T->name ] = 0;
  }
  
}
