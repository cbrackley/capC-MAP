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


#ifndef GENOME_H
#define GENOME_H

#include "targets.h"

#include<string>
#include<set>
#include<map>
#include<vector>

namespace CAPCMAIN_NS {

  // Forward Declarations
  
  struct rest_fragment;
  //struct target;
  struct samfrag;

  
  // Structures
  
  struct genome {
    // Genome is defined by its restriction fragments, stored as vectors in
    // a map where the keys are the chomosome names
    
    int Nchrom,
      Ntargets;
    
    std::map<std::string, std::set<rest_fragment> > restriction_fragments;
    typedef std::map<std::string,
      std::set<rest_fragment> >::iterator it_rest_map;
    typedef std::map<std::string,
      std::set<rest_fragment> >::const_iterator const_it_rest_map;
    typedef std::set<rest_fragment>::iterator it_rest_set;
    typedef std::set<rest_fragment>::const_iterator const_it_rest_set;
	
    std::set<target> targets;
    typedef std::set<target>::iterator it_targs;

    std::map<std::string,int> list_for_duplicates;

    bool are_targets_loaded,
      are_restfrags_loaded;

    void load_targets(const std::string &);
    void load_rest_frags(const std::string &);
    bool is_duplicate(const std::vector<samfrag> &);

    genome() : count(*this) {
      // constructor
      are_targets_loaded=0;
      are_restfrags_loaded=0;
    }

    struct counters {
      // structure of various counters
      genome &me;
      
      // global counters
      long unsigned int total_read_frags,
	total_read_sets,
	duplicates_removed,
	none_mapped,
	no_targets,
	multiple_targets,
	no_reporters,
	exclusion,
	multiple_reporters,
	total_interchrom,
	total_validPairs;

      // per target counters
      std::map<std::string, long unsigned int> validPairs,
	onlyInter;

      // constructor
      counters(genome &g) : me(g) {};
      
      // setup function
      void setup();
      
      // outputs
      void output_interchrom(const std::string &) const;
      void output_report(const std::string &) const;
      
    } count;
    
  };

  
}

#endif
