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


#include "samfragments.h"
#include "genome.h"
#include "bedfiles.h"

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>


using namespace CAPCMAIN_NS;


samfrag samfrag::samline2samfrag(const std::string &line) {

  samfrag thefragment;
  std::stringstream sline;
  std::string junk,
    cigar;
  int flag,
    mapq;
  
  sline.clear(); sline.str(line);
  sline>>thefragment.name;
  thefragment.setname = name2setname(thefragment.name);

  sline>>flag>>thefragment.chrom>>thefragment.start>>mapq>>
    cigar>>junk>>junk>>junk>>junk;

  thefragment.length=junk.size();

  // if it didn't align, store the sequence for checking duplicates
  if (thefragment.chrom=="*") {
    thefragment.sequence=junk;
  }

  // sam files store coordinates in a 1-based coordinate system
  // so subtract 1 to get 0 based
  thefragment.start--;
  
  return thefragment;
  
}




rest_fragment samfrag::expand_to_restfrag(const genome& gnm) const {
  // get the restriction enzyme fragment which this samfrag belongs to

  rest_fragment originalfragment;
  genome::const_it_rest_map achromosome;
  genome::const_it_rest_set afragment;

  originalfragment = *this;
  
  // get a pointer to the list for this chromosome
  achromosome = gnm.restriction_fragments.find(chrom);
  if ( achromosome == gnm.restriction_fragments.end() ) {
      throw std::runtime_error("a sam line was mapped to a chromosome"
			       " not present in the fragments list");
  }

  // find the fragment
  afragment = achromosome->second.upper_bound(originalfragment);

  if ( afragment == achromosome->second.begin() ) {
    // it is at the start -- double check
    if ( !(originalfragment.start >= afragment->start &&
	   originalfragment.start < afragment->end ) ) {
      throw std::runtime_error("something has gone wrong in deterimning"
			       " the fragment a sam line belongs to"
			       " (begining).");
    } else {
      return *afragment;
    }
  }

  if ( afragment == achromosome->second.end() ) {
    --afragment; // go back one
    // it is at the end -- double check
    if ( !(originalfragment.start >= afragment->start &&
	   originalfragment.start <= afragment->end ) ) {
      throw std::runtime_error("something has gone wrong in deterimning"
			       " the fragment a sam line belongs to (end).");
    } else {
      return *afragment;
    }
  }

  --afragment; // go back one
  
  if ( originalfragment.start >= afragment->start &&
       originalfragment.start < afragment->end  ) {
    return *afragment;
  } else {
    throw std::runtime_error("somethinge has gone wrong in deterimning"
			     " the fragment a sam line belongs to.");
  }
  
}
