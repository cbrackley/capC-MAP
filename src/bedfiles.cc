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


#include "bedfiles.h"
#include "samfragments.h"

#include <string>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>

using namespace CAPCMAIN_NS;


std::string CAPCMAIN_NS::name2setname(const std::string& name) {
  // Set names are the read names up to the first / character
  unsigned found;
  found = name.find("DIGEST");
  return name.substr(0,found);
}



// Member functions

void bed_feature::setstrand(const std::string &s) {
  strand = s;
  if (s=="+") { istrand = 1; }
  else if (s=="-") { istrand = -1; }
  else { istrand = 0; }
}

bed_feature bed_feature::line2bed_feature(const std::string &line) {
  // Read from a string to generate a bed line
  // Order must be: chrom start end name score strand
  // Additional information is ignored
  
  std::string c,ss,se,nn,sc,st;
  std::stringstream sline;
    
  sline.str( line );
  sline>>c>>ss>>se;
  if ( ss.find_first_not_of("0123456789") != std::string::npos ||
       se.find_first_not_of("0123456789") != std::string::npos ) {
    throw std::runtime_error("unreadable entry in bed file");
  }

  bed_feature newfeature( c, atoi( ss.c_str() ), atoi( se.c_str() ) );
    
  if ( sline.peek() == EOF ) { return newfeature; }
    
  sline>>nn;
  newfeature.name=nn;
    
  if ( sline.peek() == EOF ) { return newfeature; }

  sline>>sc;
  if ( sc.find_first_not_of("0123456789") != std::string::npos ) {
    throw std::runtime_error("unreadable entry in bed file");
  }
  newfeature.score=atoi( sc.c_str() );
    
  if ( sline.peek() == EOF ) { return newfeature; }

  sline>>st;
  if ( st.find_first_not_of("+-.") != std::string::npos ) {
    throw std::runtime_error("unreadable entry in bed file");
  }
  newfeature.setstrand(st);

  return newfeature;
  
}

bool bed_feature::is_inside(const bed_feature &second) const {
  // check if the feature this is inside second
  // must be wholly inside the region

  if (this->chrom != second.chrom) {return false;}
  
  if (this->start >= second.start &&
      this->end <= second.end ) {
    return true;
  } else {
    return false;
  }
  
}


bool rest_fragment::operator< (const rest_fragment &f2) const {
  // rest frags are sorted by start position
  return start < f2.start;
}




