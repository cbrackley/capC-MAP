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


#ifndef BEDFILES_H
#define BEDFILES_H

#include<string>

namespace CAPCMAIN_NS {

  // Forward declarations
  
  struct genome;
  struct samfrag;
  struct rest_fragment;
  
  // Function declarations

  std::string name2setname(const std::string&);

  // Structures
  
  struct bed_feature {

    std::string chrom,
      name,
      strand;
    long int start,
      end;
    int istrand;
    double score;
    
    std::string restofline; // sometimes we want everything after chr:ss-ee

    bed_feature() {};
    bed_feature(const std::string &chrom, const long int &start,
		const long int &end) : chrom(chrom), start(start), end(end) {};

    static bed_feature line2bed_feature(const std::string &);
    void setstrand(const std::string&);
    bool is_inside(const bed_feature&) const;
    
  };


  struct rest_fragment : bed_feature {

    mutable bool is_target;
    mutable std::string targetname;
    
    rest_fragment() {
      // default constructor
      is_target = 0;
    }

    rest_fragment(const bed_feature&B) : bed_feature(B.chrom,B.start,B.end) { 
      // construct from base struct
      is_target = 0;
    }

    rest_fragment(const std::string &chrom, const long int &start,
		  const long int &end) : bed_feature(chrom,start,end) {
    }

    bool operator< (const rest_fragment &) const;
      
  };

  
}


#endif
