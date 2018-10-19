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


#ifndef LOCATION2FRAGMENT_H
#define LOCATION2FRAGMENT_H

#include <string>
#include <vector>

// Forward declarations
namespace CAPCMAIN_NS {
  struct genome;
  struct rest_fragment;
}

namespace LOC2FRAG_NS {

  struct parameters {

    std::string fragfile,
      bedfile,
      outfile;

    std::string location_string,
      location_chrom;
    long int location_start,
      location_end;


    unsigned int rfileflag;

    unsigned int bedfileflag,
      locflag,
      outfileflag;

  };

  void parse_loc2frag_command_line(const int &, char **, parameters &);

  std::vector<CAPCMAIN_NS::rest_fragment
    > locs2frags(const CAPCMAIN_NS::genome &, 
		 const std::vector<CAPCMAIN_NS::rest_fragment> &);

  std::vector<CAPCMAIN_NS::rest_fragment> load_locations(const parameters &);

}

#endif
