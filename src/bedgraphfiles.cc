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


#include "bedgraphfiles.h"

#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

using namespace POSTPROCESS_NS;

bool bedgraph_entry::operator< (const bedgraph_entry &f2) const {
  // binned bed features are sorted by start position
  return start < f2.start;
}

bedgraph_entry bedgraph_entry::line2bedgraph_entry(const std::string &line) {
  // Read from a string to generate a bedgraph line
  // Order must be: chrom start end score
  // Additional information is ignored
  
  std::string c,ss,se,sc;
  std::stringstream sline;

  sline.str( line );
  sline>>c>>ss>>se>>sc;
  if ( ss.find_first_not_of("0123456789") != std::string::npos ||
       se.find_first_not_of("0123456789") != std::string::npos ) {
    throw std::runtime_error("unreadable entry in bedgraph file");
  }

  return bedgraph_entry( c, atoi( ss.c_str() ), atoi( se.c_str() ), atof( sc.c_str() ) );

}
