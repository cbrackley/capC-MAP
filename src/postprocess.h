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

#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <string>
#include <vector>
#include <set>

namespace POSTPROCESS_NS {

  struct parameters {

    unsigned int Ninfiles;
    std::string targetname;
    std::string chromsizesname;

    std::vector<std::string> infile,
      location;
    
    std::set<std::string> chrom;
    typedef std::set<std::string>::const_iterator itchrom;
    
    std::string outfile;

    unsigned int binwidth,
      windowwidth;
    
    unsigned int binflag,
      normflag,
      interchromflag;
    unsigned long int total_reads;

  };


}

#endif
