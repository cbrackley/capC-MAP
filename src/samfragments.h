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


#ifndef SAMFRAGMENTS_H
#define SAMFRAGMENTS_H

#include "bedfiles.h"
#include <string>

namespace CAPCMAIN_NS {

  // Forward Declarations

  struct genome;


  
  // Structures
  struct samfrag : bed_feature {

    std::string setname,
      sequence;
    unsigned int length;

    static samfrag samline2samfrag(const std::string &);
    
    rest_fragment  expand_to_restfrag(const genome&) const;
    
  };
    
}

#endif
