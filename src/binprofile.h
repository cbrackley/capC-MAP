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


#ifndef BINPROFILE_H
#define BINPROFILE_H

#include "bedgraphfiles.h"

#include <map>
#include <string>

namespace POSTPROCESS_NS {

  struct parameters;
  
  struct binprofile {

    int step,
      window;
    
    std::map<std::string,std::map <bedgraph_entry,double> > profile;
    typedef std::map<std::string,std::map <bedgraph_entry,double> >::iterator it_chromprofile;
    typedef std::map <bedgraph_entry,double>::iterator it_profile;
    

  binprofile(const int &step,const int &window) : step(step), window(window) {}

    void add(const bedgraph_entry &,const std::map<std::string,long int> &);
    
  };

  void do_binning(const parameters &, binprofile&,
		  const std::map<std::string, long int> &);
  std::map<std::string, long int> load_chrom_sizes(const std::string &);

}

#endif
