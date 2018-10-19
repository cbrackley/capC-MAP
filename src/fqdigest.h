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


#ifndef FQDIGEST_H
#define FQDIGEST_H

#include <string>

namespace DIGEST_NS {

  struct filenames {
    // a structure for file names
    std::string fq1,
      fq2,
      outfile;
  };

  struct parameters {
    // a structure for parameter values
    bool long_flag;
    std::string enzyme;
    unsigned int cutpoint;

    parameters();
    
  };

  struct counters {
    // variables for counting digestions
    long unsigned int pairs,
      cuts,
      outfrags;
  counters() : pairs(0), cuts(0), outfrags(0) {};
  };

  void parse_command_line(const int &, char **, filenames &, parameters &);
  void do_digestion(const filenames&, const parameters&, counters&);
  void do_digestion_long(const filenames&, const parameters&, counters&);
  void write_log(const filenames&, const parameters&,const counters&);
  
}

#endif

