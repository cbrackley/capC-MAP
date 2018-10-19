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


#ifndef BEDGRAPHFILES_H
#define BEDGRAPHFILES_H

#include<string>

namespace POSTPROCESS_NS {

  struct bedgraph_entry {

    std::string chrom;
    long int start,
      end;
    double score;

    bedgraph_entry() {};
    bedgraph_entry(const std::string &chrom, const long int &start,
		   const long int &end, const long int &score
		   ) : chrom(chrom), start(start), end(end), score(score) {};
    bedgraph_entry(const std::string &chrom, const long int &start,
		   const long int &end		   
		   ) : chrom(chrom), start(start), end(end), score(0.0) {};

    bool operator< (const bedgraph_entry &) const;

    static bedgraph_entry line2bedgraph_entry(const std::string &);

  };


}

#endif
