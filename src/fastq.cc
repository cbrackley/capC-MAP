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


#include "fastq.h"

#include <string>
#include <sstream>
#include <stdexcept>

using namespace DIGEST_NS;

bool fastq_entry::load_new(std::ifstream &inf) {
  // Read four line fastq entry from file
  // Trim /1 or /2 off name if its there
  
  if ( inf.peek() == EOF ) {
    // we are at the end of the file, return false
    return false;
  }
  
  getline(inf,firstline);
  getline(inf,seq);
  getline(inf,thirdline);
  getline(inf,phred);

  // check four lines were read ok
  if ( inf.fail() || inf.bad() ) {
    throw std::runtime_error("cannot read group of four lines from fastq file");
  }

  // get the name from firstline
  std::stringstream ss(firstline);
  ss>>name;   // first word is the name
  name = name.substr(0,name.find_last_of('/')); // remove /1 or /2 if its there
  
  return true;

}





fastq_entry fastq_entry::extract(const unsigned int &start,
				 const unsigned int &end) const {
  // Given a start and end point, make a new fastq entry
  // String counting starts at 0

  fastq_entry newfq;

  newfq.name = name;
  newfq.seq = seq.substr(start,end-start+1);
  newfq.thirdline = "+";
  newfq.phred = phred.substr(start,end-start+1);

  return newfq;

}
							     
