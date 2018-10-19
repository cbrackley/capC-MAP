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


#include "binprofile.h"
#include "postprocess.h"
#include "bedgraphfiles.h"

#include <map>
#include <stdexcept>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>

#include <iostream>

using namespace POSTPROCESS_NS;

void POSTPROCESS_NS::do_binning(const parameters &params, binprofile &binned,
				const std::map<std::string, long int> &chromsizes) {
  // Read bed file, and put into binned profile

  bedgraph_entry newbed;
  std::ifstream inf;
  std::string line;
  std::ofstream ouf;
  double normfactor;

  inf.open( params.outfile.c_str() );
  if ( inf.good() ) {
    throw std::runtime_error("file "+params.outfile+" already exists (will not "
			     "overwrite).");
  }
  inf.close();
  

  // normalize by bin to window ratio, so that sum of
  // reads should roughly be the same.
  normfactor = double(binned.step)/double(binned.window);

  if ( params.normflag == 1 && params.total_reads>0 ) {
    // Normalization to reads per million has been selected.
    // Keep factor above, but add reads-per-million factor.
    normfactor= normfactor*1.0e6/double(params.total_reads);
  }

  inf.open( params.infile[0].c_str() );
  if ( ! inf.good() ) {
    throw std::runtime_error("cannot open file "+params.infile[0]);
  }

  
  while ( getline(inf,line) )  {
    if ( line.compare(0,5,"track")!=0 ) { // ignore track definition lines
      newbed = bedgraph_entry::line2bedgraph_entry(line);
      binned.add(newbed,chromsizes);
    }
  }

  inf.close();

  // Now output
  ouf.open( params.outfile.c_str() );

  ouf<<"track type=bedGraph name=\""<<params.targetname<<"_"<<params.binwidth<<"_"<<params.windowwidth
     <<"\" description=\""<<"Target: "<<params.targetname
     <<" binstep="<<params.binwidth<<" window="<<params.windowwidth
     <<"\""<<std::endl;

  for (binprofile::it_chromprofile C = binned.profile.begin(); C != binned.profile.end() ;
       ++C ) {
    for (binprofile::it_profile P = C->second.begin() ; P != C->second.end() ; ++P ) {
      ouf<<P->first.chrom<<"\t"
	 <<P->first.start<<"\t"
	 <<P->first.end<<"\t"
	 <<P->second*normfactor<<std::endl;
    }
  }
  ouf.close();

  
}

std::map<std::string, long int> POSTPROCESS_NS::load_chrom_sizes(const std::string &filename) {

  std::map<std::string, long int> sizes;
  std::ifstream inf;
  std::stringstream sline;
  std::string line,
    chrom,
    slen;
  long int len;

  inf.open( filename.c_str() );
  if ( ! inf.good() ) {
    throw std::runtime_error("cannot open file "+filename);
  }

  while ( getline(inf,line) )  {
    sline.clear();
    sline.str( line );

    sline>>chrom>>slen;

    if ( sizes.count(chrom) > 0 ) {
      throw std::runtime_error("duplicate chromosome in chrom sizes file");
    }
    
    if ( slen.find_first_not_of("0123456789") != std::string::npos ) {
      throw std::runtime_error("unreadable entry in chrom sizes file");
    }

    len = atoi( slen.c_str() );

    sizes[chrom] = len;
    
  }

  inf.close();
  
  return sizes;
  
}


void binprofile::add(const bedgraph_entry &A,
		     const std::map<std::string,long int> &chromsizes) {
  // add a bed feature to the binned profile

  std::string chrom = A.chrom;
  it_chromprofile C;
  
  int lowbin,     // the lowest bin which will include this bed line
    highbin,      // the highest bin which will include this bed line
    midpoint;     // mid point of the bed line

  if ( chromsizes.count(chrom) == 0 ) {
    throw std::runtime_error("unknown chromosome "+chrom);
  }
  
  midpoint = int( 0.5*(A.start+A.end) );

  // get the left edge of the bins this will go into
  lowbin = int( double(midpoint-0.5*window)/double(step) );
  highbin = int( double(midpoint+0.5*window)/double(step) );

  if ( lowbin<0 ) {
    lowbin = 0;
  }
  while ( highbin+step > chromsizes.at(chrom) ) {
    highbin -= step;
  }

  C = profile.find(chrom);
  if ( C == profile.end() ) {
    // chromosome not already in the map
    profile.insert( std::pair<std::string,std::map<bedgraph_entry,double> >(chrom, std::map<bedgraph_entry,double>() ) );
  }
  C = profile.find(chrom);
  
  for (int start=lowbin;start<=highbin;start++) {
    C->second[ bedgraph_entry(chrom,start*step,start*step+step) ]+=A.score;
  }
  
}
