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


#include "messages.h"

#include <string>
#include <iostream>

using namespace COMMON_NS;

void COMMON_NS::message(const std::string &line) {
  // send a message to stdout
  std::cout<<line<<std::endl;
}

void COMMON_NS::warning_message(const std::string &line) {
  // send a warning message to stderr
  std::cerr<<" WARNING ::  "<<line<<std::endl;
}

