/**
 * Metagenomics Canopy Clustering Implementation
 *
 * Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
 *
 * This file is part of Metagenomics Canopy Clustering Implementation.
 *
 * Metagenomics Canopy Clustering Implementation is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Metagenomics Canopy Clustering Implementation is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <Log.hpp>

loglevel_e log_level = logDEBUG4;

Logger::Logger(loglevel_e _loglevel) {
    //buffer << _loglevel << " :" << std::string( 
    //        _loglevel > logDEBUG 
    //            ? (_loglevel - logDEBUG) * 4 
    //            : 1
    //            , ' ');
    //for(int i = _loglevel - logDEBUG + 1; i > 0 ; i--)
    //    buffer << "\t";

}

//template<typename T> Logger& Logger::operator<<(T const & value)

Logger::~Logger()
{
    buffer << std::endl;
    // This is atomic according to the POSIX standard
    // http://www.gnu.org/s/libc/manual/html_node/Streams-and-Threads.html
    std::cerr << buffer.str();
}
