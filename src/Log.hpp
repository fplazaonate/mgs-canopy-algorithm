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
#ifndef _LOG_HPP_
#define _LOG_HPP_

#include <iostream>
#include <sstream>

enum loglevel_e
    {logERR, logPROGRESS, logWARN, logINFO, logDEBUG, logDEBUG1, logDEBUG2, logDEBUG3, logDEBUG4};


extern loglevel_e log_level;

class Logger
{
public:
    Logger(loglevel_e _loglevel = logERR);

    template <typename T> Logger& operator<<(T const & value)
    {
        buffer << value;
        return *this;
    }

    ~Logger();

private:
    std::ostringstream buffer;
};


#define _log(level) \
if (level > log_level) ; \
else Logger(level)

#endif
