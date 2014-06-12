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
#ifndef TIME_PROFILE
#define TIME_PROFILE

#include <time.h>
#include <cmath>
#include <string>

#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

typedef boost::unordered_map<string, pair<bool, time_t> > TimerType;

class TimeProfile {
    private:
        TimerType timers;
    public:
        void start_timer(string timer_name);
        void restart_timer(string timer_name);
        void stop_timer(string timer_name);

        friend ostream& operator<<(ostream& ost, const TimeProfile& c);
};

#endif
