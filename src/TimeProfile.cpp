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
#include <TimeProfile.hpp>
#include <boost/foreach.hpp>

void TimeProfile::start_timer(string timer_name){
    if(timers.count(timer_name) && timers[timer_name].first)
        throw "Timer is already running";

    timers[timer_name] = make_pair(true, time(NULL)); 
}

void TimeProfile::restart_timer(string timer_name){
    if(!timers.count(timer_name)) 
        throw "Timer doesn't exist";
                                   
    timers[timer_name] = make_pair(true, time(NULL)); 
}

void TimeProfile::stop_timer(string timer_name){
    if((!timers.count(timer_name)) || (!timers[timer_name].first))
        throw "Timer doesn't exist or was already stopped";

    timers[timer_name].first = false;
    timers[timer_name].second = time(NULL) - timers[timer_name].second; 

}

std::ostream& operator<<(std::ostream& ost, const TimeProfile& tp)
{
    ost << ">>>>>>>>>>Time Statistics>>>>>>>>" << std::endl;
    BOOST_FOREACH(TimerType::value_type timer, tp.timers){
        if( !timer.second.first )
            ost << timer.first << ": " << timer.second.second/60 << "min (" << timer.second.second << "sec)" << endl;
    }
    ost << ">>>>>>>>>>END>>>>>>>>" << std::endl;

}

