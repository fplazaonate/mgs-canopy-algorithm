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

