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
