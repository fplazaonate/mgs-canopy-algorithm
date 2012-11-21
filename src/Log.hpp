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
