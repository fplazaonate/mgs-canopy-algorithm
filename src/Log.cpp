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
