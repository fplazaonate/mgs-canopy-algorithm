#include <signal_handlers.hpp>
#include <stdlib.h>

using namespace std;


int terminate_called = 0;

void signal_callback_gentle_handler(int signum){
    _log(logERR) << "Received signal: " << signum;
    terminate_called += 1;
}

void signal_callback_die_handler(int signum){
    _log(logERR) << "Received signal: " << signum << " Bye! Bye!";
    exit(1);
}

void die_if_true(int terminate_called){
    if(terminate_called)
        exit(1);
}
