#ifndef SIGNAL_HANDLERS
#define SIGNAL_HANDLERS

#include <signal.h>
#include <Log.hpp>

extern int terminate_called;

void signal_callback_gentle_handler(int signum);
void signal_callback_die_handler(int signum);

void die_if_true(int terminate_called);

#endif
