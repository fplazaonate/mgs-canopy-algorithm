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
