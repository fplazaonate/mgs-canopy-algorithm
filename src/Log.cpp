/**
 * Metagenomics Canopy Clustering Implementation
 *
 * Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
 * Copyright (C) 2015 Enterome
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
#include <Log.hpp>
#include <boost/assign/list_of.hpp> 

loglevel_e Logger::log_level = logDEBUG4;
const std::vector<std::string> Logger::valid_verbosities = 
boost::assign::list_of("error")("progress")("warn")("info")("debug")("debug1")("debug2")("debug3");

Logger::Logger() {
}

Logger::~Logger()
{
    buffer << std::endl;
    // This is atomic according to the POSIX standard
    // http://www.gnu.org/s/libc/manual/html_node/Streams-and-Threads.html
    std::cerr << buffer.str();
}


bool Logger::set_log_level(const std::string& verbosity_option)
{
	if(verbosity_option == "error"){
		Logger::log_level = logERR;
		return true;
	}else if(verbosity_option == "progress"){
		Logger::log_level = logPROGRESS;
		return true;
	}else if(verbosity_option == "warn"){
		Logger::log_level = logWARN;
		return true;
	}else if(verbosity_option == "info"){
		Logger::log_level = logINFO;
		return true;
	}else if(verbosity_option == "debug"){
		Logger::log_level = logDEBUG;
		return true;
	}else if(verbosity_option == "debug1"){
		Logger::log_level = logDEBUG1;
		return true;
	}else if(verbosity_option == "debug2"){
		Logger::log_level = logDEBUG2;
		return true;
	}else if(verbosity_option == "debug3"){
		Logger::log_level = logDEBUG3;
		return true;
	}

	return false;
}

