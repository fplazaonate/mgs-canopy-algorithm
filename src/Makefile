#
#Metagenomics Canopy Clustering Implementation
#
#Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
#Copyright (C) 2015 Enterome
#
#This file is part of Metagenomics Canopy Clustering Implementation.
#
#Metagenomics Canopy Clustering Implementation is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#Metagenomics Canopy Clustering Implementation is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this software.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################
##	MACROS
############################################################


GNU_COMP = g++
GNU_OPENMP = -fopenmp

CLANG_COMP = clang++ 
CLANG_OPENMP = #-fopenmp

EXECUTABLE_NAME = cc.bin

ifndef COMP
COMP = $(GNU_COMP)
COMP_OPENMP = $(GNU_OPENMP)
endif

CFLAGS = -Wall -Wextra -O3 -march=native
BOOST_INCDIR=
INCLUDES = -I./ $(BOOST_INCDIR)

LDFLAGS = 
BOOST_LIBDIR=
LIBS =  $(BOOST_LIBDIR) -lboost_program_options

ADDITIONAL_FILES_TO_BE_DELETED = 

############################################################
##	OBJECTS
############################################################

SRC = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp, %.o, $(SRC))

### EXECUTABLE

$(EXECUTABLE_NAME): $(OBJECTS) 
	$(COMP) -o $@ $(COMP_OPENMP) $(LDFLAGS) $^ $(LIBS) 

%.o: %.cpp
	$(COMP) -o $@ $(COMP_OPENMP) -c $(CFLAGS) $(INCLUDES) $^

#### Targets

all: $(EXECUTABLE_NAME) 

############################################################
##	CLEAN
############################################################

.PHONY: clean

clean:
	rm -f *.bin *.o *.d $(ADDITIONAL_FILES_TO_BE_DELETED) $(COMMON_DIR)*.o
 
