#
#Metagenomics Canopy Clustering Implementation
#
#Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
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

#ADDITIONAL_OPTS = -p -g
ADDITIONAL_OPTS = -O3 -msse4.2 

INCLUDES = -I./ -I/tools/boost/include/ #Add path to local boost here if needed

LIBS =  -L./ -L/tools/boost/lib/ -Wl,-Bstatic -lboost_program_options  -Wl,-Bdynamic

ADDITIONAL_FILES_TO_BE_DELETED = 

############################################################
##	OBJECTS
############################################################

### EXECUTABLE

$(EXECUTABLE_NAME): main.o Point.o Canopy.o CanopyClustering.o Stats.o Log.o signal_handlers.o TimeProfile.o
	$(COMP) -o $@ $(COMP_OPENMP) $(ADDITIONAL_OPTS) $(INCLUDES) $^ $(LIBS) 

main.o:	main.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^
	
Point.o: Point.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

Canopy.o: Canopy.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

CanopyClustering.o: CanopyClustering.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

Stats.o: Stats.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

TimeProfile.o: TimeProfile.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

signal_handlers.o: signal_handlers.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

Log.o: Log.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

#### Targets

all: $(EXECUTABLE_NAME) 

############################################################
##	CLEAN
############################################################

.PHONY: clean

clean:
	-rm -f *.bin *.o *.d $(ADDITIONAL_FILES_TO_BE_DELETED) $(COMMON_DIR)*.o
 
