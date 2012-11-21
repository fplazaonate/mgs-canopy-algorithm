############################################################
##	MACROS
############################################################


GNU_COMP = g++ 
GNU_OPENMP = #-fopenmp

CLANG_COMP = clang++ 
CLANG_OPENMP = #-fopenmp

EXECUTABLE_NAME = cc.bin
TEST_EXECUTABLE_NAME = 

ifndef COMP
COMP = $(GNU_COMP)
COMP_OPENMP = $(GNU_OPENMP)
endif

ADDITIONAL_OPTS = -p -g #-O3 -DAE_CPU=AE_INTEL#-Wall 

ALGLIB_PATH = ../lib/alglib/src/

INCLUDES = -I./ -I../lib/alglib/src/ -I/opt/local/include #Add path to local boost here if needed

LIBS =  -lboost_program_options-mt -lboost_iostreams-mt -L../lib/alglib/src/  -L./ -L/opt/local/lib/

ADDITIONAL_FILES_TO_BE_DELETED = 

############################################################
##	OBJECTS
############################################################

### EXECUTABLE

$(EXECUTABLE_NAME): main.o Point.o Canopy.o CanopyClustering.o Stats.o Log.o 
	$(COMP) -o $@ $(COMP_OPENMP) $(ADDITIONAL_OPTS) $(INCLUDES) $(LIBS) $^

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

Log.o: Log.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

#### Targets

all: $(EXECUTABLE_NAME) 

############################################################
##	CLEAN
############################################################

.PHONY: clean

clean:
	-rm -f *.bin *.o $(ADDITIONAL_FILES_TO_BE_DELETED) $(COMMON_DIR)*.o
 
