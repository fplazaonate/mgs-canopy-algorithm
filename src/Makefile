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

ADDITIONAL_OPTS = -g -p #-03 #-O3#-Wall 

ALGLIB_PATH = ../lib/alglib/src/

INCLUDES = -I./ -I../lib/alglib/src/ -I/opt/local/include #Add path to local boost here if needed

LIBS =  -lboost_program_options-mt -lboost_iostreams-mt -L../lib/alglib/src/  -L./ -L/opt/local/lib/

ADDITIONAL_FILES_TO_BE_DELETED = 

############################################################
##	OBJECTS
############################################################

### EXECUTABLE

$(EXECUTABLE_NAME): main.o Point.o CanopyClustering.o Log.o alglib_ap.o alglib_statistics.o alglib_internal.o alglib_specialfunctions.o alglib_misc.o alglib_linalg.o 
	$(COMP) -o $@ $(COMP_OPENMP) $(ADDITIONAL_OPTS) $(INCLUDES) $(LIBS) $^

main.o:	main.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^
	
Point.o: Point.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

CanopyClustering.o: CanopyClustering.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

Log.o: Log.cpp 
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

alglib_ap.o: $(ALGLIB_PATH)ap.cpp
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

alglib_statistics.o: $(ALGLIB_PATH)statistics.cpp
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

alglib_specialfunctions.o: $(ALGLIB_PATH)specialfunctions.cpp
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

alglib_internal.o: $(ALGLIB_PATH)alglibinternal.cpp
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

alglib_misc.o: $(ALGLIB_PATH)alglibmisc.cpp
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^

alglib_linalg.o: $(ALGLIB_PATH)linalg.cpp
	$(COMP) -o$@ $(COMP_OPENMP) -c $(ADDITIONAL_OPTS) $(INCLUDES) $^
### Targets

all: $(EXECUTABLE_NAME) 

############################################################
##	CLEAN
############################################################

.PHONY: clean

clean:
	-rm -f *.bin *.o $(ADDITIONAL_FILES_TO_BE_DELETED) $(COMMON_DIR)*.o
 
