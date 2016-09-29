DEFS = -DFP -DSELFENERGY -DSECUMUL

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o .h 

# Compiler
CC=g++ -g -std=c++11 $(DEFS)

#Compiler Options
WARNINGS = -Wall -Wextra -Wpointer-arith  -Wcast-qual -Wcast-align \
 -pedantic-errors  -Wno-long-long -pedantic -Wchar-subscripts -Wcomment \
 -Wdisabled-optimization -Wformat -Wformat=2 -Wformat-nonliteral \
 -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winline \
 -Winvalid-pch -Wunsafe-loop-optimizations -Wmissing-braces \
 -Wmissing-field-initializers -Wmissing-format-attribute \
 -Wmissing-include-dirs -Wmissing-noreturn -Wpacked -Wparentheses \
 -Wreturn-type -Wsequence-point  -Wsign-compare \
 -Wstack-protector -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch \
 -Wswitch-default -Wswitch-enum -Wtrigraphs -Wuninitialized \
 -Wunknown-pragmas -Wunreachable-code -Wunused -Wunused-function \
 -Wunused-label  -Wunused-value -Wunused-variable -Wvariadic-macros \
 -Wvolatile-register-var -Wwrite-strings -Wextra -Wfloat-equal \
 -Wredundant-decls -Wconversion -Wno-unused-parameter
 
#CFLAGS= -O0 -funroll-all-loops -ffast-math -fopenmp $(WARNINGS)
#CFLAGS = $(WARNINGS)
CFLAGS= -O3 -funroll-all-loops -ffast-math -fopenmp
#CFLAGS = -pg -fopenmp

LDFLAGS= -fopenmp
#LDFLAGS= -pg -fopenmp

#Include Path
INC=/project/theorie/h/H.Guertner/lib

#Sources
SOURCE=Adaption.cpp mystructs.cpp dvector.cpp tmap.cpp DiagMC_run.cpp DiagMC_config.cpp DiagMC_io.cpp DiagMC_measure.cpp DiagMC_updates.cpp DiagMC_estimator.cpp DiagMC_test.cpp DiagMC_secumul.cpp Diagram_test.cpp Diagram.cpp

#Objects
OBJ= $(SOURCE:.cpp=.o)

#Header
DEPS= DiagMC.h Diagram.h dvector.h mystructs.h Adaption.h tmap.h DiagMCException.h DiagramException.h RunException.h

#Executable
EXE = ../DiagMC_BEC

#----------------------------------------------

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< -I$(INC)
	
$(EXE): $(OBJ) 
	$(CC) $(LDFLAGS) $(OBJ) -o $@ 

.PHONY: clean
clean:
	rm -f $(OBJ)
	
full: $(EXE) clean
