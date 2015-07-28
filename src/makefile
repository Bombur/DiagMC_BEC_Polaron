SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o .h 

# Compiler
CC=g++ -g -std=c++11 

#Compiler Options
WARNINGS = -Wall -Wextra -Wpointer-arith  -Wcast-qual -Wcast-align \
 -pedantic-errors  -Wno-long-long -pedantic -Wchar-subscripts -Wcomment \
 -Wdisabled-optimization -Wformat -Wformat=2 -Wformat-nonliteral \
 -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winline \
 -Winvalid-pch -Wunsafe-loop-optimizations -Wmissing-braces \
 -Wmissing-field-initializers -Wmissing-format-attribute \
 -Wmissing-include-dirs -Wmissing-noreturn -Wpacked -Wparentheses \
 -Wreturn-type -Wsequence-point -Wshadow -Wsign-compare \
 -Wstack-protector -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch \
 -Wswitch-default -Wswitch-enum -Wtrigraphs -Wuninitialized \
 -Wunknown-pragmas -Wunreachable-code -Wunused -Wunused-function \
 -Wunused-label  -Wunused-value -Wunused-variable -Wvariadic-macros \
 -Wvolatile-register-var -Wwrite-strings -Wextra -Wfloat-equal \
 -Wredundant-decls -Wconversion -Wno-unused-parameter
 
#CFLAGS= -O3 -funroll-all-loops -ffast-math $(WARNINGS)
#CFLAGS = $(WARNINGS)
CFLAGS= -O3 -funroll-all-loops -ffast-math

LDFLAGS= 

#Include Path
INC=/project/theorie/h/H.Guertner/lib

#Sources
SOURCE=Ising_markov_main.cpp Ising.cpp Ising_config.cpp Ising_io.cpp Ising_updates.cpp NeighborTable.cpp 


#Objects
OBJ= $(SOURCE:.cpp=.o)

#Header
DEPS= Ising.h NeighborTable.h myexceptions.h

#Executeable
EXE= ../Ising_markov

#----------------------------------------------

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< -I$(INC)

	
$(EXE): $(OBJ) 
	$(CC) $(LDFLAGS) $(OBJ) -o $@ 

.PHONY: clean
clean:
	rm -f $(OBJ)
	
full: $(EXE) clean