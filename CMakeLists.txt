cmake_minimum_required (VERSION 2.8.12)
project(DiagMC_BEC CXX)
SET(CMAKE_BUILD_TYPE Release) # Release, Debug or RelWithDebInfo

#ALPS
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_PREFIX_PATH "/project/theorie/h/H.Guertner/lib/ALPSCore/build")

#Compiler FLAGS
set(COMP_FLAGS "-std=c++11 -fopenmp")

option(COMP_WARN
	"Use Warning Flags" OFF)
option(COMP_OPT
	"Use Optimization Flags" ON)
option(COMP_DEBUG
	"Use gdb Flag" ON)
option(COMP_TIME
	"gprof time profile" OFF)

message("Compile Options")
if (COMP_WARN)
	set(CFLAGS "${CFLAGS} -Wall -Wextra -Wpointer-arith  -Wcast-qual -Wcast-align \
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
 -Wredundant-decls -Wconversion -Wno-unused-parameter")
 message("-Warnings")
endif (COMP_WARN)

if (COMP_OPT) 
	set(CFLAGS "${CFLAGS} -O3 -funroll-all-loops -ffast-math")
	message("-Optimization")
endif (COMP_OPT)

if (COMP_DEBUG)
	set(CFLAGS "${CFLAGS} -g")
	set(LDFLAGS "${LDFLAGS} -g")
	message("-gdb")
endif (COMP_DEBUG)

if(COMP_TIME)
	set(CFLAGS "${CFLAGS} -pg")
	set(LDFLAGS "${LDFLAGS} -pg")
	message("-gprof")
endif (COMP_TIME)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMP_FLAGS} ${CFLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COMP_FLAGS} ${CFLAGS}")

#Definitions
option(FP 
	"Froehlich Polaron" ON)
option(BEC
	"BEC Polaron" OFF)
option(SELFENERGY 
	"Selfenergy Sampling" ON)
option(SECUMUL
	"Cumulative Sampling" ON)
option(NCHECK 
	"No Debug" OFF)

message("Definitions")
if(FP)
	add_definitions(-DFP)
	message("-FP")
else()
	message("bla")
endif(FP)

if(BEC)
	add_definitions(-DBEC)
	message("-BEC")
endif(BEC)

if(SELFENERGY)
	add_definitions(-DSELFENERGY)
	message("-SELFENERGY")
endif(SELFENERGY)

if(SECUMUL)
	add_definitions(-DSECUMUL)
	message("-SECUMUL")
endif(SECUMUL)

if(NCHECK)
	add_definitions(-DNDEBUG -DNCHECK)
	message("-NDEBUG")
endif(NCHECK)

#Include
include_directories(/project/theorie/h/H.Guertner/lib)
include_directories(/project/theorie/h/H.Guertner/BEC_Polaron/src)

#Source Files
file(GLOB SOURCES "*.cpp")

#Cscope File
file(GLOB HEADERS "*.h")
set(AllFiles ${SOURCES} ${HEADERS})
string(REPLACE ";" "\n" AllFiles "${AllFiles}")
file(WRITE ${CMAKE_SOURCE_DIR}/cscope.files.in "${AllFiles}")
configure_file(${CMAKE_SOURCE_DIR}/cscope.files.in ${CMAKE_SOURCE_DIR}/cscope.files COPYONLY)

#CTags
set_source_files_properties(tags PROPERTIES GENERATED true)
add_custom_target(tags COMMAND ctags --declarations -I -l c++ ./src WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} DEPENDS ${PROJECT_NAME})

#create EXE
add_executable(${PROJECT_NAME} ${SOURCES})

find_package(ALPSCore REQUIRED accumulators)
target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES})
