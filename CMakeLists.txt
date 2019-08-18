#######################################################################################
# The original CMakeList.txt file is downloaded from:
# https://github.com/jedbrown/dohp/blob/master/CMakeLists.txt
# It was modified for my own project.
#######################################################################################

SET (CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")

# These are stuffs I copied and pasted here
# From here
SET (CMAKE_C_COMPILER             "/opt/moose/gcc-8.3.0/bin/gcc")
#my own gcc
#SET (CMAKE_C_COMPILER             "/usr/local/bin/gcc")

#SET (CMAKE_C_FLAGS                "-Wall -std=c99")
#SET (CMAKE_C_FLAGS_DEBUG          "-g")
#SET (CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
#SET (CMAKE_C_FLAGS_RELEASE        "-O4 -DNDEBUG")
#SET (CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

SET (CMAKE_CXX_COMPILER             "/opt/moose/gcc-8.3.0/bin/g++")
#my own temp g++
#SET (CMAKE_CXX_COMPILER             "/usr/local/bin/g++")

#SET (CMAKE_CXX_FLAGS                "-Wall")
#SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
#SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
#SET (CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
#SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
# To here

cmake_minimum_required (VERSION 3.10)
project (PETScApp)

# Get root directory path, show it
list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
message ("The root dir is: " ${CMAKE_SOURCE_DIR})

# Normally PETSc is built with MPI, if not, use CC=mpicc, etc
find_package (PETSc REQUIRED)

include (CheckCSourceCompiles) # The name is misleading, this also tries to link
check_c_source_compiles ("
#define PragmaQuote(a) _Pragma(#a)
#define PragmaGCC(a) PragmaQuote(GCC a)
int main(int argc,char *argv[]) {
  PragmaGCC(diagnostic ignored \"-Wconversion\")
  char c = (int)argv[0][0] + argv[argc-1][0];
  return c;
}" dHAVE_PRAGMA_GCC)

# LZ: Adds flags to the compiler command line. I probably don't need it atm.
# add_definitions (-std=c99)
# additional comments: seems give error message with g++

# Essential: include our directories first otherwise we can get internal headers from some installed path
include_directories (${PETSC_INCLUDES})
add_definitions (${PETSC_DEFINITIONS})
# additional include(s)
include_directories (${CMAKE_SOURCE_DIR})
include_directories (${CMAKE_SOURCE_DIR}/include)

#link_directories (${CMAKE_SOURCE_DIR}/solver)

set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Lists of all source file subdirectories
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR} SRC_ROOT)
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR}/src SRC)
#message ("SRC_ROOT" ${SRC_ROOT})
#message ("SRC_PROBLEMS" ${SRC_PROBLEMS})


ADD_EXECUTABLE (inc_exec ${SRC_ROOT} ${SRC})
#TARGET_LINK_LIBRARIES(PETScApp ${PETSC_LIBRARIES} PETScInterface)
#TARGET_LINK_LIBRARIES(inc_exec ${PETSC_LIBRARIES})
TARGET_LINK_LIBRARIES(inc_exec ${PETSC_LIBRARIES})

# CTEST testing #########################################################
ENABLE_TESTING()

set(System_TESTS CavityMesh)

foreach(TEST ${System_TESTS})
  ADD_TEST(NAME ${TEST}
            COMMAND ${CMAKE_COMMAND}
            -DTEST_PROG=$<TARGET_FILE:inc_exec>
            -DSOURCEDIR=${CMAKE_SOURCE_DIR}
            -DTEST_NAME=${TEST}
            -DTOLERANCE=1e-6  # comment this line for an exact test
            -DNEARZERO=1e-8 # treat any number smaller than 1e-8 as 0
            -P ${CMAKE_SOURCE_DIR}/Testing/runtest.cmake)
            #  message("COMMAND= " ${CMAKE_COMMAND})
            #  message("-DSOURCEDIR= " ${CMAKE_SOURCE_DIR})
            #  message("-DTEST_PROG= " $<TARGET_FILE:PETScApp>)
            #  message("-DTEST_NAME= " ${TEST})
            #  message("-P " ${CMAKE_SOURCE_DIR}/Testing/runtest.cmake)
endforeach()
