#######################################################################################
# The original CMakeList.txt file is downloaded from:
# https://github.com/jedbrown/dohp/blob/master/CMakeLists.txt
# It was modified for my own project.
#######################################################################################

SET (CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")

# Compiler selection: if the caller hasn't already forced a compiler (via the CC/CXX
# environment variables, or -DCMAKE_C_COMPILER/-DCMAKE_CXX_COMPILER), prefer the MPI
# compiler wrappers on PATH. This matches whatever PETSc/MPI environment is currently
# active (e.g. a conda environment like 'moose' activated via 'conda activate moose'),
# since PETSc itself is typically built with mpicc/mpicxx. Falls back to CMake's normal
# compiler auto-detection if no MPI wrapper is found.
if (NOT CMAKE_C_COMPILER AND NOT DEFINED ENV{CC})
  find_program (_default_mpicc mpicc)
  if (_default_mpicc)
    SET (CMAKE_C_COMPILER "${_default_mpicc}")
  endif ()
endif ()
if (NOT CMAKE_CXX_COMPILER AND NOT DEFINED ENV{CXX})
  find_program (_default_mpicxx mpicxx)
  if (_default_mpicxx)
    SET (CMAKE_CXX_COMPILER "${_default_mpicxx}")
  endif ()
endif ()

cmake_minimum_required (VERSION 3.10)
project (PETScApp)

# No build type was ever set here, which means neither GCC nor Clang add any
# optimization flag (i.e. an implicit -O0). Default to Release (-O3 -DNDEBUG)
# unless the caller asked for something else, e.g. -DCMAKE_BUILD_TYPE=Debug.
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif ()

# Get root directory path, show it
list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
message ("The root dir is: " ${CMAKE_SOURCE_DIR})

# Normally PETSc is built with MPI, if not, use CC=mpicc, etc
find_package (PETSc REQUIRED)

# Some PETSc builds (configured with --download-exodusii) ship a real ExodusII
# library alongside libpetsc. When available, link it directly too: ExodusOutput.C
# calls the raw ExodusII C API (ex_put_time) itself, which isn't reachable through
# PETSC_LIBRARIES alone -- PETSc's own minimal link check only needs libpetsc.
get_filename_component (PETSC_LIB_DIR ${PETSC_LIBRARIES} DIRECTORY)
find_library (EXODUS_LIBRARY NAMES exodus HINTS ${PETSC_LIB_DIR})
if (EXODUS_LIBRARY)
  message (STATUS "Found ExodusII library: ${EXODUS_LIBRARY}")
endif ()

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
include_directories (/Users/lingzou/projects/Fluent2DMesh FLUENT_MESH)

#link_directories (${CMAKE_SOURCE_DIR}/solver)

set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Lists of all source file subdirectories
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR} SRC_ROOT)
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR}/src SRC)
AUX_SOURCE_DIRECTORY(/Users/lingzou/projects/Fluent2DMesh FLUENT_MESH)
#message ("SRC_ROOT" ${SRC_ROOT})
#message ("SRC_PROBLEMS" ${SRC_PROBLEMS})


ADD_EXECUTABLE (inc_exec ${SRC_ROOT} ${SRC} ${FLUENT_MESH})
#TARGET_LINK_LIBRARIES(PETScApp ${PETSC_LIBRARIES} PETScInterface)
#TARGET_LINK_LIBRARIES(inc_exec ${PETSC_LIBRARIES})
TARGET_LINK_LIBRARIES(inc_exec ${PETSC_LIBRARIES} ${EXODUS_LIBRARY})

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
