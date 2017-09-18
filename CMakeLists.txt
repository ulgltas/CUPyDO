# Main CMake file.
# Configure the compilation of the C++ core of CUPyDO and the Py wrapper.
#
# Authors : D. THOMAS.
#
# COPYRIGHT (C) University of Liege, 2017.

PROJECT(CUPyDO)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../ccupydo CACHE PATH
                        "Single output directory for building all libraries.")
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../ccupydo CACHE PATH
                        "Single output directory for building all executables.")

MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")

# Check for C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# Build for parallel run using MPI (default is OFF)
OPTION(WITH_MPI "Build for parallel run" OFF)
# Build Python wrapper (default is ON)
OPTION(WITH_PYWRAPPER "Build with the Python bindings" ON)

IF(WITH_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  IF(MPI_FOUND)
    MESSAGE(STATUS "MPI_INCLUDE_PATH=${MPI_INCLUDE_PATH}")
    MESSAGE(STATUS "MPI_LIBRARIES=${MPI_LIBRARIES}")
    ADD_DEFINITIONS(-DHAVE_MPI)
    FIND_PACKAGE(PETSc REQUIRED)
  ENDIF(MPI_FOUND)
ENDIF(WITH_MPI)

IF(WITH_PYWRAPPER)
  FIND_PACKAGE(PythonInterp REQUIRED)
  FIND_PACKAGE(PythonLibs REQUIRED)
  MESSAGE(STATUS "PYTHONINTERP_FOUND=${PYTHONINTERP_FOUND}")
  MESSAGE(STATUS "PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}")
  MESSAGE(STATUS "PYTHONLIBS_FOUND=${PYTHONLIBS_FOUND}")
  MESSAGE(STATUS "PYTHON_LIBRARIES=${PYTHON_LIBRARIES}")
  MESSAGE(STATUS "PYTHON_INCLUDE_PATH=${PYTHON_INCLUDE_PATH}")
  FIND_PACKAGE(NUMPY REQUIRED)
  FIND_PACKAGE(SWIG REQUIRED)
  SET(CMAKE_SWIG_OUTDIR "${EXECUTABLE_OUTPUT_PATH}")
  MESSAGE(STATUS "SWIG_EXECUTABLE=${SWIG_EXECUTABLE}")
  IF(WITH_MPI)
    FIND_PACKAGE(PETSc4Py REQUIRED)
    FIND_PACKAGE(Mpi4Py REQUIRED)
  ENDIF(WITH_MPI)
ENDIF(WITH_PYWRAPPER)

INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} ${PETSC_INCLUDE_PATH} )

#FILE(GLOB SRCS test.cpp)
#ADD_EXECUTABLE(cTestGeneral src/cTest_general.cpp src/cInterpolator.cpp src/cInterfaceMatrix.cpp src/cManager.cpp src/adtcore.cpp src/cAdt.cpp include/cInterpolator.h include/cManager.h include/cInterfaceMatrix.h include/adtcore.h include/cAdt.h )
#TARGET_LINK_LIBRARIES(cTestGeneral ${MPI_LIBRARIES} ${PETSC_LIBRARIES})

#-------------------------------------------------------------
INCLUDE(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include ${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_DIR} ${PETSC4PY_INCLUDE_DIR} ${MPI4PY_INCLUDE_DIR})

FILE(GLOB SRCS src/cInterpolator.cpp src/cManager.cpp src/cFlexInterfaceData.cpp src/cInterfaceMatrix.cpp src/cLinearSolver.cpp src/adtcore.cpp src/cAdt.cpp include/cInterpolator.h include/cManager.h include/cMpi.h include/cFlexInterfaceData.h include/cInterfaceMatrix.h include/cLinearSolver.h include/adtcore.h include/cAdt.h )
FILE(GLOB ISRCS include/CCupydo.i)

IF(WITH_MPI)
  SET(CMAKE_SWIG_FLAGS "-DHAVE_MPI")
ELSE(WITH_MPI)
  SET(CMAKE_SWIG_FLAGS "")
ENDIF(WITH_MPI)

SET(CMAKE_SWIG_OUTDIR ${PROJECT_SOURCE_DIR}/ccupydo)

SET_SOURCE_FILES_PROPERTIES(${ISRCS} PROPERTIES CPLUSPLUS ON)

SET(SWINCFLAGS -I${PROJECT_SOURCE_DIR} )

SET_SOURCE_FILES_PROPERTIES(${ISRCS} PROPERTIES SWIG_FLAGS "${SWINCFLAGS}")

SWIG_ADD_MODULE(CCupydo python ${ISRCS} ${SRCS})

SWIG_LINK_LIBRARIES(CCupydo ${PYTHON_LIBRARIES} ${MPI_LIBRARIES} ${PETSC_LIBRARIES})
