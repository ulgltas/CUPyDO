# Copyright 2018 University of Liège
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Main CMake file.
# Configure the compilation of the C++ core of CUPyDO and the Py wrapper.
#
# Authors : D. THOMAS.

CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
PROJECT(CUPyDO)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")

# macros/fcts
INCLUDE(cupyMacros)

# Type of build
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE )
ENDIF()
MESSAGE(STATUS "Build type : ${CMAKE_BUILD_TYPE}")

# Options
OPTION(WITH_MPI "Build for parallel run" OFF)

# -- C++11
SET(CMAKE_CXX_STANDARD 11) # newer way to set C++11 (requires cmake>=3.1)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add warnings during compilation when build type is Debug
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
ENDIF()

# Build for parallel run using MPI (default is OFF)
IF(WITH_MPI)
    FIND_PACKAGE(MPI REQUIRED)
    IF(MPI_FOUND)
        MESSAGE(STATUS "MPI_INCLUDE_PATH=${MPI_INCLUDE_PATH}")
        MESSAGE(STATUS "MPI_LIBRARIES=${MPI_LIBRARIES}")
        ADD_DEFINITIONS(-DHAVE_MPI)
        FIND_PACKAGE(PETSc REQUIRED)
    ENDIF(MPI_FOUND)
ENDIF(WITH_MPI)

# Python/SWIG
FIND_PACKAGE(PythonInterp 2.7 REQUIRED)
FIND_PACKAGE(PythonLibs 2.7 REQUIRED)
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

INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} ${PETSC_INCLUDE_PATH} )

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#-------------------------------------------------------------
INCLUDE(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES(  ${PROJECT_SOURCE_DIR}/src 
                      ${PYTHON_INCLUDE_PATH} 
                      ${NUMPY_INCLUDE_DIR} 
                      ${PETSC4PY_INCLUDE_DIR} 
                      ${MPI4PY_INCLUDE_DIR})

FILE(GLOB SRCS src/*.h src/*.cpp)

FILE(GLOB ISRCS _src/*.i)

IF(WITH_MPI)
    SET(CMAKE_SWIG_FLAGS "-DHAVE_MPI")
ELSE(WITH_MPI)
    SET(CMAKE_SWIG_FLAGS "")
ENDIF(WITH_MPI)

SET_SOURCE_FILES_PROPERTIES(${ISRCS} PROPERTIES CPLUSPLUS ON)

SET(SWINCFLAGS -I${PROJECT_SOURCE_DIR} )

SET_SOURCE_FILES_PROPERTIES(${ISRCS} PROPERTIES SWIG_FLAGS "${SWINCFLAGS}")

IF(${CMAKE_VERSION} VERSION_LESS "3.8.0")
    SWIG_ADD_MODULE(CCupydo python ${ISRCS} ${SRCS})
ELSE()
    SWIG_ADD_LIBRARY(CCupydo
                     LANGUAGE python
                     SOURCES ${ISRCS} ${SRCS})
ENDIF()

SWIG_LINK_LIBRARIES(CCupydo ${PYTHON_LIBRARIES} ${MPI_LIBRARIES} ${PETSC_LIBRARIES})

# --- Install ---
SET(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
MESSAGE(STATUS "Installation directory : ${CMAKE_INSTALL_PREFIX}")

INSTALL_TARGETS(/ccupydo _CCupydo)
INSTALL_FILES(/ccupydo FILES "${PROJECT_BINARY_DIR}/CCupydo.py")

# --- Enable CTest ---
ENABLE_TESTING()
ADD_SUBDIRECTORY(tests)


