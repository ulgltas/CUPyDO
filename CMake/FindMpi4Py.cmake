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

# Find the mpi4py Python module
#
# Authors : D. THOMAS
#
# This will define
#
#  MPI4PY_FOUND          - TRUE if mpi4py is found
#  MPI4PY_VERSION        - mpi4py version number
#  MPI4PY_INCLUDE_DIR    - where to find mpi4py/mpi4py.i, etc
#  MPI4PY_HEARDERS       - path to mpi4py.h
#  MPI4PY_SWIG_FILE      - path to mpi4py.i
#  MPI4PY_LIBRARIES      - path to MPI library

IF(MPI4PY_INCLUDE_DIR)
    SET(MPI4PY_FIND_QUIETLY TRUE)
ENDIF(MPI4PY_INCLUDE_DIR)

IF(NOT PYTHON_EXECUTABLE)
    IF(MPI4PY_FIND_REQUIRED)
        MESSAGE(SEND_ERROR
          "Python executable not found, so required Mpi4Py module not found"
          )
    ENDIF(MPI4PY_FIND_REQUIRED)
  
ELSE(NOT PYTHON_EXECUTABLE)
    EXECUTE_PROCESS(
        COMMAND ${PYTHON_EXECUTABLE} -c "import distutils.sysconfig as cg; print(cg.get_python_lib(1,0))"
        OUTPUT_VARIABLE PYTHON_SITEDIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    EXECUTE_PROCESS(
        COMMAND ${PYTHON_EXECUTABLE} -c "import mpi4py; from sys import stdout; stdout.write(mpi4py.get_include())"
        OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
        RESULT_VARIABLE MPI4PY_NOT_FOUND
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(Mpi4Py DEFAULT_MSG MPI4PY_INCLUDE_DIR  )
  
    IF(MPI4PY_FOUND)  
        EXECUTE_PROCESS( # Find name of MPI .so
            COMMAND ${PYTHON_EXECUTABLE} -c "from mpi4py import MPI; from sys import stdout; from os import path; stdout.write(path.basename(MPI.__file__))"
            OUTPUT_VARIABLE MPI4PY_LIB_NAME
            RESULT_VARIABLE MPI4PY_NOT_FOUND
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        EXECUTE_PROCESS(  
            COMMAND ${PYTHON_EXECUTABLE} -c "import mpi4py; print(mpi4py.__version__)"
            OUTPUT_VARIABLE MPI4PY_VERSION  
            RESULT_VARIABLE  MPI4PY_NOT_FOUND  
            OUTPUT_STRIP_TRAILING_WHITESPACE  
        )  
        IF(NOT MPI4PY_FIND_QUIETLY)  
            MESSAGE(STATUS "mpi4py version ${MPI4PY_VERSION} found")  
        ENDIF(NOT MPI4PY_FIND_QUIETLY)  
        FIND_FILE(MPI4PY_HEADERS mpi4py.h HINTS ${MPI4PY_INCLUDE_DIR}/mpi4py ${PYTHON_SITEDIR}/mpi4py/include/mpi4py)
        IF(NOT MPI4PY_HEADERS)  
            MESSAGE(STATUS "mpi4py.h not found !")  
        ENDIF(NOT MPI4PY_HEADERS)  
        FIND_FILE(MPI4PY_SWIG_FILE mpi4py.i HINTS ${MPI4PY_INCLUDE_DIR}/mpi4py   ${PYTHON_SITEDIR}/mpi4py/include/mpi4py)
        IF(NOT MPI4PY_SWIG_FILE)  
            MESSAGE(STATUS "mpi4py.i not found !")  
        ENDIF(NOT MPI4PY_SWIG_FILE)  
        FIND_FILE(MPI4PY_LIBRARIES NAMES ${MPI4PY_LIB_NAME} HINTS ${MPI4PY_INCLUDE_DIR}/.. ${PYTHON_SITEDIR}/mpi4py)
    ELSE(MPI4PY_FOUND)  
        IF(MPI4PY_FIND_REQUIRED)  
            MESSAGE(FATAL_ERROR "mpi4py not found !")  
        ENDIF(MPI4PY_FIND_REQUIRED)
    ENDIF(MPI4PY_FOUND)

ENDIF(NOT PYTHON_EXECUTABLE)

MESSAGE(STATUS "MPI4PY_INCLUDE_DIR=${MPI4PY_INCLUDE_DIR}")
MESSAGE(STATUS "MPI4PY_HEADERS=${MPI4PY_HEADERS}")
MESSAGE(STATUS "MPI4PY_SWIG_FILE=${MPI4PY_SWIG_FILE}")
MESSAGE(STATUS "MPI4PY_LIBRARIES=${MPI4PY_LIBRARIES}")

MARK_AS_ADVANCED(MPI4PY_INCLUDE_DIR, MPI4PY_VERSION, MPI4PY_HEARDERS, MPI4PY_SWIG_FILE, MPI4PY_LIBRARIES)
