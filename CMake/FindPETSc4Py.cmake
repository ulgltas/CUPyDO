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

# Find the petsc4py Python module
#
# Authors : D. THOMAS
#
# This will define
#
#  PETSC4PY_INCLUDE_DIR       - where to find petsc4py/petsc4py.i, etc.
#  PETSC4PY_FOUND             - TRUE if petsc4py is found
#  PETSC4PY_VERSION           - version of petsc4py
#  PETSC4PY_HEADERS           - path to petsc4py.h
#  PETSC4PY_SWIG_FILE         - path to petsc4py.i
#  PETSC4PY_LIBRARIES         - path to PETSc library

IF(PETSC4PY_INCLUDE_DIR)
    SET(PETSC4PY_FIND_QUIETLY TRUE)
ENDIF(PETSC4PY_INCLUDE_DIR)

IF(NOT PYTHON_EXECUTABLE)
    IF(PETSC4PY_FIND_REQUIRED)
        MESSAGE(SEND_ERROR
            "Python executable not found, so required PETSc4Py module not found"
        )
    ENDIF(PETSC4PY_FIND_REQUIRED)

ELSE(NOT PYTHON_EXECUTABLE)
    EXECUTE_PROCESS(
        COMMAND ${PYTHON_EXECUTABLE} -c "import distutils.sysconfig as cg; print(cg.get_python_lib(1,0))"
        OUTPUT_VARIABLE PYTHON_SITEDIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    EXECUTE_PROCESS(
        COMMAND ${PYTHON_EXECUTABLE} -c "import petsc4py; from sys import stdout; stdout.write(petsc4py.get_include())"
        OUTPUT_VARIABLE PETSC4PY_INCLUDE_DIR
        RESULT_VARIABLE PETSC4PY_NOT_FOUND
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(PETSc4Py DEFAULT_MSG PETSC4PY_INCLUDE_DIR)

    IF(PETSC4PY_FOUND)
        EXECUTE_PROCESS( # Find name of PETSc .so
            COMMAND ${PYTHON_EXECUTABLE} -c "from petsc4py import PETSc; from sys import stdout; from os import path; stdout.write(path.basename(PETSc.__file__))"
            OUTPUT_VARIABLE PETSC4PY_LIB_NAME
            RESULT_VARIABLE PETSC4PY_NOT_FOUND
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        EXECUTE_PROCESS(
            COMMAND ${PYTHON_EXECUTABLE} -c "import petsc4py; print(petsc4py.__version__)"
            OUTPUT_VARIABLE PETSC4PY_VERSION
            RESULT_VARIABLE PETSC4PY_NOT_FOUND
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        IF(NOT PETSC4PY_FIND_QUIETLY)
            MESSAGE(STATUS "petsc4py version ${PETSC4PY_VERSION} found")
        ENDIF(NOT PETSC4PY_FIND_QUIETLY)
        FIND_FILE(PETSC4PY_HEADERS petsc4py.h HINTS ${PETSC4PY_INCLUDE_DIR}/petsc4py ${PYTHON_SITEDIR}/petsc4py/include/petsc4py)
        IF(NOT PETSC4PY_HEADERS)
            MESSAGE(STATUS "petsc4py.h not found !")
        ENDIF(NOT PETSC4PY_HEADERS)
        FIND_FILE(PETSC4PY_SWIG_FILE petsc4py.i HINTS ${PETSC4PY_INCLUDE_DIR}/petsc4py ${PYTHON_SITEDIR}/petsc4py/include/petsc4py)
        IF(NOT PETSC4PY_SWIG_FILE)
            MESSAGE(STATUS "petsc4py.i not found !")
        ENDIF(NOT PETSC4PY_SWIG_FILE)
        FIND_FILE(PETSC4PY_LIBRARIES NAMES ${PETSC4PY_LIB_NAME} HINTS ${PETSC4PY_INCLUDE_DIR}/../lib ${PETSC4PY_INCLUDE_DIR}/../lib/$ENV{PETSC_ARCH} ${PYTHON_SITEDIR}/petsc4py/lib ${PYTHON_SITEDIR}/petsc4py/lib/$ENV{PETSC_ARCH})
    ELSE(PETSC4PY_FOUND)
        IF(PETSC4PY_FIND_REQUIRED)
            MESSAGE(FATAL_ERROR "petsc4py headers missing")
        ENDIF(PETSC4PY_FIND_REQUIRED)
    ENDIF(PETSC4PY_FOUND)

ENDIF(NOT PYTHON_EXECUTABLE)

MESSAGE(STATUS "PETSC4PY_INCLUDE_DIR=${PETSC4PY_INCLUDE_DIR}")
MESSAGE(STATUS "PETSC4PY_HEADERS=${PETSC4PY_HEADERS}")
MESSAGE(STATUS "PETSC4PY_SWIG_FILE=${PETSC4PY_SWIG_FILE}")
MESSAGE(STATUS "PETSC4PY_LIBRARIES=${PETSC4PY_LIBRARIES}")

MARK_AS_ADVANCED(PETSC4PY_INCLUDE_DIR, PETSC4PY_VERSION, PETSC4PY_HEADERS, PETSC4PY_SWIG_FILE, PETSC4PY_LIBRARIES)
