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

# Find the native numpy Python module
#
# Authors : D. THOMAS
#
# This will define
#
#  NUMPY_INCLUDE_DIR       - where to find numpy/arrayobject.h, etc.
#  NUMPY_FOUND             - TRUE if numpy is found

IF(NUMPY_INCLUDE_DIR)
    SET(NUMPY_FIND_QUIETLY TRUE)
ENDIF(NUMPY_INCLUDE_DIR)

IF(NOT PYTHON_EXECUTABLE)
    IF(NUMPY_FIND_REQUIRED)
        MESSAGE(SEND_ERROR
          "Python executable not found, so required NumPy module not found"
          )
    ENDIF(NUMPY_FIND_REQUIRED)

ELSE(NOT PYTHON_EXECUTABLE)
    EXECUTE_PROCESS(
        COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; from sys import stdout; stdout.write(numpy.get_include())"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        RESULT_VARIABLE NUMPY_NOT_FOUND
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(NumPy DEFAULT_MSG NUMPY_INCLUDE_DIR)

    IF(NUMPY_FOUND)
        EXECUTE_PROCESS(
            COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print(numpy.__version__)"
            OUTPUT_VARIABLE NUMPY_VERSION
            RESULT_VARIABLE  NUMPY_NOT_FOUND
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    IF(NOT NUMPY_FIND_QUIETLY)
        MESSAGE(STATUS "numpy version ${NUMPY_VERSION} found")
    ENDIF(NOT NUMPY_FIND_QUIETLY)
    ELSE(NUMPY_FOUND)
        IF(NUMPY_FIND_REQUIRED)
            MESSAGE(FATAL_ERROR "numpy not found !")
        ENDIF(NUMPY_FIND_REQUIRED)
    ENDIF(NUMPY_FOUND)
ENDIF(NOT PYTHON_EXECUTABLE)  

MESSAGE(STATUS "NUMPY_INCLUDE_DIR=${NUMPY_INCLUDE_DIR}")

MARK_AS_ADVANCED(NUMPY_INCLUDE_DIR, NUMPY_VERSION)

