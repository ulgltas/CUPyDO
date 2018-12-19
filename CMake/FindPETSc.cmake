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

# Find the PETSc library
#
# Authors : D. THOMAS & co
#
# This will define
#
#  <PETSC_DIR>                - root folder of petsc
#  PETSC_INCLUDE_PATH         - where to find petsc.h
#  PETSC_LIBRARIES            - path to libpetsc.so

FIND_PATH(PETSC_DIR include/petsc.h HINTS ENV PETSC_DIR PATHS /opt/local/lib/petsc $ENV{HOME}/petsc $ENV{PETSC_DIR})

IF (PETSC_DIR)
  FIND_LIBRARY(PETSC_LIBRARIES NAMES petsc PATHS ${PETSC_DIR}/lib  ${PETSC_DIR}/$ENV{PETSC_ARCH}/lib)

  INCLUDE(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(PETSc DEFAULT_MSG PETSC_DIR PETSC_LIBRARIES)

  IF(PETSC_FOUND)
    SET(PETSC_INCLUDE_PATH "${PETSC_DIR}/include;${PETSC_DIR}/$ENV{PETSC_ARCH}/include" CACHE PATH "" FORCE)
    MESSAGE(STATUS "PETSC_DIR=${PETSC_DIR}")
    MESSAGE(STATUS "PETSC_LIBRARIES=${PETSC_LIBRARIES}")
    MESSAGE(STATUS "PETSC_INCLUDE_PATH=${PETSC_INCLUDE_PATH}")
  ENDIF()
ELSE()
  FIND_PATH(PETSC_INCLUDE_PATH "petsc.h")
  FIND_LIBRARY(PETSC_LIBRARIES "petsc")

  INCLUDE(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(PETSc DEFAULT_MSG PETSC_INCLUDE_PATH PETSC_LIBRARIES)

  IF(PETSC_FOUND)
    MESSAGE(STATUS "PETSC_LIBRARIES=${PETSC_LIBRARIES}")
    MESSAGE(STATUS "PETSC_INCLUDE_PATH=${PETSC_INCLUDE_PATH}")
  ENDIF()
ENDIF()
