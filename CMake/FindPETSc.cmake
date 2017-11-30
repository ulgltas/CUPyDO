# Find the PETSc library
#
# Authors : D. THOMAS
#
# COPYRIGHT (C) University of Liege, 2017.
#
# This will define
#
#  PETSC_DIR                  - where to find petsc4py/petsc4py.i, etc.
#  PETSC4PY_FOUND             - TRUE if PETSc is found
#  PETSC_INCLUDE_PATH         - where to find petsc.h
#  PETSC_LIBRARIES            - path to libpetsc.so

FIND_PATH(PETSC_DIR include/petsc.h HINTS ENV PETSC_DIR PATHS /usr/local/lib/python2.7/dist-packages/petsc /opt/local/lib/petsc $ENV{HOME}/petsc $ENV{PETSC_DIR})

FIND_LIBRARY(PETSC_LIBRARIES NAMES petsc PATHS ${PETSC_DIR}/lib  ${PETSC_DIR}/$ENV{PETSC_ARCH}/lib)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PETSc DEFAULT_MSG PETSC_DIR PETSC_LIBRARIES)

IF(PETSC_FOUND)
  SET(PETSC_INCLUDE_PATH "${PETSC_DIR}/include;${PETSC_DIR}/$ENV{PETSC_ARCH}/include" CACHE PATH "" FORCE)
  MESSAGE(STATUS "PETSC_DIR=${PETSC_DIR}")
  MESSAGE(STATUS "PETSC_LIBRARIES=${PETSC_LIBRARIES}")
  MESSAGE(STATUS "PETSC_INCLUDE_PATH=${PETSC_INCLUDE_PATH}")
ELSE()
ENDIF()
