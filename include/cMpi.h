/*!
 * Header for the MPI wrapper.
 * Authors : D. THOMAS.
 *
 * COPYRIGHT (C) University of Liège, 2017.
 */

#pragma once

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifdef HAVE_MPI
  typedef MPI_Comm Cupydo_Comm;
#else //HAVE_MPI
  typedef const int* Cupydo_Comm;
#endif //HAVE_MPI