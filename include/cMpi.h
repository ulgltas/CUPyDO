/*!
 * Header for the MPI wrapper.
 * Authors : D. THOMAS.
 *
 * COPYRIHGHT (C) University of Liège, 2017.
 */

#pragma once

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifdef HAVE_MPI
  typedef MPI_Comm Cupydo_Comm;
#else //HAVE_MPI
  typedef int Cupydo_Comm;
#endif //HAVE_MPI
