/*!
 * Header for the linear solver.
 * Authors : D. THOMAS.
 *
 * COPYRIHGHT (C) University of Liège, 2017.
 */

#pragma once

#ifdef HAVE_MPI
#include "petscksp.h"
#include "petscvec.h"
#endif  //HAVE_MPI

#include "cInterfaceMatrix.h"

class CLinearSolver{
#ifdef HAVE_MPI
  KSP KSPSolver;
  PC Precond;
#endif
public:
  CLinearSolver(CInterfaceMatrix* val_matrixOperator);
  virtual ~CLinearSolver();
#ifdef HAVE_MPI
  void solve(Vec &VecX, Vec &VecB) const;
#endif  //HAVE_MPI
};
