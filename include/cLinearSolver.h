/*!
 * Header for the linear solver.
 * Authors : D. THOMAS.
 *
 * COPYRIGHT (C) University of Liège, 2017.
 */

#pragma once

#include <string>

#ifdef HAVE_MPI
#include "petscksp.h"
#include "petscvec.h"
#endif  //HAVE_MPI

#include "cInterfaceMatrix.h"

class CLinearSolver{
#ifdef HAVE_MPI
  KSP KSPSolver;
  PC Precond;
  int nInt, maxInt;
  double rNorm, relTol, absTol, divTol;

#endif
public:
  CLinearSolver(CInterfaceMatrix* val_matrixOperator);
  virtual ~CLinearSolver();
#ifdef HAVE_MPI
  void solve(Vec &VecX, Vec &VecB);
#endif  //HAVE_MPI
  void setMaxNumberIterations(const int& val_maxInt);
  void setPreconditioner(const std::string& val_precond);
  void monitor();
  void printTolerances();

};
