/*!
 * Source for the linear solver.
 * Authors : D. THOMAS.
 *
 * COPYRIGHT (C) University of Liège, 2017.
 */

#ifdef HAVE_MPI
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#endif

#include "../include/cMpi.h"
#include "../include/cLinearSolver.h"

using namespace std;

CLinearSolver::CLinearSolver(CInterfaceMatrix *val_matrixOperator){

#ifdef HAVE_MPI
  KSPCreate(MPI_COMM_WORLD, &KSPSolver);
  KSPSetType(KSPSolver, KSPFGMRES);
  KSPGetPC(KSPSolver, &Precond);
  PCSetType(Precond, PCJACOBI);
  KSPSetOperators(KSPSolver, val_matrixOperator->getMat(),val_matrixOperator->getMat());
  //KSPSetFromOptions(KSPSolver);
  KSPSetInitialGuessNonzero(KSPSolver, PETSC_TRUE);
  KSPSetUp(KSPSolver);
#endif  //HAVE_MPI
}

CLinearSolver::~CLinearSolver(){

#ifdef HAVE_MPI
  if(KSPSolver) KSPDestroy(&KSPSolver);
#endif  //HAVE_MPI

}

#ifdef HAVE_MPI
void CLinearSolver::solve(Vec &VecX, Vec &VecB) const{

  KSPSolve(KSPSolver, VecB, VecX);

}
#endif  //HAVE_MPI
