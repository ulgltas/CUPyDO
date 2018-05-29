/*
 * Copyright 2018 University of Liège
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*!
 * Source for the linear solver.
 * Authors : D. THOMAS.
 */

 #include <iostream>
 #include <cassert>

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
  PCSetUp(Precond);
  KSPSetOperators(KSPSolver, val_matrixOperator->getMat(),val_matrixOperator->getMat());
  //KSPSetFromOptions(KSPSolver);
  KSPSetInitialGuessNonzero(KSPSolver, PETSC_TRUE);
  KSPSetUp(KSPSolver);
  printTolerances();
#endif  //HAVE_MPI
}

CLinearSolver::~CLinearSolver(){

#ifndef NDEBUG
  cout << "Calling CLinearSolver::~CLinearSolver()" << endl;
#endif  //NDEBUG

#ifdef HAVE_MPI
  if(KSPSolver) KSPDestroy(&KSPSolver);
#endif  //HAVE_MPI

}

#ifdef HAVE_MPI
void CLinearSolver::solve(CFlexInterfaceData* B, CFlexInterfaceData* X){

  assert(X->getDim() == B->getDim());

  for(int i=0; i<X->getDim(); i++){
    KSPSolve(KSPSolver, B->getData(i), X->getData(i));
    monitor();
  }

}
#endif  //HAVE_MPI

void CLinearSolver::setMaxNumberIterations(const int& val_maxInt){

#ifdef HAVE_MPI
  KSPGetTolerances(KSPSolver, &relTol, &absTol, &divTol, &maxInt);

  maxInt = val_maxInt;
  KSPSetTolerances(KSPSolver, relTol, absTol, divTol, maxInt);
#endif //HAVE_MPI

}

void CLinearSolver::setRelativeTolerance(const double& val_relTol){

#ifdef HAVE_MPI
  KSPGetTolerances(KSPSolver, &relTol, &absTol, &divTol, &maxInt);
  relTol = val_relTol;
  KSPSetTolerances(KSPSolver, relTol, absTol, divTol, maxInt);

  printTolerances();

#endif //HAVE_MPI


}

void CLinearSolver::setPreconditioner(const std::string& val_precond){  //DO NOT USE !

#ifdef HAVE_MPI

  int nlocal, first;
  KSP *subksp;
  PC subpc;

  KSPGetPC(KSPSolver, &Precond);

  if(val_precond.compare("JACOBI") == 0){
    PCSetType(Precond, PCJACOBI);
    PCSetUp(Precond);
  }
  else if(val_precond.compare("BJACOBI") == 0){
    PCSetType(Precond, PCBJACOBI);
    PCBJacobiSetTotalBlocks(Precond,14,NULL);
    PCSetUp(Precond);
    PCBJacobiGetSubKSP(Precond,&nlocal,&first,&subksp);
    for(int ii=0; ii<nlocal; ii++){
      KSPGetPC(subksp[ii],&subpc);
      PCSetType(subpc,PCSOR);
      KSPSetType(subksp[ii],KSPFGMRES);
    }
  }
  else if(val_precond.compare("SCHWARZ") == 0){
    PCSetType(Precond, PCASM);
    PCASMSetType(Precond, PC_ASM_BASIC);
    PCASMSetTotalSubdomains(Precond,14, NULL,NULL);
    PCASMSetOverlap(Precond, 1);
    PCSetUp(Precond);
    PCASMGetSubKSP(Precond, &nlocal, &first, &subksp);
    for(int ii=0; ii<nlocal; ii++){
      KSPGetPC(subksp[ii],&subpc);
      PCSetType(subpc,PCJACOBI);
      KSPSetType(subksp[ii],KSPFGMRES);
    }

  }
  else if(val_precond.compare("SOR") == 0){
    PCSetType(Precond, PCSOR);
    PCSetUp(Precond);
  }
  else {
    cout << "Preconditioner not recognized, using default JACOBI" << endl;
    PCSetType(Precond, PCJACOBI);
    PCSetUp(Precond);
  }
#endif //HAVE_MPI

}

void CLinearSolver::monitor(){

#ifdef HAVE_MPI
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  KSPGetIterationNumber(KSPSolver, &nInt);
  KSPGetResidualNorm(KSPSolver, &rNorm);

  if(rank==0) cout << "KSP solved in " << nInt << " iterations, residual norm: " << rNorm << endl;
#endif //HAVE_MPI

}

void CLinearSolver::printTolerances(){

#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  KSPGetTolerances(KSPSolver, &relTol, &absTol, &divTol, &maxInt);

  if(rank==0) {
    cout << "Relative tolerance: " << relTol << endl;
    cout << "Absolute tolerance: " << absTol << endl;
    cout << "Diverging tolerance: " << divTol << endl;
    cout << "Maximum number of iterations: " << maxInt << endl;
  }
#endif //HAVE_MPI

}

