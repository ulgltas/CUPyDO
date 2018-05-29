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
 * Header for the linear solver.
 * Authors : D. THOMAS.
 */

#pragma once

#include <string>

#ifdef HAVE_MPI
#include "petscksp.h"
#include "petscvec.h"
#endif  //HAVE_MPI

#include "cInterfaceMatrix.h"
#include "cFlexInterfaceData.h"

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
  void solve(CFlexInterfaceData* B, CFlexInterfaceData* X);
#endif  //HAVE_MPI
  void setMaxNumberIterations(const int& val_maxInt);
  void setRelativeTolerance(const double& val_relTol);
  void setPreconditioner(const std::string& val_precond);
  void monitor();
  void printTolerances();

};
