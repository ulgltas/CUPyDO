/*
 * Copyright 2018 University of Liege
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
 * Header for CInterfaceMatrix.
 * Authors : D. THOMAS.
 */

#ifndef CINTERFACEMATRIX_H
#define CINTERFACEMATRIX_H

#include <vector>

#include "cMpi.h"

#ifdef HAVE_MPI
#include "petscmat.h"
#endif //HAVE_MPI

#include "cFlexInterfaceData.h"

class CInterfaceMatrix
{
#ifdef HAVE_MPI
  Mat H;
#else  //HAVE_MPI
  std::vector<double> H;
#endif //HAVE_MPI
  int M, N;

public:
  CInterfaceMatrix(int const &val_M, int const &val_N);
  virtual ~CInterfaceMatrix();
  void createDense();
  void createSparse(int val_dnz, int val_onz);
  void createSparseFullAlloc();
  void setValue(int const &iGlobalIndex, int const &jGlobalIndex, double const &value);
  void setValues(int const &m, int const iGlobalIndices[], int const &n, int const jGlobalIndices[], double const values[]);
  void assemble();
  void mult(CFlexInterfaceData *B, CFlexInterfaceData *X);
#ifdef HAVE_MPI
  Mat getMat();
#else  //HAVE_MPI
  void getMat(int *size1, int *size2, double **mat_array);
#endif //HAVE_MPI
};

#endif //CINTERFACEMATRIX_H
