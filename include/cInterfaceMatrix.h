/*!
 * Header for CInterfaceMatrix.
 * Authors : D. THOMAS.
 *
 * COPYRIGHT (C) University of Liège, 2017.
 */

#pragma once

#include <vector>

#include "cMpi.h"

#ifdef HAVE_MPI
#include "petscmat.h"
#endif  //HAVE_MPI

class CInterfaceMatrix{
#ifdef HAVE_MPI
Mat H;
#else //HAVE_MPI
std::vector<double> H;
#endif //HAVE_MPI
int M,N;
public:
  CInterfaceMatrix(int const& val_M, int const& val_N);
  virtual ~CInterfaceMatrix();
  void createDense();
  void createSparse(int val_dnz, int val_onz);
  void createSparseFullAlloc();
  void setValue(int const& iGlobalIndex, int const& jGlobalIndex, double const& value);
  void assemble();
#ifdef HAVE_MPI
  Mat getMat();
#else  //HAVE_MPI
  void getMat(int* size1, int* size2, double** mat_array);
#endif  //HAVE_MPI
};
