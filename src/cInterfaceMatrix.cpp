/*!
 * Source for CInterfaceMatrix.
 * Authors : D. THOMAS.
 *
 * COPYRIGHT (C) University of Liège, 2017.
 */

#include <iostream>
#include <vector>

#ifdef HAVE_MPI
#include "petscmat.h"
#endif  //HAVE_MPI

#include "../include/cMpi.h"
#include "../include/cInterfaceMatrix.h"

using namespace std;

CInterfaceMatrix::CInterfaceMatrix(int const& val_M, int const& val_N):M(val_M), N(val_N){ }

CInterfaceMatrix::~CInterfaceMatrix(){

#ifdef HAVE_MPI
  if(H){
    MatDestroy(&H);
  }
  else{  }
#endif  //HAVE_MPI

}

void CInterfaceMatrix::createDense(){

#ifdef HAVE_MPI
  //MatCreateDense(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, M, N, NULL, &H);
  MatCreateAIJ(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, M, N, N, NULL, N, NULL, &H);
#else  //HAVE_MPI
  H.resize(M*N);
#endif  //HAVE_MPI

}

void CInterfaceMatrix::createSparse(int val_dnz, int val_onz){

#ifdef HAVE_MPI
  MatCreateAIJ(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, M, N, val_dnz, NULL, val_onz, NULL, &H);
#else //HAVE_MPI
  H.resize(M*N);
#endif //HAVE_MPI

}

void CInterfaceMatrix::createSparseFullAlloc(){

#ifdef HAVE_MPI
  MatCreateAIJ(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, M, N, N, NULL, N, NULL, &H);
#else //HAVE_MPI
  H.resize(M*N);
#endif //HAVE_MPI

}

void CInterfaceMatrix::setValue(const int &iGlobalIndex, const int &jGlobalIndex, double const &value){

#ifdef HAVE_MPI
  MatSetValue(H, iGlobalIndex, jGlobalIndex, value, INSERT_VALUES);
#else  //HAVE_MPI
  H[iGlobalIndex*N+jGlobalIndex] = value;
#endif  //HAVE_MPI
}

void CInterfaceMatrix::setValues(int const& m, int const iGlobalIndices[], int const& n, int const jGlobalIndices[], double const values[]){

#ifdef HAVE_MPI
  MatSetValues(H, m, iGlobalIndices, n, jGlobalIndices, values, INSERT_VALUES);
#else //HAVE_MPI
  for(int ii=0; ii<m; ii++){
    for(int jj=0; jj<n; jj++){
      H[iGlobalIndices[ii]*N+jGlobalIndices[jj]] = values[ii*n+jj];
    }
  }
#endif //HAVE_MPI
}

void CInterfaceMatrix::assemble(){

#ifdef HAVE_MPI
  MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
#endif  //HAVE_MPI
}

#ifdef HAVE_MPI
Mat CInterfaceMatrix::getMat(){

  return H;
}
#else //HAVE_MPI
void CInterfaceMatrix::getMat(int* size1, int* size2, double** mat_array){

  *size1 = M;
  *size2 = N;
  *mat_array = &(H.front());;
}

#endif  //HAVE_MPI
