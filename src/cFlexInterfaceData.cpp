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
 * Source for CFlexInterfaceData.
 * Authors : D. THOMAS., M.L. CERQUAGLIA
 */


#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

#ifdef HAVE_MPI
#include "petscvec.h"
#endif  //HAVE_MPI

#include "../include/cMpi.h"
#include "../include/cFlexInterfaceData.h"

using namespace std;

CFlexInterfaceData::CFlexInterfaceData(int const& val_nPoint, int const& val_nDim, Cupydo_Comm val_comm):nPoint(val_nPoint), nDim(val_nDim), comm(val_comm){

  dataContainer.resize(nDim);

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecCreateMPI(comm, PETSC_DECIDE, nPoint, &(dataContainer[ii]));
    VecSet(dataContainer[ii], 0.0);
  }
#else //HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    dataContainer[ii] = new double[nPoint];
    fill(dataContainer[ii], dataContainer[ii]+nPoint, 0.0);
  }
#endif  //HAVE_MPI
}

CFlexInterfaceData::CFlexInterfaceData(CFlexInterfaceData &data){

  nPoint = data.nPoint;
  nDim = data.nDim;
  comm = data.comm;

  dataContainer.resize(nDim);

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecCreateMPI(comm, PETSC_DECIDE, nPoint, &(dataContainer[ii]));
    VecCopy(data.getData(ii), dataContainer[ii]);
  }
#else  //HAVE_MPI
  double* dataToCopy;
  int size;
  for(int ii=0; ii<nDim; ii++){
    dataContainer[ii] = new double[nPoint];
    data.getData(ii, &size, &dataToCopy);
    assert(size==nPoint);
    std::copy(dataToCopy, dataToCopy+size, dataContainer[ii]);
  }
#endif  //HAVE_MPI
}

CFlexInterfaceData::~CFlexInterfaceData(){

#ifndef NDEBUG
  cout << "Calling CFlexInterfaceData::~CFlexInterfaceData()" << endl;
#endif  //NDEBUG

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    if(dataContainer[ii]){
      VecDestroy(&(dataContainer[ii]));
    }
  }
#else //HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    if(dataContainer[ii] != NULL) delete [] dataContainer[ii];
  }
#endif  //HAVE_MPI
}

void CFlexInterfaceData::destroy(){

#ifndef NDEBUG
  cout << "Calling CFlexInterfaceData::destroy()" << endl;
#endif  //NDEBUG

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    if(dataContainer[ii]){
      VecDestroy(&(dataContainer[ii]));
    }
  }
#else //HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    if(dataContainer[ii] != NULL) delete [] dataContainer[ii];
  }
#endif  //HAVE_MPI

}

void CFlexInterfaceData::view(const int &iDim){

#ifdef HAVE_MPI
  VecView(dataContainer[iDim], PETSC_VIEWER_STDOUT_(comm));
#else  //HAVE_MPI
  for(int jj=0; jj<nPoint; jj++){
    cout << dataContainer[iDim][jj] << endl;
  }
#endif  //HAVE_MPI
}

int CFlexInterfaceData::getnPoint(){

  return nPoint;
}

int CFlexInterfaceData::getDim(){

  return nDim;
}

Cupydo_Comm CFlexInterfaceData::getComm(){

  return comm;
}

int CFlexInterfaceData::getLocalLength(){

  int localLength;

#ifdef HAVE_MPI
  VecGetLocalSize(dataContainer[0], &localLength);
#else  //HAVE_MPI
  localLength = nPoint;
#endif  //HAVE_MPI

  return localLength;
}

void CFlexInterfaceData::setValue(const int& iDim, const int& index, const double& value){

#ifdef HAVE_MPI
  VecSetValue(dataContainer[iDim], index, value, INSERT_VALUES);
#else
  dataContainer[iDim][index] = value;
#endif

}

void CFlexInterfaceData::setValues(const int& iDim, int size_indices, int *indices_list, int size_values, double *values_array){

  assert(size_indices == size_values);

#ifdef HAVE_MPI
  VecSetValues(dataContainer[iDim], size_indices, indices_list, values_array, INSERT_VALUES);
#else //HAVE_MPI
  for(int ii=0; ii<size_indices; ii++){
    dataContainer[iDim][indices_list[ii]] = values_array[ii];
  }
#endif //HAVE_MPI

}

void CFlexInterfaceData::setAllValues(const int& iDim, const double& value){

#ifdef HAVE_MPI
  VecSet(dataContainer[iDim], value);
#else //HAVE_MPI
  fill(dataContainer[iDim], dataContainer[iDim]+nPoint, value);
#endif //HAVE_MPI
}

#ifdef HAVE_MPI
  Vec CFlexInterfaceData::getData(const int& iDim){
    return dataContainer[iDim];
  }
#else //HAVE_MPI
  void CFlexInterfaceData::getData(const int& iDim, int* size, double** data_array){
    *size = nPoint;
    *data_array = dataContainer[iDim];
  }

  void CFlexInterfaceData::setData(const int& iDim, int size, double* data){
	assert(size == nPoint);
	std::copy(data, data+nPoint, dataContainer[iDim]);
  }

#endif //HAVE_MPI

void CFlexInterfaceData::getDataArray(const int& iDim, int* size, double** data_array){



#ifdef HAVE_MPI
  VecGetLocalSize(dataContainer[iDim],size);
  VecGetArray(dataContainer[iDim], &(*data_array));
#else //HAVE_MPI
  *size = nPoint;
  *data_array = dataContainer[iDim];
#endif

}

void CFlexInterfaceData::assemble(){

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecAssemblyBegin(dataContainer[ii]);
    VecAssemblyEnd(dataContainer[ii]);
  }
#endif //HAVE_MPI

}

vector<double> CFlexInterfaceData::norm(){

    vector<double> norm_list(nDim);
#ifdef HAVE_MPI

  for(int ii=0; ii<nDim; ii++){
    VecNorm(dataContainer[ii], NORM_2, &(norm_list[ii]));
  }
#endif //HAVE_MPI
    return norm_list;
}

vector<double> CFlexInterfaceData::sum(){

    vector<double> sum_list(nDim);
#ifdef HAVE_MPI

  for(int ii=0; ii<nDim; ii++){
    VecSum(dataContainer[ii], &(sum_list[ii]));
  }
#endif  //HAVE_MPI
    return sum_list;
}

void CFlexInterfaceData::copy(CFlexInterfaceData& target){

  assert(nPoint == target.nPoint);
  assert(nDim == target.nDim);

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecCopy(dataContainer[ii], target.getData(ii));
  }
#else //HAVE_MPI
  int dataContainerTargetSize;
  double* dataContainerTarget;
  for(int ii=0; ii<nDim; ii++){
    target.getData(ii, &dataContainerTargetSize, &dataContainerTarget);
    assert(dataContainerTargetSize == nPoint);
    std::copy(dataContainer[ii], dataContainer[ii]+nPoint, dataContainerTarget);
  }
#endif //HAVE_MPI

}

void CFlexInterfaceData::set(CFlexInterfaceData& donor){

  assert(nPoint == donor.nPoint);
  assert(nDim == donor.nDim);

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecCopy(donor.getData(ii), dataContainer[ii]);
  }
#else  //HAVE_MPI
  double* dataToCopy;
  int size;
  for(int ii=0; ii<nDim; ii++){
    donor.getData(ii, &size, &dataToCopy);
    assert(size==nPoint);
    std::copy(dataToCopy, dataToCopy+size, dataContainer[ii]);
  }
#endif  //HAVE_MPI

}

vector<double> CFlexInterfaceData::dot(CFlexInterfaceData& data){

    vector<double> dot_list(nDim);
#ifdef HAVE_MPI

  for(int ii=0; ii<nDim; ii++){
    VecDot(dataContainer[ii],data.getData(ii),&(dot_list[ii]));
  }
#endif  //HAVE_MPI
    return dot_list;
}

vector<int> CFlexInterfaceData::getOwnershipRange() const{

  vector<int> ownershipRange(2);
#ifdef HAVE_MPI
  VecGetOwnershipRange(dataContainer[0], &(ownershipRange[0]), &(ownershipRange[1]));
#else  //HAVE_MPI
  ownershipRange[0] = 0;
  ownershipRange[1] = nPoint;
#endif  //HAVE_MPI

  return ownershipRange;
}

void CFlexInterfaceData::add(CFlexInterfaceData& data){

#ifndef NDEBUG
  cout << "Calling CFlexInterfaceData::add()" << endl;
#endif  //NDEBUG

  assert(nPoint == data.nPoint);
  assert(nDim == data.nDim);


#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecAXPY(dataContainer[ii],1.0,data.getData(ii));
  }
#else  //HAVE_MPI
  double* dataToAdd;
  int size;
  for(int ii=0; ii<nDim; ii++){
    data.getData(ii, &size, &dataToAdd);
    assert(nPoint==size);
    for(int jj=0; jj<nPoint; jj++){
      dataContainer[ii][jj] += dataToAdd[jj];
    }
  }
#endif  //HAVE_MPI

}

void CFlexInterfaceData::add(const double & scalar){

#ifndef NDEBUG
  cout << "Calling CFlexInterfaceData::add()" << endl;
#endif  //NDEBUG

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecShift(dataContainer[ii], scalar);
  }
#else  //HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    for(int jj=0; jj<nPoint; jj++){
      dataContainer[ii][jj] += scalar;
    }
  }
#endif  //HAVE_MPI

}

void CFlexInterfaceData::add(const int & scalar){

#ifndef NDEBUG
  cout << "Calling CFlexInterfaceData::add()" << endl;
#endif  //NDEBUG

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecShift(dataContainer[ii], scalar);
  }
#else  //HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    for(int jj=0; jj<nPoint; jj++){
      dataContainer[ii][jj] += scalar;
    }
  }
#endif  //HAVE_MPI

}

void CFlexInterfaceData::sub(CFlexInterfaceData& data){

#ifndef NDEBUG
  cout << "Calling CFlexInterfaceData::sub()" << endl;
#endif  //NDEBUG

  assert(nPoint == data.nPoint);
  assert(nDim == data.nDim);

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecAXPY(dataContainer[ii],-1.0,data.getData(ii));
  }
#else  //HAVE_MPI
  double* dataToAdd;
  int size;
  for(int ii=0; ii<nDim; ii++){
    data.getData(ii, &size, &dataToAdd);
    assert(nPoint==size);
    for(int jj=0; jj<nPoint; jj++){
      dataContainer[ii][jj] -= dataToAdd[jj];
    }
  }
#endif  //HAVE_MPI

}

void CFlexInterfaceData::sub(const double & scalar){

#ifndef NDEBUG
  cout << "Calling CFlexInterfaceData::sub()" << endl;
#endif  //NDEBUG

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecShift(dataContainer[ii], -scalar);
  }
#else  //HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    for(int jj=0; jj<nPoint; jj++){
      dataContainer[ii][jj] -= scalar;
    }
  }
#endif  //HAVE_MPI

}

void CFlexInterfaceData::sub(const int & scalar){

#ifndef NDEBUG
  cout << "Calling CFlexInterfaceData::sub()" << endl;
#endif  //NDEBUG

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecShift(dataContainer[ii], -scalar);
  }
#else  //HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    for(int jj=0; jj<nPoint; jj++){
      dataContainer[ii][jj] -= scalar;
    }
  }
#endif  //HAVE_MPI

}

void CFlexInterfaceData::scale(const double& value){

#ifndef NDEBUG
  cout << "Calling CFlexInterfaceData::scale()" << endl;
#endif  //NDEBUG

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecScale(dataContainer[ii], value);
  }
#else  //HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    for(int jj=0; jj<nPoint; jj++){
      dataContainer[ii][jj] *= value;
    }
  }
#endif  //HAVE_MPI

}

/*CFlexInterfaceData & CFlexInterfaceData::operator=(CFlexInterfaceData& data){

  cout << "Calling CFlexInterfaceData::operator=()" << endl;
  if(this == &data) return *this;

  data.copy(*this);
}*/

/*CFlexInterfaceData & CFlexInterfaceData::operator+=(CFlexInterfaceData& data){

  cout << "Calling CFlexInterfaceData::operator+=()" << endl;
  assert(nPoint == data.nPoint);
  assert(nDim == data.nDim);

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecAXPY(dataContainer[ii],1.0,data.getData(ii));
  }
#else  //HAVE_MPI
  double* dataToAdd;
  int size;
  for(int ii=0; ii<nDim; ii+){
    data.getData(ii, &size, &dataToAdd);
    assert(nPoint==size);
    for(int jj=0; jj<nPoint; jj++){
      dataContainer[ii][jj] += dataToAdd[jj];
    }
  }
#endif  //HAVE_MPI

  return *this;

}*/
