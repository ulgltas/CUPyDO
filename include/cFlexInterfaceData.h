/*!
 * Header for CFlexInterfaceData.
 * Authors : D. THOMAS.
 *
 * COPYRIGHT (C) University of Liège, 2017.
 */

#pragma once

#include <vector>

#ifdef HAVE_MPI
#include "petscvec.h"
#endif  //HAVE_MPI

#include "cMpi.h"

class CFlexInterfaceData{
#ifdef HAVE_MPI
  std::vector<Vec> dataContainer;
#else //HAVE_MPI
  std::vector<double*> dataContainer;
#endif //HAVE_MPI
public:
  //Public members
  CFlexInterfaceData(int const& val_nPoint, int const& val_nDim, Cupydo_Comm val_comm);
  CFlexInterfaceData(CFlexInterfaceData& data);
  virtual ~CFlexInterfaceData();
  void destroy();
  void view(const int& iDim);
  int getnPoint();
  int getDim();
  Cupydo_Comm getComm();
  int getLocalLength();
  void setValue(const int& iDim, const int& index, const double& value);
  void setValues(const int& iDim, int size_indices, int *indices_list, int size_values, double *values_array);
  void setAllValues(const int& iDim, const double& value);
  //void getDataContainer();
#ifdef HAVE_MPI
  Vec getData(const int& iDim);
#else //HAVE_MPI
  void getData(const int& iDim, int* size, double** data);
#endif //HAVE_MPI
  void getDataArray(const int& iDim, int* size, double** data_array);
  void assemble();
  std::vector<double> norm();
  std::vector<double> sum();
  void copy(CFlexInterfaceData& target);
  void set(CFlexInterfaceData& donor);
  std::vector<double> dot(CFlexInterfaceData& data);
  std::vector<int> getOwnershipRange() const;
  void add(CFlexInterfaceData& data);
  void add(const double & scalar);
  void add(const int & scalar);
  void sub(CFlexInterfaceData& data);
  void sub(const double & scalar);
  void sub(const int & scalar);
  void scale(const double& value);
  //CFlexInterfaceData & operator=(CFlexInterfaceData& data);
  //CFlexInterfaceData & operator+=(CFlexInterfaceData& data);
  //Public attributes
  int nPoint, nDim;
  Cupydo_Comm comm;
};
