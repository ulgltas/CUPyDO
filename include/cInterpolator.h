/*!
 * Header for CInterpolator.
 * Authors : D. THOMAS.
 *
 * COPYRIGHT (C) University of Liège, 2017.
 */

#pragma once

#include "cManager.h"
#include "cInterfaceMatrix.h"

class CInterpolator{
CManager* manager;
public:
  CInterpolator(CManager* val_manager);
  virtual ~CInterpolator();
  void matching_fillMatrix(int size_loc_x, double* array_loc_x, int size_loc_y, double* array_loc_y, int size_loc_z, double* array_loc_z,
                       int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                       CInterfaceMatrix *H, CInterfaceMatrix *H_T,
                       int iProc) const;
  void TPS_fillMatrixA(int size_loc_x, double* array_loc_x, int size_loc_y, double* array_loc_y, int size_loc_z, double* array_loc_z,
                       int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                       CInterfaceMatrix *A, CInterfaceMatrix *A_T,
                       int iProc) const;
  void TPS_fillMatrixB(int size_loc_x, double* array_loc_x, int size_loc_y, double* array_loc_y, int size_loc_z, double* array_loc_z,
                       int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                       CInterfaceMatrix *B, CInterfaceMatrix *B_T,
                       int iProc) const;
  void RBF_fillMatrixA(int size_loc_x, double* array_loc_x, int size_loc_y, double* array_loc_y, int size_loc_z, double* array_loc_z,
                       int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                       CInterfaceMatrix *A, CInterfaceMatrix *A_T,
                       int iProc, double const& radius) const;
  void RBF_fillMatrixB(int size_loc_x, double* array_loc_x, int size_loc_y, double* array_loc_y, int size_loc_z, double* array_loc_z,
                       int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                       CInterfaceMatrix *B, CInterfaceMatrix *B_T,
                       int iProc, double const& radius) const;

  double PHI_TPS(double& distance) const;
  double PHI_RBF(double& distance, double const& radius) const;
  double distance(int val_size1, double* array1, int val_size2, double* array2) const;

  int nf_loc, ns_loc;
  int ns, nf;
  int nDim;
  int myid;
};
