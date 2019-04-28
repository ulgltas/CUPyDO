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
 * Header for CInterpolator.
 * Authors : D. THOMAS.
 */

#ifndef CINTERPOLATOR_H
#define CINTERPOLATOR_H

#include "cManager.h"
#include "cInterfaceMatrix.h"

class CInterpolator
{
    CManager *manager;
    double *minDist;
    int *jGlobalVertexSolid_array;

public:
    CInterpolator(CManager *val_manager);

    virtual ~CInterpolator();

    void matching_initSearch();

    void matching_search(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                         int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                         int iProc) const;

    void matching_fillMatrix(CInterfaceMatrix *H, CInterfaceMatrix *H_T) const;

    void TPS_fillMatrixA(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                         int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                         CInterfaceMatrix *A, CInterfaceMatrix *A_T,
                         int iProc) const;

    void TPS_fillMatrixB(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                         int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                         CInterfaceMatrix *B, CInterfaceMatrix *B_T,
                         int iProc) const;

    void consistent_TPS_fillMatrixA(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                    int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                    CInterfaceMatrix *A,
                                    int iProc) const;

    void consistent_TPS_fillMatrixBD(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                     int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                     CInterfaceMatrix *B, CInterfaceMatrix *D,
                                     int iProc) const;

    void consistent_TPS_fillMatrixC(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                    int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                    CInterfaceMatrix *C,
                                    int iProc) const;

    void RBF_fillMatrixA(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                         int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                         CInterfaceMatrix *A, CInterfaceMatrix *A_T,
                         int iProc, double const &radius) const;

    void RBF_fillMatrixB(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                         int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                         CInterfaceMatrix *B, CInterfaceMatrix *B_T,
                         int iProc, double const &radius) const;

    void consistent_RBF_fillMatrixA(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                    int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                    CInterfaceMatrix *A,
                                    int iProc, double const &radius) const;

    void consistent_RBF_fillMatrixBD(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                     int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                     CInterfaceMatrix *B, CInterfaceMatrix *D,
                                     int iProc, double const &radius) const;

    void consistent_RBF_fillMatrixC(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                    int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                    CInterfaceMatrix *C,
                                    int iProc, double const &radius) const;

    double PHI_TPS(double &distance) const;

    double PHI_RBF(double &distance, double const &radius) const;

    double distance(int val_size1, double *array1, int val_size2, double *array2) const;

    int nf_loc, ns_loc;
    int ns, nf;
    int nDim;
    int myid;
};

#endif //CINTERPOLATOR_H
