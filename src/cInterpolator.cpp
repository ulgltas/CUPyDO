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
 * Source for CInterpolator.
 * Authors : D. THOMAS.
 */

#include <iostream>
#include <cassert>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#include "cInterpolator.h"
#include "cInterfaceMatrix.h"
#include "cAdt.h"

using namespace std;

CInterpolator::CInterpolator(CManager *val_manager) : manager(val_manager)
{

    ns = 0;
    nf = 0;
    ns_loc = 0;
    nf_loc = 0;

    minDist = nullptr;
    jGlobalVertexSolid_array = nullptr;
}

CInterpolator::~CInterpolator()
{

#ifndef NDEBUG
    cout << "Calling CInterpolator::~CInterpolator()" << endl;
#endif //NDEBUG

    if (minDist != nullptr)
        delete[] minDist;
    if (jGlobalVertexSolid_array != nullptr)
        delete[] jGlobalVertexSolid_array;
}

void CInterpolator::matching_initSearch()
{

    minDist = new double[nf_loc];
    fill(minDist, minDist + nf_loc, 1E6);

    jGlobalVertexSolid_array = new int[nf_loc];
    fill(jGlobalVertexSolid_array, jGlobalVertexSolid_array + nf_loc, 0);
}

void CInterpolator::matching_search(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                    int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                    int iProc) const
{

    double fluidPoint[3] = {0.0, 0.0, 0.0};
    double dist;
    int jVertex(100);

    assert(nf_loc == size_loc_x);
    assert(nf_loc == size_loc_y);
    assert(nf_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    ADTPoint ADT(size_buff_x, buff_x, size_buff_y, buff_y, size_buff_z, buff_z);
    for (int iVertex = 0; iVertex < nf_loc; iVertex++)
    {
        fluidPoint[0] = array_loc_x[iVertex];
        fluidPoint[1] = array_loc_y[iVertex];
        fluidPoint[2] = array_loc_z[iVertex];
        ADT.queryNN(3, fluidPoint, jVertex, dist);
        if (dist < minDist[iVertex])
        {
            minDist[iVertex] = dist;
            jGlobalVertexSolid_array[iVertex] = manager->getGlobalIndex("solid", iProc, jVertex);
        }
    }
}

void CInterpolator::matching_fillMatrix(CInterfaceMatrix *H, CInterfaceMatrix *H_T) const
{

    int iGlobalVertexFluid, jGlobalVertexSolid;

    for (int iVertex = 0; iVertex < nf_loc; iVertex++)
    {
        iGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, iVertex);
        jGlobalVertexSolid = jGlobalVertexSolid_array[iVertex];
        if (minDist[iVertex] > 1e-6)
            cout << "WARNING : Tolerance for matching meshes is not matched between node F" << iGlobalVertexFluid << " and S" << jGlobalVertexSolid << " DISTANCE : " << minDist[iVertex] << " !" << endl;
        H->setValue(iGlobalVertexFluid, jGlobalVertexSolid, 1.0);
        H_T->setValue(jGlobalVertexSolid, iGlobalVertexFluid, 1.0);
    }
}

void CInterpolator::TPS_fillMatrixA(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                    int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                    CInterfaceMatrix *A, CInterfaceMatrix *A_T,
                                    int iProc) const
{

    double solidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0, 0.0, 0.0};
    double phi, dist;
    int iGlobalVertexSolid, jGlobalVertexSolid;

    assert(ns_loc == size_loc_x);
    assert(ns_loc == size_loc_y);
    assert(ns_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    vector<int> jGlobalVertexSolid_list(size_buff_x);
    vector<int> iGlobalVertexSolid_list(1);
    vector<double> phi_value(size_buff_x);

    for (int iVertex = 0; iVertex < ns_loc; iVertex++)
    {
        solidPoint[0] = array_loc_x[iVertex];
        solidPoint[1] = array_loc_y[iVertex];
        solidPoint[2] = array_loc_z[iVertex];
        iGlobalVertexSolid = manager->getGlobalIndex("solid", myid, iVertex);
        iGlobalVertexSolid_list[0] = iGlobalVertexSolid;
        for (int jVertex = 0; jVertex < size_buff_x; jVertex++)
        {
            jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
            jGlobalVertexSolid_list[jVertex] = jGlobalVertexSolid;
            solidQuery[0] = buff_x[jVertex];
            solidQuery[1] = buff_y[jVertex];
            solidQuery[2] = buff_z[jVertex];
            dist = distance(3, solidPoint, 3, solidQuery);
            phi = PHI_TPS(dist);
            phi_value[jVertex] = phi;
        }
        //Set PHI block
        A->setValues(1, &(iGlobalVertexSolid_list.front()), static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), &(phi_value.front()));   //set the entire row
        A_T->setValues(static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), 1, &(iGlobalVertexSolid_list.front()), &(phi_value.front())); //set the entire column of the transposed matrix
        //Set P block
        A->setValue(iGlobalVertexSolid, ns, 1.0);
        A->setValue(iGlobalVertexSolid, ns + 1, solidPoint[0]);
        A->setValue(iGlobalVertexSolid, ns + 2, solidPoint[1]);
        A_T->setValue(ns, iGlobalVertexSolid, 1.0);
        A_T->setValue(ns + 1, iGlobalVertexSolid, solidPoint[0]);
        A_T->setValue(ns + 2, iGlobalVertexSolid, solidPoint[1]);
        //Set P^T block
        A->setValue(ns, iGlobalVertexSolid, 1.0);
        A->setValue(ns + 1, iGlobalVertexSolid, solidPoint[0]);
        A->setValue(ns + 2, iGlobalVertexSolid, solidPoint[1]);
        A_T->setValue(iGlobalVertexSolid, ns, 1.0);
        A_T->setValue(iGlobalVertexSolid, ns + 1, solidPoint[0]);
        A_T->setValue(iGlobalVertexSolid, ns + 2, solidPoint[1]);
        if (nDim == 3)
        {
            A->setValue(iGlobalVertexSolid, ns + 3, solidPoint[2]);
            A->setValue(ns + 3, iGlobalVertexSolid, solidPoint[2]);
            A_T->setValue(ns + 3, iGlobalVertexSolid, solidPoint[2]);
            A_T->setValue(iGlobalVertexSolid, ns + 3, solidPoint[2]);
        }
    }
}

void CInterpolator::TPS_fillMatrixB(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                    int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                    CInterfaceMatrix *B, CInterfaceMatrix *B_T,
                                    int iProc) const
{

    double fluidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0, 0.0, 0.0};
    double phi, dist;
    int iGlobalVertexFluid, jGlobalVertexSolid;

    assert(nf_loc == size_loc_x);
    assert(nf_loc == size_loc_y);
    assert(nf_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    vector<int> jGlobalVertexSolid_list(size_buff_x);
    vector<int> iGlobalVertexFluid_list(1);
    vector<double> phi_value(size_buff_x);

    for (int iVertex = 0; iVertex < nf_loc; iVertex++)
    {
        fluidPoint[0] = array_loc_x[iVertex];
        fluidPoint[1] = array_loc_y[iVertex];
        fluidPoint[2] = array_loc_z[iVertex];
        iGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, iVertex);
        iGlobalVertexFluid_list[0] = iGlobalVertexFluid;
        for (int jVertex = 0; jVertex < size_buff_x; jVertex++)
        {
            jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
            jGlobalVertexSolid_list[jVertex] = jGlobalVertexSolid;
            solidQuery[0] = buff_x[jVertex];
            solidQuery[1] = buff_y[jVertex];
            solidQuery[2] = buff_z[jVertex];
            dist = distance(3, fluidPoint, 3, solidQuery);
            phi = PHI_TPS(dist);
            phi_value[jVertex] = phi;
        }
        //Set PHI block
        B->setValues(1, &(iGlobalVertexFluid_list.front()), static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), &(phi_value.front()));   //set the entire row
        B_T->setValues(static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), 1, &(iGlobalVertexFluid_list.front()), &(phi_value.front())); //set the entire column of the transposed matrix
        //Set P block
        B->setValue(iGlobalVertexFluid, ns, 1.0);
        B->setValue(iGlobalVertexFluid, ns + 1, fluidPoint[0]);
        B->setValue(iGlobalVertexFluid, ns + 2, fluidPoint[1]);
        B_T->setValue(ns, iGlobalVertexFluid, 1.0);
        B_T->setValue(ns + 1, iGlobalVertexFluid, fluidPoint[0]);
        B_T->setValue(ns + 2, iGlobalVertexFluid, fluidPoint[1]);
        if (nDim == 3)
        {
            B->setValue(iGlobalVertexFluid, ns + 3, fluidPoint[2]);
            B_T->setValue(ns + 3, iGlobalVertexFluid, fluidPoint[2]);
        }
    }
}

void CInterpolator::consistent_TPS_fillMatrixA(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                               int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                               CInterfaceMatrix *A,
                                               int iProc) const
{

    double solidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0, 0.0, 0.0};
    double phi, dist;
    int iGlobalVertexSolid, jGlobalVertexSolid;

    assert(ns_loc == size_loc_x);
    assert(ns_loc == size_loc_y);
    assert(ns_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    vector<int> jGlobalVertexSolid_list(size_buff_x);
    vector<int> iGlobalVertexSolid_list(1);
    vector<double> phi_value(size_buff_x);

    for (int iVertex = 0; iVertex < ns_loc; iVertex++)
    {
        solidPoint[0] = array_loc_x[iVertex];
        solidPoint[1] = array_loc_y[iVertex];
        solidPoint[2] = array_loc_z[iVertex];
        iGlobalVertexSolid = manager->getGlobalIndex("solid", myid, iVertex);
        iGlobalVertexSolid_list[0] = iGlobalVertexSolid;
        for (int jVertex = 0; jVertex < size_buff_x; jVertex++)
        {
            jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
            jGlobalVertexSolid_list[jVertex] = jGlobalVertexSolid;
            solidQuery[0] = buff_x[jVertex];
            solidQuery[1] = buff_y[jVertex];
            solidQuery[2] = buff_z[jVertex];
            dist = distance(3, solidPoint, 3, solidQuery);
            phi = PHI_TPS(dist);
            phi_value[jVertex] = phi;
        }
        //Set PHI block
        A->setValues(1, &(iGlobalVertexSolid_list.front()), static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), &(phi_value.front())); //set the entire row
        //Set P block
        A->setValue(iGlobalVertexSolid, ns, 1.0);
        A->setValue(iGlobalVertexSolid, ns + 1, solidPoint[0]);
        A->setValue(iGlobalVertexSolid, ns + 2, solidPoint[1]);
        //Set P^T block
        A->setValue(ns, iGlobalVertexSolid, 1.0);
        A->setValue(ns + 1, iGlobalVertexSolid, solidPoint[0]);
        A->setValue(ns + 2, iGlobalVertexSolid, solidPoint[1]);
        if (nDim == 3)
        {
            A->setValue(iGlobalVertexSolid, ns + 3, solidPoint[2]);
            A->setValue(ns + 3, iGlobalVertexSolid, solidPoint[2]);
        }
    }
}

void CInterpolator::consistent_TPS_fillMatrixBD(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                                int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                                CInterfaceMatrix *B, CInterfaceMatrix *D,
                                                int iProc) const
{

    double fluidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0, 0.0, 0.0};
    double phi, dist;
    int iGlobalVertexFluid, jGlobalVertexSolid;

    assert(nf_loc == size_loc_x);
    assert(nf_loc == size_loc_y);
    assert(nf_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    vector<int> jGlobalVertexSolid_list(size_buff_x);
    vector<int> iGlobalVertexFluid_list(1);
    vector<double> phi_value(size_buff_x);

    //Build B (donor = solid, target = fluid)
    //Build D (donor = fluid, target = solid)
    for (int iVertex = 0; iVertex < nf_loc; iVertex++)
    {
        fluidPoint[0] = array_loc_x[iVertex];
        fluidPoint[1] = array_loc_y[iVertex];
        fluidPoint[2] = array_loc_z[iVertex];
        iGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, iVertex);
        iGlobalVertexFluid_list[0] = iGlobalVertexFluid;
        for (int jVertex = 0; jVertex < size_buff_x; jVertex++)
        {
            jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
            jGlobalVertexSolid_list[jVertex] = jGlobalVertexSolid;
            solidQuery[0] = buff_x[jVertex];
            solidQuery[1] = buff_y[jVertex];
            solidQuery[2] = buff_z[jVertex];
            dist = distance(3, fluidPoint, 3, solidQuery);
            phi = PHI_TPS(dist);
            phi_value[jVertex] = phi;
            D->setValue(jGlobalVertexSolid, nf, 1.0);
            D->setValue(jGlobalVertexSolid, nf + 1, solidQuery[0]);
            D->setValue(jGlobalVertexSolid, nf + 2, solidQuery[1]);
            if (nDim == 3)
                D->setValue(jGlobalVertexSolid, nf + 3, solidQuery[2]);
        }
        B->setValues(1, &(iGlobalVertexFluid_list.front()), static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), &(phi_value.front()));
        B->setValue(iGlobalVertexFluid, ns, 1.0);
        B->setValue(iGlobalVertexFluid, ns + 1, fluidPoint[0]);
        B->setValue(iGlobalVertexFluid, ns + 2, fluidPoint[1]);
        if (nDim == 3)
        {
            B->setValue(iGlobalVertexFluid, ns + 3, fluidPoint[2]);
        }
        D->setValues(static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), 1, &(iGlobalVertexFluid_list.front()), &(phi_value.front()));
    }
}

void CInterpolator::consistent_TPS_fillMatrixC(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                               int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                               CInterfaceMatrix *C,
                                               int iProc) const
{

    double fluidPoint[3] = {0.0, 0.0, 0.0}, fluidQuery[3] = {0.0, 0.0, 0.0};
    double phi, dist;
    int iGlobalVertexFluid, jGlobalVertexFluid;

    assert(nf_loc == size_loc_x);
    assert(nf_loc == size_loc_y);
    assert(nf_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    vector<int> jGlobalVertexFluid_list(size_buff_x);
    vector<int> iGlobalVertexFluid_list(1);
    vector<double> phi_value(size_buff_x);

    for (int iVertex = 0; iVertex < nf_loc; iVertex++)
    {
        fluidPoint[0] = array_loc_x[iVertex];
        fluidPoint[1] = array_loc_y[iVertex];
        fluidPoint[2] = array_loc_z[iVertex];
        iGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, iVertex);
        iGlobalVertexFluid_list[0] = iGlobalVertexFluid;
        for (int jVertex = 0; jVertex < size_buff_x; jVertex++)
        {
            jGlobalVertexFluid = manager->getGlobalIndex("fluid", iProc, jVertex);
            jGlobalVertexFluid_list[jVertex] = jGlobalVertexFluid;
            fluidQuery[0] = buff_x[jVertex];
            fluidQuery[1] = buff_y[jVertex];
            fluidQuery[2] = buff_z[jVertex];
            dist = distance(3, fluidPoint, 3, fluidQuery);
            phi = PHI_TPS(dist);
            phi_value[jVertex] = phi;
        }
        //Set block PHI
        C->setValues(1, &(iGlobalVertexFluid_list.front()), static_cast<int>(jGlobalVertexFluid_list.size()), &(jGlobalVertexFluid_list.front()), &(phi_value.front()));
        //Set block P
        C->setValue(iGlobalVertexFluid, nf, 1.0);
        C->setValue(iGlobalVertexFluid, nf + 1, fluidPoint[0]);
        C->setValue(iGlobalVertexFluid, nf + 2, fluidPoint[1]);
        //Set block P^T
        C->setValue(nf, iGlobalVertexFluid, 1.0);
        C->setValue(nf + 1, iGlobalVertexFluid, fluidPoint[0]);
        C->setValue(nf + 2, iGlobalVertexFluid, fluidPoint[1]);
        if (nDim == 3)
        {
            C->setValue(iGlobalVertexFluid, nf + 3, fluidPoint[2]);
            C->setValue(nf + 3, iGlobalVertexFluid, fluidPoint[2]);
        }
    }
}

void CInterpolator::RBF_fillMatrixA(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                    int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                    CInterfaceMatrix *A, CInterfaceMatrix *A_T,
                                    int iProc, double const &radius) const
{

    double solidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0, 0.0, 0.0};
    double phi, dist;
    int iGlobalVertexSolid, jGlobalVertexSolid;
    vector<int> solidVertices;
    vector<int> jGlobalVertexSolid_list;
    vector<int> iGlobalVertexSolid_list(1);
    vector<double> phi_value;

    assert(ns_loc == size_loc_x);
    assert(ns_loc == size_loc_y);
    assert(ns_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    ADTPoint ADT(size_buff_x, buff_x, size_buff_y, buff_y, size_buff_z, buff_z);
    for (int iVertex = 0; iVertex < ns_loc; iVertex++)
    {
        solidPoint[0] = array_loc_x[iVertex];
        solidPoint[1] = array_loc_y[iVertex];
        solidPoint[2] = array_loc_z[iVertex];
        iGlobalVertexSolid = manager->getGlobalIndex("solid", myid, iVertex);
        iGlobalVertexSolid_list[0] = iGlobalVertexSolid;
        jGlobalVertexSolid_list.clear();
        phi_value.clear();
        ADT.queryBallNN(3, solidPoint, radius, solidVertices);
        for (int jVertex : solidVertices)
        {
            jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
            jGlobalVertexSolid_list.push_back(jGlobalVertexSolid);
            solidQuery[0] = buff_x[jVertex];
            solidQuery[1] = buff_y[jVertex];
            solidQuery[2] = buff_z[jVertex];
            dist = distance(3, solidPoint, 3, solidQuery);
            phi = PHI_RBF(dist, radius);
            phi_value.push_back(phi);
        }
        //Set block PHI
        A->setValues(1, &(iGlobalVertexSolid_list.front()), static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), &(phi_value.front()));
        A_T->setValues(static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), 1, &(iGlobalVertexSolid_list.front()), &(phi_value.front()));
        //set block P
        A->setValue(iGlobalVertexSolid, ns, 1.0);
        A->setValue(iGlobalVertexSolid, ns + 1, solidPoint[0]);
        A->setValue(iGlobalVertexSolid, ns + 2, solidPoint[1]);
        A_T->setValue(ns, iGlobalVertexSolid, 1.0);
        A_T->setValue(ns + 1, iGlobalVertexSolid, solidPoint[0]);
        A_T->setValue(ns + 2, iGlobalVertexSolid, solidPoint[1]);
        //Set block P^T
        A->setValue(ns, iGlobalVertexSolid, 1.0);
        A->setValue(ns + 1, iGlobalVertexSolid, solidPoint[0]);
        A->setValue(ns + 2, iGlobalVertexSolid, solidPoint[1]);
        A_T->setValue(iGlobalVertexSolid, ns, 1.0);
        A_T->setValue(iGlobalVertexSolid, ns + 1, solidPoint[0]);
        A_T->setValue(iGlobalVertexSolid, ns + 2, solidPoint[1]);
        if (nDim == 3)
        {
            A->setValue(iGlobalVertexSolid, ns + 3, solidPoint[2]);
            A->setValue(ns + 3, iGlobalVertexSolid, solidPoint[2]);
            A_T->setValue(ns + 3, iGlobalVertexSolid, solidPoint[2]);
            A_T->setValue(iGlobalVertexSolid, ns + 3, solidPoint[2]);
        }
    }
}

void CInterpolator::RBF_fillMatrixB(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                    int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                    CInterfaceMatrix *B, CInterfaceMatrix *B_T,
                                    int iProc, double const &radius) const
{

    double fluidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0, 0.0, 0.0};
    double phi, dist;
    int iGlobalVertexFluid, jGlobalVertexSolid;
    vector<int> solidVertices;
    vector<int> jGlobalVertexSolid_list;
    vector<double> phi_value;
    vector<int> iGlobalVertexFluid_list(1);

    assert(nf_loc == size_loc_x);
    assert(nf_loc == size_loc_y);
    assert(nf_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    ADTPoint ADT(size_buff_x, buff_x, size_buff_y, buff_y, size_buff_z, buff_z);
    for (int iVertex = 0; iVertex < nf_loc; iVertex++)
    {
        fluidPoint[0] = array_loc_x[iVertex];
        fluidPoint[1] = array_loc_y[iVertex];
        fluidPoint[2] = array_loc_z[iVertex];
        iGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, iVertex);
        iGlobalVertexFluid_list[0] = iGlobalVertexFluid;
        jGlobalVertexSolid_list.clear();
        phi_value.clear();
        ADT.queryBallNN(3, fluidPoint, radius, solidVertices);
        for (int jVertex : solidVertices)
        {
            jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
            jGlobalVertexSolid_list.push_back(jGlobalVertexSolid);
            solidQuery[0] = buff_x[jVertex];
            solidQuery[1] = buff_y[jVertex];
            solidQuery[2] = buff_z[jVertex];
            dist = distance(3, fluidPoint, 3, solidQuery);
            phi = PHI_RBF(dist, radius);
            phi_value.push_back(phi);
        }
        B->setValues(1, &(iGlobalVertexFluid_list.front()), static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), &(phi_value.front()));
        B_T->setValues(static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), 1, &(iGlobalVertexFluid_list.front()), &(phi_value.front()));
        B->setValue(iGlobalVertexFluid, ns, 1.0);
        B->setValue(iGlobalVertexFluid, ns + 1, fluidPoint[0]);
        B->setValue(iGlobalVertexFluid, ns + 2, fluidPoint[1]);
        B_T->setValue(ns, iGlobalVertexFluid, 1.0);
        B_T->setValue(ns + 1, iGlobalVertexFluid, fluidPoint[0]);
        B_T->setValue(ns + 2, iGlobalVertexFluid, fluidPoint[1]);
        if (nDim == 3)
        {
            B->setValue(iGlobalVertexFluid, ns + 3, fluidPoint[2]);
            B_T->setValue(ns + 3, iGlobalVertexFluid, fluidPoint[2]);
        }
    }
}

void CInterpolator::consistent_RBF_fillMatrixA(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                               int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                               CInterfaceMatrix *A,
                                               int iProc, double const &radius) const
{

    double solidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0, 0.0, 0.0};
    double phi, dist;
    int iGlobalVertexSolid, jGlobalVertexSolid;
    vector<int> solidVertices;
    vector<int> jGlobalVertexSolid_list;
    vector<int> iGlobalVertexSolid_list(1);
    vector<double> phi_value;

    assert(ns_loc == size_loc_x);
    assert(ns_loc == size_loc_y);
    assert(ns_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    ADTPoint ADT(size_buff_x, buff_x, size_buff_y, buff_y, size_buff_z, buff_z);
    for (int iVertex = 0; iVertex < ns_loc; iVertex++)
    {
        solidPoint[0] = array_loc_x[iVertex];
        solidPoint[1] = array_loc_y[iVertex];
        solidPoint[2] = array_loc_z[iVertex];
        iGlobalVertexSolid = manager->getGlobalIndex("solid", myid, iVertex);
        iGlobalVertexSolid_list[0] = iGlobalVertexSolid;
        jGlobalVertexSolid_list.clear();
        phi_value.clear();
        ADT.queryBallNN(3, solidPoint, radius, solidVertices);
        for (int jVertex : solidVertices)
        {
            jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
            jGlobalVertexSolid_list.push_back(jGlobalVertexSolid);
            solidQuery[0] = buff_x[jVertex];
            solidQuery[1] = buff_y[jVertex];
            solidQuery[2] = buff_z[jVertex];
            dist = distance(3, solidPoint, 3, solidQuery);
            phi = PHI_RBF(dist, radius);
            phi_value.push_back(phi);
        }
        //Set block PHI
        A->setValues(1, &(iGlobalVertexSolid_list.front()), static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), &(phi_value.front()));
        //Set block P
        A->setValue(iGlobalVertexSolid, ns, 1.0);
        A->setValue(iGlobalVertexSolid, ns + 1, solidPoint[0]);
        A->setValue(iGlobalVertexSolid, ns + 2, solidPoint[1]);
        //Set block P^T
        A->setValue(ns, iGlobalVertexSolid, 1.0);
        A->setValue(ns + 1, iGlobalVertexSolid, solidPoint[0]);
        A->setValue(ns + 2, iGlobalVertexSolid, solidPoint[1]);
        if (nDim == 3)
        {
            A->setValue(iGlobalVertexSolid, ns + 3, solidPoint[2]);
            A->setValue(ns + 3, iGlobalVertexSolid, solidPoint[2]);
        }
    }
}

void CInterpolator::consistent_RBF_fillMatrixBD(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                                int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                                CInterfaceMatrix *B, CInterfaceMatrix *D,
                                                int iProc, double const &radius) const
{

    double fluidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0, 0.0, 0.0};
    double solidPoint[3] = {0.0, 0.0, 0.0}, fluidQuery[3] = {0.0, 0.0, 0.0};
    double phi, dist;
    int iGlobalVertexFluid, jGlobalVertexSolid;
    int iGlobalVertexSolid, jGlobalVertexFluid;
    vector<int> solidVertices;
    vector<int> fluidVertices;

    vector<int> jGlobalVertexSolid_list;
    vector<int> iGlobalVertexFluid_list(1);
    vector<int> jGlobalVertexFluid_list;
    vector<int> iGlobalVertexSolid_list(1);
    vector<double> phi_value;

    assert(nf_loc == size_loc_x);
    assert(nf_loc == size_loc_y);
    assert(nf_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    //Build B (donor = solid, target = fluid)
    ADTPoint ADTDonor(size_buff_x, buff_x, size_buff_y, buff_y, size_buff_z, buff_z);
    for (int iVertex = 0; iVertex < nf_loc; iVertex++)
    {
        fluidPoint[0] = array_loc_x[iVertex];
        fluidPoint[1] = array_loc_y[iVertex];
        fluidPoint[2] = array_loc_z[iVertex];
        iGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, iVertex);
        iGlobalVertexFluid_list[0] = iGlobalVertexFluid;
        jGlobalVertexSolid_list.clear();
        phi_value.clear();
        ADTDonor.queryBallNN(3, fluidPoint, radius, solidVertices);
        for (int jVertex : solidVertices)
        {
            jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
            jGlobalVertexSolid_list.push_back(jGlobalVertexSolid);
            solidQuery[0] = buff_x[jVertex];
            solidQuery[1] = buff_y[jVertex];
            solidQuery[2] = buff_z[jVertex];
            dist = distance(3, fluidPoint, 3, solidQuery);
            phi = PHI_RBF(dist, radius);
            phi_value.push_back(phi);
        }
        B->setValues(1, &(iGlobalVertexFluid_list.front()), static_cast<int>(jGlobalVertexSolid_list.size()), &(jGlobalVertexSolid_list.front()), &(phi_value.front()));
        B->setValue(iGlobalVertexFluid, ns, 1.0);
        B->setValue(iGlobalVertexFluid, ns + 1, fluidPoint[0]);
        B->setValue(iGlobalVertexFluid, ns + 2, fluidPoint[1]);
        if (nDim == 3)
        {
            B->setValue(iGlobalVertexFluid, ns + 3, fluidPoint[2]);
        }
    }

    //Build D (donor = fluid, target = solid)
    ADTPoint ADTTarget(size_loc_x, array_loc_x, size_loc_y, array_loc_y, size_loc_z, array_loc_z);
    for (int iVertex = 0; iVertex < size_buff_x; iVertex++)
    {
        solidPoint[0] = buff_x[iVertex];
        solidPoint[1] = buff_y[iVertex];
        solidPoint[2] = buff_z[iVertex];
        iGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, iVertex);
        iGlobalVertexSolid_list[0] = iGlobalVertexSolid;
        jGlobalVertexFluid_list.clear();
        phi_value.clear();
        ADTTarget.queryBallNN(3, solidPoint, radius, fluidVertices);
        for (int jVertex : fluidVertices)
        {
            jGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, jVertex);
            jGlobalVertexFluid_list.push_back(jGlobalVertexFluid);
            fluidQuery[0] = array_loc_x[jVertex];
            fluidQuery[1] = array_loc_y[jVertex];
            fluidQuery[2] = array_loc_z[jVertex];
            dist = distance(3, solidPoint, 3, fluidQuery);
            phi = PHI_RBF(dist, radius);
            phi_value.push_back(phi);
        }
        D->setValues(1, &(iGlobalVertexSolid_list.front()), static_cast<int>(jGlobalVertexFluid_list.size()), &(jGlobalVertexFluid_list.front()), &(phi_value.front()));
        D->setValue(iGlobalVertexSolid, nf, 1.0);
        D->setValue(iGlobalVertexSolid, nf + 1, solidPoint[0]);
        D->setValue(iGlobalVertexSolid, nf + 2, solidPoint[1]);
        if (nDim == 3)
        {
            D->setValue(iGlobalVertexSolid, nf + 3, solidPoint[2]);
        }
    }
}

void CInterpolator::consistent_RBF_fillMatrixC(int size_loc_x, double *array_loc_x, int size_loc_y, double *array_loc_y, int size_loc_z, double *array_loc_z,
                                               int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                               CInterfaceMatrix *C,
                                               int iProc, double const &radius) const
{

    double fluidPoint[3] = {0.0, 0.0, 0.0}, fluidQuery[3] = {0.0, 0.0, 0.0};
    double phi, dist;
    int iGlobalVertexFluid, jGlobalVertexFluid;
    vector<int> fluidVertices;
    vector<int> iGlobalVertexFluid_list(1);
    vector<int> jGlobalVertexFluid_list;
    vector<double> phi_value;

    assert(nf_loc == size_loc_x);
    assert(nf_loc == size_loc_y);
    assert(nf_loc == size_loc_z);

    assert(size_buff_x == size_buff_y);
    assert(size_buff_y == size_buff_z);
    assert(size_buff_x == size_buff_z);

    ADTPoint ADT(size_buff_x, buff_x, size_buff_y, buff_y, size_buff_z, buff_z);
    for (int iVertex = 0; iVertex < nf_loc; iVertex++)
    {
        fluidPoint[0] = array_loc_x[iVertex];
        fluidPoint[1] = array_loc_y[iVertex];
        fluidPoint[2] = array_loc_z[iVertex];
        iGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, iVertex);
        iGlobalVertexFluid_list[0] = iGlobalVertexFluid;
        jGlobalVertexFluid_list.clear();
        phi_value.clear();
        ADT.queryBallNN(3, fluidPoint, radius, fluidVertices);
        for (int jVertex : fluidVertices)
        {
            jGlobalVertexFluid = manager->getGlobalIndex("fluid", iProc, jVertex);
            jGlobalVertexFluid_list.push_back(jGlobalVertexFluid);
            fluidQuery[0] = buff_x[jVertex];
            fluidQuery[1] = buff_y[jVertex];
            fluidQuery[2] = buff_z[jVertex];
            dist = distance(3, fluidPoint, 3, fluidQuery);
            phi = PHI_RBF(dist, radius);
            phi_value.push_back(phi);
        }
        //Set block PHI
        C->setValues(1, &(iGlobalVertexFluid_list.front()), static_cast<int>(jGlobalVertexFluid_list.size()), &(jGlobalVertexFluid_list.front()), &(phi_value.front()));
        //Set block P
        C->setValue(iGlobalVertexFluid, nf, 1.0);
        C->setValue(iGlobalVertexFluid, nf + 1, fluidPoint[0]);
        C->setValue(iGlobalVertexFluid, nf + 2, fluidPoint[1]);
        //Set block P^T
        C->setValue(nf, iGlobalVertexFluid, 1.0);
        C->setValue(nf + 1, iGlobalVertexFluid, fluidPoint[0]);
        C->setValue(nf + 2, iGlobalVertexFluid, fluidPoint[1]);
        if (nDim == 3)
        {
            C->setValue(iGlobalVertexFluid, nf + 3, fluidPoint[2]);
            C->setValue(nf + 3, iGlobalVertexFluid, fluidPoint[2]);
        }
    }
}

double CInterpolator::PHI_TPS(double &distance) const
{

    if (distance > 0.0)
        return pow(distance, 2) * log10(distance);
    else
        return 0.0;
}

double CInterpolator::PHI_RBF(double &distance, const double &radius) const
{

    double eps(distance / radius);

    if (eps < 1.0)
        return pow((1.0 - eps), 4) * (4.0 * eps + 1.0);
    else
        return 0.0;
}

double CInterpolator::distance(int val_size1, double *array1, int val_size2, double *array2) const
{

    assert(val_size1 <= 3);
    assert(val_size2 <= 3);

    return sqrt(pow(array2[0] - array1[0], 2) + pow(array2[1] - array1[1], 2) + pow(array2[2] - array1[2], 2));
}
