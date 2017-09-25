/*!
 * Source for CInterpolator.
 * Authors : D. THOMAS.
 *
 * COPYRIHGHT (C) University of Liège, 2017.
 */

#include <iostream>
#include <cassert>
#include <cmath>
#include <ctime>
#include <vector>

#include "../include/cInterpolator.h"
#include "../include/cInterfaceMatrix.h"
#include "../include/cAdt.h"

using namespace std;

CInterpolator::CInterpolator(CManager *val_manager):ns(0),nf(0),ns_loc(0),nf_loc(0),manager(val_manager){}

CInterpolator::~CInterpolator(){}

void CInterpolator::matching_fillMatrix(int size_loc_x, double* array_loc_x, int size_loc_y, double* array_loc_y, int size_loc_z, double* array_loc_z,
                       int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                       CInterfaceMatrix *H, CInterfaceMatrix *H_T,
                       int iProc) const {

  double fluidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0, 0.0, 0.0};
  double dist;
  int iGlobalVertexFluid, jGlobalVertexSolid, jVertex(100);

  assert(nf_loc == size_loc_x);
  assert(nf_loc == size_loc_y);
  assert(nf_loc == size_loc_z);

  assert(size_buff_x == size_buff_y);
  assert(size_buff_y == size_buff_z);
  assert(size_buff_x == size_buff_z);

  ADTPoint ADT(size_buff_x, buff_x, size_buff_y, buff_y, size_buff_z, buff_z);
  for(int iVertex=0; iVertex < nf_loc; iVertex++){
    fluidPoint[0] = array_loc_x[iVertex];
    fluidPoint[1] = array_loc_y[iVertex];
    fluidPoint[2] = array_loc_z[iVertex];
    ADT.queryNN(3, fluidPoint, jVertex, dist);
    iGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, iVertex);
    jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
    if(dist > 1e-6) cout << "WARNING : Tolerance for matching meshes is not matched between node F" << iGlobalVertexFluid << " and S" << jGlobalVertexSolid <<
    " : (" << fluidPoint[0] <<", " << fluidPoint[1] <<", " <<fluidPoint[2] << ")<-->(" << buff_x[jVertex] << ", " << buff_y[jVertex] << ", " << buff_z[jVertex] << ") , DISTANCE : " << dist << " !" << endl;
    H->setValue(iGlobalVertexFluid, jGlobalVertexSolid, 1.0);
    H_T->setValue(jGlobalVertexSolid, iGlobalVertexFluid, 1.0);
  }

}

void CInterpolator::TPS_fillMatrixA(int size_loc_x, double* array_loc_x, int size_loc_y, double* array_loc_y, int size_loc_z, double* array_loc_z,
                                    int size_buff_x, double* buff_x, int size_buff_y, double* buff_y, int size_buff_z, double* buff_z,
                                    CInterfaceMatrix* A, CInterfaceMatrix* A_T,
                                    int iProc) const {

  double solidPoint[3] = {0.0,0.0,0.0}, solidQuery[3] = {0.0,0.0,0.0};
  double phi, dist;
  int iGlobalVertexSolid, jGlobalVertexSolid;

  assert(ns_loc == size_loc_x);
  assert(ns_loc == size_loc_y);
  assert(ns_loc == size_loc_z);

  assert(size_buff_x == size_buff_y);
  assert(size_buff_y == size_buff_z);
  assert(size_buff_x == size_buff_z);

  for(int iVertex=0; iVertex<ns_loc; iVertex++){
    solidPoint[0] = array_loc_x[iVertex];
    solidPoint[1] = array_loc_y[iVertex];
    solidPoint[2] = array_loc_z[iVertex];
    iGlobalVertexSolid = manager->getGlobalIndex("solid", myid, iVertex);
    for(int jVertex=0; jVertex<size_buff_x; jVertex++){
      jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
      solidQuery[0] = buff_x[jVertex];
      solidQuery[1] = buff_y[jVertex];
      solidQuery[2] = buff_z[jVertex];
      dist = distance(3, solidPoint, 3, solidQuery);
      phi = PHI_TPS(dist);
      A->setValue(iGlobalVertexSolid, jGlobalVertexSolid, phi);
      A_T->setValue(jGlobalVertexSolid, iGlobalVertexSolid, phi);
    }
    A->setValue(iGlobalVertexSolid, ns, 1.0);
    A->setValue(iGlobalVertexSolid, ns+1, solidPoint[0]);
    A->setValue(iGlobalVertexSolid, ns+2, solidPoint[1]);
    A_T->setValue(ns, iGlobalVertexSolid, 1.0);
    A_T->setValue(ns+1, iGlobalVertexSolid, solidPoint[0]);
    A_T->setValue(ns+2, iGlobalVertexSolid, solidPoint[1]);
    if(nDim == 3){
      A->setValue(iGlobalVertexSolid, ns+3, solidPoint[2]);
      A_T->setValue(ns+3, iGlobalVertexSolid, solidPoint[2]);
    }
  }

}

void CInterpolator::TPS_fillMatrixB(int size_loc_x, double* array_loc_x, int size_loc_y, double* array_loc_y, int size_loc_z, double* array_loc_z,
                                    int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                                    CInterfaceMatrix *B, CInterfaceMatrix *B_T,
                                    int iProc) const {

  double fluidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0,0.0,0.0};
  double phi, dist;
  int iGlobalVertexFluid, jGlobalVertexSolid;

  assert(nf_loc == size_loc_x);
  assert(nf_loc == size_loc_y);
  assert(nf_loc == size_loc_z);

  assert(size_buff_x == size_buff_y);
  assert(size_buff_y == size_buff_z);
  assert(size_buff_x == size_buff_z);

  for(int iVertex=0; iVertex < nf_loc; iVertex++){
    fluidPoint[0] = array_loc_x[iVertex];
    fluidPoint[1] = array_loc_y[iVertex];
    fluidPoint[2] = array_loc_z[iVertex];
    iGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, iVertex);
    for(int jVertex=0; jVertex < size_buff_x; jVertex++){
      jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
      solidQuery[0] = buff_x[jVertex];
      solidQuery[1] = buff_y[jVertex];
      solidQuery[2] = buff_z[jVertex];
      dist = distance(3, fluidPoint, 3, solidQuery);
      phi = PHI_TPS(dist);
      B->setValue(iGlobalVertexFluid, jGlobalVertexSolid,phi);
      B_T->setValue(jGlobalVertexSolid, iGlobalVertexFluid, phi);
    }
    B->setValue(iGlobalVertexFluid, ns, 1.0);
    B->setValue(iGlobalVertexFluid, ns+1, fluidPoint[0]);
    B->setValue(iGlobalVertexFluid, ns+2, fluidPoint[1]);
    B_T->setValue(ns, iGlobalVertexFluid, 1.0);
    B_T->setValue(ns+1, iGlobalVertexFluid, fluidPoint[0]);
    B_T->setValue(ns+2, iGlobalVertexFluid, fluidPoint[1]);
    if(nDim == 3){
      B->setValue(iGlobalVertexFluid, ns+3, fluidPoint[2]);
      B_T->setValue(ns+3, iGlobalVertexFluid, fluidPoint[2]);
    }
  }

}

void CInterpolator::RBF_fillMatrixA(int size_loc_x, double* array_loc_x, int size_loc_y, double* array_loc_y, int size_loc_z, double* array_loc_z,
                     int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                     CInterfaceMatrix *A, CInterfaceMatrix *A_T,
                     int iProc, double const& radius) const {

  double solidPoint[3] = {0.0,0.0,0.0}, solidQuery[3] = {0.0,0.0,0.0};
  double phi, dist;
  int iGlobalVertexSolid, jGlobalVertexSolid;
  vector<int> solidVertices;

  assert(ns_loc == size_loc_x);
  assert(ns_loc == size_loc_y);
  assert(ns_loc == size_loc_z);

  assert(size_buff_x == size_buff_y);
  assert(size_buff_y == size_buff_z);
  assert(size_buff_x == size_buff_z);

  ADTPoint ADT(size_buff_x, buff_x, size_buff_y, buff_y, size_buff_z, buff_z);
  for(int iVertex=0; iVertex<ns_loc; iVertex++){
    solidPoint[0] = array_loc_x[iVertex];
    solidPoint[1] = array_loc_y[iVertex];
    solidPoint[2] = array_loc_z[iVertex];
    iGlobalVertexSolid = manager->getGlobalIndex("solid", myid, iVertex);
    ADT.queryBallNN(3, solidPoint, radius, solidVertices);
    for(int jVertex : solidVertices){
      jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
      solidQuery[0] = buff_x[jVertex];
      solidQuery[1] = buff_y[jVertex];
      solidQuery[2] = buff_z[jVertex];
      dist = distance(3, solidPoint, 3, solidQuery);
      phi = PHI_RBF(dist, radius);
      A->setValue(iGlobalVertexSolid, jGlobalVertexSolid, phi);
      A_T->setValue(jGlobalVertexSolid, iGlobalVertexSolid, phi);
    }
    A->setValue(iGlobalVertexSolid, ns, 1.0);
    A->setValue(iGlobalVertexSolid, ns+1, solidPoint[0]);
    A->setValue(iGlobalVertexSolid, ns+2, solidPoint[1]);
    A_T->setValue(ns, iGlobalVertexSolid, 1.0);
    A_T->setValue(ns+1, iGlobalVertexSolid, solidPoint[0]);
    A_T->setValue(ns+2, iGlobalVertexSolid, solidPoint[1]);
    if(nDim == 3){
      A->setValue(iGlobalVertexSolid, ns+3, solidPoint[2]);
      A_T->setValue(ns+3, iGlobalVertexSolid, solidPoint[2]);
    }
  }

}

void CInterpolator::RBF_fillMatrixB(int size_loc_x, double* array_loc_x, int size_loc_y, double* array_loc_y, int size_loc_z, double* array_loc_z,
                     int size_buff_x, double *buff_x, int size_buff_y, double *buff_y, int size_buff_z, double *buff_z,
                     CInterfaceMatrix *B, CInterfaceMatrix *B_T,
                     int iProc, double const& radius) const {

  double fluidPoint[3] = {0.0, 0.0, 0.0}, solidQuery[3] = {0.0,0.0,0.0};
  double phi, dist;
  int iGlobalVertexFluid, jGlobalVertexSolid;
  vector<int> solidVertices;

  assert(nf_loc == size_loc_x);
  assert(nf_loc == size_loc_y);
  assert(nf_loc == size_loc_z);

  assert(size_buff_x == size_buff_y);
  assert(size_buff_y == size_buff_z);
  assert(size_buff_x == size_buff_z);

  ADTPoint ADT(size_buff_x, buff_x, size_buff_y, buff_y, size_buff_z, buff_z);
  for(int iVertex=0; iVertex < nf_loc; iVertex++){
    fluidPoint[0] = array_loc_x[iVertex];
    fluidPoint[1] = array_loc_y[iVertex];
    fluidPoint[2] = array_loc_z[iVertex];
    iGlobalVertexFluid = manager->getGlobalIndex("fluid", myid, iVertex);
    ADT.queryBallNN(3, fluidPoint, radius, solidVertices);
    for(int jVertex : solidVertices){
      jGlobalVertexSolid = manager->getGlobalIndex("solid", iProc, jVertex);
      solidQuery[0] = buff_x[jVertex];
      solidQuery[1] = buff_y[jVertex];
      solidQuery[2] = buff_z[jVertex];
      dist = distance(3, fluidPoint, 3, solidQuery);
      phi = PHI_RBF(dist, radius);
      B->setValue(iGlobalVertexFluid, jGlobalVertexSolid,phi);
      B_T->setValue(jGlobalVertexSolid, iGlobalVertexFluid, phi);
    }
    B->setValue(iGlobalVertexFluid, ns, 1.0);
    B->setValue(iGlobalVertexFluid, ns+1, fluidPoint[0]);
    B->setValue(iGlobalVertexFluid, ns+2, fluidPoint[1]);
    B_T->setValue(ns, iGlobalVertexFluid, 1.0);
    B_T->setValue(ns+1, iGlobalVertexFluid, fluidPoint[0]);
    B_T->setValue(ns+2, iGlobalVertexFluid, fluidPoint[1]);
    if(nDim == 3){
      B->setValue(iGlobalVertexFluid, ns+3, fluidPoint[2]);
      B_T->setValue(ns+3, iGlobalVertexFluid, fluidPoint[2]);
    }
  }
}

double CInterpolator::PHI_TPS(double &distance) const{

  if(distance > 0.0) return pow(distance,2)*log10(distance);
  else return 0.0;

}

double CInterpolator::PHI_RBF(double &distance, const double &radius) const{

  double eps(distance/radius);

  if(eps < 1.0) return pow((1.0-eps),4)*(4.0*eps+1.0);
  else return 0.0;

}

double CInterpolator::distance(int val_size1, double *array1, int val_size2, double *array2) const{

  return sqrt(pow(array2[0] - array1[0],2) + pow(array2[1] - array1[1],2) + pow(array2[2]-array1[2],2));
}
