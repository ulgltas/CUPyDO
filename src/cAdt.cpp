/*!
 * Source for CADT.
 * Authors : D. THOMAS.
 *
 * COPYRIHGHT (C) University of Liège, 2017.
 */

#include <iostream>
#include <cassert>

#include "../include/adtcore.h"
#include "../include/cAdt.h"

using namespace std;

ADTPoint::ADTPoint(int size_x, double* data_x, int size_y, double* data_y, int size_z, double* data_z){

  assert(size_x==size_y);
  assert(size_y==size_z);

  nDim = 3;
  size = size_x;

  data = NULL;
  dataIDs = NULL;
  dataTree = NULL;

  data = new double[nDim*size];
  dataIDs = new int[size];

  //Static cast to avoid warnings at compilation
  unsigned long size_unsigned_long = static_cast<unsigned long>(size);

  for(unsigned long ii=0; ii<size_unsigned_long; ii++){
    data[nDim*ii] = data_x[ii];
    data[nDim*ii+1] = data_y[ii];
    data[nDim*ii+2] = data_z[ii];
  }

  for(unsigned long ii=0; ii<size_unsigned_long; ii++){
    dataIDs[ii] = ii;
  }

  dataTree = new ADT_PointType(nDim, size, data, dataIDs);

}

ADTPoint::~ADTPoint(){

  cout << "CALL->ADTPoint::~ADTPoint()" << endl;
  if(data != NULL) delete data;
  if(dataIDs != NULL) delete dataIDs;
  if(dataTree != NULL) delete dataTree;
}

void ADTPoint::queryNN(int size, double *coord, int &pointID, double &distance){

  int rank = 0;
  dataTree->queryNearestNeighboor(coord, distance, pointID, rank);
}

void ADTPoint::queryBallNN(int size, double* coord, double radius, std::vector<int> &allIDs){

  vector<double> allDist;
  //vector<int> allIDS;
  int rank = 0;

  dataTree->queryBallNeighboors(coord, radius, allDist, allIDs, rank);

}
