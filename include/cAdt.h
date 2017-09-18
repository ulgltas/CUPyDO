/*!
 * Header for CADT.
 * Authors : D. THOMAS.
 *
 * COPYRIHGHT (C) University of Liège, 2017.
 */

#pragma once

#include <vector>

#include "../include/adtcore.h"

class ADTPoint{
protected:
  double* data;
  int *dataIDs;
  int size, nDim;
  ADT_PointType *dataTree;
public:
  //ADTPoint(int size_x, double* data_x, int val_size_y, double* val_data_y);
  ADTPoint(int size_x, double* data_x, int val_size_y, double* val_data_y, int size_z, double* data_z);
  ~ADTPoint();
  void queryNN(int size, double* coord, int &pointID, double &distance);
  void queryBallNN(int size, double* coord, double radius, std::vector<int> &allIDs);
};
