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
 * Header for CADT.
 * Authors : D. THOMAS.
 */

#ifndef CADT_H
#define CADT_H

#include <vector>

#include "adtcore.h"

class ADTPoint
{
protected:
  double *data;
  int *dataIDs;
  int size, nDim;
  ADT_PointType *dataTree;

public:
  ADTPoint(int size_x, double *data_x, int val_size_y, double *val_data_y, int size_z, double *data_z);
  ~ADTPoint();
  void queryNN(int size, double *coord, int &pointID, double &distance);
  void queryBallNN(int size, double *coord, double radius, std::vector<int> &allIDs);
};

#endif //CADT_H
