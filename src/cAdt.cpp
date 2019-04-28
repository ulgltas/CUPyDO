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
 * Source for CADT.
 * Authors : D. THOMAS.
 */

#include <iostream>
#include <cassert>

#include "adtcore.h"
#include "cAdt.h"

using namespace std;

ADTPoint::ADTPoint(int size_x, double *data_x, int size_y, double *data_y, int size_z, double *data_z)
{

    assert(size_x == size_y);
    assert(size_y == size_z);

    nDim = 3;
    size = size_x;

    data = NULL;
    dataIDs = NULL;
    dataTree = NULL;

    data = new double[nDim * size];
    dataIDs = new int[size];

    //Static cast to avoid warnings at compilation
    unsigned long size_unsigned_long = static_cast<unsigned long>(size);

    for (unsigned long ii = 0; ii < size_unsigned_long; ii++)
    {
        data[nDim * ii] = data_x[ii];
        data[nDim * ii + 1] = data_y[ii];
        data[nDim * ii + 2] = data_z[ii];
    }

    for (unsigned long ii = 0; ii < size_unsigned_long; ii++)
    {
        dataIDs[ii] = ii;
    }

    dataTree = new ADT_PointType(nDim, size, data, dataIDs);
}

ADTPoint::~ADTPoint()
{

    if (data != NULL)
        delete data;
    if (dataIDs != NULL)
        delete dataIDs;
    if (dataTree != NULL)
        delete dataTree;
}

void ADTPoint::queryNN(int size, double *coord, int &pointID, double &distance)
{

    assert(size <= 3);

    int rank = 0;
    dataTree->queryNearestNeighboor(coord, distance, pointID, rank);
}

void ADTPoint::queryBallNN(int size, double *coord, double radius, std::vector<int> &allIDs)
{

    assert(size <= 3);

    vector<double> allDist;
    int rank = 0;

    dataTree->queryBallNeighboors(coord, radius, allDist, allIDs, rank);
}
