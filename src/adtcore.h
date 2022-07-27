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
 * Header for adtcore.
 * This implementation of the ADT is an adaption of the implementation
 * from the SU2 package, adapted for CUPyDO, and originally written
 * by E. van der Weide (University of Twente).
 *
 * Authors : D. THOMAS.
 */

#ifndef ADTCORE_H
#define ADTCORE_H

#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

const unsigned short CHILD_LEFT(0);
const unsigned short CHILD_RIGHT(1);

class ADT_Node
{
public:
    ADT_Node();                                          //Constructor
    ~ADT_Node();                                         //Destructor
    ADT_Node(const ADT_Node &Other_ADT_Node);            //Copy constructor
    ADT_Node &operator=(const ADT_Node &Other_ADT_Node); //Assignment operator

    void setDim(int dim);
    int getDim();
    void setNPoints(int val_nPoints);
    void setCentralNodeID(int val_node);
    void setChildTerminal(int numChild, bool val);
    void setChild(int numChild, int val);
    void setParent(int val);
    void setMinCoords(double *val_minCoords);
    void setMaxCoords(double *val_maxCoords);
    void setValMinCoords(int dim, double val);
    void setValMaxCoords(int dim, double val);
    void setIfMinValMinCoords(int dim, double val);
    void setIfMaxValMaxCoords(int dim, double val);

    int getNPoints(void) const;
    double getValMinCoords(int dim);
    double getValMaxCoords(int dim);
    int getChild(int numChild);
    int getParent() const;
    bool isChildTerminal(int numChil);
    int getCentralNodeID(void) const;

private:
    void copy(const ADT_Node &Other_ADT_Node); //Copy function
    double *minCoords;
    double *maxCoords;
    bool childrenAreTerminal[2]; //
    int children[2];
    int parent;
    int centralNodeID;
    int nPoints;
    int nDim;
};

class ADT_BaseType
{
protected:
    int nLeaves;
    int nDimADT;
    bool empty;
    std::vector<ADT_Node> leaves;
    std::vector<int> nodeParents;

    ADT_BaseType();                                                             //Constructor
    virtual ~ADT_BaseType();                                                    //Virutal destructor
    void buildADT(int nDim, int nPoints, double *coord);                        //Build the ADT
    double computeDistanceSquare(double *pointA, double *pointB);               //Compute the Euclidean distance between two points
    double computeDistanceSquare(double *pointA, ADT_Node &BBox);               //Compute the Euclidean distance between a point and a node (bounding box)
    bool intersectSphere(double *center, double const &radius, ADT_Node &BBox); //Determine if a sphere intersect a bounding box or not
    bool isIncluded(double *pointA, ADT_Node &BBox);                            //Determine if a poit is included into a bounding box

private:
    std::vector<double> coordMinLeaves;
    std::vector<double> coordMaxLeaves;

    ADT_BaseType(const ADT_BaseType &);            //Copy constructor
    ADT_BaseType &operator=(const ADT_BaseType &); //Assignment operator

public:
    bool isEmpty(void) const; //Determine if the tree is empty
    void display(void);       //Print
};

class ADT_PointType : public ADT_BaseType
{
private:
    std::vector<int> frontLeaves;
    std::vector<int> frontLeavesNew;
    std::vector<double> coordPoints;
    std::vector<int> localPointIDs;
    std::vector<int> ranksOfPoints;
    ADT_PointType();
    ADT_PointType(const ADT_PointType &);
    ADT_PointType &operator=(const ADT_PointType &);

public:
    ADT_PointType(int nDim, int nPoints, double *coord, int *pointID); //Constructor
    ~ADT_PointType();
    void queryNearestNeighboor(double *coord, double &dist, int &pointID, int &rankID);                                               //Query a nearest neighboor point
    void queryBallNeighboors(double *coord, double const &radius, std::vector<double> &dist, std::vector<int> &pointID, int &rankID); //Query all neighboor within a sphere
};

class ADT_Compare
{
private:
    const double *pointCoord;
    const int splitDirection;
    const int nDim;
    ADT_Compare();

public:
    ADT_Compare(const double *coord, const int splitDir, const int nDimADT);
    ~ADT_Compare();
    bool operator()(const int p0, const int p1) const;
};

#endif //ADTCORE_H
