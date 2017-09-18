/*!
 * Source for adtcore.
 * This implementation of the ADT is an adaption of the implementation
 * from the SU2 package, adapted for Cupydo, and originally written
 * by E. van der Weide (University of Twente).
 *
 * Authors : D. THOMAS.
 *
 * COPYRIHGHT (C) 2012-2017 SU2, the open-source CFD code.
 * COPYRIHGHT (C) University of Liège, 2017.
 */

#include "../include/adtcore.h"

#include <iomanip>
#include <stdexcept>

using namespace std;

/* CLASS ADT_NODE*/

ADT_Node::ADT_Node(){}

ADT_Node::~ADT_Node(){}

ADT_Node::ADT_Node(const ADT_Node &Other_ADT_Node){
    copy(Other_ADT_Node);
}

ADT_Node& ADT_Node::operator=(const ADT_Node &Other_ADT_Node){
    copy(Other_ADT_Node);
    return *this;
}

void ADT_Node::copy(const ADT_Node &Other_ADT_Node){

    nDim = Other_ADT_Node.nDim;

    childrenAreTerminal[0] = Other_ADT_Node.childrenAreTerminal[0];
    childrenAreTerminal[1] = Other_ADT_Node.childrenAreTerminal[1];

    children[0] = Other_ADT_Node.children[0];
    children[1] = Other_ADT_Node.children[1];

    centralNodeID = Other_ADT_Node.centralNodeID;

    minCoords = Other_ADT_Node.minCoords;
    maxCoords = Other_ADT_Node.maxCoords;

    parent = Other_ADT_Node.parent;

}

void ADT_Node::setDim(int dim){
    nDim = dim;
}

int ADT_Node::getDim(){
  return nDim;
}

void ADT_Node::setNPoints(int val_nPoints){
    nPoints = val_nPoints;
}

void ADT_Node::setCentralNodeID(int val_node){
    centralNodeID = val_node;
}

void ADT_Node::setMinCoords(double *val_minCoords){
    minCoords = val_minCoords;
}

void ADT_Node::setMaxCoords(double *val_maxCoords){
    maxCoords = val_maxCoords;
}

void ADT_Node::setValMinCoords(int dim, double val){
    if(dim > nDim) throw out_of_range("Exceed node dimension");
    minCoords[dim] = val;
}

void ADT_Node::setValMaxCoords(int dim, double val){
    if(dim > nDim) throw out_of_range("Exceed node dimension");
    maxCoords[dim] = val;
}

void ADT_Node::setChild(int numChild, int val){
    children[numChild] = val;
}

void ADT_Node::setParent(int val){
    parent = val;
}

void ADT_Node::setChildTerminal(int numChild, bool val){
    childrenAreTerminal[numChild] = val;
}

int ADT_Node::getNPoints(void) const{
    return nPoints;
}

void ADT_Node::setIfMinValMinCoords(int dim, double val){
    if(dim > nDim) throw out_of_range("Exceed node dimension");
    minCoords[dim] = min(minCoords[dim], val);
}

void ADT_Node::setIfMaxValMaxCoords(int dim, double val){
    if(dim > nDim) throw out_of_range("Exceed node dimension");
    maxCoords[dim] = max(maxCoords[dim], val);
}

double ADT_Node::getValMinCoords(int dim){
    if(dim > nDim) throw out_of_range("Exceed node dimension");
    return minCoords[dim];
}

double ADT_Node::getValMaxCoords(int dim){
    if(dim > nDim) throw out_of_range("Exceed node dimension");
    return maxCoords[dim];
}

int ADT_Node::getChild(int numChild){
    return children[numChild];
}

int ADT_Node::getParent() const{
    return parent;
}

bool ADT_Node::isChildTerminal(int numChil){
    return childrenAreTerminal[numChil];
}

int ADT_Node::getCentralNodeID(void) const{
    return centralNodeID;
}

/* CLASS ADT_BASETYPE */

ADT_BaseType::ADT_BaseType(){}

ADT_BaseType::~ADT_BaseType(){}

void ADT_BaseType::buildADT(int nDim, int nPoints, double *coord){

    /* Compute the total number of leaves. For an optimally balanced tree, this number is equal to nPoints-1. */
    nDimADT = nDim;
    empty = false;
    nLeaves = nPoints -1;
    if(nPoints <= 1) ++nLeaves; //If 0 then nLeaves = 0, if 1 then nLeaves = 1
    if(nLeaves == 0) {empty = true; return;}

    /* Allocate memory for the leaves container as well as for the bounding boxes (min and max coords) */
    leaves.resize(nLeaves);
    nodeParents.resize(nPoints);
    coordMinLeaves.resize(nDim*nLeaves);
    coordMaxLeaves.resize(nDim*nLeaves);

    /* Define the controllers, which follow and control the subdivision of the tree. */
    /* The New versin of the controllers are used to set them for the next round. */
    /* This is the number of couples of point we can create from the total number of points. */
    int maxCouples = (nPoints+1)/2;
    /* These will contain the indices of the points (from 0 to nPoints-1). The order can vary during the building process. */
    vector<int> pointIDs(nPoints), pointIDsNew(nPoints);
    /*These are used to store the points distribution among the leaves at each round (incremental).
      nPoints[i] is the cumulated number of points belonging to the preivous leaves on the same round. */
    vector<int> nPointIDs(maxCouples+1),   nPointIDsNew(maxCouples+1);
    /* These will store the current leaves indices that will be treated at each round. */
    vector<int> curLeaf(maxCouples),       curLeafNew(maxCouples);

    /* Start building the ADT. */

    /* Initialize the controllers such that all points belong to the root leaf. */
    nPointIDs[0] = 0; nPointIDs[1] = nPoints;
    curLeaf[0] = 0;

    for(int i=0; i<nPoints; ++i) pointIDs[i] = i;

    int nLeavesToDivide = 1, nLeavesTot = 1;       //Number of leaves to be divided during the current round and total number of leaves created
    int nLeavesToDivideNew = 0;                    //Number of leaves that will have to be divided during the next round
    int leafNPoints(nPoints);                      //Number of points in the current leaf

    //cout << "Building ADT..." << endl;
    int round(0);

    leaves[0].setParent(0);

    /*Loop to subdivide the leaves such that the tree is optimally balanced. Each round will set a tree level.*/
    while(1) {

      /* Criterion to exit the loop. Happens when all the leaves of the previous round are terminal. */
      if(nLeavesToDivide == 0) break;

      round++;

      /* Initializations for the next round of subdivisions. */
      nLeavesToDivideNew = 0;
      nPointIDsNew[0] = 0;

      /* Loop over the current number of leaves to be divided. */
      for(int i=0; i<nLeavesToDivide; ++i) {

        leafNPoints = nPointIDs[i+1] - nPointIDs[i];
        int mm = curLeaf[i];               //Get the leaf index
        leaves[mm].setNPoints(leafNPoints);
        leaves[mm].setDim(nDimADT);

        /*Set the pointers for the coordinates of the leaf to the correct locations in the vectors coorMinLeaves and coorMaxLeaves .*/
        leaves[mm].setMinCoords(coordMinLeaves.data() + nDim*mm);
        leaves[mm].setMaxCoords(coordMaxLeaves.data() + nDim*mm);

        /* Determine the bounding box coordinates of the leaf. */

        //First guess : we take the first point on the current leaf (nPointIDs[i])
        int index = nDim*pointIDs[nPointIDs[i]];
        for(int dim=0; dim<nDim; ++dim){
            leaves[mm].setValMinCoords(dim, coord[index+dim]);
            leaves[mm].setValMaxCoords(dim, coord[index+dim]);
        }

        //Loop over all the points (except the guess) on the current leaf to find the extrema !!maybe we can find something more efficient !!
        for(int j=(nPointIDs[i]+1); j<nPointIDs[i+1]; ++j) {
          index = nDim*pointIDs[j];
          for(int dim=0; dim<nDim; ++dim) {
            leaves[mm].setIfMinValMinCoords(dim, coord[index+dim]);
            leaves[mm].setIfMaxValMaxCoords(dim, coord[index+dim]);
          }
        }

        /* Compute the split direction for this leaf. Split direction = largest dimension of the leaf so that isotropy ids reached as quickly as possible */
        int splitDir= 0;
        double distMax = -1.0;
        for(int dim=0; dim<nDim; ++dim) {
          const double dist = leaves[mm].getValMaxCoords(dim) - leaves[mm].getValMinCoords(dim);
          if(dist > distMax){
              distMax = dist;
              splitDir = dim;
          }
        }

        /* Sort the points of the current leaf in increasing order along the split direction.
           The functor ADT_Compare is used as the sorting criterion.
           This will modify pointIDs. */
        sort(pointIDs.data() + nPointIDs[i], pointIDs.data() + nPointIDs[i+1],
             ADT_Compare(coord, splitDir, nDim));

        /* Associate a node ID (approximately the central node) to the leaf */
        leaves[mm].setCentralNodeID(pointIDs[nPointIDs[i] + leafNPoints/2]);

        /* Determine if this leaf is terminal or if it must be divided. */
        if(leafNPoints <= 2) {
          /* Terminal leaf. The children are the ID's of the points defining the bounding box.*/
          leaves[mm].setChild(CHILD_LEFT, pointIDs[nPointIDs[i]]);
          nodeParents[pointIDs[nPointIDs[i]]] = mm;
          leaves[mm].setChild(CHILD_RIGHT, pointIDs[nPointIDs[i+1]-1]);
          nodeParents[pointIDs[nPointIDs[i+1]-1]] = mm;
          leaves[mm].setChildTerminal(CHILD_LEFT, true);
          leaves[mm].setChildTerminal(CHILD_RIGHT, true);
        }
        else {

          /* The leaf must be divided. */
          int nPointsLeftChild = (leafNPoints+1)/2;   //Determine the number of points in the left leaf (at least equals 2)
          int kk  = nPointIDs[i] + nPointsLeftChild;  //Number of points that will be on the left child leaf + cumulated number of points on the previous leaves on this round
          int nfl = nPointIDsNew[nLeavesToDivideNew]; //Initialize the counter with the right number of points

          /* Copy the ID's of the left points into pointIDsNew. Update the
             corresponding entry in nPointIDsNew with the new number of points. */
          for(int k=nPointIDs[i]; k<kk; ++k)
            pointIDsNew[nfl++] = pointIDs[k];
          nPointIDsNew[nLeavesToDivideNew+1] = nfl;

          leaves[mm].setChild(CHILD_LEFT, nLeavesTot);     //The left child ID is actually the current value of the total number of leaf (root = leaf 0).
          leaves[mm].setChildTerminal(CHILD_LEFT, false);
          leaves[nLeavesTot].setParent(mm);

          curLeafNew[nLeavesToDivideNew] = nLeavesTot; //Updare the current leaves list for the next round
          ++nLeavesToDivideNew;                        //For the next round, a new leaf will have to be divided
          ++nLeavesTot;                                //Increment nLeavesTot since a new leaf has been created

          /* The right child leaf will only be created if we can put more than one point
             in it (= if the current leaf has more than three points).
             If the new leaf only has one point in it, this point will be used as a terminal child */
          if(leafNPoints == 3)
          {

            /* Only three points present in the current leaf : last point is stored as the second child of
               the current leaf. */
            leaves[mm].setChild(CHILD_RIGHT, pointIDs[nPointIDs[i+1]-1]);
            leaves[mm].setChildTerminal(CHILD_RIGHT, true);
            nodeParents[pointIDs[nPointIDs[i+1]-1]] = mm;
          }
          else {

            /* More than 3 points are present and thus the right leaf is created.
               Same principle as for the left leaf. */
            int nPointsRightChild = leafNPoints-nPointsLeftChild;
            int nfr = nPointIDsNew[nLeavesToDivideNew];

            for(int k=kk; k<nPointIDs[i+1]; ++k)
              pointIDsNew[nfr++] = pointIDs[k];
            nPointIDsNew[nLeavesToDivideNew+1] = nfr;

            leaves[mm].setChild(CHILD_RIGHT, nLeavesTot);     //The child ID is actually the current value of the total number of leaf (root = leaf 0).
            leaves[mm].setChildTerminal(CHILD_RIGHT, false);
            leaves[nLeavesTot].setParent(mm);

            curLeafNew[nLeavesToDivideNew] = nLeavesTot;  //Update the current leaves list for the next round
            ++nLeavesToDivideNew;                         //For the next round, a new leaf will have to be divide
            ++nLeavesTot;                                 //Increment nLeavesTot since a new leaf has been created
          }
        }
      }

      /* Set the data for the next round. */
      nLeavesToDivide = nLeavesToDivideNew;
      for(int i=0; i<=nLeavesToDivide; ++i)           nPointIDs[i] = nPointIDsNew[i];
      for(int i=0; i< nLeavesToDivide; ++i)           curLeaf[i]   = curLeafNew[i];
      for(int i=0; i<nPointIDs[nLeavesToDivide]; ++i) pointIDs[i]  = pointIDsNew[i];
    }

    cout << "ADT has been built in " << round << " rounds." << endl;
}

double ADT_BaseType::computeDistanceSquare(double *pointA, double *pointB){

    double distance(0.0);
    for(int dim=0; dim<nDimADT; dim ++){
        distance += pow(pointA[dim] - pointB[dim], 2);
    }

    return distance;
}

double ADT_BaseType::computeDistanceSquare(double *pointA, ADT_Node &BBox){

    double distance(0.0), ds(0.0);

    for(int dim=0; dim<nDimADT; ++dim) {
      ds = 0.0;
      if(     pointA[dim] < BBox.getValMinCoords(dim)) ds = pointA[dim] - BBox.getValMinCoords(dim);
      else if(pointA[dim] > BBox.getValMaxCoords(dim)) ds = pointA[dim] - BBox.getValMaxCoords(dim);

      distance += ds*ds;
    }

    return distance;
}

bool ADT_BaseType::intersectSphere(double *center, const double &radius, ADT_Node &BBox){

    double distanceSquare(0.0);

    distanceSquare = computeDistanceSquare(center, BBox);
    //cout << "Computed squared distance" << distanceSquare << endl;

    if (distanceSquare == 0){   //Sphere center is included in the box, so we have intersection
      return true;
    }
    else{
      if(pow(radius,2) < distanceSquare) return false; // The sphere is completely outside of the box
      else return true;
    }

}

bool ADT_BaseType::isIncluded(double *pointA, ADT_Node &BBox){

    double distance(1.0);

    distance = computeDistanceSquare(pointA, BBox);

    if(distance == 0) return true;
    else return false;

}

void ADT_BaseType::display(void){

    for(int ii=0; ii<leaves.size(); ii++){
        cout << "Leaf " << ii << endl;
        cout << "Central node ID : " << leaves[ii].getCentralNodeID() << endl;
        cout << "Child left : " << leaves[ii].getChild(CHILD_LEFT) << "--" << leaves[ii].isChildTerminal(CHILD_LEFT) << endl;
        cout << "Child right : " << leaves[ii].getChild(CHILD_RIGHT) << "--" << leaves[ii].isChildTerminal(CHILD_RIGHT) << endl;
        cout << "Parent : " << leaves[ii].getParent() << endl;
        //cout << "Xmin " << leaves[ii].minCoords[0] << endl;
        //cout << "Xmax " << leaves[ii].getValMaxCoords(0) << endl;
        //cout << "Ymin " << leaves[ii].minCoords[1] << endl;
        //cout << "Ymax " << leaves[ii].getValMaxCoords(1) << endl;
        //cout << "----------------------------------" << endl;
    }

    for(int ii=0; ii<nodeParents.size(); ii++){
        cout << nodeParents[ii] << endl;
    }
}

/* CLASS ADT_POINTTYPE */

ADT_PointType::ADT_PointType(int nDim, int nPoints, double *coord, int *pointID){

    //cout << "***************************************" << endl;
    //cout << "Creating ADT" << endl;
    //cout << "Type : Points" << endl;
    //cout << "Dimensions : " << nDim << endl;
    //cout << "Number of points : " << nPoints << endl << endl;

    coordPoints.assign(coord, coord + nDim*nPoints);
    localPointIDs.assign(pointID, pointID + nPoints);
    ranksOfPoints.assign(nPoints, 0);

    buildADT(nDim, nPoints, coord);
    //cout << "***************************************" << endl;
}

ADT_PointType::~ADT_PointType(){}

void ADT_PointType::queryNearestNeighboor(double *coord, double &dist, int &pointID, int &rankID){


    /* Initialize the nearest node to the central node of the
       root leaf. */
    int kk = leaves[0].getCentralNodeID(), minIndex;
    double *coordTarget = coordPoints.data() + nDimADT*kk;

    pointID  = localPointIDs[kk];
    rankID   = ranksOfPoints[kk];
    minIndex = kk;

    int round(0);

    dist = computeDistanceSquare(coord, coordTarget);
    //cout << "Initialize with leaf 0 at distance " << dist << " from central node " << kk << endl;

    /* Traverse the tree to find the nearest node and start at the root. */
    frontLeaves.clear();       // Make sure to wipe out any data from aprevious search.
    frontLeaves.push_back(0);  // Initialize frontLeaves such that it only contains the root leaf.


    //cout << "Initial distance " << dist << endl;

    /* Infinite loop of the tree traversal. */
    while(1) {

        round++;
        //cout << "********* ROUND *********** : " << round << endl;

      /* Clear the front for the next round. */
      frontLeavesNew.clear();

      /* Loop over the leaves of the current front. */
      for(int i=0; i<frontLeaves.size(); ++i) {

        /* Get the current leaf and loop over its children (always 2 children). */
        const int iLeaf = frontLeaves[i];
        //cout << "Looking in leaf : " << iLeaf << endl;
        for(int iChild=0; iChild<2; iChild++) {

          /* Determine whether this child is a point or a leaf. */
          int childID = leaves[iLeaf].getChild(iChild);
          if( leaves[iLeaf].isChildTerminal(iChild) ) {

            /*--- Child is a point. Evaluate the distance to that point. ---*/
            coordTarget = coordPoints.data() + nDimADT*childID;
            double distTarget = 0;
            distTarget = computeDistanceSquare(coord, coordTarget);

            //cout << "Child is a point " << childID <<  " with distance " << sqrt(distTarget) << endl;

            /* If the distance is smaller than the current minimum, update the search result. */
            if(distTarget < dist) {
              dist     = distTarget;
              pointID  = localPointIDs[childID];
              rankID   = ranksOfPoints[childID];
              minIndex = childID;
              //cout << "Updating new dist " << sqrt(dist) << endl;
            }
          }
          else {

            /* Child is a leaf. Evaluate the distance to that leaf. */
            double posDist = 0.0;
            posDist = computeDistanceSquare(coord, leaves[childID]);  //This will return 0.0 if the point is included in the bounding box associated to the leaf.

            //cout << "Child is a leaf " << childID <<  " with distance " << sqrt(posDist) << endl;

            /* If the distance is smaller thant the current minimum, store the leaf ID to be treated during the next round.
               Also evaluate the distance to the leaf central node, and use it to update the search result.*/
            if(posDist < dist) {
              //cout << "I will have to look at the leaf " << childID << " during next round !" << endl;
              frontLeavesNew.push_back(childID);

              const int jj = leaves[childID].getCentralNodeID();

              coordTarget = coordPoints.data() + nDimADT*jj;
              double distTarget = 0;
              distTarget = computeDistanceSquare(coord, coordTarget);

              //cout << "Distance with the central node " << jj << " : " << sqrt(distTarget) << endl;

              if(distTarget < dist) {
                dist     = distTarget;
                pointID  = localPointIDs[jj];
                rankID   = ranksOfPoints[jj];
                minIndex = jj;
                //cout << "Updating new dist " << sqrt(dist) << endl;
              }
            }
          }
        }
      }

      /* Update the data for the next round*/
      frontLeaves = frontLeavesNew;

      /* If the new front is empty, it means we have reached a terminal leaf and the search is over. */
      if(frontLeaves.size() == 0) break;
    }

    /* Recompute the distance. */
    coordTarget = coordPoints.data() + nDimADT*minIndex;
    dist = computeDistanceSquare(coord, coordTarget);

    /* Up to now we have used the distance squared. Take the sqaure root of the result. */
    dist = sqrt(dist);
}

void ADT_PointType::queryBallNeighboors(double *coord, double const& radius, std::vector<double> &dist, std::vector<int> &pointID, int &rankID){

    //int nearestNeighboor(0);
    double distanceSquare(0.0);
    double* coordTarget;

    int round(0);

    dist.clear();
    pointID.clear();
    rankID = -1;

    frontLeaves.clear();
    frontLeaves.push_back(0);

    while (1) {

      round++;

      frontLeavesNew.clear();

      /* Loop over the leaves of the current round */
      for(int i=0; i<frontLeaves.size(); i++){
        const int iLeaf = frontLeaves[i];
        //cout << "Looking in leaf " << iLeaf << endl;

        /* Loop over the children of the current leaf */
        for(int iChild=0; iChild<2; iChild++){
          int childID = leaves[iLeaf].getChild(iChild);

          /* If the child is a point, test if it is included in the ball search */
          if(leaves[iLeaf].isChildTerminal(iChild)){
            //cout << "Child is point " << iChild << endl;
            coordTarget = coordPoints.data() + nDimADT*childID;
            distanceSquare = computeDistanceSquare(coord, coordTarget);
            /* If the point is included in the ball search, add the point to the list */
            if (distanceSquare <= pow(radius,2)){
              //cout << "We have a point with distance squared " << distanceSquare << endl;
              dist.push_back(sqrt(distanceSquare));
              pointID.push_back(childID);
            }
          }
          /* If the child is a leaf, test if it intersects the ball search */
          else{
            /* If the leaf intersects the ball, it has to be investigated at the next round */
            if(intersectSphere(coord, radius, leaves[childID])){
              //cout << "We have intersection with leaf " << childID << " and it will be investigated udring next round !" << endl;
              frontLeavesNew.push_back(childID);
            }
            /* If not, the ball is outside the leaf and this part of the tree will never be traversed */
            else{
              //cout << "Close the way at leaf " << childID << endl;
            }
          }
        }
      }

      frontLeaves = frontLeavesNew;
      if(frontLeaves.size() == 0) break;
    }
}

/* FUNCTOR ADT_COMPARE*/

ADT_Compare::ADT_Compare(const double *coord, const int splitDir, const int nDimADT): pointCoord(coord), splitDirection(splitDir), nDim(nDimADT){}

ADT_Compare::~ADT_Compare(){}

bool ADT_Compare::operator ()(const int p0, const int p1) const{
    return pointCoord[nDim*p0+splitDirection] < pointCoord[nDim*p1+splitDirection];
}
