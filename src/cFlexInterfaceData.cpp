/*!
 * Source for CFlexInterfaceData.
 * Authors : D. THOMAS.
 *
 * COPYRIHGHT (C) University of Liège, 2017.
 */


#include <vector>

#ifdef HAVE_MPI
#include "petscvec.h"
#endif

#include "../include/cMpi.h"
#include "../include/cFlexInterfaceData.h"

using namespace std;

CFlexInterfaceData::CFlexInterfaceData(int const& val_nPoint, int const& val_nDim):nPoint(val_nPoint), nDim(val_nDim){

  dataContainer.resize(nDim);

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nPoint, &(dataContainer[ii]));
    VecSet(dataContainer[ii], 0.0);
  }
#else //HAVE_MPI

#endif  //HAVE_MPI
}

CFlexInterfaceData::~CFlexInterfaceData(){

#ifdef HAVE_MPI
  for(int ii=0; ii<nDim; ii++){
    if(dataContainer[ii]) VecDestroy(&(dataContainer[ii]));
  }
#endif  //HAVE_MPI
}

vector<int> CFlexInterfaceData::getOwnershipRange() const{

  vector<int> ownershipRange(2);
#ifdef HAVE_MPI
  VecGetOwnershipRange(dataContainer[0], &(ownershipRange[0]), &(ownershipRange[1]));
#else
  ownershipRange[0] = 0;
  ownershipRange[1] = nPoint;
#endif

  return ownershipRange;
}
