/*!
 * Source for CManager.
 * Authors : D. THOMAS.
 *
 * COPYRIHGHT (C) University of Liège, 2017.
 */

#include <iostream>
#include <cassert>
#include <vector>

#include "../include/cMpi.h"
#include "../include/cManager.h"

using namespace std;

CManager::CManager(){

#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#else //HAVE_MPI
  mpiSize = 1;
#endif
  nPhyscis = 2;
  nIndex = 2;

  globalIndexRange.resize(nPhyscis);

  for(int ii=0; ii<globalIndexRange.size(); ii++){
    globalIndexRange[ii].resize(mpiSize);
    for(int jj=0; jj < globalIndexRange[ii].size(); jj++){
      globalIndexRange[ii][jj].resize(nIndex);
    }
  }

}

CManager::~CManager(){}

int CManager::getGlobalIndex(string const& str_physics, int const& iProc, int const& iVertex){


  int physics;

  if(str_physics.compare("fluid") == 0) physics = 0;
  else if(str_physics.compare("solid") == 0) physics = 1;
  else physics=1000;

  return globalIndexRange[physics][iProc][0] + iVertex;

}

void CManager::setGlobalIndexing(std::string str_physics,std::vector<std::vector<double> > index_range){

  int physics;

  assert(index_range.size() == mpiSize);

  if(str_physics.compare("fluid") == 0) physics = 0;
  else if(str_physics.compare("solid") == 0) physics = 1;
  else physics=1000;

  for(int ii=0; ii < mpiSize; ii++){
    for(int jj=0; jj < 2; jj++){
      globalIndexRange[physics][ii][jj] = index_range[ii][jj];
    }
  }
}
