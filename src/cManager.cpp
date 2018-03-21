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
 * Source for CManager.
 * Authors : D. THOMAS.
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
#endif //HAVE_MPI
  nPhyscis = 2;
  nIndex = 2;

  globalIndexRange.resize(nPhyscis);

  for(unsigned int ii=0; ii<globalIndexRange.size(); ii++){
    globalIndexRange[ii].resize(mpiSize);
    for(unsigned int jj=0; jj < globalIndexRange[ii].size(); jj++){
      globalIndexRange[ii][jj].resize(nIndex);
    }
  }

}

CManager::~CManager(){

#ifndef NDEBUG
  cout << "Calling CManager::~CManager()" << endl;
#endif  //NDEBUG

}

int CManager::getGlobalIndex(string const& str_physics, int const& iProc, int const& iVertex){


  int physics;

  if(str_physics.compare("fluid") == 0) physics = 0;
  else if(str_physics.compare("solid") == 0) physics = 1;
  else physics=1000;

  return globalIndexRange[physics][iProc][0] + iVertex;

}

void CManager::setGlobalIndexing(std::string str_physics,std::vector<std::vector<int> > index_range){

  int physics;

  assert(index_range.size() == static_cast<unsigned int>(mpiSize));

  if(str_physics.compare("fluid") == 0) physics = 0;
  else if(str_physics.compare("solid") == 0) physics = 1;
  else physics=1000;

  for(int ii=0; ii < mpiSize; ii++){
    for(int jj=0; jj < 2; jj++){
      globalIndexRange[physics][ii][jj] = index_range[ii][jj];
    }
  }
}
