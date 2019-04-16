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
 * General tests for CCUPyDO.
 * Authors : D. THOMAS.
 */

#include <iostream>
#include <string>

#include "../include/cFlexInterfaceData.h"
#include "../include/cInterfaceMatrix.h"
#include "../include/cMpi.h"

using namespace std;

int main(int argc, char** argv){

  int rank, mpiSize;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif  //HAVE_MPI

  if(mpiSize != 4){
    throw string("MPI size must be equal to 4 for the test !");
  }

  CFlexInterfaceData flexInt(4,3);
  CInterfaceMatrix intMat(4,4);


  return 0;

}
