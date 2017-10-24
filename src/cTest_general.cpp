/*!
 * General tests for CCUPyDO.
 * Authors : D. THOMAS.
 *
 * COPYRIGHT (C) University of Liège, 2017.
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
