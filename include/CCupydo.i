// SWIG configuration file of CUPyDO.
//
// Authors : D. THOMAS.
//
// COPYRIGHT (C) University of Liege, 2017.

%feature("autodoc","1");

%module(docstring=
"'CCUPyDO' module",
directors="1",
threads="1"
) CCupydo
%{
#define SWIG_FILE_WITH_INIT
#include "cMpi.h"
#include "cInterpolator.h"
#include "cManager.h"
#include "cFlexInterfaceData.h"
#include "cInterfaceMatrix.h"
#include "cLinearSolver.h"
%}

// petite bidouille pour pouvoir compiler avec "threads=1" et iterer sur des std_vector
// (sinon ca explose à la destruction de l'iterateur)
%nothread swig::SwigPyIterator::~SwigPyIterator();

// ----------- MODULES UTILISES ------------
%import "cMpi.h"
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "typemaps.i"
%include "numpy.i"
#ifdef HAVE_MPI
  %include "mpi4py/mpi4py.i"
  %mpi4py_typemap(Comm, MPI_Comm)
  %include "petsc4py/petsc4py.i"
#endif

// ----------- NATIVE CLASSES ----------------
%init %{
  import_array();
%}

namespace std {
   %template(VecInt) vector<int>;
   %template(VecDouble) vector<double>;
   %template(VecVecInt) vector< vector<int> >;
   %template(VecVecDouble) vector< vector<double> >;
   %template() vector<unsigned long>;
}

//%numpy_typemaps(int,    NPY_INT   , int)
//%numpy_typemaps(unsigned long,    NPY_ULONGLONG   , int)
//%numpy_typemaps(unsigned long*,    NPY_INT   , int*)
//%numpy_typemaps(int*,    NPY_INT   , unsigned long*)
//%numpy_typemaps(int,    NPY_INT   , unsigned long)
//%numpy_typemaps(unsigned long*     , NPY_UINT64    , int*)

//%numpy_typemaps(int,    NPY_INT   , int)
//%numpy_typemaps(int,    NPY_INT   , int)
//%numpy_typemaps(int,    NPY_INT   , int)

//%apply (int DIM1, double* IN_ARRAY1) {(int size, double *array)};
//%apply (int DIM1, double* IN_ARRAY1) {(unsigned long size, double *array)};
//%apply (int DIM1, double* IN_ARRAY1) {(unsigned short size, double *array)};
//%apply (int DIM1, double* IN_ARRAY1) {(unsigned int size, double *array)};
//%apply (int DIM1, unsigned short* IN_ARRAY1) {(unsigned long size, unsigned short *array)};
//%apply (int DIM1, unsigned long* IN_ARRAY1) {(unsigned long size, unsigned long *array)};

//%apply (int DIM1, double* IN_ARRAY1) {(int val_size1, double *array1), 
//                                      (int val_size2, double *array2)}

%apply (int DIM1, double* IN_ARRAY1) {(int size_loc_x, double* array_loc_x),
                                      (int size_loc_y, double* array_loc_y),
                                      (int size_loc_z, double* array_loc_z),
                                      (int size_buff_x, double* buff_x),
                                      (int size_buff_y, double* buff_y),
                                      (int size_buff_z, double* buff_z)}

%apply (int DIM1, int* IN_ARRAY1) {(int size_indices, int* indices_list)}
%apply (int DIM1, double* IN_ARRAY1) {(int size_values, double *values_array)}
%apply (int DIM1, double* IN_ARRAY1) {(int size, double *data)}
%apply(int *DIM1, double** ARGOUTVIEW_ARRAY1) {(int* size, double** data_array)}
#ifndef HAVE_MPI
%apply(int *DIM1, int* DIM2, double** ARGOUTVIEW_ARRAY2) {(int* size1, int* size2, double** mat_array)}
#endif

//%apply (int DIM1, double* IN_ARRAY1) {(int size, double *in_array)};
//%apply (int DIM1, double* IN_ARRAY1) {(unsigned long val_size, double *coord)};
//%apply (int DIM1 , double* ARGOUT_ARRAY1) {(int size, double* out_array)}
//%apply (int *DIM1, double** ARGOUTVIEW_ARRAY1) {(int* size, double** out_array)}

//%apply int& OUTPUT {unsigned long &pointID};
//%apply double& OUTPUT {double &distance};

%feature("director") CInterpolator;
//%pythonappend CInterpolator "self.__disown__()"    // for directors
%include "cInterpolator.h"

%feature("director") CManager;
//%pythonappend CManager "self.__disown__()"    // for directors
%include "cManager.h"

%feature("director") CFlexInterfaceData;
//%pythonappend CFlexInterfaceData "self.__disown__()"    // for directors
%include "cFlexInterfaceData.h"

%feature("director") CInterfaceMatrix;
//%pythonappend CInterfaceMatrix "self.__disown__()"    // for directors
%include "cInterfaceMatrix.h"

%feature("director") CLinearSolver;
//%pythonappend CLinearSolver "self.__disown__()"    // for directors
%include "cLinearSolver.h"

%feature("director::except"){
    if ($error != NULL) {
        std::cout << "[in director:except]\n";
        throw std::runtime_error("Director problem");
    }
}