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
 * Header for the MPI wrapper.
 * Authors : D. THOMAS.
 */

#ifndef CMPI_H
#define CMPI_H

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifdef HAVE_MPI
typedef MPI_Comm Cupydo_Comm;
#else  //HAVE_MPI
typedef const int *Cupydo_Comm;
#endif //HAVE_MPI

/*!
 * Allow python to know if CUPyDO has been built with MPI
 * Authors : A. CROVATO.
 */
class CMpi
{
public:
    bool haveMPI;
    CMpi()
    {
#ifdef HAVE_MPI
        haveMPI = true;
#else  //HAVE_MPI
        haveMPI = false;
#endif //HAVE_MPI
    }
    virtual ~CMpi() {}
};

#endif //CMPI_H
