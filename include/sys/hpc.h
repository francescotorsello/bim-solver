/**
 *  @file      hpc.h
 *  @brief     High-Performance Computing (HPC) support.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _HPC_H_INCLUDED
#define _HPC_H_INCLUDED

#if _OPENMP

    #include <omp.h> // Use Open Multi-Processing (OMP)

    /** `OMP_parallel_for` indicates a parallelizable `for` loop.
     */
    #define OMP_parallel_for  _Pragma("omp parallel for") for

    // #define OMP_parallel_for  __pragma(omp parallel for) for  /* on windows */

#else
    #define OMP_parallel_for  for
#endif

#ifdef _USEMPI
    #include "mpiWorld.h"   // Use the real Message Passing Interface (MPI)
#else
    #include "mpiDummyWorld.h"
#endif

#endif // _HPC_H_INCLUDED
