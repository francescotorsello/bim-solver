/**
 *  @file      bispec.cpp
 *  @brief     The covariant BSSN evolution for spherically symmetric bimetric spacetimes,
 *             with the pseudospectral method.
 *  @authors   Francesco Torsello, Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <vector>
#include<functional>

#include "../include/numMethods/dataTypes.h"

#include "../include/sys/slog.h"
#include "../include/sys/paramsHolder.h"
#include "../include/bimetricModel.h"

#ifndef OBSERVER
    #define OBSERVER 1
#endif // OBSERVER

#ifndef _EVOLVE_DSIG
    #define _EVOLVE_DSIG 0
#endif // _EVOLVE_DSIG

#ifndef _DETECT_NAN
    #define _DETECT_NAN 1
#endif // _DETECT_NAN

using namespace std;

#include "../include/bispec/fields/fields.h"

//////////////////////////////////////////////////////////////////////////////////////////
///// Include class ChebyshevCoefficients
//////////////////////////////////////////////////////////////////////////////////////////

#include "../include/bispec/spectral-id/ChebyshevCoefficients.h"


//////////////////////////////////////////////////////////////////////////////////////////
///// Include class BispecInput
//////////////////////////////////////////////////////////////////////////////////////////

#include "../include/bispec/spectral-id/BispecInput.h"


//////////////////////////////////////////////////////////////////////////////////////////
///// Include class PrimaryFields
//////////////////////////////////////////////////////////////////////////////////////////

#include "../include/bispec/fields/PrimaryFields.h"


//////////////////////////////////////////////////////////////////////////////////////////
///// Include class GaugeVariables
//////////////////////////////////////////////////////////////////////////////////////////

#include "../include/bispec/fields/GaugeVariables.h"

//////////////////////////////////////////////////////////////////////////////////////////
/////  Include class DependentFields
//////////////////////////////////////////////////////////////////////////////////////////

#include "../include/bispec/fields/DependentFields.h"


//////////////////////////////////////////////////////////////////////////////////////////
/////  Include class BispecEvolve
//////////////////////////////////////////////////////////////////////////////////////////

#include "../include/bispec/evolution/BispecEvolve.h"


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int main()
{

    /** Import the values of the Chebyshev polynomials at the collocation points.
        The are constant in time, hence we can compute them in Mathematica with high
        precision and load them once here.
      */

    ChebyshevCoefficients chebyshevValues( "../bim-solver/include/chebyshev-values/chebyshev-values.dat" );
    if( ! chebyshevValues.isOK () ) {
        return -1;
    }

    /*std::cout << "Values of the Chebyshev polynomials in ChebyshevCoefficients," << std::endl << std::endl;

    std::cout << chebyshevValues.orders() << " x " << chebyshevValues.expOrd() + 1
        << " x " << chebyshevValues.colpoints() << std::endl;

    for( size_t der_order = 0; der_order < chebyshevValues.orders(); ++der_order )
    {
        for( size_t cheby_index = 0; cheby_index < chebyshevValues.expOrd() + 1;
            ++cheby_index )
        {
            std::cout << std::endl;
            for( size_t n = 0; n < chebyshevValues.colpoints(); ++n )
            {
                std::cout << "(" << der_order << "," << cheby_index << "," << n << ") = "
                          << chebyshevValues(der_order,cheby_index,n) << std::endl;
            }
        }
    }*/

    /*std::cout << "Values of the regularized first derivatives in  "
              << "ChebyshevCoefficients," << std::endl << std::endl;

    std::cout << chebyshevValues.orders() << " x " << chebyshevValues.expOrd() + 1
        << " x " << chebyshevValues.colpoints() << std::endl;

    for( size_t cheby_index = 0; cheby_index < chebyshevValues.expOrd() + 1;
        ++cheby_index )
    {
        std::cout << std::endl;
        for( size_t n = 0; n < chebyshevValues.colpoints(); ++n )
        {
            std::cout << "(" << cheby_index << "," << n << ") = "
                        << chebyshevValues.reg_der_coeff[ n +
                            cheby_index * chebyshevValues.colpoints() ] << std::endl;
        }
    }*/

    /*std::cout << "Values of the regularized second derivatives in  "
              << "ChebyshevCoefficients," << std::endl << std::endl;

    std::cout << chebyshevValues.orders() << " x " << chebyshevValues.expOrd() + 1
        << " x " << chebyshevValues.colpoints() << std::endl;

    for( size_t cheby_index = 0; cheby_index < chebyshevValues.expOrd() + 1;
        ++cheby_index )
    {
        std::cout << std::endl;
        for( size_t n = 0; n < chebyshevValues.colpoints(); ++n )
        {
            std::cout << "(" << cheby_index << "," << n << ") = "
                        << chebyshevValues.reg_derr_coeff[ n +
                            cheby_index * chebyshevValues.colpoints() ] << std::endl;
        }
    }*/

    //cc.exportChebyCoeffs();

    /*std::cout << "The following is the even evolution matrix," << std::endl << std::endl;

    for( size_t even_cheby_index = 0; even_cheby_index < chebyshevValues.chebys();
        ++even_cheby_index )
    {
        for( size_t n = 0; n < chebyshevValues.colpoints(); ++n )
        {
            std::cout << "(" << even_cheby_index << "," << n << ") = "
                        << chebyshevValues.ee_matrix_even[ n +
                            even_cheby_index * chebyshevValues.colpoints() ]
                            << std::endl;
        }
    }

    std::cout << std::endl << "The following is the odd evolution matrix,"
        << std::endl << std::endl;

    for( size_t odd_cheby_index = 0; odd_cheby_index < chebyshevValues.chebys();
        ++odd_cheby_index )
    {
        for( size_t n = 0; n < chebyshevValues.colpoints(); ++n )
        {
            std::cout << "(" << odd_cheby_index << "," << n << ") = "
                        << chebyshevValues.ee_matrix_odd[ n +
                            odd_cheby_index * chebyshevValues.colpoints() ]
                            << std::endl;
        }
    }*/

    /** Import the initial data, i.e., the values of the Chebyshev spectral coefficients
        in the Chebyshev series of the fields on the initial hypersurface.
      */

    BispecInput ID("../run/specInput.dat");
    if( ! ID.isOK () ) {
        return -1;
    }

    /*std::cout << "The following is the spectral initial data," << std::endl;

    std::cout << ID.n_allfields() << " x " << ID.n_chebycoeffs() << std::endl;

    for( size_t field = 0; field < ID.n_allfields(); ++field )
    {
        std::cout << std::endl << "This is the field " << field << std::endl << std::endl;
        for( size_t n = 0; n < ID.n_chebycoeffs(); ++n )
        {
            std::cout << "(" << field << "," << n << ") = "
                      << ID(field,n) << std::endl;
        }

    }

    std::cout << std::endl;*/

    /// - Read the run-time configuration parameters (taken from bim-solver)
    ///
    //Parameters params( argc >= 2 ? argv[1] : "config.ini" );
    Parameters params( "../run/config.ini" );

    /** Define the array containing the spectral coefficients as a function of the fields
        defined in the namespace fields, the time step m and the
        collocation point index n.
      */

    //ChebyshevExpansion chebySeries( ID );

    /*std::cout << "The following is the spectral initial data assigned to the
        coefficients," << std::endl;
    for( size_t field = 0; field < ID.n_fields(); ++field )
    {
        std::cout << "This is the field " << field << std::endl;
        for( size_t n = 0; n < ID.exp_order() + 1; ++n )
        {
            std::cout << "(" << field << "," << n << ") = "
                      << chebySeries.spectral_coeffs( 0, field, n ) << std::endl;
        }
    }*/

    /** Define and compute the values of the fields at the collocation points
      */

     //PrimaryFields defineFlds( ID, chebyshevValues );

    /** Define and compute everything on the initial hypersurface. The constructor of
        dependentFields calls primaryFields' one.
      */

     //DependentFields setupIH( ID, chebyshevValues, params );

    /** Evolve the spectral coefficients
      */

    BispecEvolve bispecEvolution( ID, chebyshevValues, params );

    //

    //bispecEvolution.

    /*std::cout <<
        "The time derivatives of the even fields arranged in the vector b,"
            << std::endl << std::endl;
    for( size_t field = 0; field < bispecEvolution.n_evenflds(); ++field )
    {
        std::cout << "These are the derivatives of the even field "
            << "(field,collocation_point), " << std::endl << std::endl;
        for( size_t n = 0; n < bispecEvolution.n_colpoints(); ++n )
        {
            std::cout << "(" << field << "," << n << ") = "
                      << bispecEvolution.even_evolution_eqs( 0, field, n )
                      << std::endl;
        }
        std::cout << std::endl;
    }*/

    /*std::cout <<
        "The time derivatives of the spectral coefficients on the initial hypersurface,"
            << std::endl << std::endl;
    for( size_t field = 0; field < bispecEvolution.n_evenflds(); ++field )
    {
        std::cout << "These are the derivatives of the even spectral coefficients of "
            << "(field,cheby_index), " << std::endl << std::endl;
        for( size_t cheby_index = 0; cheby_index < bispecEvolution.n_chebys(); ++cheby_index )
        {
            std::cout << "(" << field << "," << cheby_index << ") = "
                      << bispecEvolution.get_even_spec_t( 0, field, cheby_index ) << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
    for( size_t field = 0; field < bispecEvolution.n_oddflds(); ++field )
    {
        std::cout << "These are the derivatives of the odd spectral coefficients of "
            << "(field,cheby_index), " << std::endl << std::endl;
        for( size_t cheby_index = 0; cheby_index < bispecEvolution.n_chebys(); ++cheby_index )
        {
            std::cout << "(" << field << "," << cheby_index << ") = "
                      << bispecEvolution.get_odd_spec_t( 0, field, cheby_index ) << std::endl;
        }
        std::cout << std::endl;
    }*/

    //////////////////////////////////////////////////////////////////////////////////////

    /// The following is the last line in bim-solver. It should be the same here.
    /// - Evolve the equations of motion
    ///
    //return integrator.evolveEquations () ? 0 : -1;

    //////////////////////////////////////////////////////////////////////////////////////

    return 0;
}
