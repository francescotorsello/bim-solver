/**
 *  @file      wave-solver.cpp
 *  @brief     IVP solver for the wave equation that uses bimetric-solver framework.
 *  @author    Francesco Torsello, Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 *
 *  @mainpage  IVP solver for the wave equation
 *
 *  The code implements a pedagogical IVP solver for the wave equation that uses
 *  <a href="http://qft.nu/_private_bim-ss">bimetric-solver framework</a>
 *  @cite BimSolver:2018. <p>
 *
 *
 *  @par       Main file:
 *             wave-solver.cpp
 *
 *  @version   0, inception 2018-08-23, last modified **2018-08-23 14:28**
 *  @author    Francesco Torsello, Mikica Kocic
 *
 *  @par       References
 *  @copydoc   citelist
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <chrono>

#include "numMethods.h"        // The implemented numerical modules
#include "sys/paramsHolder.h"  // Holds 'key=value' pairs got from the parameter file
#include "sys/slog.h"          // For writing both to cerr and cout simultaneously
#include "sys/trackUsedTime.h" // Keep track of the elapsed time of the application
#include "sys/hpc.h"           // High-Performance Computing (HPC) support

#define PI 3.14159

/////////////////////////////////////////////////////////////////////////////////////////
// Declare our grid functions (their number must be known to the grid-driver beforehand)

#include "grid/gridFunctions.h"

/** Definitions of grid functions are here.
 */
namespace fld
{
    /** Identifiers of all known grid functions (the variables or fields on a grid).
     */
    enum Index
    {
        u = GFCNT, v,              //!< Evolved variables
        u_t,v_t,                   //!< Time derivatives of the evolved variables
        waveLast                   //!< The total number of grid functions
    };
    #undef GFCNT // Clear the old number of grid functions then
    #define GFCNT fld::waveLast     //!< Set the number of grid functions on a grid point

    /** The fields that will be read from the initial data.
     */
    static const std::vector<Int> waveInput = {
        u,v
    };

    /** The fields which will be written to the output.
     */
    static const std::vector<Int> waveOutput = {
        u,v
    };

    /** The grid functions that are evolved in time.
     */
    static const std::vector<EvolvedBy> waveEvolved = {
        { u, u_t }, { v, v_t },
    };
}

/////////////////////////////////////////////////////////////////////////////////////////
// The grid-driver and other modules

#include "grid/gridPoint.h"
#include "grid/gridDriver.h"
#include "grid/gridInitialData.h"
#include "grid/gridOutput.h"
#include "grid/integrator.h"

/////////////////////////////////////////////////////////////////////////////////////////

/** waveEvolve encapsulates the IVP solver for the wave equation.
 */
class waveEvolve
    : GridUser, public IntegFace
{
    /////////////////////////////////////////////////////////////////////////////////////
    // Fields (grid functions)
    /////////////////////////////////////////////////////////////////////////////////////

    emitField( u )  emitField(u_t)  emitDerivative_r( u )  emitDerivative_rr( u )
    emitField( v )  emitField(v_t)  emitDerivative_r( v )  emitDerivative_rr( v )

    /////////////////////////////////////////////////////////////////////////////////////
    // RHS of the time derivatives (equations of motion)
    /////////////////////////////////////////////////////////////////////////////////////

    Real force( Int m, Int n ) {
        return cos(m/PI)*exp(-pow(sin(1.*n)+cos(2.*PI*n),2.))+ sin(m/(sqrt(2.)*PI))/3.;
    }

    Real eq_u_t( Int m, Int n ) {
        return v(m, n);
    }


    Real eq_v_t( Int m, Int n ) {
        //return 1/r(m, n)*(2*u_r(m, n)+r(m, n)*u_rr(m, n));
        return u_rr(m, n)+force(m, n);
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // MoLInterface implementation
    /////////////////////////////////////////////////////////////////////////////////////

    /** Calculates the dependent variables from the prime state variables
     *  which are needed for the integration.
     */
    virtual void integStep_Prepare( Int m )
    {} // Unimplemented

    /** Calculate the right-hand side for time evolution.
     */
    virtual void integStep_CalcEvolutionRHS( Int m )
    {
        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            u_t (m,n) = eq_u_t (m,n);
            v_t (m,n) = eq_v_t (m,n);

        }
    }

    /** At each checkpoint, outputs the grid row every mSkip steps, also
     *  reporting the integration time and checking for eventual NaNs.
     */
    virtual bool integStep_CheckForNaNs( Int m, Int nFrom, Int nTo )
    {
        for( Int n = nFrom; n < nTo; ++n ) {
            if ( std::isnan( GF( fld::u, m, nGhost + n ) ) ) {
                return true;
            }
        }
        return false;
    }

    /** Nothing to do.
     */
    virtual void integStep_Finalize( Int mNext, Int mPrev )
    {}

    /** Apply the periodic boundary conditions on left.
     */
    virtual void applyLeftBoundaryCondition( Int m )
    {
        for( Int i = 0; i < nGhost; ++i )
        {
            t(m,i) = t(m,nLen+i);
            r(m,i) = r(m,nLen+i);
            u(m,i) = u(m,nLen+i);
            v(m,i) = v(m,nLen+i);
        }
    }

    /** Apply the periodic boundary conditions on right.
     */
    virtual void applyRightBoundaryCondition( Int m )
    {
        for( Int i = 0; i < nGhost; ++i )
        {
            t(m,nGhost+nLen+i) = t(m,nGhost+i);
            r(m,nGhost+nLen+i) = r(m,nGhost+i);
            u(m,nGhost+nLen+i) = u(m,nGhost+i);
            v(m,nGhost+nLen+i) = v(m,nGhost+i);
        }
    }

public:

    /** Initial data prepared directly
     */
    bool setupInitialData ()
    {
        //static const double PI = acos(-1.0);

        for( Int n = 0; n < nLen; ++n )
        {
            r(0,nGhost+n) = n;
            u(0,nGhost+n) = sin(n);
            v(0,nGhost+n) = cos(n);
        }

        return true;
    }

    /** Creates and configures the wave solver from the given parameters.
     */
    waveEvolve( Parameters& params, UniformGrid& ug, MoLIntegrator& integ )
        : GridUser( ug )
    {
        slog << "Wave Solver:" << std::endl << std::endl;

        // Sign up for the integration
        //
        integ.addToEvolution( this );

        // Add our grid functions to the evolution
        //
        integ.keepEvolved( fld::waveEvolved ); // Evolved
    }
};

/////////////////////////////////////////////////////////////////////////////////////////
// The main entry point of `wave-solver`.
//
int main( int argc, char* argv[] )
{
    TrackUsedTime timer;

    // Read the configuration
    //
    Parameters params( argc >= 2 ? argv[1] : "config.ini" );

    // Create the grid driver
    //
    UniformGrid ugrid( params );

    /** Read the initial Data
    //
     */
    /*GridInitialData input( params, ugrid );
    input.gridFunctions( fld::bimInput );
    if ( ! input.load () ) return -1;
     */

    // Create the output sync
    //
    GridOutputWriter output ( params, ugrid );
    output.gridFunctions( fld::waveOutput );
    if ( ! output.open () ) return -1;

    // Evolve the equations from the manually setup ID
    //
    MoLIntegrator integ( params, ugrid, output );
    waveEvolve wave( params, ugrid, integ );
    wave.setupInitialData ();

    return integ.evolveEquations () ? 0 : -1;
}
