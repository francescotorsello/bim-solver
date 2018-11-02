/**
 *  @file      gowdy-solver.cpp
 *  @brief     IVP solver for Gowdy spacetimes that uses bimetric-solver framework.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 *
 *  @mainpage  IVP solver for Gowdy spacetimes
 *
 *  The code implements a pedagogical IVP solver for Gowdy spacetimes that uses
 *  <a href="http://qft.nu/_private_bim-ss">bimetric-solver framework</a>
 *  @cite BimSolver:2018. <p>
 *  The equations of motion come from @cite Garfinkle:2004tu.
 *
 *  @par       Main file:
 *             gowdy-solver.cpp
 *
 *  @version   0.3, inception 2018-04-19, last modified **2018-08-16 22:54**
 *  @author    Mikica Kocic
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
        P = GFCNT, Q, R, S, L,      //!< Evolved variables
        P_t, Q_t, R_t, S_t, L_t,    //!< Time derivatives of the evolved variables
        C,                          //!< Constraint violation
        gowdyLast                   //!< The total number of grid functions
    };
    #undef GFCNT // Clear the old number of grid functions then
    #define GFCNT fld::gowdyLast //!< Set the number of grid functions on a grid point

    /** The fields that will be read from the initial data.
     */
    static const std::vector<Int> gowdyInput = {
        P, Q, R, L, S
    };

    /** The fields which will be written to the output.
     */
    static const std::vector<Int> gowdyOutput = {
        P, Q, R, L, S,
        C
    };

    /** The grid functions that are evolved in time.
     */
    static const std::vector<EvolvedBy> gowdyEvolved = {
        { P, P_t }, { Q, Q_t }, { R, R_t },
        { S, S_t }, { L, L_t }
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

/** GowdyEvolve encapsulates the IVP solver for Gowdy spacetime.
 */
class GowdyEvolve
    : GridUser, public IntegFace
{
    /////////////////////////////////////////////////////////////////////////////////////
    // Fields (grid functions)
    /////////////////////////////////////////////////////////////////////////////////////

    emitField( P )  emitField(P_t)  emitDerivative_r( P )  emitDerivative_rr( P )
    emitField( Q )  emitField(Q_t)  emitDerivative_r( Q )  emitDerivative_rr( Q )
    emitField( R )  emitField(R_t)  emitDerivative_r( R )  emitDerivative_rr( R )
    emitField( S )  emitField(S_t)  emitDerivative_r( S )  emitDerivative_rr( S )
    emitField( L )  emitField(L_t)  emitDerivative_r( L )  emitDerivative_rr( L )
    emitField( C )

    /////////////////////////////////////////////////////////////////////////////////////
    // RHS of the time derivatives (equations of motion)
    /////////////////////////////////////////////////////////////////////////////////////

    Real eq_P_t( Int m, Int n ) {
        return R(m, n);
    }

    Real eq_Q_t( Int m, Int n ) {
        return S(m, n);
    }

    Real eq_R_t( Int m, Int n ) {
        return P_rr(m, n) / exp(2 * t(m,n))
            - exp(2 * (P(m, n) - t(m,n))) * pow2(Q_r(m, n))
            + exp(2 * P(m, n)) * pow2(S(m, n));
    }

    Real eq_S_t( Int m, Int n ) {
        return Q_rr(m, n) / exp(2 * t(m,n))
            + 2 * P_r(m, n) * Q_r(m, n) / exp(2 * t(m,n))
            - 2 * R(m, n) * S(m, n);
    }

    Real eq_L_t( Int m, Int n ) {
        return -pow2(P_r(m, n)) / exp(2 * t(m,n))
            - exp(2 * (P(m, n) - t(m,n))) * pow2(Q_r(m, n))
            - exp(2 * P(m, n)) * pow2(S(m, n))
            - pow2(R(m, n));
    }

    // Constraints violation

    Real eq_C_violation( Int m, Int n ) {
        return L_r(m, n) + 2 * P_r(m, n) * R(m, n)
            + 2 * exp(2 * P(m, n)) * Q_r(m, n) * S(m, n);
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
            P_t (m,n) = eq_P_t (m,n);
            Q_t (m,n) = eq_Q_t (m,n);
            R_t (m,n) = eq_R_t (m,n);
            S_t (m,n) = eq_S_t (m,n);
            L_t (m,n) = eq_L_t (m,n);

            // Calculate the constraint violation
            C(m,n) = eq_C_violation(m,n);
        }
    }

    /** At each checkpoint, outputs the grid row every mSkip steps, also
     *  reporting the integration time and checking for eventual NaNs.
     */
    virtual bool integStep_CheckForNaNs( Int m, Int nFrom, Int nTo )
    {
        for( Int n = nFrom; n < nTo; ++n ) {
            if ( std::isnan( GF( fld::P, m, nGhost + n ) ) ) {
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
            P(m,i) = P(m,nLen+i);
            Q(m,i) = Q(m,nLen+i);
            R(m,i) = R(m,nLen+i);
            S(m,i) = S(m,nLen+i);
            L(m,i) = L(m,nLen+i);
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
            P(m,nGhost+nLen+i) = P(m,nGhost+i);
            Q(m,nGhost+nLen+i) = Q(m,nGhost+i);
            R(m,nGhost+nLen+i) = R(m,nGhost+i);
            S(m,nGhost+nLen+i) = S(m,nGhost+i);
            L(m,nGhost+nLen+i) = L(m,nGhost+i);
        }
    }

public:

    /** The initial data is prepared by hand (thus not loaded from any file).
     */
    bool setupInitialData ()
    {
        static const double PI = acos(-1.0);

        for( Int n = 0; n < nLen; ++n )
        {
            double x = 2.0 * PI * n / nLen;
            r(0,nGhost+n) = x;
            P(0,nGhost+n) = 0;
            Q(0,nGhost+n) = cosl(x);
            R(0,nGhost+n) = 5 * cosl(x);
            S(0,nGhost+n) = 0;
            L(0,nGhost+n) = 0;
        }

        return true;
    }

    /** Creates and configures the Gowdy solver from the given parameters.
     */
    GowdyEvolve( Parameters& params, UniformGrid& ug, MoLIntegrator& integ )
        : GridUser( ug )
    {
        slog << "Gowdy Solver:" << std::endl << std::endl;

        // Sign up for the integration
        //
        integ.addToEvolution( this );

        // Add our grid functions to the evolution
        //
        integ.keepEvolved( fld::gowdyEvolved ); // Evolved
    }
};

/////////////////////////////////////////////////////////////////////////////////////////
// The main entry point of `gowdy-solver`.
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

    // Create the output sync
    //
    GridOutputWriter output ( params, ugrid );
    output.gridFunctions( fld::gowdyOutput );
    if ( ! output.open () ) return -1;

    // Evolve the equations from the manually setup ID
    //
    MoLIntegrator integ( params, ugrid, output );
    GowdyEvolve gowdy( params, ugrid, integ );
    gowdy.setupInitialData ();

    return integ.evolveEquations () ? 0 : -1;
}
