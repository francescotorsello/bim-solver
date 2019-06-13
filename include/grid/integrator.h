/**
 *  @file      integrator.h
 *  @brief     Time evolution using the Method of Lines (MoL).
 *  @authors   Mikica Kocic, Francesco Torsello, Marc Compere, Thomas Treichl
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _INTEGRATOR_H_INCLUDED
#define _INTEGRATOR_H_INCLUDED

#include <csignal>

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g11 Numerical Methods                                                    */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
// The adaptive method to be used:
//
#define _DEBUG_ADPT 0

/** Encapsulates an adaptive stepper for Runge-Kutta methods.
 *
 *  The implementation is based on ode45.m that was originally written
 *  by Marc Compere (2001-07-03) under the GPL for octave, which was further adapted
 *  by Thomas Treichl (2006-08-10) to the new syntax that is used by the odepkg
 *  for octave to be compatible to LabMat's ode45.
 *
 *  References: @cite Dormand:1980rk, @cite Ascher:1998com, @cite Butcher:2008num
 *
 *  See also:
 *      T. Kimura. "On Dormand-Prince method". Retrieved April 27 (2009),
 *      <URL: http://depa.fquim.unam.mx/amyd/archivero/DormandPrince_19856.pdf>.
 */
class AdaptiveStepsizeControl : GridUser
{
    bool quit;          //!< Indicates whether the integration should be aborted
    bool stepSizeOK;    //!< Stepsize was OK (proceed with the next step)
    Int  nIterations;   //!< Current number of iterations (at the same stepsize)
    Int  nRejected;     //!< Number of rejected iteration steps
                        //   Configuration parameters
    bool enabled;       //!< True if the adaptive step method is used
    Real minStepsize;   //!< Minimum allowed `delta_t`
    Real maxStepsize;   //!< Maximum allowed `delta_t`
    Real absTolerance;  //!< A threshold below which the value of the solution
                        //!< component is unimportant. The absolute error tolerances
                        //!< determine the accuracy when the solution approaches zero.
                        //!< Note: `absTolerance = 10^{-AccuarcyGoal}`
    Real relTolerance;  //!< A measure of the error relative to the size of each solution
                        //!< component. Roughly, it controls the number of correct digits
                        //!< in all solution components, except those smaller than
                        //!< thresholds absTolerance.
                        //!< Note: `relTolerance = 10^{-PrecisionGoal}`
    Int maxIterations;  //!< Number of allowed (re)iterations at the same `Y^{m}`

    /** Floating-point relative accuracy.
     *  The distance from 1.0 to the next larger number; e.g., 2^(-52) for double.
     */
    Real eps = std::numeric_limits<double>::epsilon();

public:

    /** Enables usage of the adaptive method.
     */
    void enable () {
        enabled = true;
    }

    /** Returns true if the adaptive method is enabled.
     */
    bool isEnabled () const {
        return enabled;
    }

    /** Indicates whether the integration should be aborted.
     */
    bool shouldQuit () const {
        return quit;
    }

    /** Returns true if the last integration step should be repeated/(re)iterated.
     */
    bool repeatStep () const {
        return ! stepSizeOK;
    }

    /** Indicates weather the next-to-last MoL substep should to be propagated
     *  instead of the last MoL substep.
     */
    bool propagateNextToLastSubstep () const {
        return false;
    }

    /** Constructs the adaptive step method as specified in the parameter file.
     */
    AdaptiveStepsizeControl( Parameters& params, UniformGrid& ug )
        : GridUser( ug )
        , quit( false )
        , stepSizeOK( false )
        , nIterations( 0 )
        , nRejected( 0 )
    {
        params.get( "adpt.enabled",         enabled,        false        );
        params.get( "adpt.minStepsize",     minStepsize,    1e3 * eps    );
        params.get( "adpt.maxStepsize",     maxStepsize,    delta_r/2.   );
        params.get( "adpt.maxIterations",   maxIterations,  1000         );
        params.get( "adpt.absTolerance",    absTolerance,   1e-8         );
        params.get( "adpt.relTolerance",    relTolerance,   1e-8         );
    }

    /** Shows the configuration parameters.
     */
    void showParameters ()
    {
        slog << "Using adaptive stepper:" << std::endl << std::endl
             << "    minStepsize = "  << minStepsize
             << ",  maxStepsize = "   << maxStepsize << std::endl
             << "    absTolerance = " << absTolerance
             << ",  relTolerance = "  << relTolerance
             << ",  maxIterations = " << maxIterations
             << std::endl << std::endl;
    }

    /** Adjusts the stepsize and updates the indicators weather the integrator should
     *  continue, quit or repeat (reiterate) the current MoL step.
     *  @return New delta_t
     */
    Real getNewStepSize
    (
        Int m,        //!< Current slice
        Real delta_t, //!< The current stepsize
        Real err,     //!< The absolute local truncation error
        Real gfNorm,  //!< The infinity norm of Y^{m} at the beginning of the step
        Int  pOrder   //!< The RK order
        )
    {
        /// Procedure:
        /// - Calculate the acceptable error `tau`
        ///
        Real tau = /*std::*/MAX( absTolerance, relTolerance * /*std::*/MAX( 1.0, gfNorm ) );

        /// If the error is acceptable then update `Y[m+1]` by taking  the higher order
        /// estimation as "local extrapolation".
        /// (Acceptable means `actual error err <= allowable error tau`.)
        ///
        stepSizeOK = ( err <= tau );
        if( stepSizeOK ) {
            nIterations = 0;
        }
        else {
            ++nRejected; // Keep statistics of unacceptable iterations
        }

        /// - Update the step size for the next integration step.
        ///   See p.91 in Ascher & Petzold @cite Ascher:1998com
        ///
        Real resize_factor = 0.8 * pow( tau/err, 1./pOrder );
        Real suggested_dt = delta_t * resize_factor;

        /// Do not change delta_t if the error was acceptable and resize_factor < 90%
        ///
        if ( ! stepSizeOK || resize_factor >= 0.9 ) {
            delta_t = std::min( maxStepsize, suggested_dt );
        }

        /// @todo optionally quit (the alternative is to disable adaptive method)
        ///
        if( delta_t < minStepsize )
        {
            quit = true;
            slog << "Stepsize " << delta_t << " is smaller than allowed " << minStepsize
                 << ". Try to reduce the initial delta_t and/or maxStepsize."
                 << std::endl;
        }
        else if( ++nIterations > maxIterations )
        {
            quit = true;
            slog << "Maximum number of iterations " << maxIterations
                 << " has been reached. "
                 << "Try to reduce the initial delta_t and/or maxStepsize."
                 << std::endl;
        }

        /// - Track the statistics (as a special GF @ref fld::sysStat)
        ///
        GF( fld::sysStat, m, 0 ) = err;
        GF( fld::sysStat, m, 1 ) = gfNorm;
        GF( fld::sysStat, m, 2 ) = tau;
        GF( fld::sysStat, m, 3 ) = resize_factor;
        GF( fld::sysStat, m, 4 ) = suggested_dt;
        GF( fld::sysStat, m, 5 ) = delta_t;
        GF( fld::sysStat, m, 6 ) = nIterations;
        GF( fld::sysStat, m, 7 ) = nRejected;

        #if _DEBUG_ADPT
        slog << "err: " << err << ", gfNorm: " << gfNorm << ", tau: " << tau
             << ", rsz: " << resize_factor
             << ", sdt: " << suggested_dt << ", dt: " << delta_t
             << ", " << ( stepSizeOK ? "OK" : "reiterate" )
             << ", niter: " << nIterations << ", nrej: " << nRejected
             << std::endl;
        #endif

        return delta_t;
    }

    void showInfo( Int m )
    {
        std::cout << ",  sdt = " << GF( fld::sysStat, m, 4 ) // suggested delta_t
                  << ",  dt = " << GF( fld::sysStat, m, 5 ); // delta_t
    }
};

/** Time evolution using the Method of Lines (MoL).
 *
 *  The implemented methods are declared in MoL::knownMethods.
 *
 *  The integrator can evolve equations from several modules at the same time.
 *  All such modules should implement @ref IntegFace ans sign-in to the integrator.
 *  The modules should also specify which grid functions are either evolved or
 *  kept constant by the integrator. (Those grid functions that are unknown to
 *  the integrator have to be maintained by the the modules themselves.)
 *
 *  The system provides the coordinates `t` and `r` as global grid functions
 *  (they are defined in @ref fld::sysIndex). The integrator keeps time `t` evolved and
 *  the spatial coordinate `r` constant.
 *
 *  The integrator handles process signals (e.g., Ctrl-C), which will cause
 *  a premature end of the integration.
 */
class MoL : GridUser
{
    std::vector<IntegFace*> eomList;       //!< Equations of motion that we evolve
    std::vector<Int> constantGF;           //!< Grid functions that are kept constant
    std::vector<fld::EvolvedBy> evolvedGF; //!< Grid functions that are integrated in time
    Int n_evolved = evolvedGF.size();

    GridOutputWriter* output; //!< The integration output is directed here

    std::string methodID;     //!< The integration method identifier
    std::string updateJ_ID;   //!< The updateJacobian identifier
    Int      method;          //!< The method code (Euler, ICN, MoL, ...)
    time_hr  rt_start;        //!< The realtime when the integration started
    bool     running;         //!< Indicates if the integration should proceed
    Real     t_0;             //!< The initial `t`.
    Real     t_1;             //!< Integrate up to this `t`.
    Real     delta_t;         //!< The integration step
    Real     dissip;          //!< Kreiss-Oliger dissipation coefficient.
    Real     dissip_delta_r;  //!< `dissip / delta_r`
    Real     cur_t;           //!< Current time
    Real     tolerance;       //!< Tolerance required by the Newton method in DIRKs
    Int      updateJ;         //!< Parameter from config.ini setting the update of the
                              //!< Jacobian in DIRKs
    Int      modulo;          //!< When the Jacobian is updated after n steps, modulo = n

    /////////////////////////////////////////////////////////////////////////////////////
    // The Jacobian of the evolution equations
    /////////////////////////////////////////////////////////////////////////////////////

    MatReal Jacobian;

    enum UpdateJ         //!< Possibilities to update the Jacobian in DIRK
    {
        STEP             = 0,   // Update the Jacobian at every step
        STAGE            = 1,   // Update the Jacobian at every stage
        ITERATION        = 2,   // Update the Jacobian at every iteration
        MULTIPLE_STEPS   = 3,   // Update the Jacobian after multiple steps
    };

    /////////////////////////////////////////////////////////////////////////////////////
    //
    /////////////////////////////////////////////////////////////////////////////////////


    AdaptiveStepsizeControl adpt;  //!< Dynamically adjusts delta_t

    /** A list of the implemented integration methods.
     */
    static std::map<std::string,int> knownMethods;

    /** Keep the list of all created integrator objects (used by signalHandler).
     */
    static std::vector<MoL*> knownIntegrators;

    /** Forward the signal to all integrators (signaling them to quit).
     */
    static void signalHandler( int signum )
    {
        std::cout << std::endl;
        std::cerr << "*** Signal " << signum << " received" << std::endl;

        if( knownIntegrators.size () == 0 ) {
            exit( -1 );
        }
        for( auto i: knownIntegrators ) {
            i->quit ();
        }
    }

   /** Returns max(abs(GF)) for all evolved GF.
    */
    Real getInfNormOfEvolvedGF( Real m )
    {
        Real norm = 0;

        _Pragma("omp parallel for reduction(max : norm)")

        for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            for( auto e: evolvedGF )
            {
                Real absgf = abs( GF( e.f, m, n ) );
                if( absgf > norm ) {
                    norm = absgf;
                }
            }
        }

        return norm;
    }

    /** Find the minimum of an array of length n_evolved
    */
    Real findMinimum( VecReal _array )
    {
        Real minimum = _array[0];
        for ( Int i = 0; i < n_evolved; ++i )
        {
            if ( _array[i] <= minimum )
            {
                minimum = _array[i];
            }
        }
        return minimum;
    }

   /** Returns an overall absolute error estimate between two time slices.
    *  It is calculated as max(abs(difference)) which is norm(difference,infinity).
    */
    Real getErrorEstimate( Real m1, Real m2 )
    {
        Real max_err = 0;

        _Pragma("omp parallel for reduction(max : max_err)")

        for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            for( auto e: evolvedGF )
            {
                Real cur_err = abs( GF( e.f, m1, n ) - GF( e.f, m2, n ) );
                if( cur_err > max_err ) {
                    max_err = cur_err;
                }
            }
        }

        return max_err;
    }

    /** Propagates (copies) evolved GFs from one time slice to another.
     */
    void copyEvolvedGF( Int mTo, Int mFrom )
    {
        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            for( auto e: evolvedGF ) {
                GF( e.f, mTo, n ) = GF( e.f, mFrom, n );
            }
        }
    }

    /** Shows the configuration parameters.
     */
    void showParameters ()
    {
        // Expected output file data size (in bytes)
        //
        const Real sz = ( ( t_1 - t_0 ) / delta_t / output->get_mSkip() + 1 )
                           * output->recordSizeInBytes();

        slog << "Integrator:" << std::endl << std::endl
             << "    t_0 = " << t_0 << ",  t_1 = " << t_1 << ",  delta_t = " << delta_t
             << ",  method = " << methodID << " (#" << method << ")" << std::endl
             << "    Kreiss-Oliger dissipation = " << dissip
             << ", order = " << KO_ORDER << std::endl
             << "    Expected output data size = "
             << round( ( sz < 1e3 ? sz : sz < 1e6 ? sz/1e3 : sz/1e6 ) * 10 ) / 10
             << ( sz < 1e3 ? "" : sz < 1e6 ? " kB" : " MB" )
             << ",  GF count = " << output->GFCount ()
             << std::endl << std::endl;

        if( false && delta_r != r(0,nGhost+1) - r(0,nGhost) )
        {
            slog << "    Error: delta_r = " << delta_r << " != grid spacing = "
                 << r(0,nGhost+1) - r(0,nGhost)
                 << ", diff = " << delta_r - ( r(0,nGhost+1) - r(0,nGhost) )
                 << std::endl;
        }
    }

    /** Report the integration time and the estimated real-time to finish.
     */
    void reportIntegrationTime( Int m )
    {
        if ( mpiRank() > 1 ) {
            return; // only the first in rank can access the standard output
        }

        Real tnow = GF( fld::t, m, nGhost ); // current integration time

        auto elapsed = 1e-3 * std::chrono::duration_cast<std::chrono::milliseconds>
                ( std::chrono::high_resolution_clock::now () - rt_start ).count ();

        auto rt_factor = elapsed / ( tnow - t_0 );  // real time/integration time
        auto rt_total = round( ( t_1 - t_0 ) * rt_factor ); // projected total time
        auto rt_left = round( ( t_1 - tnow ) * rt_factor ); // estimated time to cover

        /// @todo convert to sprintf()

        std::cout << "    t = " << std::setw(6) << std::setprecision(2) << std::fixed
             << tnow << std::setw(-1) << std::setprecision(-1) << std::defaultfloat
             << ",   real = " << round( elapsed )
             << " s   (" << round( 100 * ( tnow - t_0 ) / ( t_1 - t_0 ) )
             << "% of " << rt_total
             << " s),   ETA " << rt_left << " s";

        if( adpt.isEnabled () ) {
            adpt.showInfo( m );
        }

        std::cout << "                        \r";
        std::flush( std::cout );
    }

    /** Executed at the beginning of each integration step.
     *  Eventually calculates the RHS (`f_t`) needed for the integration.
     *  Alternatively, clears `f_t` for all evolved GFs if `zeroRHS` is true.
     */
    void integStep_Begin( Int m, bool zeroRHS = false )
    {
        if( zeroRHS )
        {
            OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
            {
                for( auto e: evolvedGF ) {
                    GF( e.f_t, m, n ) = 0;
                }
            }
            return;
        }

        /// - Prepare the intermediate variables.
        ///
        for( auto eom: eomList ) {
            eom->integStep_Prepare( m );
        }

        /// - Exchange boundaries with the neighboring MPI ranks and/or fix the
        ///   left/right boundary conditions (depending on our MPI rank).
        ///   Quit if MPI has failed or we received abort message.
        ///
        if( gridDriver->mpiSize() > 1
            && ! gridDriver->exchangeBoundaries( m ) )
        {
            running = false;
            return;
        }
        if ( gridDriver->isFirstInRank () ) {
            for( auto eom: eomList ) {
                eom->applyLeftBoundaryCondition( m );
            }
        }
        if ( gridDriver->isLastInRank () ) {
            for( auto eom: eomList ) {
                eom->applyRightBoundaryCondition( m );
            }
        }

        /// - Calculate the RHS of time evolution equations.
        ///
        for( auto eom: eomList ) {
            eom->integStep_CalcEvolutionRHS( m );
        }
    }

    /** Executed at the end of each integration step.
     *  Returns `false` if we need to break the integration loop.
     */
    void integStep_End( Int m1, Int m )
    {
        /// - Copy the grid functions that are kept constant to the next slice m1.
        ///
        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            for( auto f: constantGF ) {
                GF( f, m1, n ) = GF( f, m, n );
            }
        }

        /// - Finalize the integration step (e.g., doing diagnostic and post
        ///   processing of EoM like calculating the constraint violations).
        ///
        for( auto eom: eomList ) {
            eom->integStep_Finalize( m1, m );
        }

        /// - Wait on finishing the exchange of the boundaries (when using MPI).
        ///
        if ( ! gridDriver->waitExchange () ) {
            quit( m ); // Quit if MPI has failed
        }
    }

    /** Executed at each checkpoint.
     */
    void integStep_Checkpoint( Int m )
    {
        /// - Report the elapsed (both integration and real) time.
        ///
        reportIntegrationTime( m );

        /// - Check for eventual NaNs in the central/lower region of the grid.
        ///
        for( auto eom: eomList )
        {
            if ( ! eom->integStep_Diagnostics( m,
                            output->get_nFrom (), output->get_nTo () ) ) {
                if ( running ) { // quit only once
                    quit( m );
                }
            }
        }

        /// - Send the grid slice to the output.
        ///
        output->write( m, cur_t, delta_t );
    }

    /** Evolves the state variables using the FT scheme.
     */
    void integStep_Euler( Int mNext, Int m );

    /** Evolves the sate variables using the ICN averaging.
     */
    void integStep_AvgICN( Int mNext, Int mMid, Int m );

    /** Time evolution using the iterated Crank-Nicolson (ICN).
     */
    void integrate_AvgICN( int ICN_steps );

    /** Prepares the state variables for the intermediate MoL steps.
     */
    void MoL_BetaInit( Int m1, Int m, Real beta_delta_t );

    /** Accumulates a sum of the intermediate method of lines (MoL) steps.
     */
    void MoL_AlphaSum( Int m1, Int m, Real alpha );

    /** Compute the LU decomposition of the Newton iteration matrix
     */
    void computeNewtonIterationMatrix(
                                        Int stage_i,
                                        const ButcherTable& BT,
                                        MatReal& NewItMat
                                     );

    /** Compute the stage value for a given DIRK method
     */
    void DIRK_computeStage   (
                                Int m,                     // time step
                                Int stage_i,               // stage to be computed
                                const ButcherTable& BT,    // the DIRK's Butcher table
                                LUDecomposition& Newton,   // Newton iteration matrix
                                Real tol                   // Newton method's tolerance
                             );

    /** Compute the step value for a given DIRK method
     */
    void DIRK_computeStep    (
                                Int m,                      // time step to be computed
                                const ButcherTable& BT,     // the DIRK's Butcher table
                                LUDecomposition& Newton,    // Newton iteration matrix
                                Real tol                    // Newton method's tolerance
                             );

    /** Time evolution using the method of lines (MoL).
     */
    void integrate_MoL( const MoLDescriptor& MoL, bool adaptiveStepSize = false );
    void integrate_MoL( const ButcherTable& BT, Real tol, bool adaptiveStepSize = false );

public:

    /** Constructs the MoL integrator as specified in the parameter file.
     */
    MoL( Parameters& params, UniformGrid& ug, GridOutputWriter& out )
        : GridUser( ug ), output( &out ), adpt( params, ug ),
          Jacobian( n_evolved, n_evolved, Real(0) )
    {

        //#include "jacobian-BSSN/DIRK_Jacobian_cBSSN.h"

        constantGF.reserve( 64 );
        evolvedGF.reserve( 128 );

        static std::map<std::string,int> updateJacobian =
        {
            { "step",      STEP      },     { "stage",         STAGE          },
            { "iteration", ITERATION },     { "multipleSteps", MULTIPLE_STEPS }
        };

        constantGF.push_back( fld::r ); // Radial coordinate is kept constant by default

        params.get( "integ.t_0",     t_0,      0.0 );
        params.get( "integ.t_1",     t_1,      1.0 );
        params.get( "integ.delta_t", delta_t,  0.1 );
        params.get( "integ.dissip",  dissip,   0.0 );

        methodID = params.get( "integ.method", method, 0, knownMethods );
        updateJ_ID = params.get( "DIRK.updateJ", updateJ, STEP, updateJacobian );

        if( updateJ = MULTIPLE_STEPS )
        {
            params.get( "DIRK.modulo",  modulo,   5 );
        }

        // Fix the sign of delta_t
        //
        delta_t = t_1 < t_0 ? -abs(delta_t) : abs(delta_t);

        cur_t = t_0;  // The initial `t`
        running = false; // Stopped integration by default (will be enabled later on)

        // Calculate Kreiss-Oliger dissipation corrected delta r
        //
        dissip_delta_r = dissip / delta_r;

        // Register to receive the signals from the system
        //
        knownIntegrators.push_back( this );
        if ( knownIntegrators.size () == 1 ) {
            signal( SIGINT, signalHandler );
            signal( SIGTERM, signalHandler );
        }
    }

    /** The destructor (here only responsible for cleaning up standard output).
     */
    ~MoL()
    {
        if( mpiRank() == 0 ) {
            std::cout << std::endl << std::endl;
        }
    }

    /** Moves the integration final time to the current time and marks
     *  the integration stopped.
     */
    void quit( Int m = -1 )
    {
        running = false; // Stop the integration

        t_1 = m < 0 ? cur_t : GF( fld::t, m, nGhost );

        if( mpiRank() == 0 ) std::cout << std::endl;
        slog << "*** Quitting at t = " << t_1 << std::endl;

        /// - In case of NaNs, signal the other MPI ranks that we are quitting.
        /// @todo Distribute the abort msg to all ranks instead of calling MPI_Abort.
        ///
        gridDriver->abortExchange ();
    }

    /** Adds the EoM to the list of all evolved subsystems.
     */
    void addToEvolution( IntegFace* eomInterface )
    {
        eomList.push_back( eomInterface );
    }

    /** The given fields will be kept constant by the integrator.
     */
    void keepConstant( const std::vector<Int>& gfs )
    {
        constantGF.insert( constantGF.end(), gfs.begin(), gfs.end() );
    }

    /** Add the given functions to the evolution list.
     */
    void keepEvolved( const std::vector<fld::EvolvedBy>& gfs )
    {
        evolvedGF.insert( evolvedGF.end(), gfs.begin(), gfs.end() );
    }

    /** Returns the time step.
     */
    Real dt () const {
        return delta_t;
    }

    /** Integrates the evolution equations in the given grid realm.
     */
    bool evolveEquations ()
    {
        showParameters ();

        output->writeHeader ();

        // Start timing the real-time
        rt_start = std::chrono::high_resolution_clock::now ();

        switch( method )
        {
            case  0:  integrate_AvgICN ( 0    );  break;
            case  1:  integrate_AvgICN ( 1    );  break;
            case  2:  integrate_AvgICN ( 2    );  break;
            case  3:  integrate_AvgICN ( 3    );  break;

            case  4:  integrate_MoL( ICN2 );  break;
            case  5:  integrate_MoL( ICN3 );  break;
            case  6:  integrate_MoL( RK1  );  break;
            case  7:  integrate_MoL( RK2  );  break;
            case  8:  integrate_MoL( RK3  );  break;
            case  9:  integrate_MoL( RK4  );  break;

            case 10:  integrate_MoL( RK5DP_7M, /* adaptive = */ false ); break;
            case 11:  integrate_MoL( RK5DP_7M, /* adaptive = */ true  ); break;
            case 12:  integrate_MoL( RK5DP_7S, /* adaptive = */ false ); break;
            case 13:  integrate_MoL( RK5DP_7S, /* adaptive = */ true  ); break;

            //////////////////////////////////////////////////////////////////////////////
            // Below here, DIRK
            //////////////////////////////////////////////////////////////////////////////

            case 14:  integrate_MoL( ESDIRK32, tolerance, /* adaptive = */ false ); break;
        }

        return true;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of MoL static members
/////////////////////////////////////////////////////////////////////////////////////////

std::vector<MoL*> MoL::knownIntegrators;

std::map<std::string,int> MoL::knownMethods =
{
    { "Euler",     0 },   { "avgICN2",    1 },
    { "avgICN3",   2 },   { "avgICN4",    3 },
    { "ICN2",      4 },   { "ICN3",       5 },
    { "RK1",       6 },   { "RK2",        7 },
    { "RK3",       8 },   { "RK4",        9 },
    { "RK5DP_7M", 10 },   { "RK5DP_7MA", 11 },
    { "RK5DP_7S", 12 },   { "RK5DP_7SA", 13 },
    { "ESDIRK32", 14 }
};

/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of MoL explicit methods
/////////////////////////////////////////////////////////////////////////////////////////

void MoL::integStep_Euler( Int m1, Int m )
{
    OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
    {
        GF( fld::t, m1, n ) = GF( fld::t, m, n ) + delta_t; // Evolve time

        for( auto e: evolvedGF )
        {
            GF( e.f, m1, n ) = GF( e.f, m, n ) + delta_t * GF( e.f_t, m, n )
                             - KreissOligerTerm( e.f, delta_t );
        }
    }
}

void MoL::integStep_AvgICN
(
    Int m2,  // final time
    Int m1,  // middle time
    Int m    // initial time
    )
{
    const Real half_dt = delta_t / 2;

    OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
    {
        GF( fld::t, m2, n ) = GF( fld::t, m, n ) + delta_t; // Evolve time

        for( auto e: evolvedGF )
        {
            GF( e.f, m2, n ) = GF( e.f, m, n )
                    + half_dt * ( GF( e.f_t, m1, n) + GF( e.f_t, m, n ) )
                    - KreissOligerTerm( e.f, delta_t );
        }
    }
}

/** Integrates the evolution equations using the iterated Crank-Nicolson (ICN) method.
 *  If ICN_steps == 0, then the method effectively becomes FTCS.
 *  See @cite Press:2007nrec, @cite Baumgarte:2010nr
 */
void MoL::integrate_AvgICN( int ICN_steps )
{
    slog << "Employing Method of Lines, Iterative Crank-Nicolson (avg)"
         << ", N = " << ICN_steps << std::endl << std::endl;

    cur_t = t_0;  // The initial `t`
    running = true; // Enable integration

    for( Int n = 0; n < nTotal; ++n ) {
        GF( fld::t, 0, n ) = t_0;
    }

    Int  mStep = 0; // The step counter (used to filter checkpoints)
    while( running && abs(cur_t) < abs(t_1) )
    {
        for( Int m = 0; running && m < mLen && abs(cur_t) < abs(t_1); ++m )
        {
            Int mNext = m + 1 >= mLen ? 0 : m + 1;

            /////////////////////////////////////////////////////////////////////////////
            // Predictor. Use FTCS for the first step.
            // Note that this will be the only step, if ICN_steps == 0.
            //
            Int m_1 = ICN_steps >= 1 ? mLen + 1 : mNext;

            integStep_Begin ( m      );   if( ! running ) break;
            integStep_Euler ( m_1, m );
            integStep_End   ( m_1, m );

            /////////////////////////////////////////////////////////////////////////////
            // Corrector iterations. There will be no iterations, if ICN_steps == 0.
            //
            for( int i = 0; running && i < ICN_steps; ++i )
            {
                Int m_2 = ( i + 1 >= ICN_steps ? mNext : m_1 + 1 );

                integStep_Begin  ( m_1         );   if( ! running ) break;
                integStep_AvgICN ( m_2, m_1, m );
                integStep_End    ( m_2, m_1    );

                m_1 = m_2;
            }

            if( ( mStep++ % output->get_mSkip () ) == 0 ) {
                integStep_Checkpoint( m );
            }

            cur_t = GF( fld::t, m, nGhost ); // The current time step
        }
    }
}

void MoL::MoL_AlphaSum( Int m1, Int m, Real alpha )
{
    if( alpha != 0 )
    {
        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            GF( fld::t, m1, n ) += alpha * GF( fld::t, m, n ); // Evolve time

            for( auto e: evolvedGF ) {
                GF( e.f, m1, n ) += alpha * GF( e.f, m, n );
            }
        }
    }
}

void MoL::MoL_BetaInit( Int m1, Int m, Real beta_delta_t )
{
    if( beta_delta_t == 0 )
    {
        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            GF( fld::t, m1, n ) = 0;

            for( auto e: evolvedGF ) {
                GF( e.f, m1, n ) = 0;
            }
        }
    }
    else
    {
        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            GF( fld::t, m1, n ) = beta_delta_t; // Evolve time

            for( auto e: evolvedGF ) {
                GF( e.f, m1, n ) = beta_delta_t * GF( e.f_t, m, n )
                                 - KreissOligerTerm( e.f, beta_delta_t );
            }
        }
    }
}

/** The method of lines (MoL) separates the time integration from the rest of
 *  an evolution scheme. A `N`-step MoL method evolves the equation:
 *  <pre>
 *       partial_t Y = L[Y],  </pre>
 *
 *  using the following algorithm:
 *  <pre>
 *      Y^{(0)}   = Y^{m},
 *      Y^{(i+1)} = Sum_{j=0}^{i} ( alpha_{ij} Y^{(j)} ) +            <- MoL_AlphaSum
 *                  Delta t beta_i L[Y^{(i)}],   for i = 0, ..., N-1, <- MoL_BetaInit
 *      Y^{m+1}   = Y^{(N)}.    </pre>
 *
 *  The variables `Y^{(i)}`, `i = 0, ..., N`, are intermediate.
 *
 *  The method is completely specified by `N`, and arrays alpha and beta
 *  in the structure MoLDescriptor.
 */
void MoL::integrate_MoL( const MoLDescriptor& MoL, bool adaptiveStepSize )
{
    if ( adaptiveStepSize )
    {
        adpt.enable ();
        adpt.showParameters ();
    }

    slog << "Employing Method of Lines, " << MoL.name
         << ", N = " << MoL.N << std::endl << std::endl;

    cur_t = t_0;  // The initial `t`
    running = true; // Enable integration

    for( Int n = 0; n < nTotal; ++n ) {
        GF( fld::t, 0, n ) = t_0;
    }

    // The array mk[] contains redirections to grid rows (m's):
    //   mk[0]                   points to Y^{m}   == Y^{(0)},
    //   mk[1], ..., mk[Mol.N-1] point to the intermediate Y^{(1)}, ..., Y^{(N-1)},
    //   mk[MoL.N]               points to Y^{m+1} == Y^{(N)}.
    //
    int mk[ 20 ]; // int mk[ MoL.N + 1 ]; // Note that N is at most 16
    for( int i = 0; i <= MoL.N + 1; ++i ) {
        mk[i] = mLen + i;
    }

    Int  mStep = 0; // The step counter
    while( running && abs(cur_t) < abs(t_1) )
    {
        for( Int m = 0; running && m < mLen && abs(cur_t) < abs(t_1); /* nothing */ )
        {
            mk[0] = m; // Y^{(m)}
            mk[MoL.N] = m + 1 >= mLen ? 0 : m + 1; // Y^{(m+1)}

            for( int i = 0; running && i < MoL.N; ++i )
            {
                integStep_Begin( mk[i] );   if( ! running ) break;

                MoL_BetaInit( mk[i+1], mk[i], MoL.beta[i] * delta_t );

                for( int j = 0; j <= i; ++j ) {
                    MoL_AlphaSum( mk[i+1], mk[j], MoL.alpha[i][j] );
                }

                integStep_End( mk[i+1], mk[i] );
            }

            /////////////////////////////////////////////////////////////////////////////
            //
            if( adaptiveStepSize )
            {
                // Do an alternative last step, i = N - 1, but slightly modified.
                // Note:
                //      mk[i+2], beta2[i]    in MoL_BetaInit,
                //      mk[i+1], alpha[i+1]  in MoL_AlphaSum, and
                //      mk[i+2]              in integStep_End
                //
                Int i = MoL.N - 1;
                integStep_Begin( mk[i] );   if( ! running ) break;

                MoL_BetaInit( mk[i+2], mk[i], MoL.beta2[i] * delta_t );

                for( int j = 0; j <= i; ++j ) {
                    MoL_AlphaSum( mk[i+2], mk[j], MoL.alpha[i+1][j] );
                }

                integStep_End( mk[i+2], mk[i] );

                Real error  = getErrorEstimate( mk[MoL.N+1], mk[MoL.N] );
                Real gfNorm = getInfNormOfEvolvedGF( m );

                delta_t = adpt.getNewStepSize( m, delta_t, error, gfNorm, MoL.pOrder );

                if( adpt.shouldQuit () ) {
                    quit ();
                    break; // Should quit (if, e.g., max nr of iterations was achieved)
                }
                else if( adpt.repeatStep () ) {
                    continue; // Repeat the same MoL step (with the new delta_t)
                }
                else if( adpt.propagateNextToLastSubstep () ) {
                    copyEvolvedGF( mk[MoL.N], mk[MoL.N-1] );
                    // ... then we take Y^{(m+1)} = mk[MoL.N]
                }
                else {
                    // default: Y^{(m+1)} = mk[MoL.N]
                }
            }
            /////////////////////////////////////////////////////////////////////////////

            if( ( mStep++ % output->get_mSkip () ) == 0 ) {
                integStep_Checkpoint( m );
            }

            cur_t = GF( fld::t, m, nGhost ); // Get `t` from the current slice
            ++m; // Next time slice

            /*OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
            {
                for( auto e: evolvedGF )
                {
                    if( GF( e.f, m, n ) < 1e-6 )
                    {
                        GF( e.f, m, n ) = 0;
                    }
                }
            }*/
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// Definitions of MoL DIRK methods
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/** References
 *  [1] @cite Kennedy, Christopher A., Carpenter, Mark H.,
 *            "Diagonally Implicit Runge-Kutta Methods for Ordinary Differential
 *            Equations. A Review."
 *            https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20160005923.pdf
 */

/** Overloading of integrate_MoL.
 *  An implicit Runge-Kutta method (IRK) requires a different input than an explicit one.
 *  Hence, we overload integrate_MoL in order to accept the needed input.
 *  integrate_MoL updates the Jacobian at each time step, if required, and uses
 *  MoL_DIRK_computeStep, defined below.
 */
void MoL::integrate_MoL(
                         const ButcherTable& BT,              // the DIRK's Butcher table
                         Real                tol,             // Newton method's tolerance
                         bool                adaptiveStepSize
                       )
{
    // First computation of the Newton iteration matrix and its LU decomposition
    MatReal NewtonItMat ( n_evolved, n_evolved, Real(0) );

    /*
        Note that if the Jacobian is updated after every step (or
        after many steps), we can update it considering any
        of the values on the diagonal of the Butcher table. If the
        method is a SDIRK, then this choice doesn't affect the result.
        Otherwise, it does. Nonetheless, one can choose anyone of the
        diagonal elements.
    */
    computeNewtonIterationMatrix( /* stage = */ 3, BT, NewtonItMat );
    LUDecomposition Newton( NewtonItMat );

    // Adaptive stepsize control

    if ( adaptiveStepSize )
    {
        adpt.enable ();
        adpt.showParameters ();
    }

    slog << "Employing Method of Lines, " << BT.name
         << ", s = " << BT.s << std::endl << std::endl;

    cur_t = t_0;    // The initial `t`
    running = true; // Enable integration

    for( Int n = 0; n < nTotal; ++n ) {
        GF( fld::t, 0, n ) = t_0;
    }

    Int  mStep = 0; // The step counter

    while( running && abs(cur_t) < abs(t_1) )
    {
        for( Int m = 0; running && m < mLen && abs(cur_t) < abs(t_1); /* nothing */ )
        {

            if( updateJ == STEP )
            {
                computeNewtonIterationMatrix( 3, BT, NewtonItMat );
                Newton.LUDecompose( NewtonItMat );
            }
            else if( updateJ == MULTIPLE_STEPS && mStep % modulo == 0 )
            {
                computeNewtonIterationMatrix( 3, BT, NewtonItMat );
                Newton.LUDecompose( NewtonItMat );
            }

            DIRK_computeStep( m, BT, Newton, tol );

            //////////////////////////////////////////////////////////////////////////////

            if( ( mStep++ % output->get_mSkip () ) == 0 ) {
                integStep_Checkpoint( m );
            }

            cur_t = GF( fld::t, m, nGhost ); // Get `t` from the current slice
            ++m; // Next time slice
        }
    }
}

/** MoL_DIRK_computeStep updates the Jacobian at each stage, if required, and computes all
 *  the stages within a time step, using MoL_DIRK_computeStage, defined below.
 */
void MoL::DIRK_computeStep   (
                                Int                 m,       // time step to be computed
                                const ButcherTable& BT,      // the DIRK's Butcher table
                                LUDecomposition&    Newton,  // Newton iteration matrix
                                Real                tol      // Newton method's tolerance
                             )
{

    MatReal NewtonItMat ( n_evolved, n_evolved, Real(0) );

    for( Int stage_i = 0; stage_i < BT.s; ++stage_i ) // for each stage
    {
        // If the Jacobian has to be updated at every stage
        if( updateJ == STAGE )
        {
            // update it
            computeNewtonIterationMatrix( stage_i, BT, NewtonItMat );
            // and LU decompose it
            Newton.LUDecompose( NewtonItMat );
        }

        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n ) // for each grid point
        {
            // set the initial guess for this stage
            Int field_i = 0;
            for( auto e: evolvedGF ) // for each field
            {
                if( stage_i == 0 ) // if it is the first iteration
                {
                    // Euler method using the previous time step
                    /// TODO: pass m_previous
                    GF( e.f, m, n ) = GF( e.f, m - 1, n )
                                                + delta_t * GF( e.f_t, m - 1, n );
                } else // otherwise
                {
                    // Euler method using the previous iteration
                    // (updated by MoL_DIRK_computeStage)
                    GF( e.f, m, n ) = GF( e.f, m, n )
                                      + BT.c[stage_i] * delta_t * GF( e.f_t, m, n );
                }
                ++field_i;
            }
        }

        // Compute the stage value
        DIRK_computeStage( m, stage_i, BT, Newton, tol );

    }
    // At this point, the final value of the grid functions at time step m and grid point
    // n are computed

}

/** MoL_DIRK_computeStage computes an individual stage of a DIRK.
 *  It performs the Newton iteration, which stops once the error is smaller than the
 *  user-specified tolerance (in config.ini).
 *  TODO: tolerance and updateJacobian are not read from config.ini yet. Also,
 *        some formulas should be checked.
 */
void MoL::DIRK_computeStage  (
                                Int                 m,       // time step
                                Int                 stage_i, // stage to be computed
                                const ButcherTable& BT,      // the DIRK's Butcher table
                                LUDecomposition&    Newton,  // Newton iteration matrix
                                Real                tol      // Newton method's tolerance
                             )
{
    VecReal dis( n_evolved, Real(10e+3) );
    VecReal res( n_evolved, Real(10e+3) );

    OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n ) // for each grid point
    {
        // Compute the right-hand sides of the evolution equations at this stage
        for( auto eom: eomList ) {
            eom->integStep_CalcEvolutionRHS( m );
        }

        VecReal X( n_evolved, Real(0) ); // A vector storing part of the residual vector

        Int field_i = 0;
        for( auto e: evolvedGF ) // for each field
        // note that evolvedGF must contain the evolved fields in the correct order
        {
            // Compute the sum in the residual [1, p.42, eq.(81)]
            for( Int j = 0; j < stage_i - 1; ++j )
            {
                /// TODO: check this formula
                X[field_i] += delta_t * BT.A[stage_i][j] * GF( e.f_t, m, n );
            }
            ++field_i;
        }

        //////////////////////////////////////////////////////////////////////////////////
        // The quantities above do not depend on the Newton iteration step.
        // Below is the implementation of the Newton iterative method.
        //////////////////////////////////////////////////////////////////////////////////

        /** The iteration works by overwriting the value of GF( e.f, m, n )
         *  at each iteration.
         */
        // while the desired accuracy isn't met yet, and the number of steps is less than
        // a reasonable upper bound
        Int iteration_counter = 0;
        while( MoL::findMinimum( dis ) > tol && iteration_counter < 10e+3 )
        {
            // Compute the residuals for each field
            field_i = 0;
            for( auto e: evolvedGF )
            {
                /// TODO: check this formula
                res[ field_i ] = -( GF( e.f, m, n ) - GF( e.f, m - 1, n ) ) + X[field_i]
                                 + delta_t * BT.A[stage_i - 1][stage_i - 1]
                                    * GF( e.f_t, m, n );
                ++field_i;
            }

            // The Jacobian is not updated at every iteration of the Newton method.
            // Hence, we can solve the system L * U * dis = res for dis.

            // solve the system for the displacements of the fields
            Newton.solve( res, dis );

            // displace the fields (ask Mikica if this is the best way to store the
            // iterations; that is, by overwriting the grid functions)
            field_i = 0;
            for( auto e: evolvedGF )
            {
                GF( e.f, m, n ) += dis[field_i];
                ++field_i;
            }

            ++iteration_counter;

        }
        // after this while loop, the grid functions at ( m, n ), iteration stage_i
        // have updated values
    }
}

/** computeNewtonIterationMatrix computes the Newton iteration matrix.
 */
void MoL::computeNewtonIterationMatrix(
                                        Int stage_i,
                                        const ButcherTable& BT,
                                        MatReal& NewItMat
                                      )
{
    for( Int i = 0; i < n_evolved; ++i )
    {
        for( Int j = 0; j < n_evolved; ++j )
        {
            NewItMat[i][j] = ( i == j ) ? 1 : 0
                             - delta_t * BT.A[stage_i][stage_i];
                                //* Jacobian[i][j]; // this is not defined yet
        }
    }
}

                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _INTEGRATOR_H_INCLUDED
