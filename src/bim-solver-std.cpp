/**
 *  @file      bim-solver-std.cpp
 *  @brief     Standard 3+1 evolution for spherically symmetric bimetric spacetimes.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <chrono>

#define FIXED_EXTRAPOLATION // Extrapolation is not synced with CFDS_ORDER

// #define TWEAK_KDT
// #define TWEAK_MK2
// #define TWEAK_MK3

#include "numMethods.h"        // The implemented numerical modules
#include "sys/slog.h"          // For writing both to cerr and cout simultaneously
#include "sys/trackUsedTime.h" // Keep track of the elapsed time of the application
#include "sys/paramsHolder.h"  // Holds 'key=value' pairs got from the parameter file
#include "sys/hpc.h"           // High-Performance Computing (HPC) support

/////////////////////////////////////////////////////////////////////////////////////////
// Declare our grid functions (they must be in advance known to the grid-driver)

#include "grid/gridFunctions.h"

namespace fld
{
    /////////////////////////////////////////////////////////////////////////////////////
    /// The bimetric grid functions (indices on the grid).
    enum bimIndex { bimFirst = GFCNT - 1,
    /////////////////////////////////////////////////////////////////////////////////////

        gA, gB, gK, gKD, gDA, gDB, gSig,  //!< Evolved variables in the `g`-sector
        fA, fB, fK, fKD, fDA, fDB, fSig,  //!< Evolved variables in the `f`-sector
        pfD, pfS, pftau,                  //!< State variables for the PF

        p,                         //!< Separation between two metrics (relative shift)
        q,                         //!< Overall (geometric mean) shift
        gAlp, gDAlp, fAlp, fDAlp,  //!< Lapses and their log-derivatives

        Lt,                        //!< Lorentz factor (derived from `p`)
        R,                         //!< Ratio `fB/gB`
        pfv,                       //!< The three-velocity of the PF

        /////////////////////////////////////////////////////////////////////////////////
        mpiBoundary = pfv,         //!< The last GF for which BCs must be fixed.
                                   //!< Used in MPI when calling defineGhostChunk.
        /////////////////////////////////////////////////////////////////////////////////

        g_rho, g_JK, g_JKD,        //!< Cumulative sources in th `g`-sector
        f_rho, f_JK, f_JKD,        //!< Cumulative sources in the `f`-sector

        gW, fW, cW,                //!< Lapse ratio, `( gW * gAlp + fW * fAlp ) / cW = 0`
        gW0, gW1, gW2, fW0, fW1, fW2,

        // Diagnostics
        //
        C_1, C_2, C_3, C_4,        //!< Constraints (more precise: constraint violations)
        gHorz, fHorz,              //!< Apparent horizon finders
        p_g, p_f,                  //!< Alt. expr. for the separation between two metrics

        // The RHS of the evolution equations
        //
        gA_t, gB_t, gK_t, gKD_t, gDA_t, gDB_t, gSig_t,
        fA_t, fB_t, fK_t, fKD_t, fDA_t, fDB_t, fSig_t,
        pfD_t, pfS_t, pftau_t, p_t,
        gAlp_t, gDAlp_t,

        // The cached values of the spatial derivatives
        //
        gA_r, gB_r, gK_r, gKD_r, gDA_r, gDB_r,
        fA_r, fB_r, fK_r, fKD_r, fDA_r, fDB_r,
        pfD_r, pfS_r, pftau_r, pfv_r,
        p_r, p_rr, q_r, q_rr, eq_pr_r, eq_qr_r,
        gAlp_r, gDAlp_r, fAlp_r, fAlp_rr, fDAlp_r,

        // Grid functions uses for storing temporary data
        //
        tmp, dbg, gdbg, fdbg,

    /////////////////////////////////////////////////////////////////////////////////////
    bimLast };
    #undef GFCNT // Clear the old number of grid functions then...
    #define GFCNT fld::bimLast // update the number of grid functions on a grid point.
    /////////////////////////////////////////////////////////////////////////////////////

    /** The grid functions that will be read from the initial data.
     */
    static const std::vector<Int> bimInput =
    {
        gA, gB, gK, gKD, gDA, gDB, gSig,
        fA, fB, fK, fKD, fDA, fDB, fSig,
        q, gAlp, fAlp, gDAlp,
        p,
        pfD, pfS, pftau
    };

    /** The grid functions that are involved in time.
     */
    static const std::vector<EvolvedBy> bimEvolvedGF =
    {
        { p,   p_t   },
        { gA,  gA_t  },  { gB,  gB_t  },  { gK,   gK_t   },  { gKD, gKD_t },
        { gDA, gDA_t },  { gDB, gDB_t },  { gSig, gSig_t },
        { fA,  fA_t  },  { fB,  fB_t  },  { fK,   fK_t   },  { fKD, fKD_t },
        { fDA, fDA_t },  { fDB, fDB_t },  { fSig, fSig_t },
        { pfD, pfD_t },  { pfS, pfS_t },  { pftau, pftau_t }
    };

    /** The grid functions which will be written to the output.
     */
    static const std::vector<GF_Descriptor> bimOutput =
    {
        { gA,       "gA",       "A"                               },
        { gB,       "gB",       "B"                               },
        { gK,       "gK",       "K"                               },
        { gKD,      "gKD",      "K_\\Delta"                       },
        { fA,       "fA",       "\\tilde A"                       },
        { fB,       "fB",       "\\tilde B"                       },
        { fK,       "fK",       "\\tilde K"                       },
        { fKD,      "fKD",      "\\tilde K_\\Delta"               },
    //
        { gA_t,     "gA_t",     "\\partial_t A"                   },
        { gB_t,     "gB_t",     "\\partial_t B"                   },
        { gK_t,     "gK_t",     "\\partial_t K"                   },
        { gKD_t,    "gKD_t",    "\\partial_t K_\\Delta"           },
        { fA_t,     "fA_t",     "\\partial_t \\tilde A"           },
        { fB_t,     "fB_t",     "\\partial_t \\tilde B"           },
        { fK_t,     "fK_t",     "\\partial_t \\tilde K"           },
        { fKD_t,    "fKD_t",    "\\partial_t \\tilde K_\\Delta"   },
    //
        { gDA,      "gDA",      "D_A"                             },
        { gDB,      "gDB",      "D_B"                             },
        { gSig,     "gSig",     "\\sigma_g"                       },
        { gW,       "gW",       "W_g"                             },
        { fDA,      "fDA",      "\\tilde D_A"                     },
        { fDB,      "fDB",      "\\tilde D_B"                     },
        { fSig,     "fSig",     "\\sigma_f"                       },
        { fW,       "fW",       "W_f"                             },
    //
        { gDA_t,    "gDA_t",    "\\partial_t D_A"                 },
        { gDB_t,    "gDB_t",    "\\partial_t D_B"                 },
        { gSig_t,   "gSig_t",   "\\partial_t \\sigma_g"           },
        { gAlp,     "gAlp",     "\\alpha"                         },
        { fDA_t,    "fDA_t",    "\\partial_t \\tilde D_A"         },
        { fDB_t,    "fDB_t",    "\\partial_t \\tilde D_B"         },
        { fSig_t,   "fSig_t",   "\\partial_t \\sigma_f"           },
        { fAlp,     "fAlp",     "\\tilde\\alpha"                  },
    //
        { p,        "p",        "p"                               },
        { p_t,      "p_t",      "\\partial_t p"                   },
        { q,        "q",        "q"                               },
        { cW,       "cW",       "W_c"                             },
        { C_1,      "C_1",      "C_1"                             },
        { C_2,      "C_2",      "C_2"                             },
        { C_3,      "C_3",      "C_3"                             },
        { C_4,      "C_4",      "C_4"                             },
    //
        { g_rho,    "g_rho",    "\\rho"                           },
        { g_JK,     "g_JK",     "{JK}"                            },
        { g_JKD,    "g_JKD",    "{JK}_\\Delta"                    },
        { gHorz,    "gHorz",    "z_g"                             },
        { f_rho,    "f_rho",    "\\tilde{\\rho}"                  },
        { f_JK,     "f_JK",     "\\tilde{JK}"                     },
        { f_JKD,    "f_JKD",    "\\tilde{JK}_\\Delta"             },
        { fHorz,    "fHorz",    "z_f"                             },
    //
        { pfD,      "pfD",      "\\hat D"                         },
        { pfS,      "pfS",      "\\hat S"                         },
        { pftau,    "pftau",    "\\hat \\tau"                     },
        { pfv,      "pfv",      "\\hat v"                         },
        { pfD_t,    "pfD_t",    "\\partial_t \\hat D"             },
        { pfS_t,    "pfS_t",    "\\partial_t \\hat S"             },
        { pftau_t,  "pftau_t",  "\\partial_t \\hat\\tau"          },
        { gDAlp,    "gDAlp",    "D_\\alpha"                       },
    //
        { sysStat,  "sysStat",  "\\bullet"                        }
    };

    /** Additional output diagnostics (e.g., with spatial derivatives).
     */
    static const std::vector<GF_Descriptor> bimShowDiagnostics =
    {
        { gW0,      "gW0",      "W_{g,0}"                         },
        { gW1,      "gW1",      "W_{g,1}"                         },
        { gW2,      "gW2",      "W_{g,2}"                         },
        { fW0,      "fW0",      "W_{f,0}"                         },
        { fW1,      "fW1",      "W_{f,1}"                         },
        { fW2,      "fW2",      "W_{f,2}"                         },
    //
        { p_g,      "p_g",      "p_g"                             },
        { gA_r,     "gA_r",     "\\partial_r A"                   },
        { gB_r,     "gB_r",     "\\partial_r B"                   },
        { gK_r,     "gK_r",     "\\partial_r K"                   },
        { gKD_r,    "gKD_r",    "\\partial_r K_\\Delta"           },
        { gDA_r,    "gDA_r",    "\\partial_r D_A"                 },
        { gDB_r,    "gDB_r",    "\\partial_r D_B"                 },
        { gdbg,     "gdbg",     "g\\text{dbg}"                    },
    //
        { p_f,      "p_f",      "p_f"                             },
        { fA_r,     "fA_r",     "\\partial_r \\tilde A"           },
        { fB_r,     "fB_r",     "\\partial_r \\tilde B"           },
        { fK_r,     "fK_r",     "\\partial_r \\tilde K"           },
        { fKD_r,    "fKD_r",    "\\partial_r \\tilde K_\\Delta"   },
        { fDA_r,    "fDA_r",    "\\partial_r \\tilde D_A"         },
        { fDB_r,    "fDB_r",    "\\partial_r \\tilde D_B"         },
        { fdbg,     "fdbg",     "f\\text{dbg}"                    },
    //
        { pfD_r,    "pfD_r",    "\\partial_r \\hat D"             },
        { pfS_r,    "pfS_r",    "\\partial_r \\hat S"             },
        { pftau_r,  "pftau_r",  "\\partial_r \\hat \\tau"         },
        { pfv_r,    "pfv_r",    "\\partial_r \\hat v"             },
    //
        { p_r,      "p_r",      "\\partial_r p"                   },
        { p_rr,     "p_rr",     "\\partial_{rr} p"                },
        { q_r,      "q_r",      "\\partial_r q"                   },
        { q_rr,     "q_rr",     "\\partial_{rr} q"                },
        { eq_pr_r,  "eq_pr_r",  "\\partial_r\\text{\"(p/r)\"}"    },
        { eq_qr_r,  "eq_qr_r",  "\\partial_r\\text{\"(q/r)\"}"    },
    //
        { gAlp_r,   "gAlp_r",   "\\partial_r \\alpha"             },
        { gDAlp_r,  "gDAlp_r",  "\\partial_r D_\\alpha"           },
        { fAlp_r,   "fAlp_r",   "\\partial_r \\tilde\\alpha"      },
        { fAlp_rr,  "fAlp_rr",  "\\partial_{rr} \\tilde\\alpha"   },
        { fDAlp_r,  "fDAlp_r",  "\\partial_r \\tilde D_\\alpha"   }
    };
}

/////////////////////////////////////////////////////////////////////////////////////////
// After defining the grid functions, we can use the grid-driver and the other modules

#include "grid/gridDriver.h"
#include "grid/gridInitialData.h"
#include "grid/gridOutput.h"
#include "grid/integrator.h"

/////////////////////////////////////////////////////////////////////////////////////////

#include "bimetricModel.h"

/** BimetricEvolve encapsulates a 3+1 evolution solver for bimetric spacetimes.
 *
 *  @todo EvolvedBy f, f_t --> add dissipation here?
 *  @todo Implement PIRK?
 */
class BimetricEvolve
    : BimetricModel,    // Based on the bimetric model
      GridUser,         // To access grid functions on the grid
      public IntegFace  // Implement the integration interface
{
    enum Slicing        //!< Known slicing methods
    {
        SLICE_CONSTG  = 0,  // Constant slice in g (f calculated)
        SLICE_CONSTGF = 1,  // Constant slice in both g and f
        SLICE_MS2     = 2,  // Maximal slicing, 2nd order FD
        SLICE_MS4     = 3,  // Maximal slicing, 4th order FD
        SLICE_MS2OPT  = 4,  // Maximal slicing, 2nd order FD, optimized algorithm
        SLICE_BM      = 5   // Bona-Masso slicing
    };

    Int slicing;   //!< Select the slicing: maximal, Bona-Masso, ...
    Int lin2n;     //!< Left grid-zone linear smoothing (default: `nGhost`)
    Int cub2n;     //!< Left grid-zone cubic spline smoothing (default: `5*nGhost/2+6`)
    Real delta_t;  //!< The integration step (obtained from the integrator)
    Int smooth;    //!< Smooth the fields (level of smoothness)

    /////////////////////////////////////////////////////////////////////////////////////
    /** @defgroup g5 Macros to access data in a grid                                   */
    /////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
    /////////////////////////////////////////////////////////////////////////////////////

    emitField(gA)      emitField(gA_t)      emitField(fA)      emitField(fA_t)
    emitField(gB)      emitField(gB_t)      emitField(fB)      emitField(fB_t)
    emitField(gK)      emitField(gK_t)      emitField(fK)      emitField(fK_t)
    emitField(gKD)     emitField(gKD_t)     emitField(fKD)     emitField(fKD_t)
    emitField(gDA)     emitField(gDA_t)     emitField(fDA)     emitField(fDA_t)
    emitField(gDB)     emitField(gDB_t)     emitField(fDB)     emitField(fDB_t)
    emitField(gSig)    emitField(gSig_t)    emitField(fSig)    emitField(fSig_t)

    emitField(p)       emitField(p_g)       emitField(p_f)     emitField(p_t)

    emitField(q)
    emitField(gAlp)    emitField(gAlp_t)
    emitField(gDAlp)   emitField(gDAlp_t)
    emitField(fAlp)    emitField(fDAlp)
    emitField(gW)      emitField(gW0)      emitField(gW1)      emitField(gW2)
    emitField(fW)      emitField(fW0)      emitField(fW1)      emitField(fW2)
    emitField(cW)

    emitField(g_rho)   emitField(f_rho)
    emitField(g_JK)    emitField(f_JK)
    emitField(g_JKD)   emitField(f_JKD)

    emitField(pfD)     emitField(pfD_t)
    emitField(pfS)     emitField(pfS_t)
    emitField(pftau)   emitField(pftau_t)
    emitField(pfv)     emitField(pfv_r)

    emitField(C_1)     emitField(C_2)       emitField(C_3)     emitField(C_4)
    emitField(gHorz)   emitField(fHorz)     emitField(Lt)      emitField(R)
                                                                                   /*@}*/
    /////////////////////////////////////////////////////////////////////////////////////
    /** @defgroup g6 Functions of the prime state variables                            */
    /////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
    emitDerivative_r( gA   )     emitDerivative_r( fA   )
    emitDerivative_r( gB   )     emitDerivative_r( fB   )
    emitDerivative_r( gK   )     emitDerivative_r( fK   )
    emitDerivative_r( gKD  )     emitDerivative_r( fKD  )
    emitDerivative_r( gDA  )     emitDerivative_r( fDA  )
    emitDerivative_r( gDB  )     emitDerivative_r( fDB  )

    emitDerivative_r( p     )
    emitDerivative_r( q     )
    emitDerivative_r( gAlp  )    emitDerivative_r( fAlp  )
    emitDerivative_r( gDAlp )    emitDerivative_r( fDAlp )

    emitDerivative_r( pfD   )
    emitDerivative_r( pfS   )
    emitDerivative_r( pftau )

    // These are the only 2nd spatial derivatives which are needed:

    emitDerivative_rr( p    )   // used in eq_gDA_t and eq_fDA_t
    emitDerivative_rr( q    )   // used in eq_gDA_t and eq_fDA_t
    emitDerivative_rr( fAlp )   // used in fDAlp_r

    // The extrinsic curvature relations

    #if 1
        Real gK1  ( Int m, Int n ) { return ( gK(m,n) + 2 * gKD(m,n) ) / 3; }
        Real gK2  ( Int m, Int n ) { return ( gK(m,n) - gKD(m,n) ) / 3; }
        Real fK1  ( Int m, Int n ) { return ( fK(m,n) + 2 * fKD(m,n) ) / 3; }
        Real fK2  ( Int m, Int n ) { return ( fK(m,n) - fKD(m,n) ) / 3; }
    #else
        Real gK   ( Int m, Int n ) { return gK1(m,n) + 2 * gK2(m,n); }
        Real gK_r ( Int m, Int n ) { return gK1_r(m,n) + 2 * gK2_r(m,n); }
        Real fK   ( Int m, Int n ) { return fK1(m,n) + 2 * fK2(m,n); }
        Real fK_r ( Int m, Int n ) { return fK1_r(m,n) + 2 * fK2_r(m,n); }
    #endif

    // The following macros point to equations and not to the grid functions
    /// @todo gBet, fBet to be grid functions (after calc of gAlp, fAlp)

    #define gBet(m,n)   eq_gBet(m,n)
    #define fBet(m,n)   eq_fBet(m,n)
    #define gBet_r(m,n) eq_gBet_r(m,n)
    #define fBet_r(m,n) eq_fBet_r(m,n)

    /////////////////////////////////////////////////////////////////////////////////////

    #define eq_pr(m,n) (p(m,n)/r(m,n))
    #define eq_qr(m,n) (q(m,n)/r(m,n))

    emitDerivative_r(eq_pr)
    emitDerivative_r(eq_qr)

    /////////////////////////////////////////////////////////////////////////////////////
    // The evolution equations for gAlp and gDAlp
    // Here: Bona-Masso slicing condition with f(alp) = 2/alp

    inline Real eq_BM_gAlp_t( Int m, Int n ) {
        return  -2 * gAlp(m,n) * gK(m,n);
    }

    inline Real eq_BM_gDAlp_t( Int m, Int n ) {
        return  -2 * gK_r(m,n);
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // The time derivatives (these are the evolution equations)
    /////////////////////////////////////////////////////////////////////////////////////

    Real eq_gA_t   ( Int m, Int n );        Real eq_fA_t   ( Int m, Int n );
    Real eq_gB_t   ( Int m, Int n );        Real eq_fB_t   ( Int m, Int n );
    Real eq_gDA_t  ( Int m, Int n );        Real eq_fDA_t  ( Int m, Int n );
    Real eq_gDB_t  ( Int m, Int n );        Real eq_fDB_t  ( Int m, Int n );
    Real eq_gSig_t ( Int m, Int n );        Real eq_fSig_t ( Int m, Int n );

    /////////////////////////////////////////////////////////////////////////////////////
    // The split form of operators for the evolution of gK, gKD, fK, and fKD
    /////////////////////////////////////////////////////////////////////////////////////

    Real eq_base_gK_t( Int m, Int n );
    Real eq_invr_gK_t( Int m, Int n );

    Real eq_gK_t( Int m, Int n ) {
        return eq_base_gK_t(m,n) + eq_invr_gK_t(m,n) / r(m,n);
    }

    Real eq_base_gKD_t( Int m, Int n );
    Real eq_invr_gKD_t( Int m, Int n );

    Real eq_gKD_t( Int m, Int n ) {
        return eq_base_gKD_t(m,n) + eq_invr_gKD_t(m,n) / r(m,n);
    }

    Real eq_base_fK_t( Int m, Int n );
    Real eq_invr_fK_t( Int m, Int n );

    Real eq_fK_t( Int m, Int n ) {
        return eq_base_fK_t(m,n) + eq_invr_fK_t(m,n) / r(m,n);
    }

    Real eq_base_fKD_t( Int m, Int n );
    Real eq_invr_fKD_t( Int m, Int n );

    Real eq_fKD_t( Int m, Int n ) {
        return eq_base_fKD_t(m,n) + eq_invr_fKD_t(m,n) / r(m,n);
    }

    // Here defined for completeness (though never used):
    //
    Real eq_base_gK1_t( Int m, Int n );
    Real eq_invr_gK1_t( Int m, Int n );
    Real eq_base_gK2_t( Int m, Int n );
    Real eq_invr_gK2_t( Int m, Int n );
    Real eq_base_fK1_t( Int m, Int n );
    Real eq_invr_fK1_t( Int m, Int n );
    Real eq_base_fK2_t( Int m, Int n );
    Real eq_invr_fK2_t( Int m, Int n );

    /////////////////////////////////////////////////////////////////////////////////////
    // The p-equations
    /////////////////////////////////////////////////////////////////////////////////////

    Real eq_p_t          ( Int m, Int n );
    Real eq_p1_t         ( Int m, Int n );

    Real eq_base_p_g     ( Int m, Int n );
    Real eq_invr_p_g     ( Int m, Int n );

    Real eq_p_g( Int m, Int n ) {
        return eq_base_p_g(m,n) + eq_invr_p_g(m,n) / r(m,n);
    }

    Real eq_base_p_f     ( Int m, Int n );
    Real eq_invr_p_f     ( Int m, Int n );

    Real eq_p_f( Int m, Int n ) {
        return eq_base_p_f(m,n) + eq_invr_p_f(m,n) / r(m,n);
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // The constraints (here used as the error estimators) and the equation for `p`
    /////////////////////////////////////////////////////////////////////////////////////

    Real eq_gC_rho       ( Int m, Int n );
    Real eq_fC_rho       ( Int m, Int n );
    Real eq_gC_j         ( Int m, Int n );
    Real eq_fC_j         ( Int m, Int n );
    Real eq_C_bimConsLaw ( Int m, Int n );
    Real eq_j_over_p     ( Int m, Int n );

    Real eq_C_1          ( Int m, Int n );
    Real eq_C_2          ( Int m, Int n );
    Real eq_C_3          ( Int m, Int n );
    Real eq_C_4          ( Int m, Int n );

    /////////////////////////////////////////////////////////////////////////////////////
    // The sources
    /////////////////////////////////////////////////////////////////////////////////////

    Real eq_g_rho ( Int m, Int n );       Real eq_f_rho ( Int m, Int n );
    Real eq_g_j   ( Int m, Int n );       Real eq_f_j   ( Int m, Int n );
    Real eq_g_JK1 ( Int m, Int n );       Real eq_f_JK1 ( Int m, Int n );
    Real eq_g_JK2 ( Int m, Int n );       Real eq_f_JK2 ( Int m, Int n );
    Real eq_g_JK  ( Int m, Int n );       Real eq_f_JK  ( Int m, Int n );
    Real eq_g_JKD ( Int m, Int n );       Real eq_f_JKD ( Int m, Int n );

    /////////////////////////////////////////////////////////////////////////////////////
    // The perfect fluid
    /////////////////////////////////////////////////////////////////////////////////////

    Real eq_pf_D      ( Int m, Int n );
    Real eq_pf_S      ( Int m, Int n );
    Real eq_pf_tau    ( Int m, Int n );

    Real eq_pf_v      ( Int m, Int n );
    Real eq_pf_v_r    ( Int m, Int n );

    Real eq_pf_D_t    ( Int m, Int n );
    Real eq_pf_S_t    ( Int m, Int n );
    Real eq_pf_tau_t  ( Int m, Int n );

    Real eq_pf_rho    ( Int m, Int n );
    Real eq_pf_j      ( Int m, Int n );
    Real eq_pf_J1     ( Int m, Int n );
    Real eq_pf_J2     ( Int m, Int n );

    /////////////////////////////////////////////////////////////////////////////////////
    // The gauge grid-functions
    /////////////////////////////////////////////////////////////////////////////////////

    Real eq_q    ( Int m, Int n );
    Real eq_gAlp ( Int m, Int n );
    Real eq_fAlp ( Int m, Int n );

    Real eq_gW   ( Int m, Int n );
    Real eq_gW_0 ( Int m, Int n );
    Real eq_gW_1 ( Int m, Int n );
    Real eq_gW_2 ( Int m, Int n );

    Real eq_fW   ( Int m, Int n );
    Real eq_fW_0 ( Int m, Int n );
    Real eq_fW_1 ( Int m, Int n );
    Real eq_fW_2 ( Int m, Int n );

    Real eq_cW   ( Int m, Int n );
    Real eq_cW_0 ( Int m, Int n );
    Real eq_cW_1 ( Int m, Int n );
    Real eq_cW_2 ( Int m, Int n );

    /////////////////////////////////////////////////////////////////////////////////////
    // The apparent horizon finders
    /////////////////////////////////////////////////////////////////////////////////////

    Real eq_gHorz ( Int m, Int n );
    Real eq_fHorz ( Int m, Int n );

    /////////////////////////////////////////////////////////////////////////////////////
    // The shift equations
    /////////////////////////////////////////////////////////////////////////////////////

    Real eq_gBet    ( Int m, Int n );
    Real eq_gBet_r  ( Int m, Int n );
    Real eq_gBetr_r ( Int m, Int n );
    Real eq_fBet    ( Int m, Int n );
    Real eq_fBet_r  ( Int m, Int n );
    Real eq_fBetr_r ( Int m, Int n );
                                                                                   /*@}*/
    /////////////////////////////////////////////////////////////////////////////////////
    /** @defgroup g7 Maximal slicing                                                   */
    /////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
    void maximalSlice_2( Int m, Int N, Real gAlp_at_N );
    void maximalSlice_4( Int m, Int N, Real gAlp_at_N );
    void maximalSlice_2opt( Int m, Int N, Real gAlp_at_N );
    void maximalSlice_PostSteps( Int m, Int N );
                                                                                   /*@}*/
    /////////////////////////////////////////////////////////////////////////////////////
    /** @defgroup g8 Integration                                                       */
    /////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
    /** Calculate the derived variables R, Lt, and pfv.
     *  @note TINY_Real is added to the denominator of pfv to avoid dividing by zero.
     */
    void calculateDerivedVariables( Int m, Int n )
    {
        R(m,n)   = isGR() ? 1 : fB(m,n) / gB(m,n);
        Lt(m,n)  = sqrt( 1 + p(m,n) * p(m,n) );
        pfv(m,n) = pfS(m,n) == 0 ? 0 : eq_pf_v(m,n);
    }

    /** Determine gAlp (e.g. by maximal slicing). Then find fAlp and also calculate
     *  the derivatives of the lapse related grid-functions.
     */
    void determineGaugeFunctions( Int m );

    /** Sommerfeld outgoing (radiative) wave condition on a variable.
     */
    void applySommerfeldBC
    (
        Int m,                 //!< The time slice
        Int n,                 //!< The spatial position
        Int m_old,             //!< The earlier time
        Int gf,                //!< A grid function index
        Real background = 0.0, //!< Background value; 1.0 for ALPHA
        Int fall_off    = 1    //!< Fall-off power exponent; 2 for VET, BET
    );

    /////////////////////////////////////////////////////////////////////////////////////
    // IntegFace implementation
    /////////////////////////////////////////////////////////////////////////////////////

    /** Calculate the dependent variables from the prime state variables
     *  which are needed for the integration.
     */
    virtual void integStep_Prepare( Int m );

    /** Prepare the right-hand side of the evolution equations.
     */
    virtual void integStep_CalcEvolutionRHS( Int m );

    /** Perform the additional steps after each integration step.
     */
    virtual void integStep_Finalize( Int mNext, Int mPrev );

    /** At each checkpoint calculate diagnostics variables and check for eventual NaNs.
     *  Returns 'false' as the indicator to abort the integration.
     */
    virtual bool integStep_Diagnostics( Int m, Int chkNaNs_nFrom, Int chkNaNs_nTo );

    /** Fix the state variables at the left boundary.
     */
    virtual void applyLeftBoundaryCondition( Int m );

    /** Fix the state variables at the right boundary.
     */
    virtual void applyRightBoundaryCondition( Int m );
                                                                                   /*@}*/
    /////////////////////////////////////////////////////////////////////////////////////
    /** @defgroup g9 Public methods                                                    */
    /////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
public:

    /** Creates and configures the bimetric solver from the given parameters.
     */
    BimetricEvolve( Parameters& params,
        UniformGrid& ug, GridOutputWriter& output, MoL& integ );
};

/////////////////////////////////////////////////////////////////////////////////////////
// The constructor
/////////////////////////////////////////////////////////////////////////////////////////

BimetricEvolve::BimetricEvolve( Parameters& params,
            UniformGrid& ug, GridOutputWriter& output, MoL& integ )
    : BimetricModel( params ), GridUser( ug )
{
    static std::map<std::string,int> knownSlicings =
    {
        { "const",  SLICE_CONSTG },
        { "constg", SLICE_CONSTG },  { "constgf", SLICE_CONSTGF },
        { "MS2OPT", SLICE_MS2OPT },  { "MS2",     SLICE_MS2     },
        { "MS4",    SLICE_MS4    },  { "BM",      SLICE_BM      }
    };
    std::string name = params.get( "slicing.method", slicing, 0, knownSlicings );

    params.get( "slicing.lin2n",  lin2n,  nGhost             );
    params.get( "slicing.cub2n",  cub2n,  5 * nGhost / 2 + 6 );
    params.get( "slicing.smooth", smooth, 0                  );

    slog << "Bimetric Solver:" << std::endl << std::endl
         << "    slicing = " << name << " (#" << slicing << ")"
         << ", lin2n = " << lin2n << ", cub2n = " << cub2n
         << ", smooth = " << smooth
         << std::endl << std::endl;

    if ( mpiSize() > 1 &&
        ( slicing == SLICE_MS2 || slicing == SLICE_MS2OPT || slicing == SLICE_MS4 ) )
    {
        slog << "*** Error: Maximal slicing is not compatible with MPI." << std::endl;
        gridDriver->quit( -1 );
    }

    // Cached from the integrator
    //
    delta_t = integ.dt ();

    // Sign up for the integration
    //
    integ.addToEvolution( this );

    // Add our grid functions to the evolution
    //
    integ.keepConstant( { fld::q } );  // GFs that are kept constant in time
    integ.keepEvolved( fld::bimEvolvedGF ); // GFs that are evolved by the integrator

    if( isGR () || slicing == SLICE_CONSTGF ) {
        integ.keepConstant( { fld::fAlp, fld::fDAlp } );
    }

    if ( slicing == SLICE_CONSTG || slicing == SLICE_CONSTGF ) {
        integ.keepConstant( { fld::gAlp, fld::gDAlp } );
    }
    else if ( slicing == SLICE_BM )
    {
        const static std::vector<fld::EvolvedBy> evolvedGaugeGF = {
            { fld::gAlp,  fld::gAlp_t  },
            { fld::gDAlp, fld::gDAlp_t }
        };
        integ.keepEvolved( evolvedGaugeGF );
    }

    // The list of the grid functions to be written to the output.
    //
    output.gridFunctions( fld::bimOutput );

    bool showDiagnostics = false;
    params.get( "out.diagnostics", showDiagnostics, false );

    if( showDiagnostics ) {
        output.gridFunctions( fld::bimShowDiagnostics );
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
// Implementation of the integration interface (BimetricEvolve::IntegFace methods)
/////////////////////////////////////////////////////////////////////////////////////////

void BimetricEvolve::applyLeftBoundaryCondition( Int m )
{
    for( Int i = 0; i < nGhost; ++i )
    {
        Int n  = nGhost - i - 1;
        Int nR = nGhost + i;

        /// Here we impose the parity conditions at the left boundary.

        t     (m,n) = t      (m,nR);
        r     (m,n) = -r     (m,nR);

        p     (m,n) = -p     (m,nR);

        gA    (m,n) = gA     (m,nR);      fA   (m,n) = fA    (m,nR);
        gB    (m,n) = gB     (m,nR);      fB   (m,n) = fB    (m,nR);
        gK    (m,n) = gK     (m,nR);      fK   (m,n) = fK    (m,nR);
        gKD   (m,n) = gKD    (m,nR);      fKD  (m,n) = fKD   (m,nR);
        gDA   (m,n) = -gDA   (m,nR);      fDA  (m,n) = -fDA  (m,nR);
        gDB   (m,n) = -gDB   (m,nR);      fDB  (m,n) = -fDB  (m,nR);
        gSig  (m,n) = -gSig  (m,nR);      fSig (m,n) = -fSig (m,nR);

        pfD   (m,n) = pfD    (m,nR);
        pfS   (m,n) = -pfS   (m,nR);
        pftau (m,n) = pftau  (m,nR);

        q     (m,n) = -q     (m,nR);
        gAlp  (m,n) = gAlp   (m,nR);
        fAlp  (m,n) = fAlp   (m,nR);
        gDAlp (m,n) = -gDAlp (m,nR);

        Lt    (m,n) = Lt     (m,nR);
        R     (m,n) = R      (m,nR);
        pfv   (m,n) = -pfv   (m,nR);
    }
}

void BimetricEvolve::applyRightBoundaryCondition( Int m )
{
    for( Int n = nGhost + nLen; n < nTotal; ++n )
    {
        t(m,n) = t(m,n-1);
        r(m,n) = r(m,n-1) + delta_r;

        extrapolate_R( fld::p,     m, n );   // Extrapolate [n] from [n-1], [n-2], ...

        extrapolate_R( fld::gA,    m, n );      extrapolate_R( fld::fA,    m, n );
        extrapolate_R( fld::gB,    m, n );      extrapolate_R( fld::fB,    m, n );
        extrapolate_R( fld::gK,    m, n );      extrapolate_R( fld::fK,    m, n );
        extrapolate_R( fld::gKD,   m, n );      extrapolate_R( fld::fKD,   m, n );
        extrapolate_R( fld::gDA,   m, n );      extrapolate_R( fld::fDA,   m, n );
        extrapolate_R( fld::gDB,   m, n );      extrapolate_R( fld::fDB,   m, n );
        extrapolate_R( fld::gSig,  m, n );      extrapolate_R( fld::fSig,  m, n );

        extrapolate_R( fld::pfD,   m, n );
        extrapolate_R( fld::pfS,   m, n );
        extrapolate_R( fld::pftau, m, n );

        extrapolate_R( fld::q,     m, n );
        extrapolate_R( fld::gAlp,  m, n );
        extrapolate_R( fld::gDAlp, m, n );
        extrapolate_R( fld::fAlp,  m, n );

        calculateDerivedVariables( m, n ); // Calculate R, Lt, and pfv.

        extrapolate_R( fld::pfv, m, n ); // Nevertheless, we extrapolate pfv!
    }
}

void BimetricEvolve::applySommerfeldBC
(
    Int m1, Int n, Int m, Int gf,
    const Real background,
    const Int fall_off
    )
{
    const Real wave_speed = 1.0;
    const Real Courant = wave_speed * delta_t / delta_r;

    for( Int n = nGhost + nLen + 1; n < nTotal - 1; ++n )
    {
        const Real r0     = r( m, n-2 );
        const Real r1     = r( m, n-1 );
        const Real r2     = r( m, n   );
        const Real r_last = Courant * r1 + ( 1.0 - Courant ) * r2;

        Real factor = r_last / r2;
        if( fall_off == 2 ) {
            factor *= factor;
        }

        const Real P0   = GF( gf, m, n-2 );
        const Real P1   = GF( gf, m, n-1 );
        const Real P2   = GF( gf, m, n );
        const Real P01  = ( ( r_last - r1 ) * P0  + ( r0 - r_last ) * P1  ) / ( r0 - r1 );
        const Real P12  = ( ( r_last - r2 ) * P1  + ( r1 - r_last ) * P2  ) / ( r1 - r2 );
        const Real P012 = ( ( r_last - r2 ) * P01 + ( r0 - r_last ) * P12 ) / ( r0 - r2 );

        GF( gf, m1, n ) = factor * ( P012 - background ) + background;
    }
}

void BimetricEvolve::integStep_Prepare( Int m )
{
    OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
    {
        calculateDerivedVariables( m, n );
    }
}

void BimetricEvolve::determineGaugeFunctions( Int m )
{
    /////////////////////////////////////////////////////////////////////////////////////
    /// - Maximal slicing (optional) -- this will calculate gAlp and gDAlp

    Int sliceBC = smooth >= 1 ? round(40/delta_r) : nLen/2;
    // sliceBC = round(20/delta_r);
    switch( slicing )
    {
        case SLICE_MS2OPT: maximalSlice_2opt ( m, sliceBC, 1 );  break;
        case SLICE_MS2:    maximalSlice_2    ( m, sliceBC, 1 );  break;
        case SLICE_MS4:    maximalSlice_4    ( m, sliceBC, 1 );  break;
    }

#ifdef TWEAK_KDT
    #define KDT 0*
    if( true )
    {
        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            gW   ( m, n ) = eq_gW   ( m, n );
            fW   ( m, n ) = eq_fW   ( m, n );
            cW   ( m, n ) = eq_cW   ( m, n );
            fAlp ( m, n ) = eq_fAlp ( m, n );

            gDAlp_r( m, n ) = 0;
            fAlp_r ( m, n ) = 0;
            fAlp_rr( m, n ) = 0;
            fDAlp  ( m, n ) = 0;
            fDAlp_r( m, n ) = 0;
        }
    } else
#else
    #define KDT
#endif

    /// @todo Go through the order of evaluation for the gauge grid functions
    /// @todo GF_r and GF_rr requires fixed ghost zones
    ///
    if ( isGR () || slicing == SLICE_CONSTGF )
    {
        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            gAlp_r ( m, n ) = GF_r( gAlp,  m, n );
            gDAlp_r( m, n ) = GF_r( gDAlp, m, n );
            fAlp_r ( m, n ) = GF_r( fAlp,  m, n );
            fAlp_rr( m, n ) = GF_rr( fAlp, m, n );
            fDAlp  ( m, n ) = fAlp_r(m,n) / fAlp(m,n);
            fDAlp_r( m, n ) = ( fAlp_rr(m,n) - fAlp_r(m,n) * fAlp_r(m,n) / fAlp(m,n) )
                              / fAlp(m,n);
        }
    }
    else // if not GR
    {
        // gAlp and gDAlp must the BC fixed at this point
        //
        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n ) {
            gAlp_r ( m, n ) = GF_r( gAlp,  m, n );
            gDAlp_r( m, n ) = GF_r( gDAlp, m, n );
        }

        /// @todo smoothen gAlp_r, gDAlp_r

        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            gW0( m, n ) = eq_gW_0( m, n );
            gW1( m, n ) = eq_gW_1( m, n );
            gW2( m, n ) = eq_gW_2( m, n );

            fW0( m, n ) = eq_fW_0( m, n );
            fW1( m, n ) = eq_fW_1( m, n );
            fW2( m, n ) = eq_fW_2( m, n );

            cW( m, n ) = eq_cW_2( m, n );

            gW( m, n ) = gW0(m,n) + r(m,n) * ( gW1(m,n) + r(m,n) * gW2(m,n) );
            fW( m, n ) = fW0(m,n) + r(m,n) * ( fW1(m,n) + r(m,n) * fW2(m,n) );

            fAlp( m, n ) = eq_fAlp( m, n );

            // Finally, normalize gW and fW (this will be displayed at the output)
            gW( m, n ) /= cW( m, n ) * r(m,n) * r(m,n);
            fW( m, n ) /= cW( m, n ) * r(m,n) * r(m,n);
        }

        if( smooth ) {
            smoothenGF( m, fld::fAlp, fld::tmp, fld::fAlp, 1 );
        }
        else {
            applyBoundaryConditions( m, fld::fAlp, +1 );
        }

        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n ) {
            fAlp_r( m, n ) = GF_r( fAlp, m, n );
        }

        if( smooth ) {
            smoothenGF( m, fld::fAlp_r, fld::tmp, fld::fAlp_r, -1 );
        }
        else {
            applyBoundaryConditions( m, fld::fAlp_r, -1 );
        }

        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n ) {
            fAlp_rr( m, n ) = GF_r( fAlp_r, m, n );
        }

        if( smooth ) {
            smoothenGF( m, fld::fAlp_rr, fld::tmp, fld::fAlp_rr, 1 );
        }
        else {
            applyBoundaryConditions( m, fld::fAlp_rr, 1 );
        }

        OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            fDAlp  ( m, n ) = fAlp_r(m,n) / fAlp(m,n);
            fDAlp_r( m, n ) = ( fAlp_rr(m,n) - fAlp_r(m,n) * fAlp_r(m,n) / fAlp(m,n) )
                               / fAlp(m,n);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
// Equations of Motion (generated in Mathematica)
/////////////////////////////////////////////////////////////////////////////////////////

#include "eom-std.h"

void BimetricEvolve::integStep_CalcEvolutionRHS( Int m )
{
    /////////////////////////////////////////////////////////////////////////////////////
    /// - First, calculate the values of the spatial derivatives

    OMP_parallel_for( Int n = CFDS_ORDER/2; n < nTotal - CFDS_ORDER/2; ++n )
    {
        gA_r   (m,n) = GF_r( gA,    m,n );    fA_r   (m,n) = GF_r( fA,    m,n );
        gB_r   (m,n) = GF_r( gB,    m,n );    fB_r   (m,n) = GF_r( fB,    m,n );
        gK_r   (m,n) = GF_r( gK,    m,n );    fK_r   (m,n) = GF_r( fK,    m,n );
        gKD_r  (m,n) = GF_r( gKD,   m,n );    fKD_r  (m,n) = GF_r( fKD,   m,n );
        gDA_r  (m,n) = GF_r( gDA,   m,n );    fDA_r  (m,n) = GF_r( fDA,   m,n );
        gDB_r  (m,n) = GF_r( gDB,   m,n );    fDB_r  (m,n) = GF_r( fDB,   m,n );

        p_r    (m,n) = GF_r( p,     m,n );    q_r    (m,n) = GF_r( q,     m,n );
        p_rr   (m,n) = GF_rr( p,    m,n );    q_rr   (m,n) = GF_rr( q,    m,n );
        eq_pr_r(m,n) = GF_r( eq_pr, m,n );    eq_qr_r(m,n) = GF_r( eq_qr, m,n );

        pfD_r  (m,n) = GF_r( pfD,   m,n );
        pfS_r  (m,n) = GF_r( pfS,   m,n );
        pftau_r(m,n) = GF_r( pftau, m,n );

        // pfv_r(m,n) should be last as it depends on pfD_r, pfS_r, and pftau_r
        //
        pfv_r(m,n) = pfS(m,n) == 0 ? 0 : eq_pf_v_r( m, n );
    }

    if( smooth >= 2 )
    {
        smoothenGF( m, fld::p_r,      fld::tmp, fld::p_r,      +1 );
        smoothenGF( m, fld::p_rr,     fld::tmp, fld::p_rr,     -1 );
        smoothenGF( m, fld::eq_pr_r,  fld::tmp, fld::eq_pr_r,  +1 );

        smoothenGF( m, fld::gA_r,  fld::tmp, fld::gA_r,  -1 );
        smoothenGF( m, fld::gB_r,  fld::tmp, fld::gB_r,  -1 );
        smoothenGF( m, fld::gK_r,  fld::tmp, fld::gK_r,  -1 );
        smoothenGF( m, fld::gKD_r, fld::tmp, fld::gKD_r, -1 );
        smoothenGF( m, fld::gDA_r, fld::tmp, fld::gDA_r, +1 );
        smoothenGF( m, fld::gDB_r, fld::tmp, fld::gDB_r, +1 );

        smoothenGF( m, fld::fA_r,  fld::tmp, fld::fA_r,  -1 );
        smoothenGF( m, fld::fB_r,  fld::tmp, fld::fB_r,  -1 );
        smoothenGF( m, fld::fK_r,  fld::tmp, fld::fK_r,  -1 );
        smoothenGF( m, fld::fKD_r, fld::tmp, fld::fKD_r, -1 );
        smoothenGF( m, fld::fDA_r, fld::tmp, fld::fDA_r, +1 );
        smoothenGF( m, fld::fDB_r, fld::tmp, fld::fDB_r, +1 );
    }

    /// - Calculate the gauge (e.g., do maximal slicing)
    ///
    determineGaugeFunctions( m );

    /////////////////////////////////////////////////////////////////////////////////////
    /// - Calculate the intermediate variables that do not depend on derivatives.

    OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
    {
        // The sources (dependent on both p and the primary dynamical fields)
        //
        g_rho (m,n) = eq_g_rho (m,n);      f_rho (m,n) = eq_f_rho (m,n);
        g_JK  (m,n) = eq_g_JK  (m,n);      f_JK  (m,n) = eq_f_JK  (m,n);
        g_JKD (m,n) = eq_g_JKD (m,n);      f_JKD (m,n) = eq_f_JKD (m,n);
        #ifdef TWEAK_MK3
        g_JKD (m,n) = f_JKD (m,n) = 0;
        #endif
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// - Calculate the variables that depend on the spatial derivatives

    OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
    {
        // The time evolution right-hand side
        //
        p_t     (m,n) = isGR() ? 0 : KDT eq_p_t(m,n); // + eq_C_bimConsLaw (m,n);

        gA_t    (m,n) = eq_gA_t    (m,n);       fA_t    (m,n) = eq_fA_t    (m,n);
        gB_t    (m,n) = eq_gB_t    (m,n);       fB_t    (m,n) = eq_fB_t    (m,n);
        gK_t    (m,n) = eq_gK_t    (m,n);       fK_t    (m,n) = eq_fK_t    (m,n);
        gKD_t   (m,n) = KDT eq_gKD_t   (m,n);   fKD_t   (m,n) = KDT eq_fKD_t   (m,n);
        gDA_t   (m,n) = KDT eq_gDA_t   (m,n);   fDA_t   (m,n) = KDT eq_fDA_t   (m,n);
        gDB_t   (m,n) = KDT eq_gDB_t   (m,n);   fDB_t   (m,n) = KDT eq_fDB_t   (m,n);
        gSig_t  (m,n) = KDT eq_gSig_t  (m,n);   fSig_t  (m,n) = KDT eq_fSig_t  (m,n);

        pfD_t   (m,n) = eq_pf_D_t  (m,n);
        pfS_t   (m,n) = eq_pf_S_t  (m,n);
        pftau_t (m,n) = eq_pf_tau_t(m,n);

        if( slicing == SLICE_BM )
        {
            gAlp_t  (m,n) = eq_BM_gAlp_t  (m,n);
            gDAlp_t (m,n) = eq_BM_gDAlp_t (m,n);
        }
    }

    if( smooth >= 2 )
    {
        // smoothenGF2( m, fld::p_t,   fld::tmp, fld::p_t,   -1 );
        // smoothenGF( m, fld::gA_t,  fld::tmp, fld::gA_t,  +1 );
        // smoothenGF( m, fld::gB_t,  fld::tmp, fld::gB_t,  +1 );
        smoothenGF2( m, fld::gK_t,  fld::tmp, fld::gK_t,  +1 );
        smoothenGF2( m, fld::gKD_t, fld::tmp, fld::gKD_t, +1 );
        // smoothenGF2( m, fld::gDA_t, fld::tmp, fld::gDA_t, -1 );
        // smoothenGF2( m, fld::gDB_t, fld::tmp, fld::gDB_t, -1 );
        // smoothenGF( m, fld::fA_t,  fld::tmp, fld::fA_t,  +1 );
        // smoothenGF( m, fld::fB_t,  fld::tmp, fld::fB_t,  +1 );
        smoothenGF2( m, fld::fK_t,  fld::tmp, fld::fK_t,  +1 );
        smoothenGF2( m, fld::fKD_t, fld::tmp, fld::fKD_t, +1 );
        // smoothenGF2( m, fld::fDA_t, fld::tmp, fld::fDA_t, -1 );
        // smoothenGF2( m, fld::fDB_t, fld::tmp, fld::fDB_t, -1 );
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// - Smoothen the time derivatives inside the grid zone near the outer boundary
    /// @todo Smoothen the time derivatives inside the left grid zone

    if ( gridDriver->isLastInRank() )
    {
        for( Int n = nGhost + nLen - 2*CFDS_ORDER; n < nGhost + nLen; ++n )
        {
            for( auto e: fld::bimEvolvedGF ) {
                extrapolate_D( e.f_t, m, n );
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
/// Post-evolution smoothing and calculation of diagnostics
///
void BimetricEvolve::integStep_Finalize( Int m1, Int m )
{
}

/////////////////////////////////////////////////////////////////////////////////////////
/// Find horizons, calculate constraints and similar. Also check for NaNs.
///
bool BimetricEvolve::integStep_Diagnostics( Int m, Int chkNaNs_nFrom, Int chkNaNs_nTo )
{
    OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
    {
        /// - Apparent horizon finder
        ///
        gHorz (m,n) = eq_gHorz (m,n);
        fHorz (m,n) = eq_fHorz (m,n);

        /// - The constraint violations
        ///
        C_1 (m,n) = eq_C_1 (m,n);
        C_2 (m,n) = eq_C_2 (m,n);
        C_3 (m,n) = eq_C_3 (m,n);
        C_4 (m,n) = eq_C_bimConsLaw (m,n);

        /// - Calculate p in an alternative way (from the constraints)
        ///
        p_g (m,n) = eq_p_g (m,n);
        p_f (m,n) = eq_p_f (m,n);
    }

    if ( chkNaNs_nTo < 0 ) { // CheckNaNs is disabled
        return true;
    }

    /// - Check for NaNs in gA in the given zone.
    ///
    Int nFrom = nGhost + std::min( nLen, std::max( Int(0), chkNaNs_nFrom - nOffset ) );
    Int nTo   = nGhost + std::min( nLen, std::max( Int(0), chkNaNs_nTo   - nOffset ) );

    // if( false )
    for( Int n = nFrom; n < nTo; ++n )
    {
        if( std::isnan( gA( m, n ) ) ) {
            std::cerr << "*** Detected gA NaN at t = " << t(m,n)
                      << ", r = " << r(m,n) << std::endl;
            return false;
        }
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Maximal slicing (implements a BVP solver)
/////////////////////////////////////////////////////////////////////////////////////////

#include "maximalSlice.h"

/////////////////////////////////////////////////////////////////////////////////////////
/// The main entry point of `bim-solver`.
///
int main( int argc, char* argv[] )
{
    TrackUsedTime timer;

    /// - Read the run-time configuration parameters
    ///
    Parameters params( argc >= 2 ? argv[1] : "config.ini" );

    /// - Create the grid driver
    ///
    UniformGrid ugrid( params );

    /// - Read the initial data (populating the grid slice at t_0)
    ///
    if( ! GridInitialData( params, ugrid ).addGFs( fld::bimInput ).load () ) {
        return -1;
    }

    /// - Create the output sink (for storing the results)
    ///
    GridOutputWriter output ( params, ugrid );
    if ( ! output.open () ) {
        return -1;
    }

    /// - Setup the MoL integrator on the grid (whose results will go to the output sink)
    ///
    MoL integrator( params, ugrid, output );

    /// - Create the bimetric model and define the equations of motion (to be integrated)
    ///
    BimetricEvolve bim( params, ugrid, output, integrator );

    /// - Evolve the equations of motion
    ///
    return integrator.evolveEquations () ? 0 : -1;
}
