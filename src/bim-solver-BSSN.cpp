/**
 *  @file      bim-solver-BSSN.cpp
 *  @brief     The covariant BSSN evolution for spherically symmetric bimetric spacetimes.
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

#include "numMethods.h"        // The implemented numerical modules
#include "sys/slog.h"          // For writing both to cerr and cout simultaneously
#include "sys/trackUsedTime.h" // Keep track of the elapsed time of the application
#include "sys/paramsHolder.h"  // Holds 'key=value' pairs got from the parameter file
#include "sys/hpc.h"           // High-Performance Computing (HPC) support
//#include "sys/nr3.h"

#ifndef _TEST_MODE
    #define _TEST_MODE 1
#endif // _TEST_MODE

#ifndef OBSERVER
    #define OBSERVER 1
#endif // OBSERVER

#ifndef _EVOLVE_DSIG
    #define _EVOLVE_DSIG 0
#endif // _EVOLVE_DSIG

#ifndef _DETECT_NAN
    #define _DETECT_NAN 1
#endif // _DETECT_NAN

#if _TEST_MODE == 1
    #define CFDS_ORDER 2
#endif // _TEST_MODE

/////////////////////////////////////////////////////////////////////////////////////////
// Declare our grid functions (they must be in advance known to the grid-driver)

#include "grid/gridFunctions.h"

namespace fld
{
    /////////////////////////////////////////////////////////////////////////////////////
    /// The bimetric grid functions (indices on the grid).
    enum bimIndex { bimFirst = GFCNT - 1,
    /////////////////////////////////////////////////////////////////////////////////////

        #if _TEST_MODE

            u, v, u_t, v_t,

        #else

            /// Evolved fields in the `g`-sector
            gconf, gDconf, gtrK, gA, gDA, gB, gDB, gA1, gA2, gL, gsig, gAsig,

            /// The evolution RHS for `g`
            gconf_t, gDconf_t, gtrK_t, gA_t, gDA_t, gB_t, gDB_t, gA1_t, gA2_t, gL_t,
            gsig_t, gAsig_t,

            /// Evolved fields in the `f`-sector
            fconf, fDconf, ftrK, fA, fDA, fB, fDB, fA1, fA2, fL, fsig, fAsig,

            /// The evolution RHS for `f`
            fconf_t, fDconf_t, ftrK_t, fA_t, fDA_t, fB_t, fDB_t, fA1_t, fA2_t, fL_t,
            fsig_t, fAsig_t,

            /// The traces of the shifted extrinsic curvatures
            gtrA, ftrA,

            /// The violations of the traces of the conformal extrinsic curvatures
            /// (only used when the previous ones are set to 0; they only appear in
            /// the evolution equations for the confomal metrics, see Brown)
            gtrAv, ftrAv,

            /// The Pfaffian derivatives of the traces
            gtrA_pff, ftrA_pff,

            /// The determinants of the conformal metrics
            gdet, fdet,

            /// The Pfaffian derivatives of the determinants
            gdet_pff, fdet_pff,

            /// The regularizing fields and their radial derivatives
            gBetr, gDconfr, fDconfr, gDAlpr, fBetr, fDAlpr, gLr, fLr,
            gBetr_r, gBetr_rr, gDconfr_r, fDconfr_r, gDAlpr_r, fBetr_r,
            fBetr_rr, fDAlpr_r, gLr_r, fLr_r, gBr, gBr_r,

            /// State variables for the perfect fluid (PF)
            pfD, pfS, pftau, pfv,

            /// The Lorentz factor for the PF
            pfW,

            /// The Valencia evolution RHS for the PF
            pfD_t, pfS_t, pftau_t,

            /// The radial derivative of r
            pfv_r,

            /// Sum of sources in the `g`-sector
            gJK, gJA1, gJA2, gJL,

            /// Sum of sources in the `f`-sector
            fJK, fJA1, fJA2, fJL,

            /// Physical sources computed from the Valencia variables
            grhobar,
            grho, frho,
            gj,   fj,
            gJ11, gJ22,
            fJ11, fJ22,

            /// Bimetric sources
            grhob,   gjb_u,   gJb1_ud,  gJb2_ud,
            frhob,   fjb_u,   fJb1_ud,  fJb2_ud,

            /// Separation (rel shift) between the metrics
            p, pr,

            /// Lorentz factor (derived from `p`)
            Lt, Lt2,

            /// RHS for the evolution of `p`
            p_t,

            /// Shift of the mean metric, and auxiliary variable for standard gauge
            q, qr, q_t,
            Bq, Bq_t, Bq_r,

            /// Radial derivatives
            // The radial derivatives of the evolved fields
            gconf_r,     fconf_r,
            gDconf_r,    fDconf_r,
            gtrK_r,      ftrK_r,
            gA_r,        fA_r,
            gB_r,        fB_r,
            gDA_r,       fDA_r,
            gDB_r,       fDB_r,
            gA1_r,       fA1_r,
            gA2_r,       fA2_r,
            gL_r,        fL_r,

            gAsig_r,     fAsig_r,

        // The radial derivatives of the gauge fields
            p_r,         q_r,
            pr_r,        qr_r,
            gAlp_r,      fAlp_r,
            gDAlp_r,     fDAlp_r,
            Lt_r,

        // The radial derivatives of the Valencia variables and the external sources
            pfD_r,
            pfS_r,
            pftau_r,
            gj_r,

        // The radial derivatives of the utility function R = fB/gB
            R_r,

            /// Convective derivatives
            //  The convective derivatives of the evolved fields; they arise from the \
                Pfaffians (i.e., from the Lie derivatives along the shifts)
            gAlp_convr,   fAlp_convr,
            gDAlp_convr,
            //gBet_convr,   fBet_convr,
            q_qconvr,     Bq_qconvr,

            gconf_convr,  fconf_convr,
            gDconf_convr, fDconf_convr,
            gtrK_convr,   ftrK_convr,

            gA_convr,     fA_convr,
            gB_convr,     fB_convr,
            gDA_convr,    fDA_convr,
            gDB_convr,    fDB_convr,

            gA1_convr,    fA1_convr,
            gA2_convr,    fA2_convr,

            gL_convr,     fL_convr,
            gL_qconvr,

            gsig_convr,   fsig_convr,
            gAsig_convr,  fAsig_convr,

            pfD_convr,    pfS_convr,
            pftau_convr,

            gBet_gDconf_r,
            fBet_fDconf_r,

            gBet_gDAlp_r,

            gBet_gDA_r,
            fBet_fDA_r,

            gBet_gDB_r,
            fBet_fDB_r,

            gBet_gL_r,
            fBet_fL_r,

            gsig_gL_r,
            fsig_fL_r,

            gDsig,
            fDsig, ///@fixme   these fields must be declared here because otherwise they \
                            cannot be read from the initial data. The present  \
                            implementation reads these two fields from the initial data \
                            anyways. Then, if they are not evolved, they are not  \
                            touched anymore. This is because otherwise the Mathematica \
                            file should be also changed.

            #if _EVOLVE_DSIG

                    gDsig_t,
                    fDsig_t,
                    gDsig_r,
                    fDsig_r,
                    gDsig_convr,
                    fDsig_convr,

            #else
                    gsig_r,
                    fsig_r,
                    gsig_rr,
                    fsig_rr,

            #endif // _EVOLVE_DSIG

        // These are the only spatial second derivatives that are needed:

        // The second radial derivatives of the evolved fields
            gconf_rr,
            fconf_rr,
            gA_rr,
            gB_rr,
            fA_rr,
            fB_rr,
            gA1_rr,
            fA1_rr,
            gtrK_rr,
            ftrK_rr,

        // The second radial derivatives of the gauge fields
            p_rr,   // used in eq_gDA_t and eq_fDA_t
            q_rr,   // used in eq_gDA_t and eq_fDA_t
            pr_rr,
            qr_rr,
            gAlp_rr,
            fAlp_rr,   // used in fDAlp_r

        // The second radial derivatives of the traces of the conformal extrinsic \
            curvatures (these traces are usually set to 0 in eomBSSNObserver.h)
            gtrA_rr,
            ftrA_rr,

            /// Lapses
            gAlp, fAlp,
            gDAlp, fDAlp,

            /// RHS for the evolution of the lapse
            gAlp_t, gDAlp_t,

            /// Lapse ratio (`gW * gAlp + fW * fAlp = 0`)
            gW, fW,

            /// Shifts
            gBet, fBet,

            /// Radial derivatives of the shifts
            gBet_r, fBet_r,
            gBet_rr, fBet_rr,

            /// The terms involving the Ricci scalars in the evolution equations
            gRicci, fRicci,

            /// The recurrent terms involving some radial derivatives
            gDers, fDers,

            /// Constraints (here error estimators)
            gHC, fHC,
            gMC, fMC,
            gLC, fLC,
            CL,

            /// Apparent horizon finders
            gHorz, fHorz,

            /// Ratio `fB/gB`
            R,

            mpiBoundary = R,       //!< The last grid function for which BC must be fixed.
                                   //!< Used in MPI when calling defineGhostChunk.

        #endif // _TEST_MODE

        mpiBoundary = v,       //!< The last grid function for which BC must be fixed.
                               //!< Used in MPI when calling defineGhostChunk.

        /// The error function for RKDP MoL integrator
        error,

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
        #if _TEST_MODE

            u, v

        #else

            gconf, fconf,
            gDconf, fDconf,
            gtrK, ftrK,
            gA, fA,
            gDA, fDA,
            gB, fB,
            gDB, fDB,
            gA1, fA1,

            #if OBSERVER==1

            #else

                gA2, fA2,

            #endif // OBSERVER

            gL, fL,
            gsig, fsig,
            gDsig, fDsig,
            gAsig, fAsig,
            q, Bq, p, ///@note the initial value for p is not free!
            gAlp, fAlp,
            gDAlp, fDAlp,
            pfD, pfS, pftau

         #endif // _TEST_MODE
    };

    /** The grid functions that are involved in time.
     */
    static const std::vector<EvolvedBy> bimEvolvedGF =
    {
        #if _TEST_MODE

            { u, u_t }, { v, v_t }

        #else

            { gconf,  gconf_t  },   { fconf,  fconf_t  },
            { gDconf,  gDconf_t  }, { fDconf,  fDconf_t  },
            { gtrK,  gtrK_t  },     { ftrK,  ftrK_t  },
            { gA,  gA_t  },         { gB,  gB_t  },
            { gDA,  gDA_t  },       { gDB,  gDB_t  },
            { fA,  fA_t  },         { fB,  fB_t  },
            { fDA,  fDA_t  },       { fDB,  fDB_t  },
            { gA1, gA1_t },

            #if OBSERVER==1

            #else

                { gA2, gA2_t },

            #endif // OBSERVER

            { fA1, fA1_t },

            #if OBSERVER==1

            #else

                { fA2, fA2_t },

            #endif // OBSERVER

            { gL, gL_t },           { fL, fL_t },
            { gsig, gsig_t },       { fsig, fsig_t },
            { gAsig, gAsig_t },     { fAsig, fAsig_t },
            { pfD, pfD_t }, { pfS, pfS_t }, { pftau, pftau_t }

        #endif // _TEST_MODE

        //////////////////////////////////////////////////////////////////////////////////

        /*{ gconf,  gconf_t  },   { gDconf,  gDconf_t  },
        { gtrK,  gtrK_t  },
        { gA,  gA_t  },         { gB,  gB_t  },
        { gDA,  gDA_t  },       { gDB,  gDB_t  },
        { gA1, gA1_t },
        { gL, gL_t },
        { gsig, gsig_t },

        { gAsig, gAsig_t },

        { fconf,  fconf_t  },   { fDconf,  fDconf_t  },
        { ftrK,  ftrK_t  },
        { fA,  fA_t  },         { fB,  fB_t  },
        { fDA,  fDA_t  },       { fDB,  fDB_t  },
        { fA1, fA1_t },
        { fL, fL_t },
        { fsig, fsig_t },

        { fAsig, fAsig_t },

        { p,   p_t   },

        { pfD, pfD_t }, { pfS, pfS_t }, { pftau, pftau_t }

        #if OBSERVER==1

        #else

            { gA2, gA2_t },
            { fA2, fA2_t },

        #endif // OBSERVER

        #if _EVOLVE_DSIG

           ,{ gDsig, gDsig_t }
           ,{ fDsig, fDsig_t }

        #endif // _EVOLVE_DSIG*/

    };

    /** The grid functions which will be written to the output.
     */
    static const std::vector<GF_Descriptor> bimOutput =
    {
        #if _TEST_MODE

        { u,         "u",          "u"                              },
        { v,         "v",          "v"                              },
        { u_t,       "u_t",        "\\partial _{t}u"                },
        { v_t,       "v_t",        "\\partial _{t}v"                }

        #else

        { gAlp,      "gAlp",       "\\alpha"                        },
        { gAlp_r,    "gAlp_r",     "\\partial _r \\alpha"           },
        { gDAlp,     "gDAlp",      "D_\\alpha"                      },
        { gAlp_rr,   "gAlp_rr",    "\\partial _{rr} \\alpha"        },
        { gDAlp_r,   "gDAlp_r",    "\\partial _r D_\\alpha"         },
        { gAlp_t,    "gAlp_t",     "\\partial _t \\alpha"           },
        { fAlp,      "fAlp",       "\\tilde \\alpha"                },
        { fDAlp,     "fDAlp",      "\\tilde D_\\alpha"              },
        { fAlp_rr,   "fAlp_rr",    "\\partial _{rr} \\tilde \\alpha"},
        { fDAlp_r,   "fDAlp_r",    "\\partial _r \\tilde D_\\alpha" },
        { gconf,     "gconf",      "\\phi"                          },
        { fconf,     "fconf",      "\\psi"                          },
        { gDconf,    "gDconf",     "D _\\phi"                       },
        { fDconf,    "fDconf",     "D _\\psi"                       },
        { gconf_r,   "gconf_r",   "\\partial _{r} \\phi"            },
        { fconf_r,   "fconf_r",   "\\partial _{r} \\psi"            },
        { gtrK,      "gtrK",       "K"                              },
        { ftrK,      "ftrK",       "\\tilde K"                      },

        { gA,        "gA",         "A"                              },
        { gB,        "gb",         "b"                              },
        { gDA,       "gDA",        "D_A"                            },
        { gDA_r,     "gDA_r",      "\\partial _r D_A"               },
        { gDB_r,     "gDb_r",      "\\partial _r  D_b"              },
        { gDB,       "gDb",        "D_b"                            },
        { gA1,       "gA1",        "A_1"                            },
        { gA2,       "gA2",        "A_2"                            },
        { gL,        "gL",         "\\Lambda"                       },
        { gsig,      "gsig",       "\\sigma"                        },
        { gAsig,     "gAsig",      "A_\\sigma"                      },
        { fA,        "fA",         "\\tilde A"                      },
        { fB,        "fb",         "\\tilde b"                      },
        { fDA,       "fDA",        "\\tilde D_A"                    },
        { fDA_r,     "fDA_r",      "\\partial _r \\tilde D_A"       },
        { fDB,       "fDb",        "\\tilde D_b"                    },
        { fDB_r,     "fDb_r",      "\\partial _r \\tilde D_b"       },
        { fA1,       "fA1",        "\\tilde A_1"                    },
        { fA2,       "fA2",        "\\tilde A_2"                    },
        { fL,        "fL",         "\\tilde \\Lambda"               },
        { fsig,      "fsig",       "\\tilde \\sigma"                },
        { gDconf_r,  "gDconf_r",   "\\partial_r D _\\phi"           },
        { fDconf_r,  "fDconf_r",   "\\partial_r  D_\\psi"           },
        { gconf_rr,  "gconf_rr",   "\\partial _{rr} \\phi"          },
        { fconf_rr,  "fconf_rr",   "\\partial _{rr} \\psi"          },

        { fAsig,     "fAsig",      "\\tilde A_\\sigma"              },
        { gBetr,     "gBetr",      "N r^{-1}"                       },
        { fBetr,     "fBetr",      "M r^{-1}"                       },
        { gDconfr,   "gDconfr",    "\\phi D_\\phi r^{-1}"           },
        { fDconfr,   "fDconfr",    "\\psi D_\\psi r^{-1}"           },
        { gDAlpr,    "gAlpr",      "\\alpha D_\\alpha r^{-1}"       },
        { fDAlpr,    "fDAlpr",     "\\tilde \\alpha \\tilde D_\\alpha r^{-1}"           },
        { gLr,       "gLr",        "\\Lambda r^{-1}"                },
        { fLr,       "fLr",        "\\tilde \\Lambda r^{-1}"        },
        { gBr,       "gbr",        "b\\cdot r"                      },
    //
        { gconf_t,   "gconf_t",    "\\partial_t \\phi"              },
        { fconf_t,   "fconf_t",    "\\partial_t \\psi"              },
        { gDconf_t,  "gDconf_t",   "\\partial_t D _\\phi"           },
        { fDconf_t,  "fDconf_t",   "\\partial_t  D_\\psi"           },
        { gtrK_t,    "gtrK_t",     "\\partial_t K"                  },
        { ftrK_t,    "ftrK_t",     "\\partial_t \\tilde K"          },
        { gA_t,      "gA_t",       "\\partial_t A"                  },
        { gB_t,      "gb_t",       "\\partial_t b"                  },
        { gDA_t,     "gDA_t",      "\\partial_t D_A"                },
        { gDB_t,     "gDb_t",      "\\partial_t D_b"                },
        { gA1_t,     "gA1_t",      "\\partial_t A_1"                },
        { gA2_t,     "gA2_t",      "\\partial_t A_2"                },
        { gL_t,      "gL_t",       "\\partial_t \\Lambda"           },
        { gsig_t,    "gsig_t",     "\\partial_t \\sigma"            },

        { gAsig_t,   "gAsig_t",    "\\partial_t A_\\sigma"          },
        { fA_t,      "fA_t",       "\\partial_t \\tilde A"          },
        { fB_t,      "fb_t",       "\\partial_t \\tilde b"          },
        { fDA_t,     "fDA_t",      "\\partial_t \\tilde D_A"        },
        { fDB_t,     "fDb_t",      "\\partial_t \\tilde D_b"        },
        { fA1_t,     "fA1_t",      "\\partial_t \\tilde A_1"        },
        { fA2_t,     "fA2_t",      "\\partial_t \\tilde A_2"        },
        { fL_t,      "fL_t",       "\\partial_t \\tilde \\Lambda"   },
        { fsig_t,    "fsig_t",     "\\partial_t \\tilde \\sigma"    },
        { fAsig_t,   "fAsig_t",    "\\partial_t \\tilde A_\\sigma"  },
    //
        { gdet,      "gdet",       "\\Delta"                        },
        { fdet,      "fdet",       "\\tilde\\Delta"                 },
        { gtrA,      "gtrA",       "\\text{tr}A"                    },
        { ftrA,      "ftrA",       "\\text{tr}\\tilde A"            },
        { gtrAv,     "gtrAv",      "\\text{tr}A_v"                  },
        { ftrAv,     "ftrAv",      "\\text{tr}\\tilde A_v"          },
        { gdet_pff,  "gdet_pff",   "\\partial_0 \\Delta"            },
        { fdet_pff,  "fdet_pff",   "\\partial_0 \\tilde\\Delta"     },
        { gtrA_pff,  "gtrA_pff",   "\\partial_0 \\text{tr}A"        },
        { ftrA_pff,  "ftrA_pff",   "\\partial_0 \\text{tr}\\tilde A"},
    //
        { gW,        "gW",         "W_g"                            },
        { fW,        "fW",         "W_f"                            },
    //
        { pr,        "pr",         "p r^{-1}"                       },
        { qr,        "qr",         "q r^{-1}"                       },
        { pr_r,      "pr_r",       "\\partial_r \\left(p r^{-1}\\right)"},
        { qr_r,      "qr_r",       "\\partial_r \\left(q r^{-1}\\right)"},
        { gBet,      "gBet",       "N"                              },
        { gBet_r,    "gBet_r",     "\\partial_r N"                  },
        { gBet_rr,   "gBet_rr",    "\\partial_{rr} N"               },
        { fBet,      "fBet",       "M"                              },
        { fBet_r,    "fBet_r",     "\\partial_r M"                  },
        { fBet_rr,   "fBet_rr",    "\\partial_{rr} M"               },

        { gAlp_convr,"gAlp_convr", "N^r \\partial _r \\alpha"       },
        { fAlp_convr,"fAlp_convr", "M^r \\partial _r \\tilde \\alpha" },

        { gconf_convr, "gconf_convr", "N^r \\partial _r \\phi"      },
        { fconf_convr, "fconf_convr", "M^r \\partial _r \\psi"      },
        { gDconf_convr,"gDconf_convr","N^r \\partial _r D_\\phi"    },
        { fDconf_convr,"fDconf_convr","M^r \\partial _r D_\\psi"    },
        { gtrK_convr,  "gtrK_convr",  "N^r \\partial _r K"          },
        { ftrK_convr,  "ftrK_convr",  "M^r \\partial _r \\tilde K"  },

        { gA_convr,    "gA_convr", "N^r \\partial _r A"             },
        { fA_convr,    "fA_convr", "M^r \\partial _r \\tilde A"     },
        { gB_convr,    "gB_convr", "N^r \\partial _r B"             },
        { fB_convr,    "fB_convr", "M^r \\partial _r \\tilde B"     },
        { gDA_convr,   "gDA_convr", "N^r \\partial _r D_A"          },
        { fDA_convr,   "fDA_convr", "M^r \\partial _r \\tilde D_A"  },
        { gDB_convr,   "gDB_convr", "N^r \\partial _r D_B"          },
        { fDB_convr,   "fDB_convr", "M^r \\partial _r \\tilde D_B"  },

        { gA1_convr,   "gA1_convr", "N^r \\partial _r A_1"          },
        { fA1_convr,   "fA1_convr", "M^r \\partial _r \\tilde A_1"  },
        { gA2_convr,   "gA2_convr", "N^r \\partial _r A_2"          },
        { fA2_convr,   "fA2_convr", "M^r \\partial _r \\tilde A_2"  },

        { gL_convr,    "gL_convr", "N^r \\partial _r \\Lambda"      },
        { fL_convr,    "fL_convr", "M^r \\partial _r \\tilde \\Lambda"},

        { gsig_convr,  "gsig_convr","N^r \\partial _r \\sigma"      },
        { fsig_convr,  "fsig_convr","M^r \\partial _r \\tilde \\sigma"},
        { gAsig_convr, "gAsig_convr","N^r \\partial _r A_\\sigma"      },
        { fAsig_convr, "fAsig_convr","M^r \\partial _r \\tilde A_\\sigma"},

        { pfD_convr,   "pfD_convr", "N^r \\partial _r \\hat D"      },
        { pfS_convr,   "pfS_convr", "N^r \\partial _r \\hat S"},
        { pftau_convr, "pftau_convr","N^r \\partial _r \\hat \\tau"    },

        #if _EVOLVE_DSIG

        { gDsig,     "gDsig",      "D_\\sigma"                      },
        { fDsig,     "fDsig",      "\\tilde D_\\sigma"              },
        { gDsig_convr,"gDsig_convr","N^r \\partial _r D_\\sigma"      },
        { fDsig_convr,"fDsig_convr","M^r \\partial _r \\tilde D_\\sigma"},
        { gDsig_t,   "gDsig_t",    "\\partial_t D_ \\sigma"         },
        { fDsig_t,   "fDsig_t",    "\\partial_t \\tilde D_ \\sigma" },

        #else

        { gsig_r,    "gsig_r",     "\\partial_r \\sigma"            },
        { fsig_r,    "fsig_r",     "\\partial_r \\tilde \\sigma"    },

        #endif // _EVOLVE_DSIG

        /*if( slicing == SG )
        {

            { q_qconvr,  "q_qconvr",   "q^r \\partial _r q"         },
            { Bq_qconvr, "Bq_qconvr",  "M^r \\partial _r Bq"        },
            { gL_qconvr,   "gL_qconvr","q^r \\partial _r \\Lambda"      },

        }*/
    //
        { gBetr_r,   "gBetr_r",    "\\partial_r \\left(N r^{-1}\\right)"                },
        { fBetr_r,   "fBetr_r",    "\\partial_r \\left(M r^{-1}\\right)"                },
        { gDconfr_r, "gDconfr_r",  "\\partial_r \\left(\\phi D_\\phi r^{-1}\\right)"    },
        { fDconfr_r, "fDconfr_r",  "\\partial_r \\left(\\psi D_\\psi r^{-1}\\right)"    },
        { gDAlpr_r,  "gAlpr_r",    "\\partial_r \\left(\\alpha D_\\alpha r^{-1}\\right)"},
        { fDAlpr_r,  "fDAlpr_r",
                "\\partial_r \\left(\\tilde \\alpha \\tilde D_\\alpha r^{-1}\\right)"   },
        { gLr_r,     "gLr_r",      "\\partial_r \\left(\\Lambda r^{-1}\\right)"         },
        { fLr_r,     "fLr_r",      "\\partial_r \\left(\\tilde \\Lambda r^{-1}\\right)" },
        { gBr_r,     "gbr_r",      "\\partial_r \\left(b\\cdot r\\right)"               },
    //
        { p,         "p",          "p"                              },
        { p_t,       "p_t",        "\\partial_t p"                  },
        { Lt,        "Lt",         "\\lambda"                       },
        { Lt2,       "Lt2",        "\\lambda ^2"                    },
        { q,         "q",          "q"                              },
        { p_r,       "p_r",        "\\partial_r p"                  },
        { p_rr,      "p_rr",       "\\partial_{rr} p"               },
        { q_r,       "q_r",        "\\partial_r q"                  },
        { q_rr,      "q_rr",       "\\partial_{rr} q"               },
    //
        { gHC,       "gHC",        "C^g_H"                          },
        { fHC,       "fHC",        "C^f_H"                          },
        { gMC,       "gMC",        "C^g_M"                          },
        { fMC,       "fMC",        "C^f_M"                          },
        { gLC,       "gLC",        "C^g_\\Lambda"                   },
        { fLC,       "fLC",        "C^f_\\Lambda"                   },
        { CL,        "CL",         "CL"                             },
    //
        { grhobar,   "grhobar",    "\\rho _\\mathrm{bar}"           },
        { grho,      "grho",       "\\rho"                          },
        { gj,        "gj",         "j"                              },
        { gJ11,      "gJ11",       "J^1_1"                          },
        { gJ22,      "gJ22",       "J^2_2"                          },
        { frho,      "frho",       "\\tilde \\rho"                  },
        { fj,        "fj",         "\\tilde j"                      },
        { fJ11,      "fJ11",       "\\tilde J^1_1"                  },
        { fJ22,      "fJ22",       "\\tilde J^2_2"                  },
    //
        { gJA1,      "gJA1",       "J _{A_1}"                       },
        { gJA2,      "gJA2",       "J _{A_2}"                       },
        { fJA1,      "fJA1",       "\\tilde J _{A_1}"               },
        { fJA2,      "fJA2",       "\\tilde J _{A_2}"               },
    //
        { grhob,     "grhob",      "\\rho^b"                        },
        { gjb_u,     "gjb_u",      "j^b"                            },
        { gJb1_ud,   "gJb1_ud",    "{J_b}^1_1"                      },
        { gJb2_ud,   "gJb2_ud",    "{J_b}^2_2"                      },

        { frhob,     "frhob",      "\\tilde \\rho^b"                },
        { fjb_u,     "fjb_u",      "\\tilde j^b"                    },
        { fJb1_ud,   "fJb1_ud",    "{\\tilde J_b}^1_1"              },
        { fJb2_ud,   "fJb2_ud",    "{\\tilde J_b}^2_2"              },
    //
        { gHorz,     "gHorz",      "z_g"                            },
        { fHorz,     "fHorz",      "z_f"                            },
    //
        { pfD,      "pfD",         "\\hat D"                        },
        { pfS,      "pfS",         "\\hat S"                        },
        { pftau,    "pftau",       "\\hat \\tau"                    },
        { pfv,      "pfv",         "\\hat v"                        },
        { pfv_r,    "pfv_r",       "\\partial_r \\hat v"            },
        { pfD_t,    "pfD_t",       "\\partial_t \\hat D"            },
        { pfS_t,    "pfS_t",       "\\partial_t \\hat S"            },
        { pftau_t,  "pftau_t",     "\\partial_t \\hat\\tau"         },
        { pfW,      "pfW",         "\\hat W"                        },
        { dbg,      "dbg",         "dbg"                            },
        { gdbg,     "gdbg",        "gdbg"                           }

        #endif // _TEST_MODE
    };

    /** Additional output that contains diagnostics (e.g. the spatial derivatives).
     */
    static const std::vector<GF_Descriptor> bimShowDiagnostics =
    {
        /* add GFs here that are optionally written to output */
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

/** BimetricEvolve encapsulates a BSSN evolution solver for bimetric spacetimes.
 */
class BimetricEvolve
    : BimetricModel,    // Based on the bimetric model
      GridUser,         // To access grid functions on the grid
      public IntegFace  // Implement the integration interface
{
    enum Slicing        //!< Known slicing methods
    {
        SLICE_CONSTG    = 0,  // Constant slice in g (f calculated)
        SLICE_CONSTGF   = 1,  // Constant slice in both g and f
        SLICE_MS2       = 2,  // Maximal slicing, 2nd order FD
        SLICE_MS4       = 3,  // Maximal slicing, 4th order FD
        SLICE_MS2OPT    = 4,  // Maximal slicing, 2nd order FD, optimized algorithm
        SLICE_SG        = 5,  // The standard gauge
        SLICE_MS6       = 6,  // Maximal slicing, 6th order FD
        SLICE_KD        = 7,  // The K driver parabolic gauge on the lapse
                              // (relaxing to maximal slicing) [B&S,p.111]
        SLICE_MS4D      = 8,  // The drived maximal slicing [B&S,p.111]
        SLICE_MS4G      = 9,  // Maximal slicing + Gamma-driver
    };

    Int slicing;   //!< Select the slicing: maximal, Bona-Masso, ...
    Int lin2n;     //!< Left grid-zone linear smoothing (default: `nGhost`)
    Real eta;      //!< Dissipation coefficient in the Gamma driver gauge condition.
    Real Kdiff;    //!< Effective diffusion constant for the K-driver parabolic gauge \
                        on the lapse.
    Real Kelas;    //!< 'Elastic' constant for the K-driver parabolic gauge on the lapse.
    Int cub2n;     //!< Left grid-zone cubic spline smoothing (default: `5*nGhost/2+6`)
    Real delta_t;  //!< The integration step (obtained from the integrator)
    Int smooth;    //!< Smooth the fields (level of smoothness)
    Int nSmoothFrom;    //!< Smooth from this point
    Int nSmoothUpTo;    //!< Smooth up to this point
    Int mSmoothUpTo;    //!< Smooth up to this point
    Int sgRadius;       //!< The Savitzky-Golay smoothing radius

    /////////////////////////////////////////////////////////////////////////////////////
    /** @defgroup g5 Macros to access data in a grid                                   */
    /////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
    /////////////////////////////////////////////////////////////////////////////////////

    #if _TEST_MODE

        emitField(u)    emitField(u_t)
        emitField(v)    emitField(v_t)

    #else

        // The evolved fields
        emitField(gconf)   emitField(gconf_t)
        emitField(fconf)   emitField(fconf_t)
        emitField(gDconf)  emitField(gDconf_t)
        emitField(fDconf)  emitField(fDconf_t)
        emitField(gtrK)    emitField(gtrK_t)
        emitField(ftrK)    emitField(ftrK_t)
        emitField(gA)      emitField(gA_t)
        emitField(fA)      emitField(fA_t)
        emitField(gDA)     emitField(gDA_t)
        emitField(fDA)     emitField(fDA_t)
        emitField(gB)      emitField(gB_t)
        emitField(fB)      emitField(fB_t)
        emitField(gDB)     emitField(gDB_t)
        emitField(fDB)     emitField(fDB_t)
        emitField(gA1)     emitField(gA1_t)
        emitField(fA1)     emitField(fA1_t)
        emitField(gA2)     emitField(gA2_t)
        emitField(fA2)     emitField(fA2_t)
        emitField(gL)      emitField(gL_t)
        emitField(fL)      emitField(fL_t)
        emitField(gsig)    emitField(gsig_t)

        emitField(fsig)    emitField(fsig_t)
        emitField(gAsig)   emitField(gAsig_t)
        emitField(fAsig)   emitField(fAsig_t)

        #if _EVOLVE_DSIG

            emitField(gDsig)   emitField(gDsig_t)
            emitField(fDsig)   emitField(fDsig_t)

        #endif // _EVOLVE_DSIG

        // The gauge fields
        emitField(p)       emitField(p_t)
        emitField(gAlp)    emitField(gAlp_t)
        emitField(gDAlp)   emitField(gDAlp_t)
        emitField(fAlp)
        emitField(fDAlp)
        emitField(gBet)    emitField(fBet)
        emitField(q)       emitField(q_t)
        emitField(Bq)
        emitField(Bq_t)
        emitField(gW)      emitField(fW)
        emitField(Lt)      emitField(Lt2)
        emitField(pr)      emitField(qr)

        // The regularizing functions
        emitField(gBetr)   emitField(fBetr)
        emitField(gDconfr) emitField(fDconfr)
        emitField(gDAlpr)  emitField(fDAlpr)
        emitField(gLr)     emitField(fLr)
        emitField(gBr)

        // The sources (bimetric+external) in the evolution equations
        emitField(gJK)     emitField(fJK)
        emitField(gJA1)    emitField(fJA1)
        emitField(gJA2)    emitField(fJA2)
        emitField(gJL)     emitField(fJL)

        // The Valencia variables
        emitField(pfD)     emitField(pfD_t)
        emitField(pfS)     emitField(pfS_t)
        emitField(pftau)   emitField(pftau_t)
        emitField(pfv)     emitField(pfW)
        emitField(pfv_r)

        // The external sources computed from the Valencia variables \
            (they are inside the sources defined above)
        emitField(grhobar)
        emitField(grho)    emitField(frho)
        emitField(gj)      emitField(fj)
        emitField(gJ11)    emitField(fJ11)
        emitField(gJ22)    emitField(fJ22)

        // The bimetric sources
        emitField(grhob)   emitField(frhob)
        emitField(gjb_u)   emitField(fjb_u)
        emitField(gJb1_ud) emitField(fJb1_ud)
        emitField(gJb2_ud) emitField(fJb2_ud)

        //The terms involving the Ricci tensor in the evolution equations
        emitField(gRicci)  emitField(fRicci)

        //The recurrent terms involving some radial derivatives
        emitField(gDers)  emitField(fDers)

        // The constraint violations
        emitField(gHC)     emitField(fHC)
        emitField(gMC)     emitField(fMC)
        emitField(gLC)     emitField(fLC)
        emitField(CL)

        // The apparent horizon trackers
        emitField(gHorz)   emitField(fHorz)

        // The utility function R = fB/gB
        emitField(R)

        emitField(dbg)

        /* Declarations of the determinants and traces
            (the definitions are in eomBSSNObserver.h)
         */
        Real gdet       (Int, Int);
        Real fdet       (Int, Int);
        Real gdet_pff   (Int, Int);
        Real gdet_pff_r (Int, Int);
        Real gdet_r     (Int, Int);
        Real fdet_pff   (Int, Int);
        Real fdet_pff_r (Int, Int);
        Real fdet_r     (Int, Int);
        Real gtrA       (Int, Int);
        Real gtrA_pff   (Int, Int);
        Real gtrA_r     (Int, Int);
        Real gtrAv      (Int, Int);
        Real ftrA       (Int, Int);
        Real ftrA_pff   (Int, Int);
        Real ftrA_r     (Int, Int);
        Real ftrAv      (Int, Int);

    #ifdef DEBUG_GF
        emitField(dbg1)    emitField(dbg2)      emitField(dbg3)    emitField(dbg4)
        emitField(dbg5)    emitField(dbg6)      emitField(dbg7)    emitField(dbg8)
    #endif
                                                                                    /*@}*/
        //////////////////////////////////////////////////////////////////////////////////
        /** @defgroup g6 Functions of the prime state variables                         */
        //////////////////////////////////////////////////////////////////////////////////
                                                                                    /*@{*/

        // The radial derivatives of the evolved fields
        emitDerivative_r( gconf )     emitDerivative_r( fconf )
        emitDerivative_r( gDconf )    emitDerivative_r( fDconf )
        emitDerivative_r( gtrK )      emitDerivative_r( ftrK )
        emitDerivative_r( gA   )      emitDerivative_r( fA   )
        emitDerivative_r( gB   )      emitDerivative_r( fB   )
        emitDerivative_r( gDA  )      emitDerivative_r( fDA  )
        emitDerivative_r( gDB  )      emitDerivative_r( fDB  )
        emitDerivative_r( gA1  )      emitDerivative_r( fA1  )
        emitDerivative_r( gA2  )      emitDerivative_r( fA2  )
        emitDerivative_r( gL  )       emitDerivative_r( fL   )
        emitDerivative_r( gAsig )     emitDerivative_r( fAsig )

        // The radial derivatives of the gauge fields
        emitDerivative_r( p    )      emitDerivative_r( q    )
        emitDerivative_r( pr )        emitDerivative_r( qr )
        emitDerivative_r( Bq    )
        emitDerivative_r( gAlp )      emitDerivative_r( fAlp )
        emitDerivative_r( gDAlp )     emitDerivative_r( fDAlp )
        emitDerivative_r( gBet )      emitDerivative_r( fBet )
        emitDerivative_r( Lt )

        // The radial derivatives of the regularizing functions
        emitDerivative_r( gBetr )     emitDerivative_r( fBetr )
        emitDerivative_r( gDconfr )   emitDerivative_r( fDconfr )
        emitDerivative_r( gDAlpr )    emitDerivative_r( fDAlpr )
        emitDerivative_r( gLr )       emitDerivative_r( fLr )
        emitDerivative_r( gBr )

        // The radial derivatives of the Valencia variables and the external sources
        emitDerivative_r( pfD )
        emitDerivative_r( pfS )
        emitDerivative_r( pftau )
        emitDerivative_r( gj   )

        // The radial derivatives of the utility function R = fB/gB
        emitDerivative_r( R    )

        // These are the only spatial second derivatives that are needed:

        // The second radial derivatives of the evolved fields
        emitDerivative_rr( gconf )
        emitDerivative_rr( fconf )
        emitDerivative_rr( gA    )
        emitDerivative_rr( gB    )
        emitDerivative_rr( fA    )
        emitDerivative_rr( fB    )
        emitDerivative_rr( gA1   )
        emitDerivative_rr( fA1   )
        emitDerivative_rr( gtrK  )
        emitDerivative_rr( ftrK  )

        // The second radial derivatives of the gauge fields
        emitDerivative_rr( p    )   // used in eq_gDA_t and eq_fDA_t
        emitDerivative_rr( q    )   // used in eq_gDA_t and eq_fDA_t
        emitDerivative_rr( pr )
        emitDerivative_rr( qr )
        emitDerivative_rr( gAlp )
        emitDerivative_rr( fAlp )   // used in fDAlp_r
        emitDerivative_rr( gBet )
        emitDerivative_rr( fBet )

        emitDerivative_rr( gBetr )
        emitDerivative_rr( fBetr )

        // The second radial derivatives of the traces of the conformal extrinsic \
            curvatures (these traces are usually set to 0 in eomBSSNObserver.h)
        emitDerivative_rr( gtrA  )
        emitDerivative_rr( ftrA  )

        #if _EVOLVE_DSIG

             emitDerivative_r( gDsig )
             emitDerivative_r( fDsig )

        #else

            emitDerivative_r( gsig )
            emitDerivative_r( fsig )
            emitDerivative_rr( gsig  )
            emitDerivative_rr( fsig  )

        #endif // _EVOLVE_DSIG

        // The convective derivatives of the evolved fields
        emitField(gAlp_convr)    emitField(fAlp_convr)
        emitField(gDAlp_convr)
        //emitField(gBet_convr)    emitField(fBet_convr)
        emitField(q_qconvr)      emitField(Bq_qconvr)

        emitField(gconf_convr)   emitField(fconf_convr)
        emitField(gDconf_convr)  emitField(fDconf_convr)
        emitField(gtrK_convr)    emitField(ftrK_convr)

        emitField(gA_convr)      emitField(fA_convr)
        emitField(gB_convr)      emitField(fB_convr)
        emitField(gDA_convr)     emitField(fDA_convr)
        emitField(gDB_convr)     emitField(fDB_convr)

        emitField(gA1_convr)     emitField(fA1_convr)
        emitField(gA2_convr)     emitField(fA2_convr)

        emitField(gL_convr)      emitField(fL_convr)
        emitField(gL_qconvr)

        emitField(gsig_convr)    emitField(fsig_convr)
        emitField(gAsig_convr)   emitField(fAsig_convr)

        emitField(gBet_gDconf_r) emitField(fBet_fDconf_r)
        emitField(gBet_gDAlp_r)
        emitField(gBet_gDA_r)    emitField(fBet_fDA_r)

        emitField(gBet_gDB_r)    emitField(fBet_fDB_r)
        emitField(gBet_gL_r)     emitField(fBet_fL_r)

        emitField(gsig_gL_r)     emitField(fsig_fL_r)

        #if _EVOLVE_DSIG

             emitField( gDsig_convr )
             emitField( fDsig_convr )

        #endif // _EVOLVE_DSIG

        emitField(pfD_convr)   emitField(pfS_convr)    emitField(pftau_convr)
        //////////////////////////////////////////////////////////////////////////////////
        /** @defgroup g12 Extrinsic curvatures
            @todo   maybe create a header file which converts everything to the standard \
                    ADM variables?
         */
        //////////////////////////////////////////////////////////////////////////////////
                                                                                    /*@{*/
        inline Real gK1 (Int m, Int n){

            return gA1(m, n) + 1/3 *gtrK(m, n);
        }

        inline Real gK2 (Int m, Int n){

            return gA2(m, n) + 1/3 *gtrK(m, n);
        }

        inline Real fK1 (Int m, Int n){

            return fA1(m, n) + 1/3 *ftrK(m, n);
        }

        inline Real fK2 (Int m, Int n){

            return fA2(m, n) + 1/3 *ftrK(m, n);
        }

    #endif // _TEST_MODE

                                                                                   /*@}*/

    /////////////////////////////////////////////////////////////////////////////////////
    // The time derivatives (these are the evolution equations)
    /////////////////////////////////////////////////////////////////////////////////////

    MatReal Jacobian;

    #if _TEST_MODE

        Real eq_u_t     ( Int m, Int n );
        Real eq_v_t     ( Int m, Int n );

    #else

        Real eq_gconf_t     ( Int m, Int n );
        Real eq_gDconf_t    ( Int m, Int n );
        Real eq_gtrK_t      ( Int m, Int n );
        Real eq_gA_t        ( Int m, Int n );
        Real eq_gB_t        ( Int m, Int n );
        Real eq_gDA_t       ( Int m, Int n );
        Real eq_gDB_t       ( Int m, Int n );
        Real eq_gA1_t       ( Int m, Int n );
        Real eq_gA2_t       ( Int m, Int n );
        Real eq_gL_t        ( Int m, Int n );
        Real eq_gsig_t      ( Int m, Int n );
        Real eq_gAsig_t     ( Int m, Int n );

        Real eq_fconf_t     ( Int m, Int n );
        Real eq_fDconf_t    ( Int m, Int n );
        Real eq_ftrK_t      ( Int m, Int n );
        Real eq_fA_t        ( Int m, Int n );
        Real eq_fB_t        ( Int m, Int n );
        Real eq_fDA_t       ( Int m, Int n );
        Real eq_fDB_t       ( Int m, Int n );
        Real eq_fA1_t       ( Int m, Int n );
        Real eq_fA2_t       ( Int m, Int n );
        Real eq_fL_t        ( Int m, Int n );
        Real eq_fsig_t      ( Int m, Int n );
        Real eq_fAsig_t     ( Int m, Int n );

        #if _EVOLVE_DSIG

            Real eq_gDsig_t     ( Int m, Int n );
            Real eq_fDsig_t     ( Int m, Int n );

        #endif // _EVOLVE_DSIG

        Real eq_gRicci      ( Int m, Int n );
        Real eq_fRicci      ( Int m, Int n );

        Real eq_gDers      ( Int m, Int n );
        Real eq_fDers      ( Int m, Int n );

        Real eq_SG_gAlp_t   ( Int m, Int n );
        Real eq_SG_gDAlp_t  ( Int m, Int n );
        Real eq_SG_gBet_t   ( Int m, Int n );
        Real eq_SG_gBq_t    ( Int m, Int n );

        Real eq_KD_gAlp_t   ( Int m, Int n );
        Real eq_KD_gDAlp_t  ( Int m, Int n );
        /////////////////////////////////////////////////////////////////////////////////////
        // The constraints (here used as the error estimators) and the equation for `p`
        /////////////////////////////////////////////////////////////////////////////////////

        Real eq_p_t          ( Int m, Int n );
        Real eq_p_r          ( Int m, Int n );

        Real eq_gHC          ( Int m, Int n );
        Real eq_fHC          ( Int m, Int n );
        Real eq_gMC          ( Int m, Int n );
        Real eq_fMC          ( Int m, Int n );
        Real eq_gLC          ( Int m, Int n );
        Real eq_fLC          ( Int m, Int n );
        Real eq_CL           ( Int m, Int n );
        /////////////////////////////////////////////////////////////////////////////////////
        // The terms in the maximal slicing
        /////////////////////////////////////////////////////////////////////////////////////

        Real W            ( Int m, Int n );
        Real PP           ( Int m, Int n );
        Real RR           ( Int m, Int n );
        Real QQ           ( Int m, Int n );
        Real gDAlp_at_N   ( Int m, Int N, Real h );
        /////////////////////////////////////////////////////////////////////////////////////
        // The sources
        /////////////////////////////////////////////////////////////////////////////////////

        Real eq_gJK       ( Int m, Int n );
        Real eq_gJA1      ( Int m, Int n );
        Real eq_gJA2      ( Int m, Int n );
        Real eq_gJL       ( Int m, Int n );

        Real eq_fJK       ( Int m, Int n );
        Real eq_fJA1      ( Int m, Int n );
        Real eq_fJA2      ( Int m, Int n );
        Real eq_fJL       ( Int m, Int n );

        Real eq_pf_gv      ( Int m, Int n );
        Real eq_pf_gv_r    ( Int m, Int n );
        Real eq_pf_gW      ( Int m, Int n );

        Real eq_pf_gD_t   ( Int m, Int n );
        Real eq_pf_gS_t   ( Int m, Int n );
        Real eq_pf_gtau_t ( Int m, Int n );

        Real eq_pf_grhobar( Int m, Int n );
        Real eq_pf_grho   ( Int m, Int n );
        Real eq_pf_gj     ( Int m, Int n );
        Real eq_pf_gJ11   ( Int m, Int n );
        Real eq_pf_gJ22   ( Int m, Int n );

        Real eq_pf_frho   ( Int m, Int n );
        Real eq_pf_fj     ( Int m, Int n );
        Real eq_pf_fJ11   ( Int m, Int n );
        Real eq_pf_fJ22   ( Int m, Int n );

        Real eq_grhob     ( Int m, Int n );
        Real eq_gjb_u     ( Int m, Int n );
        Real eq_gjb_d     ( Int m, Int n );
        Real eq_gJb1_dd   ( Int m, Int n );
        Real eq_gJb1_du   ( Int m, Int n );
        Real eq_gJb1_ud   ( Int m, Int n );
        Real eq_gJb1_uu   ( Int m, Int n );
        Real eq_gJb2_dd   ( Int m, Int n );
        Real eq_gJb2_du   ( Int m, Int n );
        Real eq_gJb2_ud   ( Int m, Int n );
        Real eq_gJb2_uu   ( Int m, Int n );

        Real eq_frhob     ( Int m, Int n );
        Real eq_fjb_u     ( Int m, Int n );
        Real eq_fjb_d     ( Int m, Int n );
        Real eq_fJb1_dd   ( Int m, Int n );
        Real eq_fJb1_du   ( Int m, Int n );
        Real eq_fJb1_ud   ( Int m, Int n );
        Real eq_fJb1_uu   ( Int m, Int n );
        Real eq_fJb2_dd   ( Int m, Int n );
        Real eq_fJb2_du   ( Int m, Int n );
        Real eq_fJb2_ud   ( Int m, Int n );
        Real eq_fJb2_uu   ( Int m, Int n );

        #define pfV(m,n) eq_pf_V(m,n)

        /////////////////////////////////////////////////////////////////////////////////////
        // The split form of operators
        /////////////////////////////////////////////////////////////////////////////////////

        /*Real eq_base_gK_t( Int m, Int n );
        Real eq_invr_gK_t( Int m, Int n );
        Real eq_base_gKD_t( Int m, Int n );
        Real eq_invr_gKD_t( Int m, Int n );

        Real eq_gK_t( Int m, Int n )
        {
    #ifdef DEBUG_GF
            return ( dbg1(m,n) = eq_base_gK_t(m,n) )
                 + ( dbg2(m,n) = eq_invr_gK_t(m,n) / r(m,n) );
    #else
            return eq_base_gK_t(m,n) + eq_invr_gK_t(m,n) / r(m,n);
    #endif
        }

        Real eq_gKD_t( Int m, Int n )
        {
    #ifdef DEBUG_GF
            return ( dbg3(m,n) = eq_base_gKD_t(m,n) )
                 + ( dbg4(m,n) = eq_invr_gKD_t(m,n) / r(m,n) );
    #else
            return eq_base_gKD_t(m,n) + eq_invr_gKD_t(m,n) / r(m,n);
    #endif
        }

        Real eq_base_fK_t( Int m, Int n );
        Real eq_invr_fK_t( Int m, Int n );
        Real eq_base_fKD_t( Int m, Int n );
        Real eq_invr_fKD_t( Int m, Int n );
        Real eq_invr_fSig( Int m, Int n );

        Real eq_fK_t( Int m, Int n )
        {
    #ifdef DEBUG_GF
            return ( dbg5(m,n) = eq_base_fK_t(m,n) )
                 + ( dbg6(m,n) = eq_invr_fK_t(m,n) / r(m,n) );
    #else
            return eq_base_fK_t(m,n) + eq_invr_fK_t(m,n) / r(m,n);
    #endif
        }

        Real eq_fKD_t( Int m, Int n )
        {
    #ifdef DEBUG_GF
            return ( dbg7(m,n) = eq_base_fKD_t(m,n) )
                 + ( dbg8(m,n) = eq_invr_fKD_t(m,n) / r(m,n) );
    #else
            return eq_base_fKD_t(m,n) + eq_invr_fKD_t(m,n) / r(m,n);
    #endif
        }*/

        Real eq_base_p_g( Int m, Int n );
        Real eq_invr_p_g( Int m, Int n );
        Real eq_base_p_f( Int m, Int n );
        Real eq_invr_p_f( Int m, Int n );

        /////////////////////////////////////////////////////////////////////////////////////
        // The gauge functions
        /////////////////////////////////////////////////////////////////////////////////////

        Real eq_q    ( Int m, Int n );
        Real eq_gAlp ( Int m, Int n );
        Real eq_gW   ( Int m, Int n );
        Real eq_fW   ( Int m, Int n );
        Real eq_fAlp ( Int m, Int n );

        /////////////////////////////////////////////////////////////////////////////////////
        // The apparent horizon finders
        /////////////////////////////////////////////////////////////////////////////////////

        Real eq_gHorz ( Int m, Int n );
        Real eq_fHorz ( Int m, Int n );

        /////////////////////////////////////////////////////////////////////////////////////
        // The shift equations
        /////////////////////////////////////////////////////////////////////////////////////

        Real eq_gBet         ( Int m, Int n );
        Real eq_fBet         ( Int m, Int n );
        Real eq_gBet_r       ( Int m, Int n );
        Real eq_fBet_r       ( Int m, Int n );
        Real eq_gBet_rr      ( Int m, Int n );
        Real eq_fBet_rr      ( Int m, Int n );
        Real eq_gBetr        ( Int m, Int n );
        Real eq_fBetr        ( Int m, Int n );
        Real eq_gBetr_r      ( Int m, Int n );
        Real eq_fBetr_r      ( Int m, Int n );
        Real eq_gBetr_rr     ( Int m, Int n );
        Real eq_fBetr_rr     ( Int m, Int n );

        //////////////////////////////////////////////////////////////////////////////////
        /** @defgroup g7 Maximal slicing                                                */
        //////////////////////////////////////////////////////////////////////////////////
                                                                                    /*@{*/
        void maximalSlice_2( Int m, Int N, Real gAlp_at_N );
        void maximalSlice_4( Int m, Int N, Real gAlp_at_N );
        void maximalSlice_drived_4( Int m, Int N, Real gAlp_at_N );
        void maximalSlice_4_gDAlp( Int m, Int N, Real gAlp_at_N );
        void maximalSlice_6_gDAlp( Int m, Int N, Real gAlp_at_N );
        void maximalSlice_2opt( Int m, Int N, Real gAlp_at_N );
        void maximalSlice_PostSteps( Int m, Int N );
        void maximalSlice_compute_gDAlp( Int m, Int N );

    #endif // _TEST_MODE
                                                                                   /*@}*/
    /////////////////////////////////////////////////////////////////////////////////////
    /** @defgroup g8 Integration                                                       */
    /////////////////////////////////////////////////////////////////////////////////////

    #if ! _TEST_MODE
                                                                                   /*@{*/
    /** Calculate the derived variables R, Lt, and pfv.
     *  @note TINY_Real is added to the denominator of pfv to avoid dividing by zero.
     */
    void calculateDerivedVariables( Int m, Int n )
    {

        R(m,n)     = isGR() ? 1 : (fB(m,n)*exp(2*fconf(m,n)))/(gB(m,n)*exp(2*gconf(m,n)));

        Lt(m,n)    = isGR() ? 1 : sqrt( 1 + p(m,n) * p(m,n) );

        Lt2(m,n)   = isGR() ? 1 : 1 + p(m,n) * p(m,n);

        pfv(m,n)   = eq_pf_gv(m,n);

        pfW(m,n)   = eq_pf_gW(m,n);

        GF(fld::gdet, m, n)     = gdet(m, n);

        GF(fld::gdet_pff, m, n) = gdet_pff(m, n);

        GF(fld::fdet, m, n)     = fdet(m, n);

        GF(fld::fdet_pff, m, n) = fdet_pff(m, n);

        GF(fld::gtrA, m, n)     = gtrA(m, n);

        GF(fld::gtrA_pff, m, n) = gtrA_pff(m, n);

        GF(fld::gtrAv, m, n)    = gtrAv(m, n);

        GF(fld::ftrA, m, n)     = ftrA(m, n);

        GF(fld::ftrA_pff, m, n) = ftrA_pff(m, n);

        GF(fld::ftrAv, m, n)    = ftrAv(m, n);

        gDconfr(m,n) = gDconf(m,n) / r(m,n); // TODO: check this

        gLr    (m,n )= gL(m,n) / r(m,n);

        gBr    (m,n )= gB(m,n) * r(m,n);

        fDconfr(m,n) = fDconf(m,n) / r(m,n);

        fLr    (m,n) = fL(m,n) / r(m,n);

        pr    (m,n) = p(m,n) / r(m,n);

        qr    (m,n) = q(m,n) / r(m,n);

    }

    #endif //_TEST_MODE

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

    /** computeNewtonIterationMatrix computes the Newton iteration matrix.
     */
    virtual void computeNewtonIterationMatrix(
            Int m, Int n, Int n_evolved,
            Int stage_i,
            const ButcherTable& BT,
            MatReal& NewItMat
        );
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
    : BimetricModel( params ), GridUser( ug ),
      Jacobian( fld::bimEvolvedGF.size(), fld::bimEvolvedGF.size(), Real(0) )
      ///TODO: add evolved gauge variables to bimEvolvedGF
{
    #if _TEST_MODE

    // Cached from the integrator
    //
    delta_t = integ.dt ();

    // Sign up for the integration
    //
    integ.addToEvolution( this );

    integ.keepEvolved( fld::bimEvolvedGF ); // GFs that are evolved by the integrator

    #else

    static std::map<std::string,int> knownSlicings =
    {
        { "const",  SLICE_CONSTG },
        { "constg", SLICE_CONSTG },  { "constgf", SLICE_CONSTGF },
        { "MS2OPT", SLICE_MS2OPT },  { "MS2",     SLICE_MS2     },
        { "MS4",    SLICE_MS4    },  { "SG",      SLICE_SG      },
        { "MS6",    SLICE_MS6    },  { "KD",      SLICE_KD      },
        { "MS4D",   SLICE_MS4D   },  { "MS4G",    SLICE_MS4G    }
    };
    std::string name = params.get( "slicing.method", slicing, 0, knownSlicings );

    params.get( "slicing.lin2n",        lin2n,       nGhost              );
    params.get( "slicing.cub2n",        cub2n,       5 * nGhost / 2 + 6  );
    params.get( "slicing.smooth",       smooth,      0                   );
    params.get( "slicing.dissipGauge",  eta,         0.0                 );
    params.get( "slicing.diffGauge",    Kdiff,         0.6                 );
    params.get( "slicing.elastGauge",   Kelas,         0.1                 );

    params.get( "smoothing.nSmoothFrom",  nSmoothFrom, 5*nGhost   );
    params.get( "smoothing.nSmoothUpTo",  nSmoothUpTo, output.get_nOut()   );
    params.get( "smoothing.mSmoothUpTo",  mSmoothUpTo, 50000   );
    /* get mSize out of the above */
    params.get( "smoothing.sgRadius",     sgRadius,    10.0                );
    //nSmoothUpTo = output.get_nOut();

    #if _EVOLVE_DSIG

    slog << "Equations:" << std::endl << std::endl
         << "    Evolve gDsig, fDsig: " << "yes"
         << std::endl << std::endl;

    #else

    slog << "Equations:" << std::endl << std::endl
         << "    Evolve gDsig, fDsig: " << "no"
         << std::endl << std::endl;

    #endif // _EVOLVE_DSIG

    slog << "Bimetric Solver:" << std::endl << std::endl
         << "    slicing = " << name << " (#" << slicing << ")"
         << ", lin2n = " << lin2n << ", cub2n = " << cub2n
         << ", smooth = " << smooth
         << std::endl << std::endl;

    if (smooth >= 1)
    {
        slog << "Smoothing:" << std::endl << std::endl
             << "    Smoothing from n = " << nSmoothFrom
             << ", to n = " << nSmoothUpTo
             << ", up to m = " << mSmoothUpTo << std::endl
             << "    Savitzky-Golay radius = " << sgRadius
             << std::endl << std::endl;
    }

    if( slicing == SLICE_SG ){

             slog << "The standard gauge (Bona-Masso & Gamma-driver):"
                  << std::endl << std::endl
                  << "    Gamma-driver dissipation = " << eta
                  << std::endl << std::endl;

        }

    if( slicing == SLICE_KD ){

             slog << "The K-driver parabolic gauge on the lapse:"
                  << std::endl << std::endl
                  << "    Effective diffusion constant e = " << Kdiff
                  << "    'Elastic' constant = " << Kelas
                  << std::endl << std::endl;

        }

    if( slicing == SLICE_MS4D ){

             slog << "The K-drived maximal slicing (BVP):" << std::endl << std::endl
                  << "    'Elastic' constant = " << Kelas
                  << std::endl << std::endl;

        }
    if( slicing == SLICE_MS4G ){

             slog << "The maximal slicing (BVP) & Gamma-driver:"
                  << std::endl << std::endl
                  << "    Gamma-driver dissipation = " << eta
                  << std::endl << std::endl;
        }

    if ( mpiSize() > 1 &&
        ( slicing == SLICE_MS2 || slicing == SLICE_MS2OPT || slicing == SLICE_MS4
            || slicing == SLICE_MS6 || slicing == SLICE_MS4D || slicing == SLICE_MS4G ) )
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
    /*if( slicing != SLICE_KD || slicing != SLICE_SG || slicing != SLICE_MS2
        || slicing != SLICE_MS2OPT || slicing != SLICE_MS4 || slicing != SLICE_MS4D
        || slicing != SLICE_MS6 )*/
    if( slicing == SLICE_KD || slicing == SLICE_SG || slicing == SLICE_MS2
        || slicing == SLICE_MS2OPT || slicing == SLICE_MS4 || slicing == SLICE_MS4D
        || slicing == SLICE_MS6 )
    {
        integ.keepConstant( { fld::q } );
    }  // GFs that are kept constant in time

    integ.keepEvolved( fld::bimEvolvedGF ); // GFs that are evolved by the integrator

    if( isGR () || slicing == SLICE_CONSTGF ) {
        integ.keepConstant( { fld::fAlp, fld::fDAlp } );
    }

    if ( slicing == SLICE_CONSTG || slicing == SLICE_CONSTGF ) {
        integ.keepConstant( { fld::gAlp, fld::gDAlp } );
    }
    else if ( slicing == SLICE_SG )
    {
        const static std::vector<fld::EvolvedBy> evolvedGaugeGF = {
            { fld::gAlp,  fld::gAlp_t  },
            { fld::gDAlp, fld::gDAlp_t },
            { fld::q,     fld::q_t     },
            { fld::Bq,    fld::Bq_t    }
        };
        integ.keepEvolved( evolvedGaugeGF );
    }
    else if ( slicing == SLICE_MS4G )
    {
        const static std::vector<fld::EvolvedBy> evolvedGaugeGF = {
            { fld::q,     fld::q_t     },
            { fld::Bq,    fld::Bq_t    }
        };
        integ.keepEvolved( evolvedGaugeGF );
    }
    else if ( slicing == SLICE_KD )
    {
        const static std::vector<fld::EvolvedBy> evolvedGaugeGF = {
            { fld::gAlp,  fld::gAlp_t  }
            //{ fld::gDAlp, fld::gDAlp_t },
        };
        integ.keepEvolved( evolvedGaugeGF );
    }
    else if( slicing == SLICE_MS2 || slicing == SLICE_MS2OPT || slicing == SLICE_MS4
                || slicing == SLICE_MS4D || slicing == SLICE_MS6 )
    {
        const static std::vector<fld::EvolvedBy> evolvedGaugeGF = {
            { fld::q,     fld::q_t     },
            { fld::Bq,    fld::Bq_t    }
        };
        integ.keepEvolved( evolvedGaugeGF );
    }

    #endif // _TEST_MODE

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
    #if ! _TEST_MODE

    for( Int i = 0; i < nGhost + 1; ++i )
    {

        Int n  = nGhost - i - 1;
        Int nR = nGhost + i;

        /// Here we impose the parity conditions at the left boundary.

        t     (m,n) = t      (m,nR);
        r     (m,n) = -r     (m,nR);

        p    (m,n) = -p    (m,nR);
        pr   (m,n) = pr    (m,nR);

        gconf (m,n) = gconf  (m,nR);      fconf (m,n)  = fconf  (m,nR);
        gDconf(m,n) = -gDconf(m,nR);      fDconf(m,n)  = -fDconf(m,nR);
        gtrK  (m,n) = gtrK   (m,nR);      ftrK  (m,n)  = ftrK   (m,nR);
        gA    (m,n) = gA     (m,nR);      fA    (m,n)  = fA     (m,nR);
        gB    (m,n) = gB     (m,nR);      fB    (m,n)  = fB     (m,nR);
        gDA   (m,n) = -gDA   (m,nR);      fDA   (m,n)  = -fDA   (m,nR);
        gDB   (m,n) = -gDB   (m,nR);      fDB   (m,n)  = -fDB   (m,nR);
        gA1   (m,n) = gA1    (m,nR);      fA1   (m,n)  = fA1    (m,nR);
        gA2   (m,n) = gA2    (m,nR);      fA2   (m,n)  = fA2    (m,nR);
        gL    (m,n) = -gL    (m,nR);      fL    (m,n)  = -fL    (m,nR);
        gsig  (m,n) = gsig   (m,nR);      fsig  (m,n)  = fsig   (m,nR);
        gAsig (m,n) = gAsig  (m,nR);      fAsig (m,n)  = fAsig  (m,nR);

        #if _EVOLVE_DSIG

            gDsig (m,n) = -gDsig (m,nR);
            fDsig (m,n) = -fDsig (m,nR);

        #endif // _EVOLVE_DSIG

        pfD  (m,n) = pfD   (m,nR);
        pfS  (m,n) = -pfS  (m,nR);
        pftau(m,n) = pftau (m,nR);

        q    (m,n) = -q    (m,nR);
        Bq   (m,n) = -Bq   (m,nR);
        qr   (m,n) = qr    (m,nR);
        gAlp (m,n) = gAlp  (m,nR);
        fAlp (m,n) = fAlp  (m,nR);
        gDAlp(m,n) = -gDAlp(m,nR);
        fDAlp(m,n) = -fDAlp(m,nR);

        Lt   (m,n) = Lt    (m,nR);
        Lt2  (m,n) = Lt2   (m,nR);
        R    (m,n) = R     (m,nR);
        pfv  (m,n) = -pfv  (m,nR);

        gBetr    (m,n) = gBetr   (m,nR);
        fBetr    (m,n) = fBetr   (m,nR);
        gDconfr  (m,n) = gDconfr (m,nR);
        gDAlpr   (m,n) = gDAlpr  (m,nR);
        gLr      (m,n) = gLr     (m,nR);
        gBr      (m,n) = -gBr     (m,nR);
        fDconfr  (m,n) = fDconfr (m,nR);
        fDAlpr   (m,n) = fDAlpr  (m,nR);
        fLr      (m,n) = fLr     (m,nR);

    }
    #endif // _TEST_MODE
}

void BimetricEvolve::applyRightBoundaryCondition( Int m )
{
    #if ! _TEST_MODE

    for( Int n = nGhost + nLen; n < 2*nGhost + nLen; ++n )
    {
        r(m,n) = r(m,n-1) + delta_r;
    }

    if( false )//&& BV == BVC )
    {

        for( Int n = nGhost + nLen; n < nTotal + 1; ++n )
        {
            t(m,n) = t(m,n-1);

            gconf (m,n) = - log( r(m,n) ) + log( r(m,n) + 1 );
            gDconf(m,n) = - 1 / ( r(m,n) * (r(m,n) + 1) );

            gA    (m,n) = 1;
            gB    (m,n) = 1;
            gDA   (m,n) = 0;
            gDB   (m,n) = 0;
            gtrK  (m,n) = 0;

            gA1   (m,n) = 0;
            gA2   (m,n) = 0;

            gL    (m,n) = 0;

            gsig  (m,n) = 0;
            gAsig (m,n) = 0;

            #if _EVOLVE_DSIG

                extrapolate_R( fld::gDsig, m, n );
                extrapolate_R( fld::fDsig, m, n );

            #endif // _EVOLVE_DSIG

            pfD   (m,n) = 0;
            pfS   (m,n) = 0;
            pftau (m,n) = 0;

            gDconfr (m,n) = - 1 / ( r(m,n) * r(m,n) * (r(m,n) + 1) );
            gDAlpr  (m,n) = gDAlp(m,n) / r(m,n);
            gLr     (m,n) = 0;

        //calculateDerivedVariables( m, n ); // Calculate R, Lt, and pfv.

            extrapolate_R( fld::pfv, m, n ); // Nevertheless, we extrapolate pfv!

        }

    //smoothenGF0( m, nLen + nGhost - 3, nLen + nGhost + 3, nGhost/2,  fld::gconf, \
            fld::tmp,  fld::gconf,      1 );
    //smoothenGF0( m, nLen + nGhost - 3, nLen + nGhost + 3, nGhost/2,  fld::gDconf, \
            fld::tmp,  fld::gDconf,    -1 );

    } else {

        for( Int n = nGhost + nLen; n < nTotal + 1; ++n )
        {
            t(m,n) = t(m,n-1);

            extrapolate_R( fld::p,     m, n );   // Extrapolate [n] from [n-1], [n-2], ...
            extrapolate_R( fld::pr,    m, n );

            extrapolate_R( fld::gconf, m, n );     extrapolate_R( fld::fconf, m, n );
            extrapolate_R( fld::gDconf,m, n );     extrapolate_R( fld::fDconf,m, n );
            extrapolate_R( fld::gtrK,  m, n );     extrapolate_R( fld::ftrK,  m, n );
            extrapolate_R( fld::gA,    m, n );     extrapolate_R( fld::fA,    m, n );
            extrapolate_R( fld::gB,    m, n );     extrapolate_R( fld::fB,    m, n );
            extrapolate_R( fld::gDA,   m, n );     extrapolate_R( fld::fDA,   m, n );
            extrapolate_R( fld::gDB,   m, n );     extrapolate_R( fld::fDB,   m, n );
            extrapolate_R( fld::gA1,   m, n );     extrapolate_R( fld::fA1,   m, n );
            extrapolate_R( fld::gA2,   m, n );     extrapolate_R( fld::fA2,   m, n );
            extrapolate_R( fld::gL,    m, n );     extrapolate_R( fld::fL,    m, n );
            extrapolate_R( fld::gsig,  m, n );     extrapolate_R( fld::fsig,  m, n );
            extrapolate_R( fld::gAsig, m, n );     extrapolate_R( fld::fAsig, m, n );

            #if _EVOLVE_DSIG

                extrapolate_R( fld::gDsig, m, n );
                extrapolate_R( fld::fDsig, m, n );

            #endif // _EVOLVE_DSIG

            extrapolate_R( fld::pfD,   m, n );
            extrapolate_R( fld::pfS,   m, n );
            extrapolate_R( fld::pftau, m, n );

            extrapolate_R( fld::q,     m, n );
            extrapolate_R( fld::Bq,    m, n );
            extrapolate_R( fld::qr,    m, n );
            //extrapolate_R( fld::gBet,  m, n );
            //extrapolate_R( fld::fBet,  m, n );
            //extrapolate_LIN4( fld::gAlp,  m, n );
            //extrapolate_R( fld::gAlp_t,  m, n );
            //extrapolate_R( fld::gDAlp, m, n );
            extrapolate_R( fld::fAlp,  m, n );
            extrapolate_R( fld::fDAlp, m, n );

            //extrapolate_R( fld::gBetr,   m, n );
            //extrapolate_R( fld::fBetr,   m, n );
            extrapolate_R( fld::gDconfr, m, n );
            //extrapolate_R( fld::gDAlpr,  m, n );
            extrapolate_R( fld::gLr,     m, n );
            extrapolate_R( fld::gBr,     m, n );
            extrapolate_R( fld::fDconfr, m, n );
            extrapolate_R( fld::fDAlpr,  m, n );
            extrapolate_R( fld::fLr,     m, n );

            //calculateDerivedVariables( m, n ); // Calculate R, Lt, and pfv.

            extrapolate_R( fld::pfv, m, n ); // Nevertheless, we extrapolate pfv!

        }
    }
    #endif // _TEST_MODE
}

void BimetricEvolve::applySommerfeldBC
(
    Int m1, Int n, Int m, Int gf,
    const Real background,
    const Int fall_off
    )
{
    #if ! _TEST_MODE

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

    #endif // _TEST_MODE
}

void BimetricEvolve::integStep_Prepare( Int m )
{
    #if !_TEST_MODE

    OMP_parallel_for( Int n = 0+0*nGhost; n < 2*nGhost + nLen + 1; ++n )
    {
        calculateDerivedVariables( m, n );
    }

    #endif // _TEST_MODE
}
/* this function is not used in this file */
void BimetricEvolve::determineGaugeFunctions( Int m )
{
    #if !_TEST_MODE
    /////////////////////////////////////////////////////////////////////////////////////
    /// - Maximal slicing (optional) -- this will calculate gAlp and gDAlp

    Int sliceBC = smooth >= 1 ? round(40/delta_r) : nLen/2;
    // sliceBC = round(20/delta_r);
    switch( slicing )
    {
        case SLICE_MS2OPT: maximalSlice_2opt       ( m, sliceBC, 1 );  break;
        case SLICE_MS2:    maximalSlice_2          ( m, sliceBC, 1 );  break;
        case SLICE_MS4:    maximalSlice_4          ( m, sliceBC, 1 );  break;
        case SLICE_MS6:    maximalSlice_6_gDAlp    ( m, sliceBC, 1 );  break;
        case SLICE_MS4D:   maximalSlice_drived_4   ( m, sliceBC, 1 );  break;
        case SLICE_MS4G:   maximalSlice_4_gDAlp    ( m, sliceBC, 1 );  break;
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


        OMP_parallel_for( Int n = 0 + 0*nGhost; n < nGhost + 1; ++n )
        {
            gAlp_r ( m, n ) = GF_right_r( gAlp,  m, n );
            gDAlp_r( m, n ) = GF_right_r( gDAlp, m, n );

            fAlp   ( m, n ) = gAlp( m, n );
            fAlp_r ( m, n ) = GF_right_r ( fAlp,  m, n );
            fAlp_rr( m, n ) = GF_right_rr( fAlp, m, n );

            fDAlp  ( m, n ) = fAlp_r(m,n) / (TINY_Real + fAlp(m,n));
            fDAlp_r( m, n ) = ( fAlp_rr(m,n) - fAlp_r(m,n) * fAlp_r(m,n)
                               / (TINY_Real + fAlp(m,n)) ) / (TINY_Real + fAlp(m,n));

            /*Real dbg1 = fAlp(m,n);
            Real dbg2 = fDAlp(m,n);
            Real dbg3 = fDAlp_r(m,n);*/

        }

        OMP_parallel_for( Int n = nGhost +1 ; n < nGhost + nLen + 1; ++n )
        {
            gAlp_r ( m, n ) = GF_r( gAlp,  m, n );
            gDAlp_r( m, n ) = GF_r( gDAlp, m, n );

            fAlp   ( m, n ) = gAlp( m, n );
            fAlp_r ( m, n ) = GF_r ( fAlp,  m, n );
            fAlp_rr( m, n ) = GF_rr( fAlp, m, n );

            fDAlp  ( m, n ) = fAlp_r(m,n) / (TINY_Real + fAlp(m,n));
            fDAlp_r( m, n ) = ( fAlp_rr(m,n) - fAlp_r(m,n) * fAlp_r(m,n)
                               / (TINY_Real + fAlp(m,n)) )/ (TINY_Real + fAlp(m,n));

            /*Real dbg1 = fAlp(m,n);
            Real dbg2 = fDAlp(m,n);
            Real dbg3 = fDAlp_r(m,n);*/

        }

        OMP_parallel_for( Int n = nGhost + nLen + 1; n < 2*nGhost + nLen + 1; ++n )
        {
            gAlp_r ( m, n ) = GF_left_r( gAlp,  m, n );
            gDAlp_r( m, n ) = GF_left_r( gDAlp, m, n );

            fAlp   ( m, n ) = gAlp( m, n );
            fAlp_r ( m, n ) = GF_left_r ( fAlp,  m, n );
            fAlp_rr( m, n ) = GF_left_rr( fAlp, m, n );

            fDAlp  ( m, n ) = fAlp_r(m,n) / (TINY_Real + fAlp(m,n));
            fDAlp_r( m, n ) = ( fAlp_rr(m,n) - fAlp_r(m,n) * fAlp_r(m,n)
                                / (TINY_Real + fAlp(m,n)) ) / (TINY_Real + fAlp(m,n));

            /*Real dbg1 = fAlp(m,n);
            Real dbg2 = fDAlp(m,n);
            Real dbg3 = fDAlp_r(m,n);*/

        }

        cubicSplineSmooth( m, fld::gDAlp_r, lin2n, cub2n );
        cubicSplineSmooth( m, fld::fDAlp_r, lin2n, cub2n );
    }
    else // if not GR and not geodesic slicing (constant lapse)
    {
        // gAlp and gDAlp must the BC fixed at this point
        //

        OMP_parallel_for( Int n = 0 + 0*nGhost; n < nGhost + 1; ++n ) {
            gAlp_r ( m, n ) = GF_right_r( gAlp,  m, n );
            gDAlp_r( m, n ) = GF_right_r( gDAlp, m, n );
        }

        OMP_parallel_for( Int n = nGhost + 1; n < nGhost + nLen + 1; ++n ) {
            gAlp_r ( m, n ) = GF_r( gAlp,  m, n );
            gDAlp_r( m, n ) = GF_r( gDAlp, m, n );
        }

        OMP_parallel_for( Int n = nGhost + nLen + 1; n < 2*nGhost + nLen + 1; ++n ) {
            gAlp_r ( m, n ) = GF_left_r( gAlp,  m, n );
            gDAlp_r( m, n ) = GF_left_r( gDAlp, m, n );
        }

        /// @fixme smoothen gAlp_r, gDAlp_r

        OMP_parallel_for( Int n = nGhost; n < 2*nGhost + nLen; ++n )
        // FT: I am including also the ghostcells
        {
            gW( m, n )   = eq_gW( m, n );
            fW( m, n )   = eq_fW( m, n );
            fAlp( m, n ) = eq_fAlp( m, n );
        }

        if( smooth ) {
            smoothenGF( m, fld::fAlp, fld::tmp, fld::fAlp, 1 );
        }
        else {
            applyBoundaryConditions( m, fld::fAlp, +1 );
        } // FT to MK: why do you set these up as alternatives?

        OMP_parallel_for( Int n = nGhost; n < 2*nGhost + nLen; ++n ) {
            fAlp_r( m, n ) = GF_r( fAlp, m, n );
        }

        if( smooth ) {
            smoothenGF( m, fld::fAlp_r, fld::tmp, fld::fAlp_r, -1 );
        }
        else {
            applyBoundaryConditions( m, fld::fAlp_r, -1 );
        }

        OMP_parallel_for( Int n = nGhost; n < 2*nGhost + nLen; ++n ) {
            fAlp_rr( m, n ) = GF_r( fAlp_r, m, n );
        }

        if( smooth ) {
            smoothenGF( m, fld::fAlp_rr, fld::tmp, fld::fAlp_rr, 1 );
        }
        else {
            applyBoundaryConditions( m, fld::fAlp_rr, 1 );
        }

        OMP_parallel_for( Int n = nGhost; n < 2*nGhost + nLen; ++n )
        {
            fDAlp  ( m, n ) = fAlp_r(m,n) / (TINY_Real + fAlp(m,n));
            fDAlp_r( m, n ) = ( fAlp_rr(m,n) - fAlp_r(m,n) * fAlp_r(m,n)
                               / (TINY_Real + fAlp(m,n)) )
                               / (TINY_Real + fAlp(m,n));

        }

        //cubicSplineSmooth( m, fld::gDAlp_r, lin2n, cub2n );
        //cubicSplineSmooth( m, fld::fDAlp_r, lin2n, cub2n );
    }

    #endif // _TEST_MODE
}

/////////////////////////////////////////////////////////////////////////////////////////
// Equations of Motion (generated in Mathematica)
/////////////////////////////////////////////////////////////////////////////////////////

#if _TEST_MODE

    #include "eom-test/eomTest.h"

#else

    #include "eom-BSSN.h"

#endif // _TEST_MODE

void BimetricEvolve::integStep_CalcEvolutionRHS( Int m )
{

    #if ! _TEST_MODE

        if( m < mSmoothUpTo && smooth >= 2 )
        {
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::q,
                         fld::tmp,  fld::q,      -1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::qr,
                         fld::tmp,  fld::qr,      1 );
            //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::p,  \
                           fld::tmp,  fld::p,      1 );
            //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::pr,  \
                           fld::tmp,  fld::pr,      1 );
            //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gBet, \
                           fld::tmp,  fld::gBet, -1 );
            //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gBetr, \
                           fld::tmp,  fld::gBetr, 1 );

            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gA,
                         fld::tmp,  fld::gA,      1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gB,
                         fld::tmp,  fld::gB,      1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gDA,
                         fld::tmp,  fld::gDA,     -1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gDB,
                         fld::tmp,  fld::gDB,     -1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gA1,
                         fld::tmp,  fld::gA1,     1 );

            #if OBSERVER == 1

            #else

                smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gA2,     \
                             fld::tmp,  fld::gA2,     1 );

            #endif // OBSERVER

            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gtrK,
                         fld::tmp,  fld::gtrK,    1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gconf,
                         fld::tmp,  fld::gconf,   1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gDconf,
                         fld::tmp,  fld::gDconf,  -1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gL,
                         fld::tmp,  fld::gL,      -1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gLr,
                         fld::tmp,  fld::gLr,      1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gDAlpr,
                         fld::tmp,  fld::gDAlpr,   1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gDconfr,
                         fld::tmp,  fld::gDconfr, 1 );
            //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::grhobar,  \
                           fld::tmp,  fld::pfD,     1 );
            //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::pfD,    \
                            fld::tmp,  fld::pfD,     1 );
            //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gsig,    \
                           fld::tmp,  fld::gsig,     1 );

            //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gAsig,   \
                           fld::tmp,  fld::gAsig,    1 );

        }
        /////////////////////////////////////////////////////////////////////////////////////
        /// - First, calculate the values of the spatial derivatives

        #if OBSERVER == 1 // if the Eulerian equations with traceless A are used, \
                                impose gA2 = -gA1/2

            OMP_parallel_for( Int n = 0 + 0*nGhost; n < 2*nGhost + nLen + 1; ++n )
            {

                gA2     (m,n) = - gA1( m, n ) / 2;
                fA2     (m,n) = - fA1( m, n ) / 2;

            }

        #endif // OBSERVER

        OMP_parallel_for( Int n = 0 + 0*nGhost; n < nGhost +1; ++n )
        {
           // The radial derivatives of the fields inside the ratio of the lapses \
                (which is inside the maximal slicing's equation) over the left ghosts

            p_r      (m,n) = GF_right_r (p, m, n);
            q_r      (m,n) = GF_right_r (q, m, n);
            Bq_r     (m,n) = GF_right_r (Bq, m, n);
            gA2_r    (m,n) = GF_right_r (gA2, m, n);
            fA1_r    (m,n) = GF_right_r (fA1, m, n);
            gtrK_r   (m,n) = GF_right_r (gtrK, m, n);
            ftrK_r   (m,n) = GF_right_r (ftrK, m, n);
            gBr_r    (m,n) = GF_right_r (gBr, m, n);

        }

        OMP_parallel_for( Int n = nGhost +1; n < 1*nGhost + nLen + 1; ++n )
        {
           // The radial derivatives of the fields inside the ratio of the lapses \
                (which is inside the maximal slicing's equation) over the entire grid \
                but not the ghosts

            p_r      (m,n) = GF_r (p, m, n);
            q_r      (m,n) = GF_r (q, m, n);
            Bq_r     (m,n) = GF_r (Bq, m, n);
            gA2_r    (m,n) = GF_r (gA2, m, n);
            fA1_r    (m,n) = GF_r (fA1, m, n);
            gtrK_r   (m,n) = GF_r (gtrK, m, n);
            ftrK_r   (m,n) = GF_r (ftrK, m, n);
            gBr_r    (m,n) = GF_r (gBr, m, n);

        }

        OMP_parallel_for( Int n = nGhost + nLen + 1; n < 2*nGhost + nLen + 1; ++n )
        {
           // The radial derivatives of the fields inside the ratio of the lapses \
                (which is inside the maximal slicing's equation) over the right ghosts

            p_r      (m,n) = GF_left_r (p, m, n);
            q_r      (m,n) = GF_left_r (q, m, n);
            Bq_r     (m,n) = GF_left_r (Bq, m, n);
            gA2_r    (m,n) = GF_left_r (gA2, m, n);
            fA1_r    (m,n) = GF_left_r (fA1, m, n);
            gtrK_r   (m,n) = GF_left_r (gtrK, m, n);
            ftrK_r   (m,n) = GF_left_r (ftrK, m, n);
            gBr_r    (m,n) = GF_left_r (gBr, m, n);

        }

        if( smooth >= 2 )
        {
            smoothenGF( m, fld::p_r,      fld::tmp, fld::p_r,     +1 );
            smoothenGF( m, fld::q_r,      fld::tmp, fld::q_r,     +1 );
            smoothenGF( m, fld::Bq_r,     fld::tmp, fld::Bq_r,    +1 );
            smoothenGF( m, fld::gA2_r,    fld::tmp, fld::gA2_r,   -1 );
            smoothenGF( m, fld::fA1_r,    fld::tmp, fld::fA1_r,   -1 );
            smoothenGF( m, fld::gtrK_r,   fld::tmp, fld::gtrK_r,  -1 );
            smoothenGF( m, fld::ftrK_r,   fld::tmp, fld::ftrK_r,  -1 );
        }


        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        /// - Calculate the gauge (e.g., do maximal slicing)
        ///
        switch( slicing )
        {
            case SLICE_MS2OPT: maximalSlice_2opt       ( m, nLen/2, 1 );  break;
            case SLICE_MS2:    maximalSlice_2          ( m, nLen/2, 1 );  break;
            case SLICE_MS4:    maximalSlice_4_gDAlp    ( m, nLen/2, 1 );  break;
            case SLICE_MS6:    maximalSlice_6_gDAlp    ( m, nLen/2, 1 );  break;
            case SLICE_MS4D:   maximalSlice_drived_4   ( m, nLen/2, 1 );  break;
            case SLICE_MS4G:   maximalSlice_4_gDAlp    ( m, nLen/2, 1 );  break;
        }

        //if( m < mSmoothUpTo ){
         //   smoothenGF0 ( m, 0, nLen + 2*nGhost + 1, 32,  fld::gAlp,  \
                           fld::tmp,  fld::gAlp,    1 );
        //}
        //smoothenGF0 ( m, 0, nLen + 2*nGhost + 1, 32,  fld::gDAlp,   \
                       fld::tmp,  fld::gDAlp,  -1 );

        /// @todo fixme:  Determine fAlp right after maximal slicing!?

        if ( isGR () /*|| slicing == SLICE_CONSTGF*/ )
        {

            if( slicing == SLICE_KD )
            {
                const Real h = delta_r;

                if( t( m, nGhost ) == 0 )
                {
                    maximalSlice_4_gDAlp    ( m, nGhost + nLen, 1 );
                }

                for( Int n = nGhost + nLen - 1; n < nTotal + 1; ++n )
                {
                    //extrapolate_LIN4( fld::gAlp,  m, n );
                    gAlp( m, n+1 ) =
                    ( ( 4  + 2 * h * h * QQ( m,n) ) * gAlp(m,n)
                    - ( 2 + h * PP( m,n) ) * gAlp(m,n-1) )
                    / ( 2 - h * PP( m,n) );
                }

                smoothenGF0 ( m, nGhost, nGhost + nLen, 15,  fld::gAlp,
                             fld::tmp,  fld::gAlp,     1 );
                //smoothenGF ( m,  fld::gAlp,    fld::tmp,  fld::gAlp,     1 );

                //for( Int n = 0; n < nTotal; ++n )
                //{

                 //   if( gAlp (m,n) < 0 ){

                   //     gAlp (m,n) = 0;

                    //}

                //}

                for( Int i = 0; i < nGhost + 1; ++i )
                {

                    Int n  = nGhost - i - 1;
                    Int nR = nGhost + i;

                    /// Here we impose the parity conditions at the left boundary.

                    gAlp (m,n) = gAlp  (m,nR);

                }

                for( Int n = nGhost + nLen - 1; n < nTotal + 1; ++n )
                {
                    //extrapolate_LIN4( fld::gAlp,  m, n );
                    gAlp( m, n+1 ) =
                    ( ( 4  + 2 * h * h * QQ( m,n) ) * gAlp(m,n)
                    - ( 2 + h * PP( m,n) ) * gAlp(m,n-1) )
                    / ( 2 - h * PP( m,n) );
                }

                //cubicSplineSmooth( m, fld::gAlp, lin2n, cub2n );
                //cubicSplineSmooth( m, fld::fAlp, lin2n, cub2n );

                OMP_parallel_for( Int n = 0 + 0*nGhost; n < nGhost +1; ++n )
                {
                    gDAlp( m, n ) = GF_right_r( gAlp, m, n );
                }

                OMP_parallel_for( Int n = nGhost +1 ; n < nGhost + nLen + 1; ++n )
                {
                    gDAlp( m, n ) = GF_r( gAlp, m, n );
                }

                OMP_parallel_for( Int n = nGhost + nLen + 1; n < nTotal + 1; ++n )
                {
                    gDAlp( m, n ) = GF_left_r( gAlp, m, n );
                    //extrapolate_LIN4( fld::gDAlp,  m, n );
                }

                for( Int i = 0; i < nGhost + 1; ++i )
                {

                    Int n  = nGhost - i - 1;
                    Int nR = nGhost + i;

                    /// Here we impose the parity conditions at the left boundary.

                    gDAlp(m,n) = -gDAlp(m,nR);

                }

                smoothenGF0 ( m, nGhost, nGhost + nLen, 10,  fld::gDAlp,
                                fld::tmp,  fld::gDAlp,     -1 );
                //smoothenGF ( m,  fld::gDAlp,    fld::tmp,  fld::gDAlp,     -1 );

                for( Int i = 0; i < nGhost + 1; ++i )
                {

                    Int n  = nGhost - i - 1;
                    Int nR = nGhost + i;

                    /// Here we impose the parity conditions at the left boundary.

                    gDAlp(m,n) = -gDAlp(m,nR);

                }

                OMP_parallel_for( Int n = nGhost + nLen + 1; n < nTotal + 1; ++n )
                {
                    gDAlp( m, n ) = GF_left_r( gAlp, m, n );
                    //extrapolate_LIN4( fld::gDAlp,  m, n );
                }

                //for( Int n = nGhost + nLen; n < nTotal + 1; ++n )
                //{
                //    extrapolate_R( fld::gDAlp,  m, n );
                //}

            }

            // At this point, gAlp and gDAlp must be known, and possibly splined

                OMP_parallel_for( Int n = 0 + 0*nGhost; n < nGhost + 1; ++n )
                {
                    gAlp_r ( m, n ) = GF_right_r( gAlp,  m, n );
                    gDAlp_r( m, n ) = GF_right_r( gDAlp, m, n );

                /*fAlp   ( m, n ) = gAlp( m, n );
                fAlp_r ( m, n ) = GF_right_r ( fAlp,  m, n );
                fAlp_rr( m, n ) = GF_right_rr( fAlp, m, n );

                fDAlp  ( m, n ) = fAlp_r(m,n) / (TINY_Real + fAlp(m,n));
                fDAlp_r( m, n ) = ( fAlp_rr(m,n) - fAlp_r(m,n) * fAlp_r(m,n)
                                   / (TINY_Real + fAlp(m,n)) ) / (TINY_Real + fAlp(m,n));*/

                }

                OMP_parallel_for( Int n = nGhost +1 ; n < nGhost + nLen + 1; ++n )
                {
                    gAlp_r ( m, n ) = GF_r( gAlp,  m, n );
                    gDAlp_r( m, n ) = GF_r( gDAlp, m, n );

                /*fAlp   ( m, n ) = gAlp( m, n );
                fAlp_r ( m, n ) = GF_r ( fAlp,  m, n );
                fAlp_rr( m, n ) = GF_rr( fAlp, m, n );

                fDAlp  ( m, n ) = fAlp_r(m,n) / (TINY_Real + fAlp(m,n));
                fDAlp_r( m, n ) = ( fAlp_rr(m,n) - fAlp_r(m,n) * fAlp_r(m,n)
                                   / (TINY_Real + fAlp(m,n)) )/ (TINY_Real + fAlp(m,n));*/

                }

                OMP_parallel_for( Int n = nGhost + nLen + 1; n <  nTotal + 1; ++n )
                {
                    gAlp_r ( m, n ) = GF_left_r( gAlp,  m, n );
                    gDAlp_r( m, n ) = GF_left_rr (gAlp, m, n);
                    //extrapolate_LIN4( fld::gDAlp_r,  m, n );
                    //gDAlp_r( m, n ) = GF_left_r( gDAlp, m, n );

                /*fAlp   ( m, n ) = gAlp( m, n );
                fAlp_r ( m, n ) = GF_left_r ( fAlp,  m, n );
                fAlp_rr( m, n ) = GF_left_rr( fAlp, m, n );

                fDAlp  ( m, n ) = fAlp_r(m,n) / (TINY_Real + fAlp(m,n));
                fDAlp_r( m, n ) = ( fAlp_rr(m,n) - fAlp_r(m,n) * fAlp_r(m,n)
                                   / (TINY_Real + fAlp(m,n)) ) / (TINY_Real +fAlp(m,n));*/

                }

                //smoothenGF ( m,  fld::gDAlp_r,    fld::tmp,  fld::gDAlp_r,     1 );
                smoothenGF0 ( m, nGhost, nGhost + nLen, 10,  fld::gDAlp_r,
                             fld::tmp,  fld::gDAlp_r,     1 );

                for( Int i = 0; i < nGhost + 1; ++i )
                {

                    Int n  = nGhost - i - 1;
                    Int nR = nGhost + i;

                    /// Here we impose the parity conditions at the left boundary.

                    gDAlp_r(m,n) = gDAlp_r(m,nR);

                }

                OMP_parallel_for( Int n = nGhost + nLen + 1; n < nTotal + 1; ++n )
                {
                    gDAlp_r( m, n ) = GF_left_r (gDAlp, m, n);
                    //extrapolate_LIN4( fld::gDAlp_r,  m, n );
                    //gDAlp_r( m, n ) = GF_left_r( gDAlp, m, n );

                }

            //cubicSplineSmooth( m, fld::gDAlp_r, lin2n, cub2n );
            //cubicSplineSmooth( m, fld::fDAlp_r, lin2n, cub2n );

            OMP_parallel_for( Int n = 0; n < 2*nGhost + nLen + 1; ++n )
            {
                gDAlpr ( m, n ) = gDAlp ( m, n ) / r( m, n );
                fDAlpr ( m, n ) = fDAlp ( m, n ) / r( m, n );

            }
        }
        else // if not GR and not geodesic slicing (constant lapse)
             // SLICE_KD is not implemented here yet
        {
            // gAlp and gDAlp must the BC fixed at this point
            //

            OMP_parallel_for( Int n = 0 + 0*nGhost; n < nGhost + 1; ++n ) {
                gAlp_r ( m, n ) = GF_right_r( gAlp,  m, n );
                gDAlp_r( m, n ) = GF_right_r( gDAlp, m, n );
            }

            OMP_parallel_for( Int n = nGhost + 1; n < nGhost + nLen + 1; ++n ) {
                gAlp_r ( m, n ) = GF_r( gAlp,  m, n );
                gDAlp_r( m, n ) = GF_r( gDAlp, m, n );
            }

            OMP_parallel_for( Int n = nGhost + nLen + 1; n < 2*nGhost + nLen + 1; ++n ) {
                gAlp_r ( m, n ) = GF_left_r( gAlp,  m, n );
                gDAlp_r( m, n ) = GF_left_r( gDAlp, m, n );
            }

            /// @fixme smoothen gAlp_r, gDAlp_r

            OMP_parallel_for( Int n = 0; n < 2*nGhost + nLen; ++n )
            /// FT: I am including also the ghostcells
            {
                gW( m, n )   = eq_gW( m, n );
                fW( m, n )   = eq_fW( m, n );
                fAlp( m, n ) = eq_fAlp( m, n );

                /*Real dbg1 = gW( m, n );
                Real dbg2 = eq_gW( m, n );
                Real dbg3 = fW( m, n );
                Real dbg4 = eq_fW( m, n );
                Real dbg5 = fAlp( m, n );
                Real dbg6 = eq_fAlp( m, n );*/

           /* Real dbg1  = eq_gW(m,n);
            Real dbg2  = eq_fW(m,n);
            Real dbg51 = eq_gW(m,n)/eq_fW(m,n);

            Real dbg3  = gA1(m,n);
            Real dbg4  = fA1(m,n);
            Real dbg5  = gconf(m,n);
            Real dbg6  = fconf(m,n);
            Real dbg7  = gDB(m,n);
            Real dbg8  = fDB(m,n);
            Real dbg9  = gDconf(m,n);
            Real dbg10 = fDconf(m,n);
            Real dbg11 = eq_pf_gJ11(m,n);
            Real dbg12 = eq_pf_fJ11(m,n);
            Real dbg13 = eq_pf_gJ22(m,n);
            Real dbg14 = eq_pf_fJ22(m,n);
            Real dbg15 = eq_pf_gj(m,n);
            Real dbg16 = eq_pf_fj(m,n);
            Real dbg17 = eq_pf_grho(m,n);
            Real dbg18 = eq_pf_frho(m,n);
            Real dbg19 = gA(m,n);
            Real dbg20 = fA(m,n);
            Real dbg21 = gB(m,n);
            Real dbg22 = fB(m,n);
            Real dbg23 = Lt(m,n);
            Real dbg24 = Lt2(m,n);
            Real dbg25 = P_0_2(R(m,n));
            Real dbg26 = P_0_3(R(m,n));
            Real dbg27 = P_1_0(R(m,n));
            Real dbg28 = P_1_1(R(m,n));
            Real dbg29 = P_1_2(R(m,n));
            Real dbg30 = P_1_3(R(m,n));
            Real dbg31 = P_2_0(R(m,n));
            Real dbg32 = P_2_1(R(m,n));
            Real dbg33 = P_2_2(R(m,n));
            Real dbg34 = p(m,n);
            Real dbg35 = R(m,n);
            Real dbg36 = gtrK(m,n);
            Real dbg37 = ftrK(m,n);

            Real dbg38 = gA2_r(m,n);
            Real dbg39 = fA1_r(m,n);
            Real dbg40 = p_r(m,n);
            Real dbg41 = gtrK_r(m,n);
            Real dbg42 = ftrK_r(m,n);

            Real dbg49 = 3 * gA1(m,n) * (1 + gDB(m,n)
          * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n));

            Real dbg50 = - r(m,n) * (3 * gA2_r(m,n)
          + gtrK_r(m,n));*/

            }

            if( smooth ) {
                smoothenGF( m, fld::fAlp, fld::tmp, fld::fAlp, 1 );
            }
            else {
                applyBoundaryConditions( m, fld::fAlp, +1 );
            } /// FT to MK: why do you set these up as alternatives?

            ////////////////////////////////////////////////////////////////

            OMP_parallel_for( Int n = 0 + 0*nGhost; n < nGhost + 1; ++n ) {
                fAlp_r ( m, n ) = GF_right_r( fAlp,  m, n );
            }

            OMP_parallel_for( Int n = nGhost + 1; n < nGhost + nLen + 1; ++n ) {
                fAlp_r( m, n ) = GF_r( fAlp, m, n );
            }

            OMP_parallel_for( Int n = nGhost + nLen + 1; n < 2*nGhost + nLen + 1; ++n ) {
                fAlp_r ( m, n ) = GF_left_r( fAlp,  m, n );
            }

            if( smooth ) {
                smoothenGF( m, fld::fAlp_r, fld::tmp, fld::fAlp_r, -1 );
            }
            else {
                applyBoundaryConditions( m, fld::fAlp_r, -1 );
            }

            ////////////////////////////////////////////////////////////////

            OMP_parallel_for( Int n = 0 + 0*nGhost; n < nGhost + 1; ++n ) {
                fAlp_rr ( m, n ) = GF_right_r( fAlp_r,  m, n );
            }

            OMP_parallel_for( Int n = nGhost + 1; n < nGhost + nLen + 1; ++n ) {
                fAlp_rr( m, n ) = GF_r( fAlp_r, m, n );
            }

            OMP_parallel_for( Int n = nGhost + nLen + 1; n < 2*nGhost + nLen + 1; ++n ) {
                fAlp_rr ( m, n ) = GF_left_r( fAlp_r,  m, n );
            }

            if( smooth ) {
                smoothenGF( m, fld::fAlp_rr, fld::tmp, fld::fAlp_rr, 1 );
            }
            else {
                applyBoundaryConditions( m, fld::fAlp_rr, 1 );
            }

            ////////////////////////////////////////////////////////////////

            OMP_parallel_for( Int n = 0; n < 2*nGhost + nLen; ++n )
            {
                fDAlp  ( m, n ) = fAlp_r(m,n) /*/ (TINY_Real + fAlp(m,n))*/;

            }
             OMP_parallel_for( Int n = 0; n < 2*nGhost + nLen; ++n )
            {
                fDAlp_r( m, n ) = GF_r( gDAlp, m, n );
                /*( fAlp_rr(m,n) - fAlp_r(m,n) * fAlp_r(m,n) / (TINY_Real + fAlp(m,n)) )
                                   / (TINY_Real + fAlp(m,n));*/

            }

            ////////////////////////////////////////////////////////////////

            cubicSplineSmooth( m, fld::gDAlp_r, lin2n, cub2n );
            cubicSplineSmooth( m, fld::fDAlp_r, lin2n, cub2n );

            OMP_parallel_for( Int n = 0; n < 2*nGhost + nLen + 1; ++n )
            {
                gDAlpr ( m, n ) = gDAlp ( m, n ) / r( m, n );
                fDAlpr ( m, n ) = fDAlp ( m, n ) / r( m, n );

            }
        }

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////


        if( m < mSmoothUpTo && smooth >= 2 )
        {
            /*smoothenGF( m, fld::gAlp_r,      fld::tmp, fld::gAlp_r,     -1 );
            smoothenGF( m, fld::fAlp_r,      fld::tmp, fld::fAlp_r,     -1 );
            smoothenGF( m, fld::gDAlp_r,     fld::tmp, fld::gDAlp_r,    +1 );
            smoothenGF( m, fld::fDAlp_r,     fld::tmp, fld::fDAlp_r,    +1 );*/
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo + CFDS_ORDER / 2, 32,  fld::gAlp_r,
                         fld::tmp,  fld::gAlp_r,     -1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo + CFDS_ORDER / 2, 32,  fld::fAlp_r,
                         fld::tmp,  fld::fAlp_r,     -1 );
            //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo + CFDS_ORDER / 2, 32, fld::gDAlp,\
                           fld::tmp,  fld::gDAlp,     -1 );
            //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo + CFDS_ORDER / 2, 32, fld::fDAlp,\
                           fld::tmp,  fld::fDAlp,     -1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo + CFDS_ORDER / 2, 32,  fld::gDAlp_r,
                         fld::tmp,  fld::gDAlp_r,     +1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo + CFDS_ORDER / 2, 32,  fld::fDAlp_r,
                         fld::tmp,  fld::fDAlp_r,     +1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo + CFDS_ORDER / 2, 32,  fld::gDAlpr,
                         fld::tmp,  fld::gDAlpr, 1 );
            smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo + CFDS_ORDER / 2, 32,  fld::fDAlpr,
                         fld::tmp,  fld::fDAlpr, 1 );
        }

        OMP_parallel_for( Int n = 0 + 0*nGhost; n < nGhost + 1; ++n )
        {

        // The radial derivatives of the evolved fields
            gconf_r  (m,n) = GF_right_r (gconf, m, n),
            fconf_r  (m,n) = GF_right_r (fconf, m, n),
            gDconf_r (m,n) = GF_right_r (gDconf, m, n),
            fDconf_r (m,n) = GF_right_r (fDconf, m, n),
            //gtrK_r   (m,n) = GF_r (gtrK, m, n),      ftrK_r   (m,n) = GF_r (ftrK, m, n),
            gA_r     (m,n) = GF_right_r (gA, m, n),
            fA_r     (m,n) = GF_right_r (fA, m, n),
            gB_r     (m,n) = GF_right_r (gB, m, n),
            fB_r     (m,n) = GF_right_r (fB, m, n),
            gDA_r    (m,n) = GF_right_r (gDA, m, n),
            fDA_r    (m,n) = GF_right_r (fDA, m, n),
            gDB_r    (m,n) = GF_right_r (gDB, m, n),
            fDB_r    (m,n) = GF_right_r (fDB, m, n),
            //gA2_r    (m,n) = GF_r (gA2, m, n),       fA1_r    (m,n) = GF_r (fA1, m, n),
            gA1_r    (m,n) = GF_right_r (gA1, m, n),
            fA2_r    (m,n) = GF_right_r (fA2, m, n),
            gL_r     (m,n) = GF_right_r (gL, m, n),
            fL_r     (m,n) = GF_right_r (fL, m, n),
            gAsig_r  (m,n) = GF_right_r (gAsig, m, n),
            fAsig_r  (m,n) = GF_right_r (fAsig, m, n),

            //Lt_r     (m,n) = GF_r (Lt, m, n),

            // The radial derivatives of the regularizing functions

            gDconfr_r(m,n) = GF_right_r (gDconfr, m, n),
            gDAlpr_r (m,n) = GF_right_r (gDAlpr, m, n),
            gLr_r    (m,n) = GF_right_r (gLr, m, n),
            fDconfr_r(m,n) = GF_right_r (fDconfr, m, n),
            fDAlpr_r (m,n) = GF_right_r (fDAlpr, m, n),
            fLr_r    (m,n) = GF_right_r (fLr, m, n);

            pr_r      (m,n) = GF_right_r (pr, m, n);
            qr_r      (m,n) = GF_right_r (qr, m, n);

        // The radial derivatives of the Valencia variables and the external sources
            pfD_r    (m,n) = GF_right_r (pfD, m, n),
            pfS_r    (m,n) = GF_right_r (pfS, m, n),
            pftau_r  (m,n) = GF_right_r (pftau, m, n),
            gj_r     (m,n) = GF_right_r (gj, m, n),
            pfv_r    (m,n) = eq_pf_gv_r(m,n);

        // The radial derivatives of the utility function R = fB/gB
            R_r      (m,n) = GF_right_r (R, m, n),

        // These are the only spatial second derivatives that are needed:

        // The second radial derivatives of the evolved fields
            gconf_rr (m,n) = GF_right_rr (gconf, m, n),
            fconf_rr (m,n) = GF_right_rr (fconf, m, n),
            gA_rr    (m,n) = GF_right_rr (gA, m, n),
            gB_rr    (m,n) = GF_right_rr (gB, m, n),
            fA_rr    (m,n) = GF_right_rr (fA, m, n),
            fB_rr    (m,n) = GF_right_rr (fB, m, n),
            gA1_rr   (m,n) = GF_right_rr (gA1, m, n),
            fA1_rr   (m,n) = GF_right_rr (fA1, m, n),
            gtrK_rr  (m,n) = GF_right_rr (gtrK, m, n),
            ftrK_rr  (m,n) = GF_right_rr (ftrK, m, n),

        // The second radial derivatives of the gauge fields
            p_rr     (m,n) = GF_right_rr (p, m, n),   // used in eq_gDA_t and eq_fDA_t
            q_rr     (m,n) = GF_right_rr (q, m, n),   // used in eq_gDA_t and eq_fDA_t
            gAlp_rr  (m,n) = GF_right_rr (gAlp, m, n),
            fAlp_rr  (m,n) = GF_right_rr (fAlp, m, n),   // used in fDAlp_r

        // The second radial derivatives of the traces of the conformal extrinsic \
            curvatures (these traces are usually set to 0 in eomBSSNObserver.h)
            gtrA_rr  (m,n) = GF_right_rr (gtrA, m, n),
            ftrA_rr  (m,n) = GF_right_rr (ftrA, m, n),

            #if _EVOLVE_DSIG

                gDsig_r  (m,n) = GF_right_r (gDsig, m, n),
                fDsig_r  (m,n) = GF_right_r (fDsig, m, n),

            #else

                gsig_r   (m,n) = GF_right_r (gsig, m, n),
                fsig_r   (m,n) = GF_right_r (fsig, m, n),
                gsig_rr  (m,n) = GF_right_rr (gsig, m, n),
                fsig_rr  (m,n) = GF_right_rr (fsig, m, n),

            #endif // _EVOLVE_DSIG

        // The second radial derivatives of the regularizing functions

            pr_rr      (m,n) = GF_right_rr (pr, m, n);
            qr_rr      (m,n) = GF_right_rr (qr, m, n);

        // The terms involving the Ricci tensors in the evolution equations
            gRicci     (m,n) = eq_gRicci    (m,n);
            fRicci     (m,n) = eq_fRicci    (m,n);

        // The recurrent terms involving some radial derivatives
            gDers      (m,n) = eq_gDers    (m,n);
            fDers      (m,n) = eq_fDers    (m,n);
        }

        OMP_parallel_for( Int n = nGhost+1; n < nGhost + nLen + 1; ++n )
        {

        // The radial derivatives of the evolved fields
            gconf_r  (m,n) = GF_r (gconf, m, n),     fconf_r  (m,n) = GF_r (fconf, m, n),
            gDconf_r (m,n) = GF_r (gDconf, m, n),    fDconf_r (m,n) = GF_r (fDconf, m, n),
            //gtrK_r   (m,n) = GF_r (gtrK, m, n),      ftrK_r   (m,n) = GF_r (ftrK, m, n),
            gA_r     (m,n) = GF_r (gA, m, n),        fA_r     (m,n) = GF_r (fA, m, n),
            gB_r     (m,n) = GF_r (gB, m, n),        fB_r     (m,n) = GF_r (fB, m, n),
            gDA_r    (m,n) = GF_r (gDA, m, n),       fDA_r    (m,n) = GF_r (fDA, m, n),
            gDB_r    (m,n) = GF_r (gDB, m, n),       fDB_r    (m,n) = GF_r (fDB, m, n),
            //gA2_r    (m,n) = GF_r (gA2, m, n),       fA1_r    (m,n) = GF_r (fA1, m, n),
            gA1_r    (m,n) = GF_r (gA1, m, n),       fA2_r    (m,n) = GF_r (fA2, m, n),
            gL_r     (m,n) = GF_r (gL, m, n),        fL_r     (m,n) = GF_r (fL, m, n),
            gAsig_r  (m,n) = GF_r (gAsig, m, n),     fAsig_r  (m,n) = GF_r (fAsig, m, n),

            //Lt_r     (m,n) = GF_r (Lt, m, n),

            // The radial derivatives of the regularizing functions

            gDconfr_r(m,n) = GF_r (gDconfr, m, n),
            gDAlpr_r (m,n) = GF_r (gDAlpr, m, n),
            gLr_r    (m,n) = GF_r (gLr, m, n),
            fDconfr_r(m,n) = GF_r (fDconfr, m, n),
            fDAlpr_r (m,n) = GF_r (fDAlpr, m, n),
            fLr_r    (m,n) = GF_r (fLr, m, n);

            pr_r      (m,n) = GF_r (pr, m, n);
            qr_r      (m,n) = GF_r (qr, m, n);

        // The radial derivatives of the Valencia variables and the external sources
            pfD_r    (m,n) = GF_r (pfD, m, n),
            pfS_r    (m,n) = GF_r (pfS, m, n),
            pftau_r  (m,n) = GF_r (pftau, m, n),
            gj_r     (m,n) = GF_r (gj, m, n),
            pfv_r    (m,n) = eq_pf_gv_r(m,n);

        // The radial derivatives of the utility function R = fB/gB
            R_r      (m,n) = GF_r (R, m, n),

        // These are the only spatial second derivatives that are needed:

        // The second radial derivatives of the evolved fields
            gconf_rr (m,n) = GF_rr (gconf, m, n),
            fconf_rr (m,n) = GF_rr (fconf, m, n),
            gA_rr    (m,n) = GF_rr (gA, m, n),
            gB_rr    (m,n) = GF_rr (gB, m, n),
            fA_rr    (m,n) = GF_rr (fA, m, n),
            fB_rr    (m,n) = GF_rr (fB, m, n),
            gA1_rr   (m,n) = GF_rr (gA1, m, n),
            fA1_rr   (m,n) = GF_rr (fA1, m, n),
            gtrK_rr  (m,n) = GF_rr (gtrK, m, n),
            ftrK_rr  (m,n) = GF_rr (ftrK, m, n),

        // The second radial derivatives of the gauge fields
            p_rr     (m,n) = GF_rr (p, m, n),   // used in eq_gDA_t and eq_fDA_t
            q_rr     (m,n) = GF_rr (q, m, n),   // used in eq_gDA_t and eq_fDA_t
            gAlp_rr  (m,n) = GF_rr (gAlp, m, n),
            fAlp_rr  (m,n) = GF_rr (fAlp, m, n),   // used in fDAlp_r

        // The second radial derivatives of the traces of the conformal extrinsic \
            curvatures (these traces are usually set to 0 in eomBSSNObserver.h)
            gtrA_rr  (m,n) = GF_rr (gtrA, m, n),
            ftrA_rr  (m,n) = GF_rr (ftrA, m, n),

            #if _EVOLVE_DSIG

                gDsig_r  (m,n) = GF_r (gDsig, m, n),
                fDsig_r  (m,n) = GF_r (fDsig, m, n),

            #else

                gsig_r   (m,n) = GF_r (gsig, m, n),
                fsig_r   (m,n) = GF_r (fsig, m, n),
                gsig_rr  (m,n) = GF_rr (gsig, m, n),
                fsig_rr  (m,n) = GF_rr (fsig, m, n),

            #endif // _EVOLVE_DSIG

        // The second radial derivatives of the regularizing functions

            pr_rr      (m,n) = GF_rr (pr, m, n);
            qr_rr      (m,n) = GF_rr (qr, m, n);

        // The terms involving the Ricci tensors in the evolution equations
            gRicci     (m,n) = eq_gRicci    (m,n);
            fRicci     (m,n) = eq_fRicci    (m,n);

        // The recurrent terms involving some radial derivatives
            gDers      (m,n) = eq_gDers    (m,n);
            fDers      (m,n) = eq_fDers    (m,n);
        }

        OMP_parallel_for( Int n = nGhost + nLen + 1; n < 2* nGhost + nLen + 1; ++n )
        {

        // The radial derivatives of the evolved fields
            gconf_r  (m,n) = GF_left_r (gconf, m, n),
            fconf_r  (m,n) = GF_left_r (fconf, m, n),
            gDconf_r (m,n) = GF_left_r (gDconf, m, n),
            fDconf_r (m,n) = GF_left_r (fDconf, m, n),
            //gtrK_r   (m,n) = GF_left_r (gtrK, m, n), ftrK_r   (m,n) = GF_r (ftrK, m, n),
            gA_r     (m,n) = GF_left_r (gA, m, n),
            fA_r     (m,n) = GF_left_r (fA, m, n),
            gB_r     (m,n) = GF_left_r (gB, m, n),
            fB_r     (m,n) = GF_left_r (fB, m, n),
            gDA_r    (m,n) = GF_left_r (gDA, m, n),
            fDA_r    (m,n) = GF_left_r (fDA, m, n),
            gDB_r    (m,n) = GF_left_r (gDB, m, n),
            fDB_r    (m,n) = GF_left_r (fDB, m, n),
            //gA2_r    (m,n) = GF_left_r (gA2, m, n),  fA1_r    (m,n) = GF_r (fA1, m, n),
            gA1_r    (m,n) = GF_left_r (gA1, m, n),
            fA2_r    (m,n) = GF_left_r (fA2, m, n),
            gL_r     (m,n) = GF_left_r (gL, m, n),
            fL_r     (m,n) = GF_left_r (fL, m, n),
            gAsig_r  (m,n) = GF_left_r (gAsig, m, n),
            fAsig_r  (m,n) = GF_left_r (fAsig, m, n),

            //Lt_r     (m,n) = GF_r (Lt, m, n),

            // The radial derivatives of the regularizing functions

            gDconfr_r(m,n) = GF_left_r (gDconfr, m, n),
            gDAlpr_r (m,n) = GF_left_r (gDAlpr, m, n),
            gLr_r    (m,n) = GF_left_r (gLr, m, n),
            fDconfr_r(m,n) = GF_left_r (fDconfr, m, n),
            fDAlpr_r (m,n) = GF_left_r (fDAlpr, m, n),
            fLr_r    (m,n) = GF_left_r (fLr, m, n);

            pr_r      (m,n) = GF_left_r (pr, m, n);
            qr_r      (m,n) = GF_left_r (qr, m, n);

        // The radial derivatives of the Valencia variables and the external sources
            pfD_r    (m,n) = GF_left_r (pfD, m, n),
            pfS_r    (m,n) = GF_left_r (pfS, m, n),
            pftau_r  (m,n) = GF_left_r (pftau, m, n),
            gj_r     (m,n) = GF_left_r (gj, m, n),
            pfv_r    (m,n) = eq_pf_gv_r(m,n);

        // The radial derivatives of the utility function R = fB/gB
            R_r      (m,n) = GF_left_r (R, m, n),

        // These are the only spatial second derivatives that are needed:

        // The second radial derivatives of the evolved fields
            gconf_rr (m,n) = GF_left_rr (gconf, m, n),
            fconf_rr (m,n) = GF_left_rr (fconf, m, n),
            gA_rr    (m,n) = GF_left_rr (gA, m, n),
            gB_rr    (m,n) = GF_left_rr (gB, m, n),
            fA_rr    (m,n) = GF_left_rr (fA, m, n),
            fB_rr    (m,n) = GF_left_rr (fB, m, n),
            gA1_rr   (m,n) = GF_left_rr (gA1, m, n),
            fA1_rr   (m,n) = GF_left_rr (fA1, m, n),
            gtrK_rr  (m,n) = GF_left_rr (gtrK, m, n),
            ftrK_rr  (m,n) = GF_left_rr (ftrK, m, n),

        // The second radial derivatives of the gauge fields
            p_rr     (m,n) = GF_left_rr (p, m, n),   // used in eq_gDA_t and eq_fDA_t
            q_rr     (m,n) = GF_left_rr (q, m, n),   // used in eq_gDA_t and eq_fDA_t
            gAlp_rr  (m,n) = GF_left_rr (gAlp, m, n),
            fAlp_rr  (m,n) = GF_left_rr (fAlp, m, n),   // used in fDAlp_r

        // The second radial derivatives of the traces of the conformal extrinsic curvatures\
            (these traces are usually set to 0 in eomBSSNObserver.h)
            gtrA_rr  (m,n) = GF_left_rr (gtrA, m, n),
            ftrA_rr  (m,n) = GF_left_rr (ftrA, m, n),

            #if _EVOLVE_DSIG

                gDsig_r  (m,n) = GF_left_r (gDsig, m, n),
                fDsig_r  (m,n) = GF_left_r (fDsig, m, n),

            #else

                gsig_r   (m,n) = GF_left_r (gsig, m, n),
                fsig_r   (m,n) = GF_left_r (fsig, m, n),
                gsig_rr  (m,n) = GF_left_rr (gsig, m, n),
                fsig_rr  (m,n) = GF_left_rr (fsig, m, n),

            #endif // _EVOLVE_DSIG

        // The second radial derivatives of the regularizing functions

            pr_rr      (m,n) = GF_left_rr (pr, m, n);
            qr_rr      (m,n) = GF_left_rr (qr, m, n);

        // The terms involving the Ricci tensors in the evolution equations
            gRicci     (m,n) = eq_gRicci    (m,n);
            fRicci     (m,n) = eq_fRicci    (m,n);

        // The recurrent terms involving some radial derivaties
            gDers      (m,n) = eq_gDers    (m,n);
            fDers      (m,n) = eq_fDers    (m,n);
        }

    /////////////////////////////////////////////////////////////////////////////////////
        /// - Calculate the intermediate variables that do not depend on derivatives.

        OMP_parallel_for( Int n = 0 + 0*nGhost; n <  2*nGhost + nLen + 1; ++n )
        {
            // The shifts and their radial derivatives
            gBet     (m,n) = eq_gBet    (m,n);
            gBet_r   (m,n) = eq_gBet_r  (m,n);
            gBet_rr  (m,n) = eq_gBet_rr (m,n);

            gBetr    (m,n) = eq_gBetr    (m,n);
            gBetr_r  (m,n) = eq_gBetr_r  (m,n);
            gBetr_rr (m,n) = eq_gBetr_rr (m,n);

            if( isGR() )
            {
                fBet     (m,n) = 0;
                fBet_r   (m,n) = 0;
                fBet_rr  (m,n) = 0;

                fBetr    (m,n) = 0;
                fBetr_r  (m,n) = 0;
                fBetr_rr (m,n) = 0;

            } else {

                fBet     (m,n) = eq_fBet    (m,n);
                fBet_r   (m,n) = eq_fBet_r  (m,n);
                fBet_rr  (m,n) = eq_fBet_rr (m,n);

                fBetr    (m,n) = eq_fBetr    (m,n);
                fBetr_r  (m,n) = eq_fBetr_r  (m,n);
                fBetr_rr (m,n) = eq_fBetr_rr (m,n);

            }
        }

        OMP_parallel_for( Int n = 0 + 0*nGhost; n <  nGhost + 1; ++n )
        {

        // The convective derivatives of the evolved fields
           gAlp_convr     (m,n) = gBet(m,n) * GF_right_r (gAlp, m ,n);
           gDAlp_convr    (m,n) = gBet(m,n) * GF_right_r (gDAlp, m ,n);
           fAlp_convr     (m,n) = fBet(m,n) * GF_right_r (fAlp, m ,n);
           //gBet_convr     (m,n) = GF_convr (gBet, gBet, m ,n);
           //fBet_convr     (m,n) = GF_convr (fBet, fBet, m ,n);
           q_qconvr       (m,n) = q(m,n)    * GF_right_r (q, m ,n);
           Bq_qconvr      (m,n) = q(m,n)    * GF_right_r (Bq, m ,n);

           gconf_convr    (m,n) = gBet(m,n) * GF_right_r (gconf, m ,n);
           fconf_convr    (m,n) = fBet(m,n) * GF_right_r (fconf, m ,n);
           gDconf_convr   (m,n) = gBet(m,n) * GF_right_r (gDconf, m ,n);
           fDconf_convr   (m,n) = fBet(m,n) * GF_right_r (fDconf, m ,n);
           gtrK_convr     (m,n) = gBet(m,n) * GF_right_r (gtrK, m ,n);
           ftrK_convr     (m,n) = fBet(m,n) * GF_right_r (ftrK, m ,n);

           gA_convr       (m,n) = gBet(m,n) * GF_right_r (gA, m ,n);
           fA_convr       (m,n) = fBet(m,n) * GF_right_r (fA, m ,n);
           gB_convr       (m,n) = gBet(m,n) * GF_right_r (gB, m ,n);
           fB_convr       (m,n) = fBet(m,n) * GF_right_r (fB, m ,n);
           gDA_convr      (m,n) = gBet(m,n) * GF_right_r (gDA, m ,n);
           fDA_convr      (m,n) = fBet(m,n) * GF_right_r (fDA, m ,n);
           gDB_convr      (m,n) = gBet(m,n) * GF_right_r (gDB, m ,n);
           fDB_convr      (m,n) = fBet(m,n) * GF_right_r (fDB, m ,n);

           gA1_convr      (m,n) = gBet(m,n) * GF_right_r (gA1, m ,n);
           fA1_convr      (m,n) = fBet(m,n) * GF_right_r (fA1, m ,n);
           gA2_convr      (m,n) = gBet(m,n) * GF_right_r (gA2, m ,n);
           fA2_convr      (m,n) = fBet(m,n) * GF_right_r (fA2, m ,n);

           gL_convr       (m,n) = gBet(m,n) * GF_right_r (gL, m ,n);
           gL_qconvr      (m,n) = q(m,n)    * GF_right_r (gL, m ,n);
           fL_convr       (m,n) = fBet(m,n) * GF_right_r (fL, m ,n);

           gsig_convr     (m,n) = gBet(m,n) * GF_right_r (gsig, m ,n);
           fsig_convr     (m,n) = fBet(m,n) * GF_right_r (fsig, m ,n);
           gAsig_convr    (m,n) = gBet(m,n) * GF_right_r (gAsig, m ,n);
           fAsig_convr    (m,n) = fBet(m,n) * GF_right_r (fAsig, m ,n);

           pfD_convr      (m,n) = gBet(m,n) * GF_right_r (pfD, m ,n);
           pfS_convr      (m,n) = gBet(m,n) * GF_right_r (pfS, m ,n);
           pftau_convr    (m,n) = gBet(m,n) * GF_right_r (pftau, m ,n);

           gBet_gDconf_r  (m,n) = gDconf(m,n) * eq_gBet_r (m, n);
           fBet_fDconf_r  (m,n) = fDconf(m,n) * eq_fBet_r (m, n);
           gBet_gDAlp_r   (m,n) = gDAlp(m,n)  * eq_gBet_r (m, n);
           gBet_gDA_r     (m,n) = gDA(m,n)    * eq_gBet_r (m, n);
           fBet_fDA_r     (m,n) = fDA(m,n)    * eq_fBet_r (m, n);
           gBet_gDB_r     (m,n) = gDB(m,n)    * eq_gBet_r (m, n);
           fBet_fDB_r     (m,n) = fDB(m,n)    * eq_fBet_r (m, n);
           gBet_gL_r      (m,n) = gL(m,n)     * eq_gBet_r (m, n);
           fBet_fL_r      (m,n) = fL(m,n)     * eq_fBet_r (m ,n);
           gsig_gL_r      (m,n) = gL(m,n)     * GF_right_r (gsig, m ,n);
           fsig_fL_r      (m,n) = fL(m,n)     * GF_right_r (fsig, m ,n);

           #if _EVOLVE_DSIG

             gDsig_convr   (m,n) = GF_convr (gBet, gDsig, m ,n);
             fDsig_convr   (m,n) = GF_convr (fBet, fDsig, m ,n);

           #endif // _EVOLVE_DSIG

        }

        OMP_parallel_for( Int n = nGhost+1; n <  nGhost + nLen + 1; ++n )
        {

        // The convective derivatives of the evolved fields
           gAlp_convr     (m,n) = GF_convr (gBet, gAlp, m ,n);
           gDAlp_convr    (m,n) = GF_convr (gBet, gDAlp, m ,n);
           fAlp_convr     (m,n) = GF_convr (fBet, fAlp, m ,n);
           //gBet_convr     (m,n) = GF_convr (gBet, gBet, m ,n);
           //fBet_convr     (m,n) = GF_convr (fBet, fBet, m ,n);
           q_qconvr       (m,n) = GF_convr (q,    q, m ,n);
           Bq_qconvr      (m,n) = GF_convr (q,    Bq, m ,n);

           gconf_convr    (m,n) = GF_convr (gBet, gconf, m ,n);
           fconf_convr    (m,n) = GF_convr (fBet, fconf, m ,n);
           gDconf_convr   (m,n) = GF_convr (gBet, gDconf, m ,n);
           fDconf_convr   (m,n) = GF_convr (fBet, fDconf, m ,n);
           gtrK_convr     (m,n) = GF_convr (gBet, gtrK, m ,n);
           ftrK_convr     (m,n) = GF_convr (fBet, ftrK, m ,n);

           gA_convr       (m,n) = GF_convr (gBet, gA, m ,n);
           fA_convr       (m,n) = GF_convr (fBet, fA, m ,n);
           gB_convr       (m,n) = GF_convr (gBet, gB, m ,n);
           fB_convr       (m,n) = GF_convr (fBet, fB, m ,n);
           gDA_convr      (m,n) = GF_convr (gBet, gDA, m ,n);
           fDA_convr      (m,n) = GF_convr (fBet, fDA, m ,n);
           gDB_convr      (m,n) = GF_convr (gBet, gDB, m ,n);
           fDB_convr      (m,n) = GF_convr (fBet, fDB, m ,n);

           gA1_convr      (m,n) = GF_convr (gBet, gA1, m ,n);
           fA1_convr      (m,n) = GF_convr (fBet, fA1, m ,n);
           gA2_convr      (m,n) = GF_convr (gBet, gA2, m ,n);
           fA2_convr      (m,n) = GF_convr (fBet, fA2, m ,n);

           gL_convr       (m,n) = GF_convr (gBet, gL, m ,n);
           gL_qconvr      (m,n) = GF_convr (q,    gL, m ,n);
           fL_convr       (m,n) = GF_convr (fBet, fL, m ,n);

           gsig_convr     (m,n) = GF_convr (gBet, gsig, m ,n);
           fsig_convr     (m,n) = GF_convr (fBet, fsig, m ,n);
           gAsig_convr    (m,n) = GF_convr (gBet, gAsig, m ,n);
           fAsig_convr    (m,n) = GF_convr (fBet, fAsig, m ,n);

           pfD_convr      (m,n) = GF_convr (gBet, pfD, m ,n);
           pfS_convr      (m,n) = GF_convr (gBet, pfS, m ,n);
           pftau_convr    (m,n) = GF_convr (gBet, pftau, m ,n);

           gBet_gDconf_r  (m,n) = gDconf(m,n) * eq_gBet_r (m, n);
           fBet_fDconf_r  (m,n) = fDconf(m,n) * eq_fBet_r (m, n);
           gBet_gDAlp_r   (m,n) = gDAlp(m,n)  * eq_gBet_r (m, n);
           gBet_gDA_r     (m,n) = gDA(m,n)    * eq_gBet_r (m, n);
           fBet_fDA_r     (m,n) = fDA(m,n)    * eq_fBet_r (m, n);
           gBet_gDB_r     (m,n) = gDB(m,n)    * eq_gBet_r (m, n);
           fBet_fDB_r     (m,n) = fDB(m,n)    * eq_fBet_r (m, n);
           gBet_gL_r      (m,n) = gL(m,n)     * eq_gBet_r (m, n);
           fBet_fL_r      (m,n) = fL(m,n)     * eq_fBet_r (m, n);
           gsig_gL_r      (m,n) = GF_convr (gL, gsig, m ,n);
           fsig_fL_r      (m,n) = GF_convr (fL, fsig, m ,n);

           #if _EVOLVE_DSIG

             gDsig_convr   (m,n) = GF_convr (gBet, gDsig, m ,n);
             fDsig_convr   (m,n) = GF_convr (fBet, fDsig, m ,n);

           #endif // _EVOLVE_DSIG

        }

        OMP_parallel_for( Int n = nGhost + nLen + 1; n <  2* nGhost + nLen + 1; ++n )
        {

        // The convective derivatives of the evolved fields
           gAlp_convr     (m,n) = gBet(m,n) * GF_left_r (gAlp, m ,n);
           gDAlp_convr    (m,n) = gBet(m,n) * GF_left_r (gDAlp, m ,n);
           fAlp_convr     (m,n) = fBet(m,n) * GF_left_r (fAlp, m ,n);
           //gBet_convr     (m,n) = GF_convr (gBet, gBet, m ,n);
           //fBet_convr     (m,n) = GF_convr (fBet, fBet, m ,n);
           q_qconvr       (m,n) = q(m,n)    * GF_left_r (q, m ,n);
           Bq_qconvr      (m,n) = q(m,n)    * GF_left_r (Bq, m ,n);

           gconf_convr    (m,n) = gBet(m,n) * GF_left_r (gconf, m ,n);
           fconf_convr    (m,n) = fBet(m,n) * GF_left_r (fconf, m ,n);
           gDconf_convr   (m,n) = gBet(m,n) * GF_left_r (gDconf, m ,n);
           fDconf_convr   (m,n) = fBet(m,n) * GF_left_r (fDconf, m ,n);
           gtrK_convr     (m,n) = gBet(m,n) * GF_left_r (gtrK, m ,n);
           ftrK_convr     (m,n) = fBet(m,n) * GF_left_r (ftrK, m ,n);

           gA_convr       (m,n) = gBet(m,n) * GF_left_r (gA, m ,n);
           fA_convr       (m,n) = fBet(m,n) * GF_left_r (fA, m ,n);
           gB_convr       (m,n) = gBet(m,n) * GF_left_r (gB, m ,n);
           fB_convr       (m,n) = fBet(m,n) * GF_left_r (fB, m ,n);
           gDA_convr      (m,n) = gBet(m,n) * GF_left_r (gDA, m ,n);
           fDA_convr      (m,n) = fBet(m,n) * GF_left_r (fDA, m ,n);
           gDB_convr      (m,n) = gBet(m,n) * GF_left_r (gDB, m ,n);
           fDB_convr      (m,n) = fBet(m,n) * GF_left_r (fDB, m ,n);

           gA1_convr      (m,n) = gBet(m,n) * GF_left_r (gA1, m ,n);
           fA1_convr      (m,n) = fBet(m,n) * GF_left_r (fA1, m ,n);
           gA2_convr      (m,n) = gBet(m,n) * GF_left_r (gA2, m ,n);
           fA2_convr      (m,n) = fBet(m,n) * GF_left_r (fA2, m ,n);

           gL_convr       (m,n) = gBet(m,n) * GF_left_r (gL, m ,n);
           gL_qconvr      (m,n) = q(m,n)    * GF_left_r (gL, m ,n);
           fL_convr       (m,n) = fBet(m,n) * GF_left_r (fL, m ,n);

           gsig_convr     (m,n) = gBet(m,n) * GF_left_r (gsig, m ,n);
           fsig_convr     (m,n) = fBet(m,n) * GF_left_r (fsig, m ,n);
           gAsig_convr    (m,n) = gBet(m,n) * GF_left_r (gAsig, m ,n);
           fAsig_convr    (m,n) = fBet(m,n) * GF_left_r (fAsig, m ,n);

           pfD_convr      (m,n) = gBet(m,n) * GF_left_r (pfD, m ,n);
           pfS_convr      (m,n) = gBet(m,n) * GF_left_r (pfS, m ,n);
           pftau_convr    (m,n) = gBet(m,n) * GF_left_r (pftau, m ,n);

           gBet_gDconf_r  (m,n) = gDconf(m,n) * eq_gBet_r (m, n);
           fBet_fDconf_r  (m,n) = fDconf(m,n) * eq_fBet_r (m, n);
           gBet_gDAlp_r   (m,n) = gDAlp(m,n)  * eq_gBet_r (m, n);
           gBet_gDA_r     (m,n) = gDA(m,n)    * eq_gBet_r (m, n);
           fBet_fDA_r     (m,n) = fDA(m,n)    * eq_fBet_r (m, n);
           gBet_gDB_r     (m,n) = gDB(m,n)    * eq_gBet_r (m, n);
           fBet_fDB_r     (m,n) = fDB(m,n)    * eq_fBet_r (m, n);
           gBet_gL_r      (m,n) = gL(m,n)     * eq_gBet_r (m, n);
           fBet_fL_r      (m,n) = fL(m,n)     * eq_fBet_r (m, n);
           gsig_gL_r      (m,n) = gL(m,n)     * GF_left_r (gsig, m ,n);
           fsig_fL_r      (m,n) = fL(m,n)     * GF_left_r (fsig, m ,n);

           #if _EVOLVE_DSIG

             gDsig_convr   (m,n) = GF_convr (gBet, gDsig, m ,n);
             fDsig_convr   (m,n) = GF_convr (fBet, fDsig, m ,n);

           #endif // _EVOLVE_DSIG

        }

        OMP_parallel_for( Int n = 0; n <  2* nGhost + nLen + 1; ++n )
        {

            // The sources (dependent on both p and the primary dynamical fields)
            //
            grhobar(m,n)= eq_pf_grhobar(m,n);
            grho (m,n)  = eq_pf_grho(m,n);    frho (m,n)  = eq_pf_frho(m,n);
            gj   (m,n)  = eq_pf_gj  (m,n);    fj   (m,n)  = eq_pf_fj  (m,n);
            gJ11 (m,n)  = eq_pf_gJ11(m,n);    fJ11 (m,n)  = eq_pf_fJ11(m,n);
            gJ22 (m,n)  = eq_pf_gJ22(m,n);    fJ22 (m,n)  = eq_pf_fJ22(m,n);

            gJK  (m,n)  = eq_gJK    (m,n);        fJK  (m,n) = eq_fJK  (m,n);
            gJA1 (m,n)  = eq_gJA1   (m,n);        fJA1 (m,n) = eq_fJA1 (m,n);
            gJA2 (m,n)  = - eq_gJA1   (m,n)/2;    fJA2 (m,n) = - eq_fJA1 (m,n)/2;

            gJL  (m,n)  = eq_gJL    (m,n);    fJL  (m,n) = eq_fJL  (m,n);

            // The bimetric sources
            //
            grhob  (m,n) = eq_grhob  (m,n);   frhob  (m,n) = eq_frhob  (m,n);
            gjb_u  (m,n) = eq_gjb_u  (m,n);   fjb_u  (m,n) = eq_fjb_u  (m,n);
            gJb1_ud(m,n) = eq_gJb1_ud(m,n);   fJb1_ud(m,n) = eq_fJb1_ud(m,n);
            gJb2_ud(m,n) = eq_gJb2_ud(m,n);   fJb2_ud(m,n) = eq_fJb2_ud(m,n);

        }

        if( m < mSmoothUpTo && smooth >= 1 )
        {
        //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::grho,   \
                       fld::tmp,  fld::grho,     1 );
        //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::grhobar, \
                       fld::tmp,  fld::grhobar,     1 );
        //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::pfD,    \
                       fld::tmp,  fld::pfD,     1 );
        //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gsig,   \
                       fld::tmp,  fld::gsig,     1 );
        //smoothenGF0 ( m, nSmoothFrom, nSmoothUpTo, sgRadius,  fld::gAsig,   \
                       fld::tmp,  fld::gAsig,    1 );

    }

    #endif // _TEST_MODE

    /////////////////////////////////////////////////////////////////////////////////////
    /// - Calculate the variables that depend on the spatial derivatives

    OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
    {
        #if _TEST_MODE

            u_t     (m,n) = eq_u_t  (m,n);
            v_t     (m,n) = eq_v_t  (m,n);

            std::cout << "\n\n m = " << m;

            std::cout << "\n\n u_t = " << u_t(m,n);
            std::cout << "\n\n v_t = " << v_t(m,n);
            std::cout << "\n\n eq_u_t = " << eq_u_t(m,n);
            std::cout << "\n\n eq_v_t = " << eq_v_t(m,n);

        #else

            // The time evolution right-hand side
            //
            p_t      (m,n) = isGR() ? 0 : eq_p_t(m,n);

            gconf_t  (m,n) = eq_gconf_t  (m,n);
            gDconf_t (m,n) = eq_gDconf_t (m,n);
            gtrK_t   (m,n) = eq_gtrK_t   (m,n);

            gA_t     (m,n) = eq_gA_t     (m,n);
            gDA_t    (m,n) = eq_gDA_t    (m,n);
            gB_t     (m,n) = eq_gB_t     (m,n);
            gDB_t    (m,n) = eq_gDB_t    (m,n);

            gA1_t    (m,n) = eq_gA1_t    (m,n);

            #if OBSERVER == 1

            #else

                gA2_t    (m,n) = eq_gA2_t    (m,n);

            #endif // OBSERVER

            gL_t     (m,n) = eq_gL_t     (m,n);

            gsig_t   (m,n) = eq_gsig_t   (m,n);
            gAsig_t  (m,n) = eq_gAsig_t  (m,n);

            pfD_t    (m,n) = eq_pf_gD_t  (m,n);
            pfS_t    (m,n) = eq_pf_gS_t  (m,n);
            pftau_t  (m,n) = eq_pf_gtau_t(m,n);

            if( isGR() )
            {
                fconf_t  (m,n) = 0;
                fDconf_t (m,n) = 0;
                ftrK_t   (m,n) = 0;

                fA_t     (m,n) = 0;
                fDA_t    (m,n) = 0;
                fB_t     (m,n) = 0;
                fDB_t    (m,n) = 0;

                fA1_t    (m,n) = 0;
                fA2_t    (m,n) = 0;

                fL_t     (m,n) = 0;

                fsig_t   (m,n) = 0;
                fAsig_t  (m,n) = 0;

            } else {

                fconf_t  (m,n) = eq_fconf_t  (m,n);
                fDconf_t (m,n) = eq_fDconf_t (m,n);
                ftrK_t   (m,n) = eq_ftrK_t   (m,n);

                fA_t     (m,n) = eq_fA_t     (m,n);
                fDA_t    (m,n) = eq_fDA_t    (m,n);
                fB_t     (m,n) = eq_fB_t     (m,n);
                fDB_t    (m,n) = eq_fDB_t    (m,n);

                fA1_t    (m,n) = eq_fA1_t    (m,n);

                #if OBSERVER == 1

                #else

                    fA2_t    (m,n) = eq_fA2_t    (m,n);

                #endif // OBSERVER

                fL_t     (m,n) = eq_fL_t     (m,n);

                fsig_t   (m,n) = eq_fsig_t   (m,n);
                fAsig_t  (m,n) = eq_fAsig_t  (m,n);

            }

            if( slicing == SLICE_MS2 || slicing == SLICE_MS2OPT || slicing == SLICE_MS4
               || slicing == SLICE_MS4D || slicing == SLICE_MS6 || slicing == SLICE_MS4G )
            {
                q_t     (m,n) = eq_SG_gBet_t  (m,n);
                Bq_t    (m,n) = eq_SG_gBq_t   (m,n);
            }

            if( slicing == SLICE_SG )
            {
                gAlp_t  (m,n) = eq_SG_gAlp_t  (m,n);
                gDAlp_t (m,n) = eq_SG_gDAlp_t (m,n);
                q_t     (m,n) = eq_SG_gBet_t  (m,n);
                Bq_t    (m,n) = eq_SG_gBq_t   (m,n);

            } else if( slicing == SLICE_KD )
            {
                gAlp_t  (m,n) = eq_KD_gAlp_t  (m,n);
                //gDAlp_t (m,n) = eq_SG_gDAlp_t (m,n);
            }

            #if _EVOLVE_DSIG

                gDsig_t  (m,n) = eq_gDsig_t  (m,n);
                fDsig_t  (m,n) = eq_fDsig_t  (m,n);

            #endif // _EVOLVE_DSIG

           /*Real dbg1 = eq_gconf_t  (m,n);
           Real dbg2 = eq_gDconf_t (m,n);

           Real dbg3 = eq_gtrK_t   (m,n);

           Real dbg5 = eq_gA_t     (m,n);
           Real dbg6 = eq_gDA_t     (m,n);
           Real dbg7 = eq_gB_t     (m,n);
           Real dbg8 = eq_gDB_t     (m,n);

           Real dbg9 = eq_gA1_t    (m,n);

           Real dbg12 = eq_gL_t     (m,n);

           Real dbg13 = eq_gsig_t   (m,n);
           Real dbg14 = eq_gAsig_t  (m,n);

           Real dbg15 = eq_pf_gD_t  (m,n);
           Real dbg16 = eq_pf_gS_t  (m,n);
           Real dbg17 = eq_pf_gtau_t(m,n);

           Real dbg18 = eq_KD_gAlp_t(m,n);
           Real dbg19 = gDAlp        (m,n);*/

           //Real dbg1 = gAlp  (m,n);
           //Real dbg2 = gDAlp (m,n);

           /*Real dbg1 = fL_convr(m,n);
           Real dbg3 = fBet_fL_r(m,n);
           Real dbg2 = k_f * fJL(m,n);
           Real dbg4 = fL(m,n);
           Real dbg5 = eq_fBet_r (m ,n);*/

           /*Real dbg1 = gAlp_convr     (m,n);
           Real dbg2 = fAlp_convr     (m,n);
           Real dbg3 = q_qconvr        (m,n);
           Real dbg5 = Bq_qconvr       (m,n);

           Real dbg6 = gconf_convr    (m,n);
           Real dbg7 = fconf_convr    (m,n);
           Real dbg8 = gDconf_convr   (m,n);
           Real dbg9 = fDconf_convr   (m,n);
           Real dbg10 = gtrK_convr     (m,n);
           Real dbg11 = ftrK_convr     (m,n);

           Real dbg12 = gA_convr       (m,n);
           Real dbg13 =  fA_convr       (m,n);
           Real dbg14 = gB_convr       (m,n);
           Real dbg15 = fB_convr       (m,n);
           Real dbg16 = gDA_convr      (m,n);
           Real dbg17 = fDA_convr      (m,n);
           Real dbg18 = gDB_convr      (m,n);
           Real dbg19 = fDB_convr      (m,n);

           Real dbg20 = gA1_convr      (m,n);
           Real dbg21 = fA1_convr      (m,n);
           Real dbg22 = gA2_convr      (m,n);
           Real dbg23 = fA2_convr      (m,n);

           Real dbg24 = gL_convr       (m,n);
           Real dbg25 = fL_convr       (m,n);

           Real dbg26 = gsig_convr     (m,n);
           Real dbg27 = fsig_convr     (m,n);
           Real dbg28 = gAsig_convr    (m,n);
           Real dbg29 = fAsig_convr    (m,n);*/

           /*Real dbg30 = gAlp    (m,n);
           Real dbg31 = gDAlp    (m,n);
           Real dbg32 = gDAlpr    (m,n);
           Real dbg33 = gAlp_r    (m,n);
           Real dbg34 = gDAlp_r    (m,n);
           Real dbg35 = gDAlpr_r    (m,n);

           Real dbg36 = PP    (m,n);
           Real dbg37 = QQ    (m,n);*/

           /*Real dbg51 = gRicci(m,n);
            Real dbg57 = gA1_convr(m,n);
            Real dbg58 = k_g * gJA1(m,n);

            Real dbg55 = exp(-4 * gconf(m,n));
            Real dbg59 = gDers(m,n);
            Real dbg60 = eq_gDers(m,n);
            Real dbg61 = (8
          * pow2(gDconf(m,n))) / (3. * pow2(gA(m,n)) + TINY_Real);

            Real dbg52 = (2 * gDA(m,n) * gDAlp(m,n)) / (3. * pow2(gA(m,n)) + TINY_Real);
            Real dbg62 = (2 * gDB(m,n)
         * gDAlp(m,n)) / (3. * pow2(gA(m,n)) + TINY_Real);
            Real dbg63 = (8 * gDconf(m,n) * gDAlp(m,n)) / (3.
        * pow2(gA(m,n)) + TINY_Real);*/

          /* Real dbg31 = gA1_convr(m,n) + k_g * gJA1(m,n) + exp(-4 * gconf(m,n))
         * (gDers(m,n) + gAlp(m,n) * ((2 * gDconf(m,n) * gDers(m,n))
         / (gDAlp(m,n)+TINY_Real) - (8
         * pow2(gDconf(m,n))) / (3. * pow2(gA(m,n)))));

           Real dbg32 = exp(-4 * gconf(m,n)) * ((-2
          * gDAlpr_r(m,n)) / (3. * pow2(gA(m,n))) - (4 * gDconfr_r(m,n) * gAlp(m,n))
        / (3. * pow2(gA(m,n)))) * r(m,n) + gtrA_pff(m,n) / 3. + (gAlp(m,n) * (3
        * gA1(m,n) - gtrA(m,n)) * (gtrA(m,n) + gtrK(m,n))) / 3.;

           Real dbg1 = (4 * gAsig(m,n) * gsig(m,n)
          * pow3(r(m,n))) / (3. * pow2(gA(m,n)) + TINY_Real);

           Real dbg3 = (4 * (gtrA_r(m,n) + gtrK_r(m,n)))
          / (3. * pow2(gA(m,n)) + TINY_Real);

           Real dbg2 =
         (4 * gDA(m,n) * pow2(r(m,n))) / (3. * pow2(gA(m,n)) + TINY_Real) + (4
       * gDB(m,n) * pow2(r(m,n))) / (3. * pow2(gA(m,n)) + TINY_Real);

           Real dgb4 = (8 * gconf(m,n)
       * gDconf(m,n) * pow2(r(m,n))) / (pow2(gA(m,n)) + TINY_Real) - (4 * gDAlp(m,n)
       * pow2(r(m,n))) / (3. * pow2(gA(m,n)) + TINY_Real);

           Real dbg5 = gAsig(m,n);

           Real dbg6 = gAsig_convr(m,n) + 2 * gAsig(m,n) * gBetr(m,n)
            + exp(-4 * gconf(m,n))
         * gAlp(m,n) * (((gsig_rr(m,n) / 2. + pow2(gsig(m,n))) * pow2(gB(m,n)))
        / Power(gA(m,n),4) - (gsig_gL_r(m,n) * pow2(gB(m,n))) / (2. * pow2(gA(m,n))));

           Real dbg7 = (-3 * (5 * gDB(m,n) + 11 * gDconf(m,n)) * gDers(m,n))
        / (gDAlp(m,n) + TINY_Real);

           Real dbg9 = (40 * gDB(m,n) * gDconf(m,n) + 7 * pow2(gDB(m,n)) + 44
        * pow2(gDconf(m,n))) / pow2(gA(m,n)) + (27 * pow2(gDers(m,n))
      * pow2(gA(m,n))) / (4. * pow2(gDAlp(m,n)) + TINY_Real);

           Real dbg8 = gA_convr(m,n) + (gdet_pff(m,n) * gA(m,n)) / (6. * gdet(m,n));

           Real dbg11 = gBet_r(m,n)
        * gA(m,n) + (gAlp(m,n) * gA(m,n) * (-3 * gA1(m,n) + gtrA(m,n)
      + gtrAv(m,n))) / 3.;*/

    #endif // _TEST_MODE

    }

    #if ! _TEST_MODE

    /*for( Int n = nGhost + nLen; n < nTotal + 1; ++n )
    {
        extrapolate_R( fld::gAlp_t,  m, n );
    }*/

    smoothenGF0 ( m, nGhost, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::gAlp_t,
                 fld::tmp,  fld::gAlp_t,  1 );
    //smoothenGF0 ( m, 0, nLen + 2 *nGhost + 1, 32,  fld::gtrK_t,    \
                   fld::tmp,  fld::gtrK_t,    1 );
    smoothenGF0 ( m, 0, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::gtrK,
                 fld::tmp,  fld::gtrK,  1 );
    smoothenGF0 ( m, 0, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::grhobar,
                 fld::tmp,  fld::grhobar,  1 );
    smoothenGF0 ( m, 0, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::grho,
                 fld::tmp,  fld::grho,  1 );
    //smoothenGF0 ( m, 0, nLen/*nLen + 2 *nGhost + 1*/, 10,  fld::pftau, \
                   fld::tmp,  fld::pftau,  1 );
    smoothenGF0 ( m, 0, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 20,  fld::pfW,
                 fld::tmp,  fld::pfW,  1 );
    //smoothenGF0 ( m, 0, nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::pfS, \
                   fld::tmp,  fld::pfS,  -1 );
    //smoothenGF0 ( m, 0, nLen/*nLen + 2 *nGhost + 1*/, 20,  fld::gA1, \
                   fld::tmp,  fld::gA1,  1 );
    smoothenGF0 ( m, nGhost, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::gconf,
                 fld::tmp,  fld::gconf,  1 );
    smoothenGF0 ( m, nGhost, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::gDconf,
                 fld::tmp,  fld::gDconf,  -1 );
    smoothenGF0 ( m, nGhost, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::gA,
                 fld::tmp,  fld::gA,  1 );
    smoothenGF0 ( m, nGhost, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::gB,
                 fld::tmp,  fld::gB,  1 );
    smoothenGF0 ( m, nGhost, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::gDA,
                 fld::tmp,  fld::gDA,  -1 );
    smoothenGF0 ( m, nGhost, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::gDB,
                 fld::tmp,  fld::gDB,  -1 );
    //smoothenGF0 ( m, nGhost, nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::gBr, \
                   fld::tmp,  fld::gBr,  -1 );
    //smoothenGF ( m,  fld::gtrK,      fld::tmp,  fld::gtrK,      1 );
    smoothenGF0 ( m, nGhost, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::q,
                 fld::tmp,  fld::q,  -1 );
    //smoothenGF0 ( m, nGhost, nGhost + nLen/*nLen + 2 *nGhost + 1*/, 15,  fld::gAlp,\
                 fld::tmp,  fld::gAlp,  1 );

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

    #endif // _TEST_MODE
}

/////////////////////////////////////////////////////////////////////////////////////////
/// Post-evolution smoothing and calculation of diagnostics
///
void BimetricEvolve::integStep_Finalize( Int m1, Int m )
{
    #if ! _TEST_MODE

    OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
    {

        gAlp( m, 0 ) = gAlp ( m , 1 );

        if( slicing == SLICE_MS2 || slicing == SLICE_MS2OPT || slicing == SLICE_MS4
           || slicing == SLICE_MS4D || slicing == SLICE_MS6 )
        {
           gAlp    ( m1, n )  =  gAlp   ( m, n );
           gDAlp   ( m1, n )  =  gDAlp  ( m, n );
        }

    }
    #endif //_TEST_MODE
}

/////////////////////////////////////////////////////////////////////////////////////////
/// Find horizons, calculate constraints and similar. Also check for NaNs.
///
bool BimetricEvolve::integStep_Diagnostics( Int m, Int chkNaNs_nFrom, Int chkNaNs_nTo )
{
    #if ! _TEST_MODE

    OMP_parallel_for( Int n = nGhost; n < nGhost + nLen; ++n )
    {
        /// - Apparent horizon finder (dependent only on the primary dynamical fields)
        ///
        gHorz(m,n)  = eq_gHorz  (m,n);
        fHorz(m,n)  = eq_fHorz  (m,n);

        /// - The constraint violations
        ///
        gHC (m,n) = eq_gHC (m,n);
        fHC (m,n) = eq_fHC (m,n);
        gMC (m,n) = eq_gMC (m,n);
        fMC (m,n) = eq_fMC (m,n);
        gLC (m,n) = eq_gLC (m,n);
        fLC (m,n) = eq_fLC (m,n);
        CL  (m,n) = eq_CL  (m,n);

        /// - The ratio of the lapses
        ///
        gW (m,n) = eq_gW (m,n);
        fW (m,n) = eq_fW (m,n);
    }

    #endif // _TEST_MODE

    if ( chkNaNs_nTo < 0 ) { // CheckNaNs is disabled
        return true;
    }

    /// - Check for NaNs in gA in the given zone.
    ///
    Int nFrom = nGhost + std::min( nLen, std::max( Int(0), chkNaNs_nFrom - nOffset ) );
    Int nTo   = nGhost + std::min( nLen, std::max( Int(0), chkNaNs_nTo   - nOffset ) );

    //if( false )
    for( Int n = nFrom; n < nTo; ++n )
    {

        #if _TEST_MODE

            if( /*std::*/ISNAN( u( m, n ) ) ) {
                std::cerr << "*** Detected u NaN at t = " << t(m,n)
                          << ", r = " << r(m,n) << std::endl;
                return false;
            } else if ( ISNAN( v( m, n ) ) ) {
                    std::cerr << "*** Detected gB NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
            }

        #else

            if( /*std::*/ISNAN( gA( m, n ) ) ) {
                std::cerr << "*** Detected gA NaN at t = " << t(m,n)
                          << ", r = " << r(m,n) << std::endl;
                return false;
            }

            #if _DETECT_NAN /// Track more fields to have more information in the output

                /*for( auto e: evolvedGF ) {

                    if( std::isnan( GF( e.f, m, n ) ) )
                    {
                        std::cerr << "*** Detected " <<  << " NaN at t = " << t(m,n)
                                << ", r = " << r(m,n) << std::endl;
                        return false;
                    }

                }*/ ///@fixme how to print the field names? Use bimOutput?

                /// The fields in the g-sector

                else if( ISNAN( gB( m, n ) ) ) {
                    std::cerr << "*** Detected gB NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                } else if( ISNAN( gDA( m, n ) ) ) {
                    std::cerr << "*** Detected gDA NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                } else if( ISNAN( gDB( m, n ) ) ) {
                    std::cerr << "*** Detected gDB NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                } else if( ISNAN( gA1( m, n ) ) ) {
                    std::cerr << "*** Detected gA1 NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                } /*else if( ISNAN( gA2( m, n ) ) ) {
                    std::cerr << "*** Detected gA2 NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                }*/ else if( ISNAN( gconf( m, n ) ) ) {
                    std::cerr << "*** Detected gconf NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                } else if( ISNAN( gtrK( m, n ) ) ) {
                    std::cerr << "*** Detected gtrK NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                }

                /// The fields in the f-sector

                else if( ISNAN( fA( m, n ) ) ) {
                    std::cerr << "*** Detected fA NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                } else if( ISNAN( fB( m, n ) ) ) {
                    std::cerr << "*** Detected fB NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                } else if( ISNAN( fDA( m, n ) ) ) {
                    std::cerr << "*** Detected fDA NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                } else if( ISNAN( fDB( m, n ) ) ) {
                    std::cerr << "*** Detected fDB NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                } else if( ISNAN( fA1( m, n ) ) ) {
                    std::cerr << "*** Detected fA1 NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                }/* else if( ISNAN( fA2( m, n ) ) ) {
                    std::cerr << "*** Detected fA2 NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                }*/ else if( ISNAN( fconf( m, n ) ) ) {
                    std::cerr << "*** Detected fconf NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                } else if( ISNAN( ftrK( m, n ) ) ) {
                    std::cerr << "*** Detected ftrK NaN at t = " << t(m,n)
                            << ", r = " << r(m,n) << std::endl;
                    return false;
                }

            #endif // _DETECT_NAN

        #endif // _TEST_MODE
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////

/** computeNewtonIterationMatrix computes the Newton iteration matrix.
 *  TODO: this function can probably be optimized; for example it can run over the grid
 */
void BimetricEvolve::computeNewtonIterationMatrix(
        Int m, Int n, Int n_evolved,
        Int stage_i,
        const ButcherTable& BT,
        MatReal& NewItMat
    )
{
    #if _TEST_MODE

        #include "jacobian-test/DIRK_Jacobian_test.h"

    #else

        #include "jacobian-BSSN/DIRK_Jacobian_cBSSN.h"

    #endif // _TEST_MODE

    for( Int i = 0; i < n_evolved; ++i )
    {
        for( Int j = 0; j < n_evolved; ++j )
        {
            NewItMat[i][j] = ( i == j ) ? 1 : 0
                             - delta_t * BT.A[stage_i][stage_i] * Jacobian[i][j];
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
// Maximal slicing (implements a BVP solver)
/////////////////////////////////////////////////////////////////////////////////////////

#if ! _TEST_MODE

    #include "cBSSNmaximalSlice.h"

#endif _TEST_MODE

/////////////////////////////////////////////////////////////////////////////////////////
/// The main entry point of `bim-solver`.
///
int main( int argc, char* argv[] )
{

    //LUDecomposition::sanityCheck(); return 0;

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
