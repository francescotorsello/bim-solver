/**
 *  @file      fields.h
 *  @brief     Defines the macro containing the name of the fields.
 *  @authors   Mikica Kocic, Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _INTEGRATOR_H_INCLUDED
#define _INTEGRATOR_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g?? Macro including the fields                                           */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/

#define FIELDS \
\
            /* Evolved fields in the `g`-sector */\
            gconf, gDconf, gtrK, gA, gDA, gB, gDB, gA1, gA2, gL, gsig, gAsig, \
\
            /* The evolution RHS for `g` */\
            gconf_t, gDconf_t, gtrK_t, gA_t, gDA_t, gB_t, gDB_t, gA1_t, gA2_t, gL_t,\
            gsig_t, gAsig_t,\
\
            /* Evolved fields in the `f`-sector */\
            fconf, fDconf, ftrK, fA, fDA, fB, fDB, fA1, fA2, fL, fsig, fAsig,\
\
            /* The evolution RHS for `f` */\
            fconf_t, fDconf_t, ftrK_t, fA_t, fDA_t, fB_t, fDB_t, fA1_t, fA2_t, fL_t,\
            fsig_t, fAsig_t,\
\
            /* The traces of the shifted extrinsic curvatures */\
            gtrA, ftrA,\
\
            /* The violations of the traces of the conformal extrinsic curvatures
              (only used when the previous ones are set to 0; they only appear in
              the evolution equations for the confomal metrics, see Brown) */\
            gtrAv, ftrAv,\
\
            /* The Pfaffian derivatives of the traces */\
            gtrA_pff, ftrA_pff,\
\
            /* The determinants of the conformal metrics */\
            gdet, fdet,\
\
            /* The Pfaffian derivatives of the determinants */\
            gdet_pff, fdet_pff,\
\
            /* The regularizing fields and their radial derivatives */\
            gBetr, gDconfr, fDconfr, gDAlpr, fBetr, fDAlpr, gLr, fLr,\
            gBetr_r, gBetr_rr, gDconfr_r, fDconfr_r, gDAlpr_r, fBetr_r,\
            fBetr_rr, fDAlpr_r, gLr_r, fLr_r, gBr, gBr_r,\
\
            /* State variables for the perfect fluid (PF) */\
            pfD, pfS, pftau, pfv,\
\
            /* The Lorentz factor for the PF */\
            pfW,\
\
            /* The Valencia evolution RHS for the PF */\
            pfD_t, pfS_t, pftau_t,\
\
            /* The radial derivative of r */\
            pfv_r,\
\
            /* Sum of sources in the `g`-sector */\
            gJK, gJA1, gJA2, gJL,\
\
            /* Sum of sources in the `f`-sector */\
            fJK, fJA1, fJA2, fJL,\
\
            /* Physical sources computed from the Valencia variables */\
            grhobar,\
            grho, frho,\
            gj,   fj,\
            gJ11, gJ22,\
            fJ11, fJ22,\
\
            /* Bimetric sources */\
            grhob,   gjb_u,   gJb1_ud,  gJb2_ud,\
            frhob,   fjb_u,   fJb1_ud,  fJb2_ud,\
\
            /* Separation (rel shift) between the metrics */\
            p, pr,\
\
            /* Lorentz factor (derived from `p`) */\
            Lt, Lt2,\
\
            /* RHS for the evolution of `p` */\
            p_t,\
\
            /* Shift of the mean metric, and auxiliary variable for standard gauge */\
            q, qr, q_t,\
            Bq, Bq_t, Bq_r,\
\
            /* Radial derivatives */\
            /* The radial derivatives of the evolved fields */\
            gconf_r,     fconf_r,\
            gDconf_r,    fDconf_r,\
            gtrK_r,      ftrK_r,\
            gA_r,        fA_r,\
            gB_r,        fB_r,\
            gDA_r,       fDA_r,\
            gDB_r,       fDB_r,\
            gA1_r,       fA1_r,\
            gA2_r,       fA2_r,\
            gL_r,        fL_r,\
\
            gAsig_r,     fAsig_r,\
\
        /* The radial derivatives of the gauge fields */\
            p_r,         q_r,\
            pr_r,        qr_r,\
            gAlp_r,      fAlp_r,\
            gDAlp_r,     fDAlp_r,\
            Lt_r,\
\
        /* The radial derivatives of the Valencia variables and the external sources */\
            pfD_r,\
            pfS_r,\
            pftau_r,\
            gj_r,\
\
        /* The radial derivatives of the utility function R = fB/gB */\
            R_r,\
\
            /* Convective derivatives\
               The convective derivatives of the evolved fields; they arise from the
               Pfaffians (i.e., from the Lie derivatives along the shifts) */\
            gAlp_convr,   fAlp_convr,\
            gDAlp_convr,\
            /*gBet_convr,   fBet_convr, */\
            q_qconvr,     Bq_qconvr,\
\
            gconf_convr,  fconf_convr,\
            gDconf_convr, fDconf_convr,\
            gtrK_convr,   ftrK_convr,\
\
            gA_convr,     fA_convr,\
            gB_convr,     fB_convr,\
            gDA_convr,    fDA_convr,\
            gDB_convr,    fDB_convr,\
\
            gA1_convr,    fA1_convr,\
            gA2_convr,    fA2_convr,\
\
            gL_convr,     fL_convr,\
            gL_qconvr,\
\
            gsig_convr,   fsig_convr,\
            gAsig_convr,  fAsig_convr,\
\
            pfD_convr,    pfS_convr,\
            pftau_convr,\
\
            gBet_gDconf_r,\
            fBet_fDconf_r,\
\
            gBet_gDAlp_r,\
\
            gBet_gDA_r,\
            fBet_fDA_r,\
\
            gBet_gDB_r,\
            fBet_fDB_r,\
\
            gBet_gL_r,\
            fBet_fL_r,\
\
            gsig_gL_r,\
            fsig_fL_r,\
\
            gDsig,\
            fDsig, /*@fixme   these fields must be declared here because otherwise they
                            cannot be read from the initial data. The present
                            implementation reads these two fields from the initial data
                            anyways. Then, if they are not evolved, they are not
                            touched anymore. This is because otherwise the Mathematica
                            file should be also changed. */\
\
            gsig_r,\
            fsig_r,\
            gsig_rr,\
            fsig_rr,\
\
        /* These are the only spatial second derivatives that are needed: */\
\
        /* The second radial derivatives of the evolved fields */\
            gconf_rr,\
            fconf_rr,\
            gA_rr,\
            gB_rr,\
            fA_rr,\
            fB_rr,\
            gA1_rr,\
            fA1_rr,\
            gtrK_rr,\
            ftrK_rr,\
\
        /* The second radial derivatives of the gauge fields */\
            p_rr,   /* used in eq_gDA_t and eq_fDA_t */\
            q_rr,   /* used in eq_gDA_t and eq_fDA_t */\
            pr_rr,\
            qr_rr,\
            gAlp_rr,\
            fAlp_rr,   /* used in fDAlp_r */\
\
        /* The second radial derivatives of the traces of the conformal extrinsic
            curvatures (these traces are usually set to 0 in eomBSSNObserver.h) */\
            gtrA_rr,\
            ftrA_rr,\
\
            /* Lapses */\
            gAlp, fAlp,\
            gDAlp, fDAlp,\
\
            /* RHS for the evolution of the lapse */\
            gAlp_t, gDAlp_t,\
\
            /* Lapse ratio (`gW * gAlp + fW * fAlp = 0`) */\
            gW, fW,\
\
            /* Shifts */\
            gBet, fBet,\
\
            /* Radial derivatives of the shifts */\
            gBet_r, fBet_r,\
            gBet_rr, fBet_rr,\
\
            /* The terms involving the Ricci scalars in the evolution equations */\
            gRicci, fRicci,\
\
            /* The recurrent terms involving some radial derivatives */\
            gDers, fDers,\
\
            /* Constraints (here error estimators) */\
            gHC, fHC,\
            gMC, fMC,\
            gLC, fLC,\
            CL,\
\
            /* Apparent horizon finders */\
            gHorz, fHorz,\
\
            /* Ratio `fB/gB` */\
            R,\
\
            mpiBoundary = R,       /* The last grid function for which BC must be fixed.
                                      Used in MPI when calling defineGhostChunk. */\

                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _INTEGRATOR_H_INCLUDED
