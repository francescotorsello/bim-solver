/**
 *  @file      macros.h
 *  @brief     Contains the macros defined in bispec.
 *  @authors   Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

/**
 *  The following vector contains the ORDERED fields whose spectral coefficients
 *  have to be set equal to the values read by the class bispecInput.
 *  This order MUST match the order of exportation in the Mathematica file.
 */

#define EVOLVED_FIELDS  gconf, fconf, gtrK, ftrK, gA, fA, gB, fB, gA1, fA1, gL, fL, \
                        gsig, fsig, gAsig, fAsig, pfD, pfS, pftau

#define EVEN_FIELDS  gconf, fconf, gtrK, ftrK, gA, fA, gB, fB, gA1, fA1, \
                     gsig, fsig, gAsig, fAsig, pfD, pftau

#define ODD_FIELDS  gL, fL, pfS

#define EVEN_FIELDS_T   gconf_t, fconf_t, gtrK_t, ftrK_t, gA_t, fA_t, gB_t, \
                        fB_t, gA1_t, fA1_t, gsig_t, fsig_t, gAsig_t, fAsig_t, gD_t, gtau_t

#define ODD_FIELDS_T    gL_t, fL_t, gS_t

#define DERS    gconf_r, fconf_r, gtrK_r, ftrK_r, gA_r, fA_r, gB_r, fB_r, gA1_r, fA1_r, \
                gL_r, fL_r, gsig_r, fsig_r, gAsig_r, fAsig_r, pfD_r, pfS_r, pftau_r

#define REG_DERS    gBetr_r, fBetr_r, gLr_r, fLr_r, gderAlpr_r, fderAlpr_r, \
                    gderconfr_r, fderconfr_r

#define EVEN_DERS   gconf_r, fconf_r, gtrK_r, ftrK_r, gA_r, fA_r, gB_r, fB_r, gA1_r, \
                    fA1_r, gL_r, fL_r, gsig_r, fsig_r, gAsig_r, fAsig_r, pfD_r, pfS_r, \
                    pftau_r

#define ODD_DERS  gL_r, fL_r, pfS_r

#define DERRS   gconf_rr, fconf_rr, gtrK_rr, ftrK_rr, gA_rr, fA_rr, gB_rr, fB_rr, \
                gA1_rr, fA1_rr, gL_rr, fL_rr, gsig_rr, fsig_rr, gAsig_rr, fAsig_rr, \
                pfD_rr, pfS_rr, pftau_rr

#define EVEN_DERRS  gconf_rr, fconf_rr, gtrK_rr, ftrK_rr, gA_rr, fA_rr, gB_rr, fB_rr, \
                    gA1_rr, fA1_rr, gL_rr, fL_rr, gsig_rr, fsig_rr, gAsig_rr, fAsig_rr, \
                    pfD_rr, pfS_rr, pftau_rr

#define ODD_DERRS   gL_rr, fL_rr, pfS_rr

#define GAUGE       gAlp, fAlp, hAlp, gBet, fBet, hBet, p

#define GAUGE_R     gAlp_r, fAlp_r, hAlp_r, gBet_r, fBet_r, hBet_r, p_r

#define GAUGE_RR    gAlp_rr, fAlp_rr, hAlp_rr, gBet_rr, fBet_rr, hBet_rr, p_rr

#define VALENCIA    pfD, pfS, pftau

#define VALENCIA_R  pfD_r, pfS_r, pftau_r

#define VALENCIA_RR pfD_rr, pfS_rr, pftau_rr

#define ALL_FIELDS  gconf, fconf, gtrK, ftrK, gA, fA, gB, fB, gA1, fA1, gL, fL, \
                    gsig, fsig, gAsig, fAsig, gAlp, fAlp, hAlp, gBet, fBet, hBet, p, \
                    pfD, pfS, pftau

#define ALL_DERS    gconf_r, fconf_r, gtrK_r, ftrK_r, gA_r, fA_r, gB_r, fB_r, gA1_r, \
                    fA1_r, gL_r, fL_r, gsig_r, fsig_r, gAsig_r, fAsig_r, gAlp_r, fAlp_r, \
                    hAlp_r, gBet_r, fBet_r, hBet_r, p_r, pfD_r, pfS_r, pftau_r

#define ALL_DERRS   gconf_rr, fconf_rr, gtrK_rr, ftrK_rr, gA_rr, fA_rr, gB_rr, fB_rr, \
                    gA1_rr, fA1_rr, gL_rr, fL_rr, gsig_rr, fsig_rr, gAsig_rr, fAsig_rr, \
                    gAlp_rr, fAlp_rr, hAlp_rr, gBet_rr, fBet_rr, hBet_rr, p_rr, \
                    pfD_rr, pfS_rr, pftau_rr


/** Note that the macro 'setfield' below, at present, only works inside bispecEvolve.
    It defines the functions that are included in the evolution equations exported by
    Mathematica
  */

#define setfield( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_fields[ ( n_flds * n_collocs ) * m \
            + n_collocs * fields::field + n ]; \
        }

#define setder( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_ders[ ( n_flds * n_collocs ) * m \
            + n_collocs * fields::field + n ]; \
        }

#define setderr( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_derrs[ ( n_flds * n_collocs ) * m \
            + n_collocs * fields::field + n ]; \
        }

#define setregder( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_reg_ders[ ( 8 * ( exp_ord + 1 ) ) * m \
            + ( exp_ord + 1 ) * fields::field + n ]; \
        }
