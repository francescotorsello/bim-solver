/**
 *  @file      macros.h
 *  @brief     Contains the definitions of the fields within the namespace 'fields'.
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

/** Namespace 'fields' contains the indexing of the fields and their derivatives,
    to be used in the other classes
  */
namespace fields
{
    /*/// 'fldCheby' contains all the fields that need to be expanded in Chebyshev series
    enum fldCheby { FIELDS };
    static const fldCheby flds[] = { FIELDS }; // this is needed to be able to iterate over the enumeration ( taken from https://stackoverflow.com/questions/261963/how-can-i-iterate-over-an-enum )

    /// 'derCheby' contains all the first radial derivatives that need to be expanded in Chebyshev series
    enum derCheby { DERS };
    static const derCheby ders[] = { DERS };

    /// 'derrCheby' contains all the second radial derivatives that need to be expanded in Chebyshev series
    enum derrCheby { DERRS };
    static const derrCheby derrs[] = { DERRS };

    /// 'gaugeVars' contains all the gauge variables
    enum gaugeVars { GAUGE };
    static const gaugeVars gauges[] = { GAUGE };

    /// 'gaugeDers' contains the first radial derivatives of the gauge variables
    enum gaugeDers { GAUGE_R };
    static const gaugeDers gauges_r[] = { GAUGE_R };

    /// 'gaugeDerrs' contains the second radial derivatives of the gauge variables
    enum gaugeDerrs { GAUGE_RR };
    static const gaugeDerrs gauges_rr[] = { GAUGE_RR };

    /// 'valenciaVars' contains all the gauge variables
    enum valenciaVars { VALENCIA };
    static const valenciaVars valencia[] = { VALENCIA };

    /// 'valenciaDers' contains the first radial derivatives of the gauge variables
    enum valenciaDers { VALENCIA_R };
    static const valenciaDers valencia_r[] = { VALENCIA_R };

    /// 'valenciaDerrs' contains the second radial derivatives of the gauge variables
    enum valenciaDerrs { VALENCIA_RR };
    static const valenciaDerrs valencia_rr[] = { VALENCIA_RR };*/

    /// TODO: instead of defining different enumerations, define different vectors over which you can run over.

    /// 'allFields' contains all the gauge variables
    enum allFields { ALL_FIELDS };
    static const allFields all_flds[]       = { ALL_FIELDS };
    static const allFields evolved_flds[]   = { EVOLVED_FIELDS };
    static const allFields even_flds[]      = { EVEN_FIELDS };
    static const allFields odd_flds[]       = { ODD_FIELDS };

    /// 'allDers' contains the first radial derivatives of the gauge variables
    enum allDers { ALL_DERS };
    static const allDers all_ders[]  = { ALL_DERS };
    static const allDers even_ders[] = { EVEN_DERS };
    static const allDers odd_ders[]  = { ODD_DERS };

    /// 'regDers' contains the regularizing radial derivatives
    enum regDers { REG_DERS };
    static const regDers reg_ders[] = { REG_DERS };

    /// 'allDerrs' contains the second radial derivatives of the gauge variables
    enum allDerrs { ALL_DERRS };
    static const allDerrs all_derrs[]   = { ALL_DERRS };
    static const allDerrs even_derrs[]  = { EVEN_DERRS };
    static const allDerrs odd_derrs[]   = { ODD_DERRS };

    /// perhaps the vector below is useless
    //static const std::vector<Int> bispecInput_fields = { FIELDS };

}
