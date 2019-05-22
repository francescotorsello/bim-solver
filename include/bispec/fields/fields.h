/**
 *  @file      fields.h
 *  @brief     Contains the definitions of the fields within the namespace 'fields'.
 *  @authors   Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

/**
 *  The following macros contains the ORDERED fields whose spectral coefficients
 *  have to be set equal to the values read by the class BispecInput.
 *  This order MUST match the order of exportation in the Mathematica file.
 */

//////////////////////////////////////////////////////////////////////////////////////////
///// Fields
//////////////////////////////////////////////////////////////////////////////////////////

#define ALL_FIELDS  gconf, fconf, gtrK, ftrK, gA, fA, gB, fB, gA1, fA1, gL, fL, \
                    gsig, fsig, gAsig, fAsig, gAlp, fAlp, hAlp, gBet, fBet, hBet, p, \
                    pfD, pfS, pftau

#define EVEN_FIELDS gconf, fconf, gtrK, ftrK, gA, fA, gB, fB, gA1, fA1, \
                    gsig, fsig, gAsig, fAsig, pfD, pftau

#define ODD_FIELDS  gL, fL, p, pfS

#define LAPSES      gAlp, fAlp, hAlp

#define SHIFTS      gBet, fBet, hBet

//////////////////////////////////////////////////////////////////////////////////////////
///// First derivatives
//////////////////////////////////////////////////////////////////////////////////////////

#define ALL_DERS    gconf_r, fconf_r, gtrK_r, ftrK_r, gA_r, fA_r, gB_r, fB_r, gA1_r, \
                    fA1_r, gL_r, fL_r, gsig_r, fsig_r, gAsig_r, fAsig_r, gAlp_r, fAlp_r, \
                    hAlp_r, gBet_r, fBet_r, hBet_r, p_r, pfD_r, pfS_r, pftau_r

#define EVEN_DERS   gconf_r, fconf_r, gtrK_r, ftrK_r, gA_r, fA_r, gB_r, fB_r, gA1_r, \
                    fA1_r, gL_r, fL_r, gsig_r, fsig_r, gAsig_r, fAsig_r, pfD_r, pftau_r

#define ODD_DERS    gL_r, fL_r, pfS_r

#define LAPSES_R    gAlp_r, fAlp_r, hAlp_r

#define SHIFTS_R    gBet_r, fBet_r, hBet_r

//////////////////////////////////////////////////////////////////////////////////////////
///// Second derivatives
//////////////////////////////////////////////////////////////////////////////////////////

#define ALL_DERRS   gconf_rr, fconf_rr, gtrK_rr, ftrK_rr, gA_rr, fA_rr, gB_rr, fB_rr, \
                    gA1_rr, fA1_rr, gL_rr, fL_rr, gsig_rr, fsig_rr, gAsig_rr, fAsig_rr, \
                    gAlp_rr, fAlp_rr, hAlp_rr, gBet_rr, fBet_rr, hBet_rr, p_rr, \
                    pfD_rr, pfS_rr, pftau_rr

#define EVEN_DERRS  gconf_rr, fconf_rr, gtrK_rr, ftrK_rr, gA_rr, fA_rr, gB_rr, fB_rr, \
                    gA1_rr, fA1_rr, gL_rr, fL_rr, gsig_rr, fsig_rr, gAsig_rr, fAsig_rr, \
                    pfD_rr, pfS_rr, pftau_rr

#define ODD_DERRS   gL_rr, fL_rr, pfS_rr

#define LAPSES_RR   gAlp_rr, fAlp_rr, hAlp_rr

#define SHIFTS_RR   gBet_rr, fBet_rr, hBet_rr

//////////////////////////////////////////////////////////////////////////////////////////
///// Regularized derivatives
//////////////////////////////////////////////////////////////////////////////////////////

#define ALL_REG_DERS            gBetr_r, fBetr_r, gLr_r, fLr_r, gderAlpr_r, fderAlpr_r, \
                                gderconfr_r, fderconfr_r

#define REG_DERS                gBetr_r, fBetr_r, gLr_r, fLr_r

#define REG_DERRS               gderAlpr_r, fderAlpr_r, gderconfr_r, fderconfr_r

#define FIELDS_REG_DERS         gBet, fBet, gL, fL, gAlp, fAlp, gconf, fconf

#define PARITY_FIELDS_REG_DERS  false, false, false, false, true, true, true, true
                    //////////////////////////////////////////////////////////////////////////////////////////
///// Time derivatives of the fields
//////////////////////////////////////////////////////////////////////////////////////////

#define EVEN_FIELDS_T   gconf_t, fconf_t, gtrK_t, ftrK_t, gA_t, fA_t, gB_t, fB_t, \
                        gA1_t, fA1_t, gsig_t, fsig_t, gAsig_t, fAsig_t, gD_t, gtau_t

#define ODD_FIELDS_T    gL_t, fL_t, gS_t

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/**
 *  Namespace 'fields' contains the indexing of the fields and their derivatives,
 *  to be used in the other classes
 */

namespace fields
{

    /// 'allFields' contains all the fields and gauge variables

    enum allFields { ALL_FIELDS };

    static const allFields all_flds []   = { ALL_FIELDS };
    static const allFields even_flds[]   = { EVEN_FIELDS };
    static const allFields odd_flds []   = { ODD_FIELDS };
    static const allFields lapses   []   = { LAPSES };
    static const allFields shifts   []   = { SHIFTS };

    //static const std::vector<allFields> foo = { EVEN_FIELDS };
    //vector<allFields>::iterator itr;


    /// 'allDers' contains the first radial derivatives of the fields and gauge variables

    enum allDers { ALL_DERS };

    static const allDers all_ders   []   = { ALL_DERS };
    static const allDers even_ders  []   = { EVEN_DERS };
    static const allDers odd_ders   []   = { ODD_DERS };
    static const allDers lapses_r   []   = { LAPSES_R };
    static const allDers shifts_r   []   = { SHIFTS_R };

    /// 'allDerrs' contains the second radial ders of the fields and gauge variables

    enum allDerrs { ALL_DERRS };

    static const allDerrs all_derrs []   = { ALL_DERRS };
    static const allDerrs even_derrs[]   = { EVEN_DERRS };
    static const allDerrs odd_derrs []   = { ODD_DERRS };
    static const allDerrs lapses_rr []   = { LAPSES_RR };
    static const allDerrs shifts_rr []   = { SHIFTS_RR };

    /// 'regDers' contains the regularized radial derivatives

    enum regDers { ALL_REG_DERS };

    static const regDers   reg_ders     []   = { REG_DERS };
    static const regDers   reg_derrs    []   = { REG_DERRS };
    static const allFields flds_reg_ders[]   = { FIELDS_REG_DERS };
    static const std::vector<bool> parity_fields_reg_ders   = { PARITY_FIELDS_REG_DERS };

}
