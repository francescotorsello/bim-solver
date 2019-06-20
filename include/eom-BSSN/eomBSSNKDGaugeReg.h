/** @file  eomBSSNKDGaugeReg.h
 *  @author Francesco Torsello
 *  @brief The K-driver parabolic gauge condition on the lapse (p. 111 of B&S).
 *  @version 2019-04-05T08:46:33
 *  @image html BSSNKDGaugeReg.png
 */

Real BimetricEvolve::eq_KD_gAlp_t( Int m, Int n )
{
    return Kdiff * ((k_g * gAlp(m,n) * (-gJ11(m,n) - 2 * gJ22(m,n) - grho(m,n))) / 2.
    /*1*/  - (GF_convr (q, gtrK, m ,n) + GF_convr (q, gtrA, m ,n)) + gAlp(m,n) * (-pow2(gK1(m,n)) - 2 * pow2(gK2(m,n)))
    /*1*/  + (exp(-4 * gconf(m,n)) * (-(gDA(m,n) * gDAlp(m,n)) + 2 * gDB(m,n)
    /*3*/  * gDAlp(m,n) + 2 * gDconf(m,n) * gDAlp(m,n) + 3 * gDAlp_r(m,n)))
    /*1*/  / ( pow2(gA(m,n)) + TINY_Real ) - (2 * exp(-4 * gconf(m,n)) * gDAlpr_r(m,n) * r(m,n))
    /*1*/  / ( pow2(gA(m,n)) + TINY_Real ) - Kelas * ( gtrK(m,n) + gtrA(m,n) ) );
}
