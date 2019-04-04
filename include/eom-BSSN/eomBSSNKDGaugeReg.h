/** @file  eomBSSNKDGaugeReg.h
 *  @author Francesco Torsello
 *  @brief The K-driver parabolic gauge condition on the lapse (p. 111 of B&S).
 *  @version 2019-04-04T13:50:01
 *  @image html BSSNKDGaugeReg.png
 */

Real BimetricEvolve::eq_KD_gAlp_t( Int m, Int n )
{
    return Kdiff*(gAlp_rr(m,n) - (k_g * exp(4 * gconf(m,n)) * gAlp(m,n) * (gJ11(m,n) + 2
    /*2*/  * gJ22(m,n) + grho(m,n)) * pow2(gA(m,n))) / 2. - (exp(4 * gconf(m,n))
    /*1*/  * gAlp(m,n) * (pow2(0) + 3 * pow2(gK1(m,n)) + 6 * pow2(gK2(m,n)))
    /*1*/  * pow2(gA(m,n))) / 3. + gAlp_r(m,n) * (2 * gconf_r(m,n) - gA_r(m,n) / ( gA(m,n) + TINY_Real )
    /*1*/  + (2 * gB_r(m,n)) / ( gB(m,n) + TINY_Real ) + 2 / r(m,n)) - Kelas * (gtrA(m,n) + gtrK(m,n)));
}
