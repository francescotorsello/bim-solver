/** @file  eomBSSNKDGaugeComp.h
 *  @author Francesco Torsello
 *  @brief The compactified K-driver parabolic gauge condition on the lapse (p. 111 of B&S).
 *  @version 2019-05-14T16:26:59
 *  @image html BSSNKDGaugeComp.png
 */

Real KD_gAlp_t( Int m, Int n )
{
    return Kdiff * (4 * Power(r_minus(m,n),4) * gAlp_rr(m,n) - (exp(4 * gconf(m,n)) 
    /*2*/  * gAlp(m,n) * (pow2(0) + 3 * pow2(gK1(m,n)) + 6 * pow2(gK2(m,n))) 
    /*2*/  * pow2(gA(m,n))) / 3. + k_g * (-(exp(2 * (fconf(m,n) + gconf(m,n))) 
    /*3*/  * fAlp(m,n) * fA(m,n) * gA(m,n) * P_1_2(R(m,n))) - (exp(4 * gconf(m,n)) 
    /*3*/  * fAlp(m,n) * pow2(gA(m,n)) * (2 * P_1_1(R(m,n)) + P_2_1(R(m,n)))) / (2. 
    /*3*/  * Lt(m,n)) + gAlp(m,n) * (-(exp(4 * gconf(m,n)) * pow2(gA(m,n)) * (gJ11(m,n)
    /*5*/  + 2 * gJ22(m,n) + grho(m,n) + 2 * P_1_0(R(m,n)))) / 2. - (exp(2 
    /*5*/  * (fconf(m,n) + gconf(m,n))) * fA(m,n) * gA(m,n) * (2 * P_1_1(R(m,n)) + (1 
    /*6*/  - 2 * pow2(Lt(m,n))) * P_2_1(R(m,n)))) / (2. * Lt(m,n)))) + gAlp_r(m,n) * (8
    /*2*/  * Power(r_minus(m,n),4) * gconf_r(m,n) - (4 * Power(r_minus(m,n),4) * gA_r(m,n)) 
    /*2*/  / gA(m,n) + (8 * Power(r_minus(m,n),4) * gB_r(m,n)) / gB(m,n) + (8 * pow3(r_minus(m,n)) 
    /*3*/  * (-1 + r(m,n))) / r(m,n)) - Kelas * gtrA(m,n) - Kelas * gtrK(m,n) - 2 
    /*1*/  * gBet(m,n) * pow2(r_minus(m,n)) * (gtrA_r(m,n) + gtrK_r(m,n)));
}
