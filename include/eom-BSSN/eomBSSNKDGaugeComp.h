/** @file  eomBSSNKDGaugeComp.h
 *  @author Francesco Torsello
 *  @brief The compactified K-driver parabolic gauge condition on the lapse (p. 111 of B&S).
 *  @version 2019-05-15T15:51:12
 *  @image html BSSNKDGaugeComp.png
 */

Real KD_gAlp_t( Int m, Int n )
{
    return Kdiff * ((Power(r_minus(m,n),4) * gAlp_rr(m,n)) / pow2(r_plus(m,n)) - (exp(4 
    /*3*/  * gconf(m,n)) * gAlp(m,n) * (pow2(0) + 3 * pow2(gK1(m,n)) + 6 
    /*3*/  * pow2(gK2(m,n))) * pow2(gA(m,n))) / 3. + k_g * (-(exp(2 * (fconf(m,n) 
    /*5*/  + gconf(m,n))) * fAlp(m,n) * fA(m,n) * gA(m,n) * P_1_2(R(m,n))) - (exp(4 
    /*4*/  * gconf(m,n)) * fAlp(m,n) * pow2(gA(m,n)) * (2 * P_1_1(R(m,n)) 
    /*4*/  + P_2_1(R(m,n)))) / (2. * Lt(m,n)) + gAlp(m,n) * (-(exp(4 * gconf(m,n)) 
    /*4*/  * pow2(gA(m,n)) * (gJ11(m,n) + 2 * gJ22(m,n) + grho(m,n) + 2 
    /*5*/  * P_1_0(R(m,n)))) / 2. - (exp(2 * (fconf(m,n) + gconf(m,n))) * fA(m,n) 
    /*4*/  * gA(m,n) * (2 * P_1_1(R(m,n)) + (1 - 2 * pow2(Lt(m,n))) * P_2_1(R(m,n)))) 
    /*3*/  / (2. * Lt(m,n)))) + gAlp_r(m,n) * ((2 * Power(r_minus(m,n),4) * gconf_r(m,n)) 
    /*2*/  / pow2(r_plus(m,n)) - (Power(r_minus(m,n),4) * gA_r(m,n)) / (gA(m,n) * pow2(r_plus(m,n))) + (2
    /*3*/  * Power(r_minus(m,n),4) * gB_r(m,n)) / (gB(m,n) * pow2(r_plus(m,n))) + (2 
    /*3*/  * Power(r_minus(m,n),4)) / (pow3(r_plus(m,n)) * r(m,n))) + (-6 * Kelas * (gtrA(m,n) 
    /*3*/  + gtrK(m,n)) - (6 * gBet(m,n) * pow2(r_minus(m,n)) * (gtrA_r(m,n) + gtrK_r(m,n))) 
    /*2*/  / r_plus(m,n)) / 6.);
}
