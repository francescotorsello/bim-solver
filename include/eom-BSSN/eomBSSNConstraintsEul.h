/** @file  eomBSSNConstraintsEul.h
 *  @author Francesco Torsello
 *  @brief The cBSSN regularized Eulerian constraint equations for metrics.
 *  @version 2019-05-15T15:50:29
 *  @image html BSSNconstraintsEul.png
 */

Real eq_gHC( Int m, Int n )
{
    return -pow2(gA1(m,n)) - 2 * pow2(gA1(m,n) - (gAsig(m,n) * pow2(r(m,n))) 
    /*1*/  / pow2(r_minus(m,n))) + (2 * pow2(gtrK(m,n))) / 3. + exp(-4 * gconf(m,n)) * (((4 
    /*3*/  * Power(r_minus(m,n),4) * gA_r(m,n) * gB_r(m,n)) / (gB(m,n) * pow2(r_plus(m,n))) + 4 
    /*2*/  * gA_r(m,n) * pow3(r_minus(m,n)) * ((2 * r_minus(m,n) * gconf_r(m,n)) / pow2(r_plus(m,n)) - 1 
    /*3*/  / (r(m,n) + Power(r(m,n),3)))) / pow3(gA(m,n)) + ((-2 * Power(r_minus(m,n),4) 
    /*3*/  * pow2(gB_r(m,n))) / (pow2(r_plus(m,n)) * pow2(gB(m,n))) - (2 * (gsig(m,n) 
    /*4*/  * pow3(r_plus(m,n)) * r(m,n) + 4 * Power(r_minus(m,n),4) * (2 * gconf_r(m,n) + r_plus(m,n) 
    /*5*/  * gconf_rr(m,n) * r(m,n) + r_plus(m,n) * pow2(gconf_r(m,n)) * r(m,n)))) 
    /*2*/  / (pow3(r_plus(m,n)) * r(m,n)) + (4 * pow3(r_minus(m,n)) * (-(r_minus(m,n) * r_plus(m,n) * gB_rr(m,n)
    /*5*/  * r(m,n)) + gB_r(m,n) * (3 - 4 * r_minus(m,n) * r_plus(m,n) * gconf_r(m,n) * r(m,n) 
    /*5*/  + Power(r(m,n),4)))) / (gB(m,n) * pow3(r_plus(m,n)) * r(m,n))) / pow2(gA(m,n))) 
    /*0*/  + k_g * (-2 * grho(m,n) + 2 * P_2_0(R(m,n)) + (2 * exp(2 * fconf(m,n) - 2 
    /*3*/  * gconf(m,n)) * fA(m,n) * P_2_1(R(m,n)) * Lt(m,n)) / gA(m,n));
}
Real eq_fHC( Int m, Int n )
{
    return -3 * pow2(fA1(m,n)) + (4 * fA1(m,n) * fAsig(m,n) * pow2(r(m,n))) 
    /*0*/  / pow2(r_minus(m,n)) + (2 * pow2(ftrK(m,n))) / 3. - (2 * pow2(fAsig(m,n)) 
    /*1*/  * Power(r(m,n),4)) / Power(r_minus(m,n),4) + exp(-4 * fconf(m,n)) * (((4 
    /*3*/  * Power(r_minus(m,n),4) * fA_r(m,n) * fB_r(m,n)) / (fB(m,n) * pow2(r_plus(m,n))) + 4 
    /*2*/  * fA_r(m,n) * pow3(r_minus(m,n)) * ((2 * r_minus(m,n) * fconf_r(m,n)) / pow2(r_plus(m,n)) - 1 
    /*3*/  / (r(m,n) + Power(r(m,n),3)))) / pow3(fA(m,n)) + ((-2 * pow2(fB_r(m,n)) * (-1
    /*4*/  + pow2(r(m,n))) * pow3(r_minus(m,n))) / (pow2(r_plus(m,n)) * pow2(fB(m,n))) - (4 
    /*3*/  * (Power(r_minus(m,n),4) * r_plus(m,n) * fB_rr(m,n) * r(m,n) + fB_r(m,n) * (-3 
    /*5*/  * pow3(r_minus(m,n)) - 6 * pow2(r(m,n)) * (1 + pow3(r_minus(m,n))) + 4 * Power(r_minus(m,n),4)
    /*5*/  * r_plus(m,n) * fconf_r(m,n) * r(m,n) + (16 - 3 * pow3(r_minus(m,n))) * Power(r(m,n),4)
    /*5*/  - 12 * Power(r(m,n),6) + 2 * Power(r(m,n),10)))) / (fB(m,n) * pow3(r_plus(m,n)) 
    /*3*/  * r(m,n)) - (2 * (fsig(m,n) * pow3(r_plus(m,n)) * r(m,n) + 4 * (Power(r_minus(m,n),4) 
    /*5*/  * r_plus(m,n) * fconf_rr(m,n) * r(m,n) + Power(r_minus(m,n),4) * r_plus(m,n) 
    /*5*/  * pow2(fconf_r(m,n)) * r(m,n) + 2 * fconf_r(m,n) * (-pow3(r_minus(m,n)) 
    /*6*/  - pow2(r(m,n)) * (3 + 2 * pow3(r_minus(m,n))) - (-8 + pow3(r_minus(m,n))) 
    /*6*/  * Power(r(m,n),4) - 6 * Power(r(m,n),6) + Power(r(m,n),10))))) / (pow3(r_plus(m,n))
    /*3*/  * r(m,n))) / pow2(fA(m,n))) - k_f * (2 * frho(m,n) - (2 * exp(-4 
    /*3*/  * fconf(m,n) + 4 * gconf(m,n)) * pow2(gB(m,n)) * P_2_2(R(m,n))) 
    /*1*/  / pow2(fB(m,n)) - (2 * exp(-6 * fconf(m,n) + 6 * gconf(m,n)) * gA(m,n) 
    /*2*/  * pow2(gB(m,n)) * P_2_1(R(m,n)) * Lt(m,n)) / (fA(m,n) * pow2(fB(m,n))));
}
Real eq_gMC( Int m, Int n )
{
    return k_g * (3 * exp(4 * gconf(m,n)) * gj(m,n) * gA(m,n) - (3 * exp(2 
    /*3*/  * fconf(m,n)) * fA(m,n) * p(m,n) * P_2_1(R(m,n))) / gA(m,n)) + ((-6 
    /*2*/  * gAsig(m,n) * gB_r(m,n) * pow2(r(m,n))) / (r_plus(m,n) * gB(m,n)) + (6 * r_plus(m,n) 
    /*2*/  * gAsig(m,n) * r(m,n) + pow3(r_minus(m,n)) * (-3 * gA1_r(m,n) + 2 * (-9 * gA1(m,n)
    /*4*/  * gconf_r(m,n) + gtrK_r(m,n)))) / (r_minus(m,n) * r_plus(m,n))) / gA(m,n);
}
Real eq_fMC( Int m, Int n )
{
    return k_f * (-3 * exp(4 * fconf(m,n)) * fj(m,n) * fA(m,n) - (3 * exp(-4 
    /*3*/  * fconf(m,n) + 6 * gconf(m,n)) * gA(m,n) * p(m,n) * pow2(gB(m,n)) 
    /*2*/  * P_2_1(R(m,n))) / (fA(m,n) * pow2(fB(m,n)))) + ((6 * fAsig(m,n) * fB_r(m,n)
    /*2*/  * pow2(r(m,n))) / (r_plus(m,n) * fB(m,n)) + (6 * r_minus(m,n) * fA1(m,n) * (r_plus(m,n) * (1
    /*4*/  + r_minus(m,n) - pow2(r(m,n))) + 3 * fconf_r(m,n) * pow2(r_minus(m,n)) * r(m,n)) 
    /*2*/  + r(m,n) * (-6 * r_plus(m,n) * fAsig(m,n) * r(m,n) + pow3(r_minus(m,n)) * (3 * fA1_r(m,n)
    /*4*/  - 2 * ftrK_r(m,n)))) / (r_minus(m,n) * r_plus(m,n) * r(m,n))) / fA(m,n);
}
Real eq_gLC( Int m, Int n )
{
    return -gL(m,n) + (gA_r(m,n) * pow3(r_minus(m,n))) / (pow3(gA(m,n)) * (-1 
    /*2*/  + Power(r(m,n),4))) + ((2 * gsig(m,n) * r(m,n)) / r_minus(m,n) - (2 * gB_r(m,n) 
    /*2*/  * pow3(r_minus(m,n))) / (gB(m,n) * (-1 + Power(r(m,n),4)))) / pow2(gA(m,n));
}
Real eq_fLC( Int m, Int n )
{
    return fL(m,n) - (fA_r(m,n) * pow3(r_minus(m,n))) / (pow3(fA(m,n)) * (-1 
    /*2*/  + Power(r(m,n),4))) - ((2 * fsig(m,n) * r(m,n)) / r_minus(m,n) - (2 * fB_r(m,n) 
    /*2*/  * pow3(r_minus(m,n))) / (fB(m,n) * (-1 + Power(r(m,n),4)))) / pow2(fA(m,n));
}
Real eq_CL( Int m, Int n )
{
    return p(m,n) * ((exp(2 * fconf(m,n) - 2 * gconf(m,n)) * ((2 * fA(m,n) * gB_r(m,n)
    /*4*/  * pow2(r_minus(m,n)) * P_1_2(R(m,n))) / (r_plus(m,n) * gB(m,n)) - (2 * fA(m,n) 
    /*4*/  * P_1_2(R(m,n)) * (-1 - 2 * gconf_r(m,n) * pow2(r_minus(m,n)) * r(m,n) 
    /*5*/  + Power(r(m,n),4))) / (r_plus(m,n) * r(m,n)))) / gA(m,n) + exp(-2 * fconf(m,n) + 2
    /*2*/  * gconf(m,n)) * gA(m,n) * ((-2 * fB(m,n) * gB_r(m,n) * pow2(r_minus(m,n)) 
    /*3*/  * P_1_1(R(m,n))) / (r_plus(m,n) * fA(m,n) * pow2(gB(m,n))) + (2 * P_1_1(R(m,n)) 
    /*3*/  * (r_plus(m,n) - r_plus(m,n) * pow2(r(m,n)) + 2 * gconf_r(m,n) * pow2(r_minus(m,n)) * r(m,n)) 
    /*3*/  * R(m,n)) / (r_plus(m,n) * fA(m,n) * r(m,n)) + (2 * pow2(r_minus(m,n)) * P_1_1(R(m,n)) 
    /*3*/  * (fB_r(m,n) + gB_r(m,n) * R(m,n))) / (r_plus(m,n) * fA(m,n) * gB(m,n)))) 
    /*0*/  + (p_r(m,n) * pow2(r_minus(m,n)) * P_2_1(R(m,n))) / (r_plus(m,n) * Lt(m,n)) + (exp(2 
    /*2*/  * gconf(m,n)) * gA(m,n) * ((6 * gAsig(m,n) * pow2(r(m,n)) * P_1_1(R(m,n))) 
    /*2*/  / pow2(r_minus(m,n)) - 3 * gA1(m,n) * (2 * P_1_1(R(m,n)) + P_2_1(R(m,n))) - 2 
    /*2*/  * P_1_1(R(m,n)) * gtrK(m,n) - P_2_1(R(m,n)) * gtrK(m,n) + 6 * fA1(m,n) 
    /*2*/  * P_1_1(R(m,n)) * R(m,n) * Lt(m,n) - (6 * fAsig(m,n) * pow2(r(m,n)) 
    /*3*/  * P_1_1(R(m,n)) * R(m,n) * Lt(m,n)) / pow2(r_minus(m,n)) + 2 * P_1_1(R(m,n)) 
    /*2*/  * R(m,n) * ftrK(m,n) * Lt(m,n))) / 3. + exp(2 * fconf(m,n)) * fA(m,n) * ((-2
    /*2*/  * fAsig(m,n) * pow2(r(m,n)) * P_1_2(R(m,n)) * R(m,n)) / pow2(r_minus(m,n)) 
    /*1*/  + fA1(m,n) * (P_2_1(R(m,n)) + 2 * P_1_2(R(m,n)) * R(m,n)) - (2 
    /*2*/  * P_1_1(R(m,n)) * ftrK(m,n)) / 3. + P_2_1(R(m,n)) * ftrK(m,n) - 2 * gA1(m,n)
    /*1*/  * P_1_2(R(m,n)) * Lt(m,n) + (2 * gAsig(m,n) * pow2(r(m,n)) * P_1_2(R(m,n))
    /*2*/  * Lt(m,n)) / pow2(r_minus(m,n)) - (2 * P_1_2(R(m,n)) * gtrK(m,n) * Lt(m,n)) / 3.);
}
