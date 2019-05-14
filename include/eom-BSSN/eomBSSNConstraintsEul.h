/** @file  eomBSSNConstraintsEul.h
 *  @author Francesco Torsello
 *  @brief The cBSSN regularized Eulerian constraint equations for metrics.
 *  @version 2019-05-14T16:23:36
 *  @image html BSSNconstraintsEul.png
 */

Real eq_gHC( Int m, Int n )
{
    return -pow2(gA1(m,n)) - 2 * pow2(gA1(m,n) - (gAsig(m,n) * pow2(r(m,n))) / (4. 
    /*2*/  * pow2(r_minus(m,n)))) + (2 * pow2(gtrK(m,n))) / 3. + exp(-4 * gconf(m,n)) * (((16
    /*3*/  * Power(r_minus(m,n),4) * gA_r(m,n) * gB_r(m,n)) / gB(m,n) + (16 * gA_r(m,n) 
    /*3*/  * pow3(r_minus(m,n)) * (-1 + 2 * r_minus(m,n) * gconf_r(m,n) * r(m,n))) / r(m,n)) 
    /*1*/  / pow3(gA(m,n)) + ((-8 * Power(r_minus(m,n),4) * pow2(gB_r(m,n))) / pow2(gB(m,n)) 
    /*2*/  - (16 * pow3(r_minus(m,n)) * (r_minus(m,n) * gB_rr(m,n) * r(m,n) + gB_r(m,n) * (-1 + 2 
    /*5*/  * r_minus(m,n) + 4 * r_minus(m,n) * gconf_r(m,n) * r(m,n)))) / (gB(m,n) * r(m,n)) - (2 
    /*3*/  * (gsig(m,n) * r(m,n) + 16 * Power(r_minus(m,n),4) * (2 * gconf_r(m,n) 
    /*5*/  + gconf_rr(m,n) * r(m,n) + pow2(gconf_r(m,n)) * r(m,n)))) / r(m,n)) 
    /*1*/  / pow2(gA(m,n))) + k_g * (-2 * grho(m,n) + 2 * P_2_0(R(m,n)) + (2 * exp(2 
    /*3*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * P_2_1(R(m,n)) * Lt(m,n)) 
    /*1*/  / gA(m,n));
}
Real eq_fHC( Int m, Int n )
{
    return -3 * pow2(fA1(m,n)) + (fA1(m,n) * fAsig(m,n) * pow2(r(m,n))) / pow2(r_minus(m,n))
    /*0*/  + (2 * pow2(ftrK(m,n))) / 3. - (pow2(fAsig(m,n)) * Power(r(m,n),4)) / (8. 
    /*1*/  * Power(r_minus(m,n),4)) + exp(-4 * fconf(m,n)) * (((16 * Power(r_minus(m,n),4) 
    /*3*/  * fA_r(m,n) * fB_r(m,n)) / fB(m,n) + (16 * fA_r(m,n) * pow3(r_minus(m,n)) * (-1 + 2
    /*4*/  * r_minus(m,n) * fconf_r(m,n) * r(m,n))) / r(m,n)) / pow3(fA(m,n)) + ((-8 
    /*3*/  * Power(r_minus(m,n),4) * pow2(fB_r(m,n))) / pow2(fB(m,n)) - (16 * pow3(r_minus(m,n)) 
    /*3*/  * (r_minus(m,n) * fB_rr(m,n) * r(m,n) + fB_r(m,n) * (-1 + 2 * r_minus(m,n) + 4 * r_minus(m,n) 
    /*5*/  * fconf_r(m,n) * r(m,n)))) / (fB(m,n) * r(m,n)) - (2 * (fsig(m,n) * r(m,n) 
    /*4*/  + 16 * Power(r_minus(m,n),4) * (2 * fconf_r(m,n) + fconf_rr(m,n) * r(m,n) 
    /*5*/  + pow2(fconf_r(m,n)) * r(m,n)))) / r(m,n)) / pow2(fA(m,n))) - k_f * (2 
    /*1*/  * frho(m,n) - (2 * P_2_2(R(m,n))) / pow2(R(m,n)) - (2 * exp(-2 * fconf(m,n) 
    /*3*/  + 2 * gconf(m,n)) * gA(m,n) * P_2_1(R(m,n)) * Lt(m,n)) / (fA(m,n) 
    /*2*/  * pow2(R(m,n))));
}
Real eq_gMC( Int m, Int n )
{
    return k_g * (-3 * exp(4 * gconf(m,n)) * gj(m,n) * gA(m,n) + (3 * exp(2 
    /*3*/  * fconf(m,n)) * fA(m,n) * p(m,n) * P_2_1(R(m,n))) / gA(m,n)) + ((3 
    /*2*/  * gAsig(m,n) * gB_r(m,n) * pow2(r(m,n))) / gB(m,n) + (-3 * gAsig(m,n) 
    /*2*/  * r(m,n) + 2 * pow3(r_minus(m,n)) * (3 * gA1_r(m,n) + 18 * gA1(m,n) * gconf_r(m,n)
    /*3*/  - 2 * gtrK_r(m,n))) / r_minus(m,n)) / gA(m,n);
}
Real eq_fMC( Int m, Int n )
{
    return k_f * (-3 * exp(4 * fconf(m,n)) * fj(m,n) * fA(m,n) - (3 * exp(-4 
    /*3*/  * fconf(m,n) + 6 * gconf(m,n)) * gA(m,n) * p(m,n) * pow2(gB(m,n)) 
    /*2*/  * P_2_1(R(m,n))) / (fA(m,n) * pow2(fB(m,n)))) + ((3 * fAsig(m,n) * fB_r(m,n)
    /*2*/  * pow2(r(m,n))) / fB(m,n) + (-3 * fAsig(m,n) * r(m,n) + 2 * pow3(r_minus(m,n)) 
    /*2*/  * (3 * fA1_r(m,n) + 18 * fA1(m,n) * fconf_r(m,n) - 2 * ftrK_r(m,n))) 
    /*1*/  / r_minus(m,n)) / fA(m,n);
}
Real eq_gLC( Int m, Int n )
{
    return -gL(m,n) + (((1 + r_minus(m,n)) * gsig(m,n)) / r_minus(m,n) - (4 * gB_r(m,n) 
    /*2*/  * pow2(r_minus(m,n))) / gB(m,n)) / pow2(gA(m,n)) + (2 * gA_r(m,n) * pow2(r_minus(m,n))) 
    /*0*/  / pow3(gA(m,n));
}
Real eq_fLC( Int m, Int n )
{
    return fL(m,n) - (((1 + r_minus(m,n)) * fsig(m,n)) / r_minus(m,n) - (4 * fB_r(m,n) 
    /*2*/  * pow2(r_minus(m,n))) / fB(m,n)) / pow2(fA(m,n)) - (2 * fA_r(m,n) * pow2(r_minus(m,n))) 
    /*0*/  / pow3(fA(m,n));
}
Real eq_CL( Int m, Int n )
{
    return (4 * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * fB_r(m,n) * gA(m,n) * p(m,n) 
    /*1*/  * pow2(r_minus(m,n)) * P_1_1(R(m,n))) / (fA(m,n) * gB(m,n)) - (4 * r_minus(m,n) * exp(2 
    /*2*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * p(m,n) * P_1_2(R(m,n))) / (gA(m,n)
    /*1*/  + r_minus(m,n) * gA(m,n)) + gconf_r(m,n) * ((8 * exp(2 * fconf(m,n) - 2 
    /*3*/  * gconf(m,n)) * fA(m,n) * p(m,n) * pow2(r_minus(m,n)) * P_1_2(R(m,n))) / gA(m,n) 
    /*1*/  + (8 * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * p(m,n) 
    /*2*/  * pow2(r_minus(m,n)) * P_1_1(R(m,n)) * R(m,n)) / fA(m,n)) + gB_r(m,n) * ((4 * exp(2
    /*3*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * p(m,n) * pow2(r_minus(m,n)) 
    /*2*/  * P_1_2(R(m,n))) / (gA(m,n) * gB(m,n)) + (4 * exp(-2 * fconf(m,n) + 2 
    /*3*/  * gconf(m,n)) * gA(m,n) * p(m,n) * pow2(r_minus(m,n)) * P_1_1(R(m,n)) * (-fB(m,n) 
    /*3*/  + gB(m,n) * R(m,n))) / (fA(m,n) * pow2(gB(m,n)))) + (2 * p_r(m,n) 
    /*1*/  * pow2(r_minus(m,n)) * P_2_1(R(m,n))) / Lt(m,n) + gAsig(m,n) * ((exp(2 
    /*3*/  * gconf(m,n)) * gA(m,n) * pow2(1 + r_minus(m,n)) * P_1_1(R(m,n))) / (2. 
    /*2*/  * pow2(r_minus(m,n))) + (exp(2 * fconf(m,n)) * fA(m,n) * pow2(1 + r_minus(m,n)) 
    /*2*/  * P_1_2(R(m,n)) * Lt(m,n)) / (2. * pow2(r_minus(m,n)))) + fAsig(m,n) * (-(exp(2 
    /*3*/  * fconf(m,n)) * fA(m,n) * pow2(1 + r_minus(m,n)) * P_1_2(R(m,n)) * R(m,n)) / (2. 
    /*2*/  * pow2(r_minus(m,n))) - (exp(2 * gconf(m,n)) * gA(m,n) * pow2(1 + r_minus(m,n)) 
    /*2*/  * P_1_1(R(m,n)) * R(m,n) * Lt(m,n)) / (2. * pow2(r_minus(m,n)))) + exp(2 
    /*1*/  * fconf(m,n)) * (fA(m,n) * ((-2 * P_1_1(R(m,n))) / 3. + P_2_1(R(m,n))) 
    /*1*/  * ftrK(m,n) - (2 * fA(m,n) * P_1_2(R(m,n)) * gtrK(m,n) * Lt(m,n)) / 3. 
    /*1*/  + fA(m,n) * (fA1(m,n) * (P_2_1(R(m,n)) + 2 * P_1_2(R(m,n)) * R(m,n)) - 2 
    /*2*/  * gA1(m,n) * P_1_2(R(m,n)) * Lt(m,n))) + exp(2 * gconf(m,n)) * ((-4 * r_minus(m,n)
    /*2*/  * exp(-2 * fconf(m,n)) * gA(m,n) * p(m,n) * P_1_1(R(m,n)) * R(m,n)) 
    /*1*/  / (fA(m,n) + r_minus(m,n) * fA(m,n)) - (gA(m,n) * (2 * P_1_1(R(m,n)) 
    /*3*/  + P_2_1(R(m,n))) * gtrK(m,n)) / 3. + (2 * gA(m,n) * P_1_1(R(m,n)) * R(m,n) 
    /*2*/  * ftrK(m,n) * Lt(m,n)) / 3. - gA(m,n) * (gA1(m,n) * (2 * P_1_1(R(m,n)) 
    /*3*/  + P_2_1(R(m,n))) - 2 * fA1(m,n) * P_1_1(R(m,n)) * R(m,n) * Lt(m,n)));
}
