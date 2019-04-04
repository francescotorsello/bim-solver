/** @file  eomBSSNConstraintsEul.h
 *  @author Francesco Torsello
 *  @brief The cBSSN regularized Eulerian constraint equations for metrics.
 *  @version 2019-04-04T09:57:34
 *  @image html BSSNconstraintsEul.png
 */

Real BimetricEvolve::eq_gHC( Int m, Int n )
{
    return -3 * pow2(gA1(m,n)) + (exp(-4 * gconf(m,n)) * (4 * gDA(m,n) * (gDB(m,n) + 2
    /*3*/  * gDconf(m,n)) - 2 * (2 * gDB_r(m,n) + 8 * gDB(m,n) * gDconf(m,n) + 4 
    /*3*/  * gDconf_r(m,n) + gsig(m,n) + 3 * pow2(gDB(m,n)) + 4 * pow2(gDconf(m,n))))) 
    /*0*/  / pow2(gA(m,n)) + 4 * gA1(m,n) * gAsig(m,n) * pow2(r(m,n)) + (2 
    /*1*/  * pow2(gtrK(m,n))) / 3. + k_g * (-2 * grho(m,n) + 2 * P_2_0(R(m,n)) + (2 
    /*2*/  * exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * Lt(m,n) * P_2_1(R(m,n))) 
    /*1*/  / gA(m,n)) + (4 * exp(-4 * gconf(m,n)) * (gDA(m,n) - 3 * gDB(m,n) - 4 
    /*2*/  * gDconf(m,n))) / (pow2(gA(m,n)) * r(m,n)) - 2 * pow2(gAsig(m,n)) 
    /*0*/  * Power(r(m,n),4);
}
Real BimetricEvolve::eq_fHC( Int m, Int n )
{
    return -3 * pow2(fA1(m,n)) + (exp(-4 * fconf(m,n)) * (4 * fDA(m,n) * (fDB(m,n) + 2
    /*3*/  * fDconf(m,n)) - 2 * (2 * fDB_r(m,n) + 8 * fDB(m,n) * fDconf(m,n) + 4 
    /*3*/  * fDconf_r(m,n) + fsig(m,n) + 3 * pow2(fDB(m,n)) + 4 * pow2(fDconf(m,n))))) 
    /*0*/  / pow2(fA(m,n)) + 4 * fA1(m,n) * fAsig(m,n) * pow2(r(m,n)) + (2 
    /*1*/  * pow2(ftrK(m,n))) / 3. + k_f * (-2 * frho(m,n) + (2 * exp(-2 * fconf(m,n) 
    /*3*/  + 2 * gconf(m,n)) * gA(m,n) * Lt(m,n) * P_2_1(R(m,n))) / (fA(m,n) 
    /*2*/  * pow2(R(m,n))) + (2 * P_2_2(R(m,n))) / pow2(R(m,n))) + (4 * exp(-4 
    /*2*/  * fconf(m,n)) * (fDA(m,n) - 3 * fDB(m,n) - 4 * fDconf(m,n))) / (pow2(fA(m,n))
    /*1*/  * r(m,n)) - 2 * pow2(fAsig(m,n)) * Power(r(m,n),4);
}
Real BimetricEvolve::eq_gMC( Int m, Int n )
{
    return -3 * gA1_r(m,n) - 18 * gA1(m,n) * gDconf(m,n) - 6 * gAsig(m,n) * gDB(m,n) 
    /*0*/  * pow2(r(m,n)) + k_g * (3 * exp(4 * gconf(m,n)) * gj(m,n) * pow2(gA(m,n)) - 3
    /*1*/  * exp(2 * fconf(m,n)) * fA(m,n) * p(m,n) * P_2_1(R(m,n))) - 6 * gAsig(m,n)
    /*0*/  * r(m,n) + 2 * gtrK_r(m,n);
}
Real BimetricEvolve::eq_fMC( Int m, Int n )
{
    return -3 * fA1_r(m,n) - 18 * fA1(m,n) * fDconf(m,n) - 6 * fAsig(m,n) * fDB(m,n) 
    /*0*/  * pow2(r(m,n)) + k_f * (3 * exp(4 * fconf(m,n)) * fj(m,n) * pow2(fA(m,n)) 
    /*1*/  + (3 * exp(-4 * fconf(m,n) + 6 * gconf(m,n)) * gA(m,n) * p(m,n) 
    /*2*/  * pow2(gB(m,n)) * P_2_1(R(m,n))) / pow2(fB(m,n))) - 6 * fAsig(m,n) * r(m,n) 
    /*0*/  + 2 * ftrK_r(m,n);
}
Real BimetricEvolve::eq_gLC( Int m, Int n )
{
    return gL(m,n) - (gDA(m,n) - 2 * gDB(m,n)) / pow2(gA(m,n)) + (2 * gsig(m,n) 
    /*1*/  * r(m,n)) / pow2(gA(m,n));
}
Real BimetricEvolve::eq_fLC( Int m, Int n )
{
    return fL(m,n) - (fDA(m,n) - 2 * fDB(m,n)) / pow2(fA(m,n)) + (2 * fsig(m,n) 
    /*1*/  * r(m,n)) / pow2(fA(m,n));
}
Real BimetricEvolve::eq_CL( Int m, Int n )
{
    return Lt(m,n) * p(m,n) * (6 * exp(4 * fconf(m,n)) * pow2(fA(m,n)) * P_1_2(R(m,n))
    /*1*/  + 6 * exp(4 * gconf(m,n)) * pow2(gA(m,n)) * P_1_1(R(m,n)) * R(m,n)) 
    /*0*/  + r(m,n) * (3 * exp(2 * fconf(m,n) + 2 * gconf(m,n)) * fA(m,n) * gA(m,n) 
    /*1*/  * p_r(m,n) * P_2_1(R(m,n)) + Lt(m,n) * p(m,n) * (6 * exp(4 * fconf(m,n)) 
    /*2*/  * (gDB(m,n) + 2 * gDconf(m,n)) * pow2(fA(m,n)) * P_1_2(R(m,n)) + exp(4 
    /*3*/  * gconf(m,n)) * pow2(gA(m,n)) * ((6 * fB(m,n) * (fDB(m,n) - gDB(m,n)) 
    /*4*/  * P_1_1(R(m,n))) / gB(m,n) + 6 * (gDB(m,n) + 2 * gDconf(m,n)) * P_1_1(R(m,n))
    /*3*/  * R(m,n))) + exp(2 * fconf(m,n) + 4 * gconf(m,n)) * fA(m,n) * Lt2(m,n) 
    /*1*/  * pow2(gA(m,n)) * P_1_1(R(m,n)) * R(m,n) * (-3 * fA1(m,n) + 2 * ftrK(m,n)) 
    /*1*/  + exp(4 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * pow2(fA(m,n)) 
    /*1*/  * pow2(Lt(m,n)) * P_1_2(R(m,n)) * (3 * gA1(m,n) - 2 * gtrK(m,n)) + Lt(m,n) 
    /*1*/  * (exp(4 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * pow2(fA(m,n)) * (3 
    /*3*/  * fA1(m,n) * (P_2_1(R(m,n)) - P_1_2(R(m,n)) * R(m,n)) + (-2 * P_1_1(R(m,n)) 
    /*4*/  + 3 * P_2_1(R(m,n))) * ftrK(m,n)) + exp(2 * fconf(m,n) + 4 * gconf(m,n)) 
    /*2*/  * fA(m,n) * pow2(gA(m,n)) * (3 * gA1(m,n) * (P_1_1(R(m,n)) - P_2_1(R(m,n))) 
    /*3*/  - (2 * P_1_1(R(m,n)) + P_2_1(R(m,n))) * gtrK(m,n))));
}
