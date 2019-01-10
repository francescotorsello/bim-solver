/** @file  eomBSSNConstraintsReg.h
 *  @author Francesco Torsello
 *  @brief The cBSSN regularized constraint equations for metrics.
 *  @version 2019-01-09T12:00:48
 *  @image html BSSNconstraintsReg.png
 */

Real BimetricEvolve::eq_gHC( Int m, Int n )
{
    return 4 * gA1(m,n) * gAsig(m,n) * pow2(r(m,n)) + 2 * k_g * (-grho(m,n) 
    /*1*/  + P_2_0(R(m,n)) + (exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * Lt(m,n) 
    /*2*/  * P_2_1(R(m,n))) / gA(m,n)) + (4 * exp(-4 * gconf(m,n)) * (gDA(m,n) - 3 
    /*2*/  * gDB(m,n) - 4 * gDconf(m,n))) / (pow2(gA(m,n)) * r(m,n)) - 2 
    /*0*/  * pow2(gAsig(m,n)) * Power(r(m,n),4) + (exp(-4 * gconf(m,n)) * (-12 
    /*2*/  * gDB_r(m,n) - 48 * gDB(m,n) * gDconf(m,n) + 12 * gDA(m,n) * (gDB(m,n) + 2 
    /*3*/  * gDconf(m,n)) - 24 * gDconf_r(m,n) - 6 * gsig(m,n) - 18 * pow2(gDB(m,n)) 
    /*2*/  - 24 * pow2(gDconf(m,n)) - 9 * exp(4 * gconf(m,n)) * pow2(gA1(m,n)) 
    /*2*/  * pow2(gA(m,n)) + 3 * exp(4 * gconf(m,n)) * pow2(gA(m,n)) * pow2(gtrA(m,n)) 
    /*2*/  + 2 * exp(4 * gconf(m,n)) * pow2(gA(m,n)) * pow2(gtrK(m,n)) + 4 * exp(4 
    /*3*/  * gconf(m,n)) * pow2(gA(m,n)) * gtrA(m,n) * gtrK(m,n))) / (3. 
    /*1*/  * pow2(gA(m,n)));
}
Real BimetricEvolve::eq_fHC( Int m, Int n )
{
    return 4 * fA1(m,n) * fAsig(m,n) * pow2(r(m,n)) + k_f * (-2 * frho(m,n) + (2 
    /*2*/  * ((exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * Lt(m,n) 
    /*4*/  * P_2_1(R(m,n))) / fA(m,n) + P_2_2(R(m,n)))) / pow2(R(m,n))) + (4 * exp(-4 
    /*2*/  * fconf(m,n)) * (fDA(m,n) - 3 * fDB(m,n) - 4 * fDconf(m,n))) / (pow2(fA(m,n))
    /*1*/  * r(m,n)) - 2 * pow2(fAsig(m,n)) * Power(r(m,n),4) + (exp(-4 * fconf(m,n))
    /*1*/  * (-12 * fDB_r(m,n) - 48 * fDB(m,n) * fDconf(m,n) + 12 * fDA(m,n) 
    /*2*/  * (fDB(m,n) + 2 * fDconf(m,n)) - 24 * fDconf_r(m,n) - 6 * fsig(m,n) - 18 
    /*2*/  * pow2(fDB(m,n)) - 24 * pow2(fDconf(m,n)) - 9 * exp(4 * fconf(m,n)) 
    /*2*/  * pow2(fA1(m,n)) * pow2(fA(m,n)) + 3 * exp(4 * fconf(m,n)) * pow2(fA(m,n)) 
    /*2*/  * pow2(ftrA(m,n)) + 2 * exp(4 * fconf(m,n)) * pow2(fA(m,n)) * pow2(ftrK(m,n))
    /*2*/  + 4 * exp(4 * fconf(m,n)) * pow2(fA(m,n)) * ftrA(m,n) * ftrK(m,n))) / (3. 
    /*1*/  * pow2(fA(m,n)));
}
Real BimetricEvolve::eq_gMC( Int m, Int n )
{
    return -3 * gA1_r(m,n) - 18 * gA1(m,n) * gDconf(m,n) - 6 * gAsig(m,n) * gDB(m,n) 
    /*0*/  * pow2(r(m,n)) + k_g * (3 * exp(4 * gconf(m,n)) * gj(m,n) * pow2(gA(m,n)) - 3
    /*1*/  * exp(2 * fconf(m,n)) * fA(m,n) * p(m,n) * P_2_1(R(m,n))) - 6 * gAsig(m,n)
    /*0*/  * r(m,n) + 6 * gDconf(m,n) * gtrA(m,n) + 3 * gtrA_r(m,n) + 2 * gtrK_r(m,n);
}
Real BimetricEvolve::eq_fMC( Int m, Int n )
{
    return -3 * fA1_r(m,n) - 18 * fA1(m,n) * fDconf(m,n) - 6 * fAsig(m,n) * fDB(m,n) 
    /*0*/  * pow2(r(m,n)) + k_f * (3 * exp(4 * fconf(m,n)) * fj(m,n) * pow2(fA(m,n)) 
    /*1*/  + (3 * exp(-4 * fconf(m,n) + 6 * gconf(m,n)) * gA(m,n) * p(m,n) 
    /*2*/  * pow2(gB(m,n)) * P_2_1(R(m,n))) / pow2(fB(m,n))) - 6 * fAsig(m,n) * r(m,n) 
    /*0*/  + 6 * fDconf(m,n) * ftrA(m,n) + 3 * ftrA_r(m,n) + 2 * ftrK_r(m,n);
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
    return 6 * Lt(m,n) * p(m,n) * (exp(4 * fconf(m,n)) * pow2(fA(m,n)) * P_1_2(R(m,n))
    /*1*/  + exp(4 * gconf(m,n)) * pow2(gA(m,n)) * P_1_1(R(m,n)) * R(m,n)) + (r(m,n) 
    /*1*/  * (6 * exp(4 * gconf(m,n)) * Lt(m,n) * p(m,n) * pow2(gA(m,n)) * P_1_1(R(m,n))
    /*2*/  * (fDB(m,n) * fB(m,n) - fB(m,n) * gDB(m,n) + (gDB(m,n) + 2 * gDconf(m,n)) 
    /*3*/  * gB(m,n) * R(m,n)) + exp(4 * fconf(m,n)) * gB(m,n) * Lt(m,n) * pow2(fA(m,n))
    /*2*/  * (-6 * exp(2 * gconf(m,n)) * gA2(m,n) * gA(m,n) * Lt(m,n) * P_1_2(R(m,n))
    /*3*/  + 6 * gDB(m,n) * p(m,n) * P_1_2(R(m,n)) + 12 * gDconf(m,n) * p(m,n) 
    /*3*/  * P_1_2(R(m,n)) + 3 * exp(2 * gconf(m,n)) * fA1(m,n) * gA(m,n) 
    /*3*/  * P_2_1(R(m,n)) + 6 * exp(2 * gconf(m,n)) * fA2(m,n) * gA(m,n) 
    /*3*/  * P_1_2(R(m,n)) * R(m,n) - 2 * exp(2 * gconf(m,n)) * gA(m,n) * P_1_1(R(m,n))
    /*3*/  * ftrK(m,n) + 3 * exp(2 * gconf(m,n)) * gA(m,n) * P_2_1(R(m,n)) * ftrK(m,n)
    /*3*/  - 2 * exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n) * P_1_2(R(m,n)) * gtrK(m,n))
    /*2*/  - exp(2 * (fconf(m,n) + gconf(m,n))) * fA(m,n) * gA(m,n) * gB(m,n) * (6 
    /*3*/  * exp(2 * gconf(m,n)) * gA2(m,n) * gA(m,n) * Lt(m,n) * P_1_1(R(m,n)) + 3 
    /*3*/  * exp(2 * gconf(m,n)) * gA1(m,n) * gA(m,n) * Lt(m,n) * P_2_1(R(m,n)) - 3 
    /*3*/  * p_r(m,n) * P_2_1(R(m,n)) - 6 * exp(2 * gconf(m,n)) * fA2(m,n) * gA(m,n) 
    /*3*/  * Lt2(m,n) * P_1_1(R(m,n)) * R(m,n) - 2 * exp(2 * gconf(m,n)) * gA(m,n) 
    /*3*/  * Lt2(m,n) * P_1_1(R(m,n)) * R(m,n) * ftrK(m,n) + 2 * exp(2 * gconf(m,n)) 
    /*3*/  * gA(m,n) * Lt(m,n) * P_1_1(R(m,n)) * gtrK(m,n) + exp(2 * gconf(m,n)) 
    /*3*/  * gA(m,n) * Lt(m,n) * P_2_1(R(m,n)) * gtrK(m,n)))) / gB(m,n);
}
