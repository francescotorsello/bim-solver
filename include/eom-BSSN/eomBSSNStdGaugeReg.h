/** @file  eomBSSNStdGaugeReg.h
 *  @author Francesco Torsello
 *  @brief The regularized standard gauge (Bona-Masso & Gamma-driver).
 *  @version 2018-11-07T15:06:39
 *  @image html BSSNStdGaugeReg.png
 */

Real BimetricEvolve::eq_BM_gAlp_t( Int m, Int n )
{
    return gAlp_convr(m,n) - 2 * gAlp(m,n) * gtrK(m,n);
}
Real BimetricEvolve::eq_BM_gBet_t( Int m, Int n )
{
    return (3 * Bq(m,n)) / (4 + TINY_Real) + q_qconvr(m,n);
}
Real BimetricEvolve::eq_BM_gBq_t( Int m, Int n )
{
    return -(eta * Bq(m,n)) + (4 * gBetr(m,n) * gL(m,n)) / (3 + TINY_Real) + (2 * (2 
    /*2*/  * gBet_rr(m,n) + gBetr_r(m,n))) / (TINY_Real + 3 * Power(gA(m,n),2)) + (2 
    /*1*/  * gBetr_r(m,n)) / (TINY_Real + Power(gB(m,n),2)) + Bq_qconvr(m,n) 
    /*0*/  - gL_qconvr(m,n) + k_g * gJL(m,n) + gBet_r(m,n) * (-gL(m,n) / 3. + gDA(m,n) 
    /*1*/  / (3. * pow2(gA(m,n))) + (2 * gDB(m,n)) / (3. * pow2(gA(m,n)))) + gBet(m,n) 
    /*0*/  * (((2 * gDA(m,n)) / 3. + (4 * gDB(m,n)) / 3.) * gL(m,n) + gL_r(m,n) 
    /*1*/  + (gDA_r(m,n) + 2 * gDB_r(m,n)) / (3. * pow2(gA(m,n)))) + gAlp(m,n) 
    /*0*/  * (gAsig(m,n) * ((4 * gDA(m,n) * pow2(r(m,n))) / (3. * pow2(gA(m,n))) + (4 
    /*3*/  * gDB(m,n) * pow2(r(m,n))) / (3. * pow2(gA(m,n))) + (8 * gconf(m,n) 
    /*3*/  * gDconf(m,n) * pow2(r(m,n))) / pow2(gA(m,n)) - (4 * gDAlp(m,n) 
    /*3*/  * pow2(r(m,n))) / (3. * pow2(gA(m,n)))) + (4 * gAsig(m,n) * gsig(m,n) 
    /*2*/  * pow3(r(m,n))) / (3. * pow2(gA(m,n))) - (4 * (gtrA_r(m,n) + gtrK_r(m,n))) 
    /*1*/  / (3. * pow2(gA(m,n))));
}
Real BimetricEvolve::eq_BM_gDAlp_t( Int m, Int n )
{
    return gDAlp(m,n) * gBet_r(m,n) - (gDAlp(m,n) * gAlp_convr(m,n)) / (TINY_Real 
    /*1*/  + gAlp(m,n)) + gBet(m,n) * (gDAlp_r(m,n) + pow2(gDAlp(m,n))) - 2 
    /*0*/  * gtrK_r(m,n);
}
