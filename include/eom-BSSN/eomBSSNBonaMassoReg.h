/** @file  eomBSSNBonaMassoReg.h
 *  @author Francesco Torsello
 *  @brief The regularized Bona-Masso slicing condition.
 *  @version 2018-10-26T15:57:28
 *  @image html BSSNBonaMassoReg.png
 */

Real BimetricEvolve::eq_BM_gAlp_t( Int m, Int n )
{
    return gAlp(m,n) * (gDAlp(m,n) * gBet(m,n) - 2 * gtrK(m,n));
}
Real BimetricEvolve::eq_BM_gBet_t( Int m, Int n )
{
    return (3 * Bq(m,n)) / (4 + TINY_Real) + gBet(m,n) * gBet_r(m,n);
}
Real BimetricEvolve::eq_BM_gBq_t( Int m, Int n )
{
    return -(eta * Bq(m,n)) + Bq_r(m,n) * gBet(m,n) + (4 * gBetr(m,n) * gL(m,n)) / (3
    /*1*/  + 2 * TINY_Real) + (2 * (2 * gBet_rr(m,n) + gBetr_r(m,n))) / (TINY_Real 
    /*1*/  + (3 + TINY_Real) * Power(gA(m,n),2)) + (2 * gBetr_r(m,n)) / (TINY_Real 
    /*1*/  + Power(gB(m,n),2)) + k_g * gJL(m,n) + gBet_r(m,n) * (-gL(m,n) / 3. 
    /*1*/  + gDA(m,n) / (3. * pow2(gA(m,n))) + (2 * gDB(m,n)) / ((3 + TINY_Real) 
    /*2*/  * pow2(gA(m,n)))) + gBet(m,n) * (((2 * gDA(m,n)) / 3. + (4 * gDB(m,n)) / (3 
    /*3*/  + TINY_Real)) * gL(m,n) + (gDA_r(m,n) + 2 * gDB_r(m,n)) / ((3 + TINY_Real) 
    /*2*/  * pow2(gA(m,n)))) + gAlp(m,n) * (gAsig(m,n) * ((4 * gDA(m,n) * pow2(r(m,n)))
    /*2*/  / (3. * pow2(gA(m,n))) + (4 * gDB(m,n) * pow2(r(m,n))) / ((3 + TINY_Real) 
    /*3*/  * pow2(gA(m,n))) + (8 * gconf(m,n) * gDconf(m,n) * pow2(r(m,n))) 
    /*2*/  / pow2(gA(m,n)) - (4 * gDAlp(m,n) * pow2(r(m,n))) / (3. * pow2(gA(m,n)))) 
    /*1*/  + (4 * gAsig(m,n) * gsig(m,n) * pow3(r(m,n))) / (3. * pow2(gA(m,n))) - (4 
    /*2*/  * (gtrA_r(m,n) + gtrK_r(m,n))) / ((3 + TINY_Real) * pow2(gA(m,n))));
}
Real BimetricEvolve::eq_BM_gDAlp_t( Int m, Int n )
{
    return gDAlp_r(m,n) * gBet(m,n) + gDAlp(m,n) * gBet_r(m,n) - 2 * gtrK_r(m,n);
}
