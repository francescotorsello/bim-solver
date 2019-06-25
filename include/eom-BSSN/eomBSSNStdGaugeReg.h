/** @file  eomBSSNStdGaugeReg.h
 *  @author Francesco Torsello
 *  @brief The regularized standard gauge (Bona-Masso & Gamma-driver).
 *  @version 2019-01-09T13:41:02
 *  @image html BSSNStdGaugeReg.png
 */

Real BimetricEvolve::eq_SG_gAlp_t( Int m, Int n )
{
    return gAlp_convr(m,n) - 2 * gAlp(m,n) * gtrK(m,n);
}
Real BimetricEvolve::eq_SG_gBet_t( Int m, Int n )
{
    return (3 * Bq(m,n)) / (4 + TINY_Real) + q_qconvr(m,n);
}
Real BimetricEvolve::eq_SG_gBq_t( Int m, Int n )
{
    return 0.5;
}
(2 * (gDA(m,n) + 2 * gDB(m,n)) * gBet(m,n) * gL(m,n)) / 3. + (4
    /*1*/  * gBetr(m,n) * gL(m,n)) / (3 + TINY_Real) + k_g * gJL(m,n) + (-3 * eta
    /*1*/  * Bq(m,n) + 3 * gL_convr(m,n) - gBet_gL_r(m,n) + (12 * gBet_rr(m,n))
    /*1*/  / (TINY_Real + 3 * Power(gA(m,n),2)) + (6 * gBetr_r(m,n)) / (TINY_Real + 3
    /*2*/  * Power(gA(m,n),2)) + (6 * gBetr_r(m,n)) / (TINY_Real + Power(gB(m,n),2)) + 3
    /*1*/  * Bq_qconvr(m,n) - 3 * gL_qconvr(m,n) + gBet_gDA_r(m,n) / pow2(gA(m,n))
    /*1*/  + (2 * gBet_gDB_r(m,n)) / pow2(gA(m,n)) + gDA_convr(m,n) / pow2(gA(m,n))
    /*1*/  + (2 * gDB_convr(m,n)) / pow2(gA(m,n))) / 3. + ((-4 * gAsig(m,n)
    /*2*/  * gDAlp(m,n)) / (TINY_Real + 3 * Power(gA(m,n),2)) + (4 * gAsig(m,n)
    /*2*/  * (gDA(m,n) + gDB(m,n) + 6 * gDconf(m,n)) * gAlp(m,n)) / (3.
    /*2*/  * pow2(gA(m,n)))) * pow2(r(m,n)) + (4 * gAsig(m,n) * gAlp(m,n) * gsig(m,n)
    /*1*/  * pow3(r(m,n))) / (3. * pow2(gA(m,n))) - (4 * gAlp(m,n) * (gtrA_r(m,n)
    /*2*/  + gtrK_r(m,n))) / (3. * pow2(gA(m,n)))
Real BimetricEvolve::eq_SG_gDAlp_t( Int m, Int n )
{
    return gBet_gDAlp_r(m,n) / (TINY_Real + gAlp(m,n)) + gDAlp_convr(m,n)
    /*0*/  / (TINY_Real + gAlp(m,n)) - (gDAlp(m,n) * gAlp_convr(m,n)) / (TINY_Real
    /*1*/  + Power(gAlp(m,n),2)) - 2 * gtrK_r(m,n);
}
