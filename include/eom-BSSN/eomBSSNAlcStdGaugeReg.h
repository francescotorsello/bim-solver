/** @file  eomBSSNAlcStdGaugeReg.h
 *  @author Francesco Torsello
 *  @brief Alcubierre's standard gauge [Alcubierre, Introduction to 3+1 numerical
           relativity (2008) p.148].
 *  @version 2019-05-28T13:15:12
 *  @image html BSSNStdGaugeReg.png
 */

Real BimetricEvolve::eq_SG_gAlp_t( Int m, Int n )
{
    return gAlp_convr(m,n) - 2 * gAlp(m,n) * gtrK(m,n);
}
Real BimetricEvolve::eq_SG_gBet_t( Int m, Int n )
{
    return (3 * Bq(m,n)) / 4. + q_qconvr(m,n);
}
Real BimetricEvolve::eq_SG_gBq_t( Int m, Int n )
{
    return -(eta * Bq(m,n)) + Bq_qconvr(m,n) + k_g * gJL(m,n) * pow2(gAlp(m,n))
    /*0*/  + pow2(gAlp(m,n)) * ((2 * (gDA(m,n) + 2 * gDB(m,n)) * gBet(m,n) * gL(m,n))
    /*1*/  / 3. + (4 * gBetr(m,n) * gL(m,n)) / 3. + (gBet_gDA_r(m,n) + 2
    /*2*/  * gBet_gDB_r(m,n) + 4 * gBet_rr(m,n) + gDA_convr(m,n) + 2 * gDB_convr(m,n)
    /*2*/  + 3 * gL_convr(m,n) * pow2(gA(m,n)) - gBet_gL_r(m,n) * pow2(gA(m,n)) - 3
    /*2*/  * gL_qconvr(m,n) * pow2(gA(m,n))) / (3. * pow2(gA(m,n))) + gBetr_r(m,n) * (2
    /*2*/  / (3. * pow2(gA(m,n))) + 2 / pow2(gB(m,n)))) + pow2(r(m,n)) * ((-4
    /*2*/  * gAsig(m,n) * gDAlp(m,n) * pow2(gAlp(m,n))) / (3. * pow2(gA(m,n))) + (4
    /*2*/  * gAsig(m,n) * (gDA(m,n) + gDB(m,n) + 6 * gDconf(m,n)) * pow3(gAlp(m,n)))
    /*1*/  / (3. * pow2(gA(m,n)))) + (4 * gAsig(m,n) * gsig(m,n) * pow3(gAlp(m,n))
    /*1*/  * pow3(r(m,n))) / (3. * pow2(gA(m,n))) - (4 * pow3(gAlp(m,n)) * (gtrA_r(m,n)
    /*2*/  + gtrK_r(m,n))) / (3. * pow2(gA(m,n)));
}
Real BimetricEvolve::eq_SG_gDAlp_t( Int m, Int n )
{
    return (gBet_gDAlp_r(m,n) + gDAlp_convr(m,n)) / gAlp(m,n) - (gDAlp(m,n)
    /*1*/  * gAlp_convr(m,n)) / pow2(gAlp(m,n)) - 2 * gtrK_r(m,n);
}
