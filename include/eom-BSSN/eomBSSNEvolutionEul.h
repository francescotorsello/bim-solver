/** @file  eomBSSNEvolutionEul.h
 *  @author Francesco Torsello
 *  @brief The regularized and simplified cBSSN Eulerian evolution equations.
 *  @version 2019-04-04T10:43:39
 *  @image html BSSNevolutionEul.png
 */

Real BimetricEvolve::eq_gconf_t( Int m, Int n )
{
    return gconf_convr(m,n) - (gAlp(m,n) * gtrK(m,n)) / (6 + TINY_Real);
}
Real BimetricEvolve::eq_fconf_t( Int m, Int n )
{
    return fconf_convr(m,n) - (fAlp(m,n) * ftrK(m,n)) / (6 + TINY_Real);
}
Real BimetricEvolve::eq_gDconf_t( Int m, Int n )
{
    return gBet_gDconf_r(m,n) + gDconf_convr(m,n) - (gDAlp(m,n) * gtrK(m,n)) / (6
    /*1*/  + TINY_Real) - (gAlp(m,n) * gtrK_r(m,n)) / (6 + TINY_Real);
}
Real BimetricEvolve::eq_fDconf_t( Int m, Int n )
{
    return fBet_fDconf_r(m,n) + fDconf_convr(m,n) - (fDAlp(m,n) * ftrK(m,n)) / (6
    /*1*/  + TINY_Real) - (fAlp(m,n) * ftrK_r(m,n)) / (6 + TINY_Real);
}
Real BimetricEvolve::eq_gtrK_t( Int m, Int n )
{
    return (3 * gDers(m,n)) / (TINY_Real + 2 * exp(4 * gconf(m,n))) + gtrK_convr(m,n)
    /*0*/  - (3 * (gDB(m,n) * gDAlp(m,n) + 2 * gDconf(m,n) * gDAlp(m,n)
    /*2*/  + gDAlp_r(m,n))) / (TINY_Real + exp(4 * gconf(m,n)) * Power(gA(m,n),2)) + k_g
    /*0*/  * gJK(m,n) + gAlp(m,n) * ((3 * pow2(gA1(m,n))) / 2. + pow2(gtrK(m,n)) / 3.)
    /*0*/  + (2 * gDAlpr_r(m,n) * r(m,n)) / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2));
}
Real BimetricEvolve::eq_ftrK_t( Int m, Int n )
{
    return (3 * fDers(m,n)) / (TINY_Real + 2 * exp(4 * fconf(m,n))) + ftrK_convr(m,n)
    /*0*/  - (3 * (fDB(m,n) * fDAlp(m,n) + 2 * fDconf(m,n) * fDAlp(m,n)
    /*2*/  + fDAlp_r(m,n))) / (TINY_Real + exp(4 * fconf(m,n)) * Power(fA(m,n),2)) + k_f
    /*0*/  * fJK(m,n) + fAlp(m,n) * ((3 * pow2(fA1(m,n))) / 2. + pow2(ftrK(m,n)) / 3.)
    /*0*/  + (2 * fDAlpr_r(m,n) * r(m,n)) / (TINY_Real + exp(4 * fconf(m,n))
    /*1*/  * Power(fA(m,n),2));
}
Real BimetricEvolve::eq_gA_t( Int m, Int n )
{
    return gA_convr(m,n) - gA1(m,n) * gAlp(m,n) * gA(m,n) + gBet_r(m,n) * gA(m,n);
}
Real BimetricEvolve::eq_gB_t( Int m, Int n )
{
    return gB_convr(m,n) + (gA1(m,n) * gAlp(m,n) * gB(m,n)) / (2 + TINY_Real)
    /*0*/  + gBetr(m,n) * gB(m,n);
}
Real BimetricEvolve::eq_gDA_t( Int m, Int n )
{
    return gBet_gDA_r(m,n) - gA1(m,n) * gDAlp(m,n) - gA1_r(m,n) * gAlp(m,n)
    /*0*/  + gBet_rr(m,n) + gDA_convr(m,n);
}
Real BimetricEvolve::eq_gDB_t( Int m, Int n )
{
    return gBet_gDB_r(m,n) + (gA1(m,n) * gDAlp(m,n)) / (2 + TINY_Real) + (gA1_r(m,n)
    /*1*/  * gAlp(m,n)) / (2 + TINY_Real) + gDB_convr(m,n) + gBetr_r(m,n);
}
Real BimetricEvolve::eq_fA_t( Int m, Int n )
{
    return fA_convr(m,n) - fA1(m,n) * fAlp(m,n) * fA(m,n) + fBet_r(m,n) * fA(m,n);
}
Real BimetricEvolve::eq_fB_t( Int m, Int n )
{
    return fB_convr(m,n) + (fA1(m,n) * fAlp(m,n) * fB(m,n)) / (2 + TINY_Real)
    /*0*/  + fBetr(m,n) * fB(m,n);
}
Real BimetricEvolve::eq_fDA_t( Int m, Int n )
{
    return fBet_fDA_r(m,n) - fA1(m,n) * fDAlp(m,n) - fA1_r(m,n) * fAlp(m,n)
    /*0*/  + fBet_rr(m,n) + fDA_convr(m,n);
}
Real BimetricEvolve::eq_fDB_t( Int m, Int n )
{
    return fBet_fDB_r(m,n) + (fA1(m,n) * fDAlp(m,n)) / (2 + TINY_Real) + (fA1_r(m,n)
    /*1*/  * fAlp(m,n)) / (2 + TINY_Real) + fDB_convr(m,n) + fBetr_r(m,n);
}
Real BimetricEvolve::eq_gA1_t( Int m, Int n )
{
    return gDers(m,n) / (TINY_Real + exp(4 * gconf(m,n))) + gRicci(m,n)
    /*0*/  + gA1_convr(m,n) + k_g * gJA1(m,n) - (2 * gDAlpr_r(m,n) * r(m,n))
    /*0*/  / (TINY_Real + 3 * exp(4 * gconf(m,n)) * Power(gA(m,n),2)) + gAlp(m,n) * ((2
    /*2*/  * exp(-4 * gconf(m,n)) * gDconf(m,n) * gDers(m,n)) / ( gDAlp(m,n) + TINY_Real ) - (8
    /*2*/  * exp(-4 * gconf(m,n)) * pow2(gDconf(m,n))) / (3. * pow2(gA(m,n)) + TINY_Real ) - (4
    /*2*/  * exp(-4 * gconf(m,n)) * gDconfr_r(m,n) * r(m,n)) / (3. * pow2(gA(m,n)) + TINY_Real )
    /*1*/  + gA1(m,n) * gtrK(m,n));
}
Real BimetricEvolve::eq_fA1_t( Int m, Int n )
{
    return fDers(m,n) / (TINY_Real + exp(4 * fconf(m,n))) + fRicci(m,n)
    /*0*/  + fA1_convr(m,n) + k_f * fJA1(m,n) - (2 * fDAlpr_r(m,n) * r(m,n))
    /*0*/  / (TINY_Real + 3 * exp(4 * fconf(m,n)) * Power(fA(m,n),2)) + fAlp(m,n) * ((2
    /*2*/  * exp(-4 * fconf(m,n)) * fDconf(m,n) * fDers(m,n)) / ( fDAlp(m,n) + TINY_Real ) - (8
    /*2*/  * exp(-4 * fconf(m,n)) * pow2(fDconf(m,n))) / (3. * pow2(fA(m,n)) + TINY_Real ) - (4
    /*2*/  * exp(-4 * fconf(m,n)) * fDconfr_r(m,n) * r(m,n)) / (3. * pow2(fA(m,n)) + TINY_Real )
    /*1*/  + fA1(m,n) * ftrK(m,n));
}
Real BimetricEvolve::eq_gL_t( Int m, Int n )
{
    return gL_convr(m,n) - gBet_gL_r(m,n) - (2 * gA1(m,n) * gDAlp(m,n)) / (TINY_Real
    /*1*/  + Power(gA(m,n),2)) + gBet_rr(m,n) / (TINY_Real + Power(gA(m,n),2)) + (2
    /*1*/  * gBetr_r(m,n)) / (TINY_Real + Power(gB(m,n),2)) + k_g * gJL(m,n) + gAlp(m,n)
    /*0*/  * ((4 * gAsig(m,n) * (gDA(m,n) + gDB(m,n) + 6 * gDconf(m,n))
    /*2*/  * pow2(r(m,n))) / (3. * pow2(gA(m,n)) + TINY_Real ) + (4 * gAsig(m,n) * gsig(m,n)
    /*2*/  * pow3(r(m,n))) / (3. * pow2(gA(m,n)) + TINY_Real ) - (4 * gtrK_r(m,n)) / (3.
    /*2*/  * pow2(gA(m,n)) + TINY_Real ));
}
Real BimetricEvolve::eq_fL_t( Int m, Int n )
{
    return fL_convr(m,n) - fBet_fL_r(m,n) - (2 * fA1(m,n) * fDAlp(m,n)) / (TINY_Real
    /*1*/  + Power(fA(m,n),2)) + fBet_rr(m,n) / (TINY_Real + Power(fA(m,n),2)) + (2
    /*1*/  * fBetr_r(m,n)) / (TINY_Real + Power(fB(m,n),2)) + k_f * fJL(m,n) + fAlp(m,n)
    /*0*/  * ((4 * fAsig(m,n) * (fDA(m,n) + fDB(m,n) + 6 * fDconf(m,n))
    /*2*/  * pow2(r(m,n))) / (3. * pow2(fA(m,n)) + TINY_Real ) + (4 * fAsig(m,n) * fsig(m,n)
    /*2*/  * pow3(r(m,n))) / (3. * pow2(fA(m,n)) + TINY_Real ) - (4 * ftrK_r(m,n)) / (3.
    /*2*/  * pow2(fA(m,n)) + TINY_Real ));
}
Real BimetricEvolve::eq_gsig_t( Int m, Int n )
{
    return gsig_convr(m,n) + 2 * gBetr(m,n) * gsig(m,n) + (2 * gAsig(m,n) * gAlp(m,n)
    /*1*/  * pow2(gA(m,n))) / (TINY_Real + Power(gB(m,n),2)) - (2 * gBetr_r(m,n)
    /*1*/  * pow2(gA(m,n))) / (TINY_Real + Power(gB(m,n),2) * r(m,n));
}
Real BimetricEvolve::eq_fsig_t( Int m, Int n )
{
    return fsig_convr(m,n) + 2 * fBetr(m,n) * fsig(m,n) + (2 * fAsig(m,n) * fAlp(m,n)
    /*1*/  * pow2(fA(m,n))) / (TINY_Real + Power(fB(m,n),2)) - (2 * fBetr_r(m,n)
    /*1*/  * pow2(fA(m,n))) / (TINY_Real + Power(fB(m,n),2) * r(m,n));
}
Real BimetricEvolve::eq_gAsig_t( Int m, Int n )
{
    return gAsig_convr(m,n) + 2 * gAsig(m,n) * gBetr(m,n) - gDAlpr_r(m,n)
    /*0*/  / (TINY_Real + exp(4 * gconf(m,n)) * Power(gA(m,n),2) * r(m,n)) + (3 * k_g
    /*1*/  * gJA1(m,n)) / (TINY_Real + 2 * Power(r(m,n),2)) + (3 * gDers(m,n))
    /*0*/  / (TINY_Real + 2 * exp(4 * gconf(m,n)) * Power(r(m,n),2)) + gAlp(m,n)
    /*0*/  * (exp(-4 * gconf(m,n)) * ((gsig_rr(m,n) * pow2(gB(m,n))) / (2.
    /*3*/  * Power(gA(m,n),4) + TINY_Real ) + (pow2(gsig(m,n)) * pow2(gB(m,n))) / ( Power(gA(m,n),4) + TINY_Real )
    /*2*/  - (gsig_gL_r(m,n) * pow2(gB(m,n))) / (2. * pow2(gA(m,n)) + TINY_Real )) - (3 * exp(-4
    /*3*/  * gconf(m,n)) * (5 * gDB(m,n) + 11 * gDconf(m,n)) * gDers(m,n)) / (gDAlp(m,n)
    /*2*/  * pow2(r(m,n)) + TINY_Real ) + (exp(-4 * gconf(m,n)) * (40 * gDB(m,n) * gDconf(m,n) + 7
    /*3*/  * pow2(gDB(m,n)) + 44 * pow2(gDconf(m,n)))) / (pow2(gA(m,n))
    /*2*/  * pow2(r(m,n)) + TINY_Real ) + (27 * exp(-4 * gconf(m,n)) * pow2(gDers(m,n))
    /*2*/  * pow2(gA(m,n))) / (4. * pow2(gDAlp(m,n)) * pow2(r(m,n)) + TINY_Real ) + (exp(-4
    /*3*/  * gconf(m,n)) * (gLr_r(m,n) + (-2 * gDconfr_r(m,n) + gsig_r(m,n))
    /*3*/  / ( pow2(gA(m,n)) + TINY_Real ) + (2 * gsig_r(m,n) * pow2(gB(m,n))) / ( Power(gA(m,n),4) + TINY_Real )
    /*3*/  + (gsig(m,n) * (-4 * gDB(m,n) - gL(m,n) * pow2(gB(m,n)))) / ( pow2(gA(m,n)) + TINY_Real )   ))
    /*1*/  / ( r(m,n) + TINY_Real ) + gAsig(m,n) * gtrK(m,n));
}
Real BimetricEvolve::eq_fAsig_t( Int m, Int n )
{
    return fAsig_convr(m,n) + 2 * fAsig(m,n) * fBetr(m,n) - fDAlpr_r(m,n)
    /*0*/  / (TINY_Real + exp(4 * fconf(m,n)) * Power(fA(m,n),2) * r(m,n)) + (3 * k_f
    /*1*/  * fJA1(m,n)) / (TINY_Real + 2 * Power(r(m,n),2)) + (3 * fDers(m,n))
    /*0*/  / (TINY_Real + 2 * exp(4 * fconf(m,n)) * Power(r(m,n),2)) + fAlp(m,n)
    /*0*/  * (exp(-4 * fconf(m,n)) * ((fsig_rr(m,n) * pow2(fB(m,n))) / (2.
    /*3*/  * Power(fA(m,n),4) + TINY_Real ) + (pow2(fsig(m,n)) * pow2(fB(m,n))) / ( Power(fA(m,n),4)  + TINY_Real )
    /*2*/  - (fsig_fL_r(m,n) * pow2(fB(m,n))) / (2. * pow2(fA(m,n)) + TINY_Real )) - (3 * exp(-4
    /*3*/  * fconf(m,n)) * (5 * fDB(m,n) + 11 * fDconf(m,n)) * fDers(m,n)) / (fDAlp(m,n)
    /*2*/  * pow2(r(m,n)) + TINY_Real ) + (exp(-4 * fconf(m,n)) * (40 * fDB(m,n) * fDconf(m,n) + 7
    /*3*/  * pow2(fDB(m,n)) + 44 * pow2(fDconf(m,n)))) / (pow2(fA(m,n))
    /*2*/  * pow2(r(m,n)) + TINY_Real ) + (27 * exp(-4 * fconf(m,n)) * pow2(fDers(m,n))
    /*2*/  * pow2(fA(m,n))) / (4. * pow2(fDAlp(m,n)) * pow2(r(m,n)) + TINY_Real ) + (exp(-4
    /*3*/  * fconf(m,n)) * (fLr_r(m,n) + (-2 * fDconfr_r(m,n) + fsig_r(m,n))
    /*3*/  / ( pow2(fA(m,n)) + TINY_Real ) + (2 * fsig_r(m,n) * pow2(fB(m,n))) / ( Power(fA(m,n),4) + TINY_Real )
    /*3*/  + (fsig(m,n) * (-4 * fDB(m,n) - fL(m,n) * pow2(fB(m,n)))) / ( pow2(fA(m,n))  + TINY_Real )     ))
    /*1*/  / ( r(m,n) + TINY_Real ) + fAsig(m,n) * ftrK(m,n));
}
