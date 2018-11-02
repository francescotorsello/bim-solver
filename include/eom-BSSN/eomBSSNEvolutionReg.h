/** @file  eomBSSNEvolutionReg.h
 *  @author Francesco Torsello
 *  @brief The double-regularized cBSSN evolution equations.
 *  @version 2018-10-16T11:10:58
 *  @image html BSSNevolutionReg.png
 */

Real BimetricEvolve::eq_gconf_t( Int m, Int n )
{
    return -gdet_pff(m,n) / (12. * gdet(m,n)) + gconf(m,n) * gDconf(m,n) * gBet(m,n)
    /*0*/  - (gAlp(m,n) * (gtrA(m,n) + gtrK(m,n))) / 6.;
}
Real BimetricEvolve::eq_fconf_t( Int m, Int n )
{
    return -fdet_pff(m,n) / (12. * fdet(m,n)) + fconf(m,n) * fDconf(m,n) * fBet(m,n)
    /*0*/  - (fAlp(m,n) * (ftrA(m,n) + ftrK(m,n))) / 6.;
}
Real BimetricEvolve::eq_gDconf_t( Int m, Int n )
{
    return -(gdet_pff_r(m,n) / (gdet(m,n) * (TINY_Real + 12 * gconf(m,n))))
    /*0*/  + (gdet_pff(m,n) * gdet_r(m,n)) / ((TINY_Real + 12 * gconf(m,n))
    /*1*/  * pow2(gdet(m,n))) - (gDAlp(m,n) * gAlp(m,n) * (gtrA(m,n) + gtrK(m,n)))
    /*0*/  / (TINY_Real + 6 * gconf(m,n)) + gDconf(m,n) * (gdet_pff(m,n) / (gdet(m,n)
    /*2*/  * (TINY_Real + 12 * gconf(m,n))) + (6 * gconf(m,n) * gBet_r(m,n) + gAlp(m,n)
    /*2*/  * (gtrA(m,n) + gtrK(m,n))) / (TINY_Real + 6 * gconf(m,n))) + (12
    /*1*/  * gconf(m,n) * gDconf_r(m,n) * gBet(m,n) - 2 * gAlp(m,n) * (gtrA_r(m,n)
    /*2*/  + gtrK_r(m,n))) / (TINY_Real + 12 * gconf(m,n));
}
Real BimetricEvolve::eq_fDconf_t( Int m, Int n )
{
    return -(fdet_pff_r(m,n) / (fdet(m,n) * (TINY_Real + 12 * fconf(m,n))))
    /*0*/  + (fdet_pff(m,n) * fdet_r(m,n)) / ((TINY_Real + 12 * fconf(m,n))
    /*1*/  * pow2(fdet(m,n))) - (fDAlp(m,n) * fAlp(m,n) * (ftrA(m,n) + ftrK(m,n)))
    /*0*/  / (TINY_Real + 6 * fconf(m,n)) + fDconf(m,n) * (fdet_pff(m,n) / (fdet(m,n)
    /*2*/  * (TINY_Real + 12 * fconf(m,n))) + (6 * fconf(m,n) * fBet_r(m,n) + fAlp(m,n)
    /*2*/  * (ftrA(m,n) + ftrK(m,n))) / (TINY_Real + 6 * fconf(m,n))) + (12
    /*1*/  * fconf(m,n) * fDconf_r(m,n) * fBet(m,n) - 2 * fAlp(m,n) * (ftrA_r(m,n)
    /*2*/  + ftrK_r(m,n))) / (TINY_Real + 12 * fconf(m,n));
}
Real BimetricEvolve::eq_gtrK_t( Int m, Int n )
{
    return k_g * gJK(m,n) - gtrA_pff(m,n) + gAlp(m,n) * (exp(-4 * gconf(m,n))
    /*1*/  * ((gDA(m,n) * gDAlp(m,n)) / pow2(gA(m,n)) - (2 * gDB(m,n) * gDAlp(m,n))
    /*2*/  / pow2(gA(m,n)) - gDAlp_r(m,n) / pow2(gA(m,n)) - pow2(gDAlp(m,n))
    /*2*/  / pow2(gA(m,n)) + gDAlp(m,n) * ((-2 * gconf(m,n) * gDconf(m,n))
    /*3*/  / pow2(gA(m,n)) - 2 / (pow2(gA(m,n)) * r(m,n)))) + (3 * pow2(gA1(m,n)) + 6
    /*2*/  * pow2(gA2(m,n)) + gtrK(m,n) * (2 * gtrA(m,n) + gtrK(m,n))) / 3.) + gBet(m,n)
    /*0*/  * gtrK_r(m,n);
}
Real BimetricEvolve::eq_ftrK_t( Int m, Int n )
{
    return k_f * fJK(m,n) - ftrA_pff(m,n) + fAlp(m,n) * (exp(-4 * fconf(m,n))
    /*1*/  * ((fDA(m,n) * fDAlp(m,n)) / pow2(fA(m,n)) - (2 * fDB(m,n) * fDAlp(m,n))
    /*2*/  / pow2(fA(m,n)) - fDAlp_r(m,n) / pow2(fA(m,n)) - pow2(fDAlp(m,n))
    /*2*/  / pow2(fA(m,n)) + fDAlp(m,n) * ((-2 * fconf(m,n) * fDconf(m,n))
    /*3*/  / pow2(fA(m,n)) - 2 / (pow2(fA(m,n)) * r(m,n)))) + (3 * pow2(fA1(m,n)) + 6
    /*2*/  * pow2(fA2(m,n)) + ftrK(m,n) * (2 * ftrA(m,n) + ftrK(m,n))) / 3.) + fBet(m,n)
    /*0*/  * ftrK_r(m,n);
}
Real BimetricEvolve::eq_gA_t( Int m, Int n )
{
    return (gdet_pff(m,n) * gA(m,n)) / (6. * gdet(m,n)) + (gDA(m,n) * gBet(m,n)
    /*1*/  + gBet_r(m,n)) * gA(m,n) + (gAlp(m,n) * gA(m,n) * (-3 * gA1(m,n) + gtrA(m,n)
    /*2*/  + gtrAv(m,n))) / 3.;
}
Real BimetricEvolve::eq_gB_t( Int m, Int n )
{
    return (gdet_pff(m,n) * gB(m,n)) / (6. * gdet(m,n)) + (gDB(m,n) * gBet(m,n)
    /*1*/  + gBetr(m,n)) * gB(m,n) + (gAlp(m,n) * gB(m,n) * (-3 * gA2(m,n) + gtrA(m,n)
    /*2*/  + gtrAv(m,n))) / 3.;
}
Real BimetricEvolve::eq_gDA_t( Int m, Int n )
{
    return gdet_pff_r(m,n) / (6. * gdet(m,n)) + gDA_r(m,n) * gBet(m,n) + gDA(m,n)
    /*0*/  * gBet_r(m,n) + gBet_rr(m,n) - (gdet_pff(m,n) * gdet_r(m,n)) / (6.
    /*1*/  * pow2(gdet(m,n))) + (gDAlp(m,n) * gAlp(m,n) * (-3 * gA1(m,n) + gtrA(m,n)))
    /*0*/  / 3. + (gAlp(m,n) * (-3 * gA1_r(m,n) + gtrA_r(m,n))) / 3.;
}
Real BimetricEvolve::eq_gDB_t( Int m, Int n )
{
    return gdet_pff_r(m,n) / (6. * gdet(m,n)) + gDB_r(m,n) * gBet(m,n) + gDB(m,n)
    /*0*/  * gBet_r(m,n) + gBetr_r(m,n) - (gdet_pff(m,n) * gdet_r(m,n)) / (6.
    /*1*/  * pow2(gdet(m,n))) + (gDAlp(m,n) * gAlp(m,n) * (-3 * gA2(m,n) + gtrA(m,n)))
    /*0*/  / 3. + (gAlp(m,n) * (-3 * gA2_r(m,n) + gtrA_r(m,n))) / 3.;
}
Real BimetricEvolve::eq_fA_t( Int m, Int n )
{
    return (fdet_pff(m,n) * fA(m,n)) / (6. * fdet(m,n)) + (fDA(m,n) * fBet(m,n)
    /*1*/  + fBet_r(m,n)) * fA(m,n) + (fAlp(m,n) * fA(m,n) * (-3 * fA1(m,n) + ftrA(m,n)
    /*2*/  + ftrAv(m,n))) / 3.;
}
Real BimetricEvolve::eq_fB_t( Int m, Int n )
{
    return (fdet_pff(m,n) * fB(m,n)) / (6. * fdet(m,n)) + (fDB(m,n) * fBet(m,n)
    /*1*/  + fBetr(m,n)) * fB(m,n) + (fAlp(m,n) * fB(m,n) * (-3 * fA2(m,n) + ftrA(m,n)
    /*2*/  + ftrAv(m,n))) / 3.;
}
Real BimetricEvolve::eq_fDA_t( Int m, Int n )
{
    return fdet_pff_r(m,n) / (6. * fdet(m,n)) + fDA_r(m,n) * fBet(m,n) + fDA(m,n)
    /*0*/  * fBet_r(m,n) + fBet_rr(m,n) - (fdet_pff(m,n) * fdet_r(m,n)) / (6.
    /*1*/  * pow2(fdet(m,n))) + (fDAlp(m,n) * fAlp(m,n) * (-3 * fA1(m,n) + ftrA(m,n)))
    /*0*/  / 3. + (fAlp(m,n) * (-3 * fA1_r(m,n) + ftrA_r(m,n))) / 3.;
}
Real BimetricEvolve::eq_fDB_t( Int m, Int n )
{
    return fdet_pff_r(m,n) / (6. * fdet(m,n)) + fDB_r(m,n) * fBet(m,n) + fDB(m,n)
    /*0*/  * fBet_r(m,n) + fBetr_r(m,n) - (fdet_pff(m,n) * fdet_r(m,n)) / (6.
    /*1*/  * pow2(fdet(m,n))) + (fDAlp(m,n) * fAlp(m,n) * (-3 * fA2(m,n) + ftrA(m,n)))
    /*0*/  / 3. + (fAlp(m,n) * (-3 * fA2_r(m,n) + ftrA_r(m,n))) / 3.;
}
Real BimetricEvolve::eq_gA1_t( Int m, Int n )
{
    return gRicci(m,n) + gA1_r(m,n) * gBet(m,n) + k_g * gJA1(m,n) + gtrA_pff(m,n) / 3.
    /*0*/  + gAlp(m,n) * (exp(-4 * gconf(m,n)) * (gDA(m,n) * ((4 * gconf(m,n)
    /*4*/  * gDconf(m,n)) / (3. * pow2(gA(m,n))) + (2 * gDAlp(m,n)) / (3.
    /*4*/  * pow2(gA(m,n)))) + gDB(m,n) * ((4 * gconf(m,n) * gDconf(m,n)) / (3.
    /*4*/  * pow2(gA(m,n))) + (2 * gDAlp(m,n)) / (3. * pow2(gA(m,n)))) - (2 * (2
    /*4*/  * gconf(m,n) * gDconf_r(m,n) + gDAlp_r(m,n))) / (3. * pow2(gA(m,n))) + (4
    /*3*/  * gconf(m,n) * (-1 + 2 * gconf(m,n)) * pow2(gDconf(m,n))) / (3.
    /*3*/  * pow2(gA(m,n))) - (2 * pow2(gDAlp(m,n))) / (3. * pow2(gA(m,n))) + gDAlp(m,n)
    /*2*/  * ((8 * gconf(m,n) * gDconf(m,n)) / (3. * pow2(gA(m,n))) + 2 / (3.
    /*4*/  * pow2(gA(m,n)) * r(m,n))) + (4 * gconf(m,n) * gDconf(m,n)) / (3.
    /*3*/  * pow2(gA(m,n)) * r(m,n))) + ((3 * gA1(m,n) - gtrA(m,n)) * (gtrA(m,n)
    /*3*/  + gtrK(m,n))) / 3.);
}
Real BimetricEvolve::eq_gA2_t( Int m, Int n )
{
    return -gRicci(m,n) / 2. + gA2_r(m,n) * gBet(m,n) + k_g * gJA2(m,n)
    /*0*/  + gtrA_pff(m,n) / 3. + gAlp(m,n) * (exp(-4 * gconf(m,n)) * (gDA(m,n) * ((-2
    /*4*/  * gconf(m,n) * gDconf(m,n)) / (3. * pow2(gA(m,n))) - gDAlp(m,n) / (3.
    /*4*/  * pow2(gA(m,n)))) + gDB(m,n) * ((-2 * gconf(m,n) * gDconf(m,n)) / (3.
    /*4*/  * pow2(gA(m,n))) - gDAlp(m,n) / (3. * pow2(gA(m,n)))) + (2 * gconf(m,n)
    /*3*/  * gDconf_r(m,n) + gDAlp_r(m,n)) / (3. * pow2(gA(m,n))) + (2 * (1 - 2
    /*4*/  * gconf(m,n)) * gconf(m,n) * pow2(gDconf(m,n))) / (3. * pow2(gA(m,n)))
    /*2*/  + pow2(gDAlp(m,n)) / (3. * pow2(gA(m,n))) + gDAlp(m,n) * ((-4 * gconf(m,n)
    /*4*/  * gDconf(m,n)) / (3. * pow2(gA(m,n))) - 1 / (3. * pow2(gA(m,n)) * r(m,n)))
    /*2*/  - (2 * gconf(m,n) * gDconf(m,n)) / (3. * pow2(gA(m,n)) * r(m,n))) + ((3
    /*3*/  * gA2(m,n) - gtrA(m,n)) * (gtrA(m,n) + gtrK(m,n))) / 3.);
}
Real BimetricEvolve::eq_fA1_t( Int m, Int n )
{
    return fRicci(m,n) + fA1_r(m,n) * fBet(m,n) + k_f * fJA1(m,n) + ftrA_pff(m,n) / 3.
    /*0*/  + fAlp(m,n) * (exp(-4 * fconf(m,n)) * (fDA(m,n) * ((4 * fconf(m,n)
    /*4*/  * fDconf(m,n)) / (3. * pow2(fA(m,n))) + (2 * fDAlp(m,n)) / (3.
    /*4*/  * pow2(fA(m,n)))) + fDB(m,n) * ((4 * fconf(m,n) * fDconf(m,n)) / (3.
    /*4*/  * pow2(fA(m,n))) + (2 * fDAlp(m,n)) / (3. * pow2(fA(m,n)))) - (2 * (2
    /*4*/  * fconf(m,n) * fDconf_r(m,n) + fDAlp_r(m,n))) / (3. * pow2(fA(m,n))) + (4
    /*3*/  * fconf(m,n) * (-1 + 2 * fconf(m,n)) * pow2(fDconf(m,n))) / (3.
    /*3*/  * pow2(fA(m,n))) - (2 * pow2(fDAlp(m,n))) / (3. * pow2(fA(m,n))) + fDAlp(m,n)
    /*2*/  * ((8 * fconf(m,n) * fDconf(m,n)) / (3. * pow2(fA(m,n))) + 2 / (3.
    /*4*/  * pow2(fA(m,n)) * r(m,n))) + (4 * fconf(m,n) * fDconf(m,n)) / (3.
    /*3*/  * pow2(fA(m,n)) * r(m,n))) + ((3 * fA1(m,n) - ftrA(m,n)) * (ftrA(m,n)
    /*3*/  + ftrK(m,n))) / 3.);
}
Real BimetricEvolve::eq_fA2_t( Int m, Int n )
{
    return -fRicci(m,n) / 2. + fA2_r(m,n) * fBet(m,n) + k_f * fJA2(m,n)
    /*0*/  + ftrA_pff(m,n) / 3. + fAlp(m,n) * (exp(-4 * fconf(m,n)) * (fDA(m,n) * ((-2
    /*4*/  * fconf(m,n) * fDconf(m,n)) / (3. * pow2(fA(m,n))) - fDAlp(m,n) / (3.
    /*4*/  * pow2(fA(m,n)))) + fDB(m,n) * ((-2 * fconf(m,n) * fDconf(m,n)) / (3.
    /*4*/  * pow2(fA(m,n))) - fDAlp(m,n) / (3. * pow2(fA(m,n)))) + (2 * fconf(m,n)
    /*3*/  * fDconf_r(m,n) + fDAlp_r(m,n)) / (3. * pow2(fA(m,n))) + (2 * (1 - 2
    /*4*/  * fconf(m,n)) * fconf(m,n) * pow2(fDconf(m,n))) / (3. * pow2(fA(m,n)))
    /*2*/  + pow2(fDAlp(m,n)) / (3. * pow2(fA(m,n))) + fDAlp(m,n) * ((-4 * fconf(m,n)
    /*4*/  * fDconf(m,n)) / (3. * pow2(fA(m,n))) - 1 / (3. * pow2(fA(m,n)) * r(m,n)))
    /*2*/  - (2 * fconf(m,n) * fDconf(m,n)) / (3. * pow2(fA(m,n)) * r(m,n))) + ((3
    /*3*/  * fA2(m,n) - ftrA(m,n)) * (ftrA(m,n) + ftrK(m,n))) / 3.);
}
Real BimetricEvolve::eq_gL_t( Int m, Int n )
{
    return -(gBet_r(m,n) * gL(m,n)) + gBet(m,n) * gL_r(m,n) + k_g * gJL(m,n)
    /*0*/  + (-(gdet_pff(m,n) * gL(m,n)) / 3. - gdet_pff_r(m,n) / (6. * pow2(gA(m,n))))
    /*0*/  / gdet(m,n) + gBet_rr(m,n) / pow2(gA(m,n)) + (gdet_pff(m,n) * gdet_r(m,n))
    /*0*/  / (6. * pow2(gdet(m,n)) * pow2(gA(m,n))) + (2 * gBetr_r(m,n))
    /*0*/  / pow2(gB(m,n)) + (4 * gAsig(m,n) * gDA(m,n) * gAlp(m,n) * pow2(r(m,n)))
    /*0*/  / (3. * pow2(gA(m,n))) + (4 * gAsig(m,n) * gDB(m,n) * gAlp(m,n)
    /*1*/  * pow2(r(m,n))) / (3. * pow2(gA(m,n))) + (8 * gAsig(m,n) * gconf(m,n)
    /*1*/  * gDconf(m,n) * gAlp(m,n) * pow2(r(m,n))) / pow2(gA(m,n)) + (4 * gAsig(m,n)
    /*1*/  * gAlp(m,n) * gsig(m,n) * pow3(r(m,n))) / (3. * pow2(gA(m,n))) - (2
    /*1*/  * gDAlp(m,n) * gAlp(m,n) * (3 * gA1(m,n) - gtrA(m,n))) / (3. * pow2(gA(m,n)))
    /*0*/  - (4 * gAlp(m,n) * (gtrA_r(m,n) + gtrK_r(m,n))) / (3. * pow2(gA(m,n)));
}
Real BimetricEvolve::eq_fL_t( Int m, Int n )
{
    return -(fBet_r(m,n) * fL(m,n)) + fBet(m,n) * fL_r(m,n) + k_f * fJL(m,n)
    /*0*/  + (-(fdet_pff(m,n) * fL(m,n)) / 3. - fdet_pff_r(m,n) / (6. * pow2(fA(m,n))))
    /*0*/  / fdet(m,n) + fBet_rr(m,n) / pow2(fA(m,n)) + (fdet_pff(m,n) * fdet_r(m,n))
    /*0*/  / (6. * pow2(fdet(m,n)) * pow2(fA(m,n))) + (2 * fBetr_r(m,n))
    /*0*/  / pow2(fB(m,n)) + (4 * fAsig(m,n) * fDA(m,n) * fAlp(m,n) * pow2(r(m,n)))
    /*0*/  / (3. * pow2(fA(m,n))) + (4 * fAsig(m,n) * fDB(m,n) * fAlp(m,n)
    /*1*/  * pow2(r(m,n))) / (3. * pow2(fA(m,n))) + (8 * fAsig(m,n) * fconf(m,n)
    /*1*/  * fDconf(m,n) * fAlp(m,n) * pow2(r(m,n))) / pow2(fA(m,n)) + (4 * fAsig(m,n)
    /*1*/  * fAlp(m,n) * fsig(m,n) * pow3(r(m,n))) / (3. * pow2(fA(m,n))) - (2
    /*1*/  * fDAlp(m,n) * fAlp(m,n) * (3 * fA1(m,n) - ftrA(m,n))) / (3. * pow2(fA(m,n)))
    /*0*/  - (4 * fAlp(m,n) * (ftrA_r(m,n) + ftrK_r(m,n))) / (3. * pow2(fA(m,n)));
}
Real BimetricEvolve::eq_gsig_t( Int m, Int n )
{
    return 2 * gBetr(m,n) * gsig(m,n) + gBet(m,n) * gsig_r(m,n) + (pow2(gA(m,n)) * (2
    /*2*/  * gAsig(m,n) * gAlp(m,n) - (2 * gBetr_r(m,n)) / r(m,n))) / pow2(gB(m,n));
}
Real BimetricEvolve::eq_fsig_t( Int m, Int n )
{
    return 2 * fBetr(m,n) * fsig(m,n) + fBet(m,n) * fsig_r(m,n) + (pow2(fA(m,n)) * (2
    /*2*/  * fAsig(m,n) * fAlp(m,n) - (2 * fBetr_r(m,n)) / r(m,n))) / pow2(fB(m,n));
}
Real BimetricEvolve::eq_gAsig_t( Int m, Int n )
{
    return gAsig_r(m,n) * gBet(m,n) + 2 * gAsig(m,n) * gBetr(m,n) + (k_g * (gJA1(m,n)
    /*2*/  - gJA2(m,n))) / pow2(r(m,n)) - (exp(-4 * gconf(m,n)) * gDAlpr_r(m,n))
    /*0*/  / (pow2(gA(m,n)) * r(m,n)) + gAlp(m,n) * (exp(-4 * gconf(m,n))
    /*1*/  * ((gsig_rr(m,n) * pow2(gB(m,n))) / (2. * Power(gA(m,n),4))
    /*2*/  + (pow2(gsig(m,n)) * pow2(gB(m,n))) / Power(gA(m,n),4) + gDB(m,n) * ((2
    /*4*/  * gconf(m,n) * gDconf(m,n)) / (pow2(gA(m,n)) * pow2(r(m,n))) + gDAlp(m,n)
    /*3*/  / (pow2(gA(m,n)) * pow2(r(m,n)))) + gDA(m,n) * ((-4 * gDB(m,n))
    /*3*/  / (pow2(gA(m,n)) * pow2(r(m,n))) + (2 * gconf(m,n) * gDconf(m,n))
    /*3*/  / (pow2(gA(m,n)) * pow2(r(m,n))) + gDAlp(m,n) / (pow2(gA(m,n))
    /*4*/  * pow2(r(m,n)))) + (4 * gconf(m,n) * gDconf(m,n) * gDAlp(m,n))
    /*2*/  / (pow2(gA(m,n)) * pow2(r(m,n))) + (3 * pow2(gDA(m,n))) / (pow2(gA(m,n))
    /*3*/  * pow2(r(m,n))) + (4 * pow2(gconf(m,n)) * pow2(gDconf(m,n))) / (pow2(gA(m,n))
    /*3*/  * pow2(r(m,n))) + gsig(m,n) * ((-4 * gDB(m,n)) / (pow2(gA(m,n)) * r(m,n))
    /*3*/  - (gL(m,n) * pow2(gB(m,n))) / (pow2(gA(m,n)) * r(m,n))) + gsig_r(m,n)
    /*2*/  * (-(gL(m,n) * pow2(gB(m,n))) / (2. * pow2(gA(m,n))) + (pow2(gA(m,n)) + 2
    /*4*/  * pow2(gB(m,n))) / (Power(gA(m,n),4) * r(m,n))) + (gLr_r(m,n) - (2
    /*4*/  * gDconfr_r(m,n)) / pow2(gA(m,n))) / r(m,n)) + gAsig(m,n) * (gtrA(m,n)
    /*2*/  + gtrK(m,n)));
}
Real BimetricEvolve::eq_fAsig_t( Int m, Int n )
{
    return fAsig_r(m,n) * fBet(m,n) + 2 * fAsig(m,n) * fBetr(m,n) + (k_f * (fJA1(m,n)
    /*2*/  - fJA2(m,n))) / pow2(r(m,n)) - (exp(-4 * fconf(m,n)) * fDAlpr_r(m,n))
    /*0*/  / (pow2(fA(m,n)) * r(m,n)) + fAlp(m,n) * (exp(-4 * fconf(m,n))
    /*1*/  * ((fsig_rr(m,n) * pow2(fB(m,n))) / (2. * Power(fA(m,n),4))
    /*2*/  + (pow2(fsig(m,n)) * pow2(fB(m,n))) / Power(fA(m,n),4) + fDB(m,n) * ((2
    /*4*/  * fconf(m,n) * fDconf(m,n)) / (pow2(fA(m,n)) * pow2(r(m,n))) + fDAlp(m,n)
    /*3*/  / (pow2(fA(m,n)) * pow2(r(m,n)))) + fDA(m,n) * ((-4 * fDB(m,n))
    /*3*/  / (pow2(fA(m,n)) * pow2(r(m,n))) + (2 * fconf(m,n) * fDconf(m,n))
    /*3*/  / (pow2(fA(m,n)) * pow2(r(m,n))) + fDAlp(m,n) / (pow2(fA(m,n))
    /*4*/  * pow2(r(m,n)))) + (4 * fconf(m,n) * fDconf(m,n) * fDAlp(m,n))
    /*2*/  / (pow2(fA(m,n)) * pow2(r(m,n))) + (3 * pow2(fDA(m,n))) / (pow2(fA(m,n))
    /*3*/  * pow2(r(m,n))) + (4 * pow2(fconf(m,n)) * pow2(fDconf(m,n))) / (pow2(fA(m,n))
    /*3*/  * pow2(r(m,n))) + fsig(m,n) * ((-4 * fDB(m,n)) / (pow2(fA(m,n)) * r(m,n))
    /*3*/  - (fL(m,n) * pow2(fB(m,n))) / (pow2(fA(m,n)) * r(m,n))) + fsig_r(m,n)
    /*2*/  * (-(fL(m,n) * pow2(fB(m,n))) / (2. * pow2(fA(m,n))) + (pow2(fA(m,n)) + 2
    /*4*/  * pow2(fB(m,n))) / (Power(fA(m,n),4) * r(m,n))) + (fLr_r(m,n) - (2
    /*4*/  * fDconfr_r(m,n)) / pow2(fA(m,n))) / r(m,n)) + fAsig(m,n) * (ftrA(m,n)
    /*2*/  + ftrK(m,n)));
}
