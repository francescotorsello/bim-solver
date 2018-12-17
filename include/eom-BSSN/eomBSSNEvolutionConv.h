/** @file  eomBSSNEvolutionConv.h
 *  @author Francesco Torsello
 *  @brief The regularized cBSSN evolution equations.
 *  @version 2018-12-17T16:06:07
 *  @image html BSSNevolutionConv.png
 */

Real BimetricEvolve::eq_gconf_t( Int m, Int n )
{
    return -gdet_pff(m,n) / (12. * gdet(m,n)) + gconf_convr(m,n) + (gAlp(m,n) 
    /*1*/  * (-gtrA(m,n) - gtrK(m,n))) / 6.;
}
Real BimetricEvolve::eq_fconf_t( Int m, Int n )
{
    return -fdet_pff(m,n) / (12. * fdet(m,n)) + fconf_convr(m,n) + (fAlp(m,n) 
    /*1*/  * (-ftrA(m,n) - ftrK(m,n))) / 6.;
}
Real BimetricEvolve::eq_gDconf_t( Int m, Int n )
{
    return -gdet_pff_r(m,n) / (12. * gdet(m,n)) + gBet_gDconf_r(m,n) 
    /*0*/  + gDconf_convr(m,n) + (gdet_pff(m,n) * gdet_r(m,n)) / (12. 
    /*1*/  * pow2(gdet(m,n))) + (gDAlp(m,n) * (-gtrA(m,n) - gtrK(m,n))) / 6. 
    /*0*/  + (gAlp(m,n) * (-gtrA_r(m,n) - gtrK_r(m,n))) / 6.;
}
Real BimetricEvolve::eq_fDconf_t( Int m, Int n )
{
    return -fdet_pff_r(m,n) / (12. * fdet(m,n)) + fBet_fDconf_r(m,n) 
    /*0*/  + fDconf_convr(m,n) + (fdet_pff(m,n) * fdet_r(m,n)) / (12. 
    /*1*/  * pow2(fdet(m,n))) + (fDAlp(m,n) * (-ftrA(m,n) - ftrK(m,n))) / 6. 
    /*0*/  + (fAlp(m,n) * (-ftrA_r(m,n) - ftrK_r(m,n))) / 6.;
}
Real BimetricEvolve::eq_gtrK_t( Int m, Int n )
{
    return gtrK_convr(m,n) + k_g * gJK(m,n) + exp(-4 * gconf(m,n)) * ((3 
    /*2*/  * gDers(m,n)) / 2. - (3 * (gDB(m,n) + 2 * gDconf(m,n)) * gDAlp(m,n)) 
    /*1*/  / pow2(gA(m,n)) - (3 * gDAlp_r(m,n)) / pow2(gA(m,n))) + (2 * exp(-4 
    /*2*/  * gconf(m,n)) * gDAlpr_r(m,n) * r(m,n)) / pow2(gA(m,n)) - gtrA_pff(m,n) 
    /*0*/  + (gAlp(m,n) * (3 * pow2(gA1(m,n)) + 6 * pow2(gA2(m,n)) + gtrK(m,n) * (2 
    /*3*/  * gtrA(m,n) + gtrK(m,n)))) / 3.;
}
Real BimetricEvolve::eq_ftrK_t( Int m, Int n )
{
    return ftrK_convr(m,n) + k_f * fJK(m,n) + exp(-4 * fconf(m,n)) * ((3 
    /*2*/  * fDers(m,n)) / 2. - (3 * (fDB(m,n) + 2 * fDconf(m,n)) * fDAlp(m,n)) 
    /*1*/  / pow2(fA(m,n)) - (3 * fDAlp_r(m,n)) / pow2(fA(m,n))) + (2 * exp(-4 
    /*2*/  * fconf(m,n)) * fDAlpr_r(m,n) * r(m,n)) / pow2(fA(m,n)) - ftrA_pff(m,n) 
    /*0*/  + (fAlp(m,n) * (3 * pow2(fA1(m,n)) + 6 * pow2(fA2(m,n)) + ftrK(m,n) * (2 
    /*3*/  * ftrA(m,n) + ftrK(m,n)))) / 3.;
}
Real BimetricEvolve::eq_gA_t( Int m, Int n )
{
    return gA_convr(m,n) + (gdet_pff(m,n) * gA(m,n)) / (6. * gdet(m,n)) + gBet_r(m,n)
    /*0*/  * gA(m,n) + (gAlp(m,n) * gA(m,n) * (-3 * gA1(m,n) + gtrA(m,n) 
    /*2*/  + gtrAv(m,n))) / 3.;
}
Real BimetricEvolve::eq_gB_t( Int m, Int n )
{
    return gB_convr(m,n) + (gdet_pff(m,n) * gB(m,n)) / (6. * gdet(m,n)) + gBetr(m,n)
    /*0*/  * gB(m,n) + (gAlp(m,n) * gB(m,n) * (-3 * gA2(m,n) + gtrA(m,n) 
    /*2*/  + gtrAv(m,n))) / 3.;
}
Real BimetricEvolve::eq_gDA_t( Int m, Int n )
{
    return gdet_pff_r(m,n) / (6. * gdet(m,n)) + gBet_gDA_r(m,n) + gBet_rr(m,n) 
    /*0*/  + gDA_convr(m,n) - (gdet_pff(m,n) * gdet_r(m,n)) / (6. * pow2(gdet(m,n))) 
    /*0*/  + (gDAlp(m,n) * (-3 * gA1(m,n) + gtrA(m,n))) / 3. + (gAlp(m,n) * (-3 
    /*2*/  * gA1_r(m,n) + gtrA_r(m,n))) / 3.;
}
Real BimetricEvolve::eq_gDB_t( Int m, Int n )
{
    return gdet_pff_r(m,n) / (6. * gdet(m,n)) + gBet_gDB_r(m,n) + gDB_convr(m,n) 
    /*0*/  + gBetr_r(m,n) - (gdet_pff(m,n) * gdet_r(m,n)) / (6. * pow2(gdet(m,n))) 
    /*0*/  + (gDAlp(m,n) * (-3 * gA2(m,n) + gtrA(m,n))) / 3. + (gAlp(m,n) * (-3 
    /*2*/  * gA2_r(m,n) + gtrA_r(m,n))) / 3.;
}
Real BimetricEvolve::eq_fA_t( Int m, Int n )
{
    return fA_convr(m,n) + (fdet_pff(m,n) * fA(m,n)) / (6. * fdet(m,n)) + fBet_r(m,n)
    /*0*/  * fA(m,n) + (fAlp(m,n) * fA(m,n) * (-3 * fA1(m,n) + ftrA(m,n) 
    /*2*/  + ftrAv(m,n))) / 3.;
}
Real BimetricEvolve::eq_fB_t( Int m, Int n )
{
    return fB_convr(m,n) + (fdet_pff(m,n) * fB(m,n)) / (6. * fdet(m,n)) + fBetr(m,n)
    /*0*/  * fB(m,n) + (fAlp(m,n) * fB(m,n) * (-3 * fA2(m,n) + ftrA(m,n) 
    /*2*/  + ftrAv(m,n))) / 3.;
}
Real BimetricEvolve::eq_fDA_t( Int m, Int n )
{
    return fdet_pff_r(m,n) / (6. * fdet(m,n)) + fBet_fDA_r(m,n) + fBet_rr(m,n) 
    /*0*/  + fDA_convr(m,n) - (fdet_pff(m,n) * fdet_r(m,n)) / (6. * pow2(fdet(m,n))) 
    /*0*/  + (fDAlp(m,n) * (-3 * fA1(m,n) + ftrA(m,n))) / 3. + (fAlp(m,n) * (-3 
    /*2*/  * fA1_r(m,n) + ftrA_r(m,n))) / 3.;
}
Real BimetricEvolve::eq_fDB_t( Int m, Int n )
{
    return fdet_pff_r(m,n) / (6. * fdet(m,n)) + fBet_fDB_r(m,n) + fDB_convr(m,n) 
    /*0*/  + fBetr_r(m,n) - (fdet_pff(m,n) * fdet_r(m,n)) / (6. * pow2(fdet(m,n))) 
    /*0*/  + (fDAlp(m,n) * (-3 * fA2(m,n) + ftrA(m,n))) / 3. + (fAlp(m,n) * (-3 
    /*2*/  * fA2_r(m,n) + ftrA_r(m,n))) / 3.;
}
Real BimetricEvolve::eq_gA1_t( Int m, Int n )
{
    return gRicci(m,n) + gA1_convr(m,n) + k_g * gJA1(m,n) + exp(-4 * gconf(m,n)) 
    /*0*/  * (gDers(m,n) + gAlp(m,n) * ((2 * gDconf(m,n) * gDers(m,n)) / gDAlp(m,n) - (8
    /*3*/  * pow2(gDconf(m,n))) / (3. * pow2(gA(m,n))))) + exp(-4 * gconf(m,n)) * ((-2
    /*2*/  * gDAlpr_r(m,n)) / (3. * pow2(gA(m,n))) - (4 * gDconfr_r(m,n) * gAlp(m,n))
    /*1*/  / (3. * pow2(gA(m,n)))) * r(m,n) + gtrA_pff(m,n) / 3. + (gAlp(m,n) * (3 
    /*2*/  * gA1(m,n) - gtrA(m,n)) * (gtrA(m,n) + gtrK(m,n))) / 3.;
}
Real BimetricEvolve::eq_gA2_t( Int m, Int n )
{
    return -gRicci(m,n) / 2. + gA2_convr(m,n) + k_g * gJA2(m,n) + exp(-4 
    /*1*/  * gconf(m,n)) * (-gDers(m,n) / 2. + gAlp(m,n) * (-((gDconf(m,n) * gDers(m,n))
    /*3*/  / gDAlp(m,n)) + (4 * pow2(gDconf(m,n))) / (3. * pow2(gA(m,n))))) + exp(-4 
    /*1*/  * gconf(m,n)) * (gDAlpr_r(m,n) / (3. * pow2(gA(m,n))) + (2 * gDconfr_r(m,n) 
    /*2*/  * gAlp(m,n)) / (3. * pow2(gA(m,n)))) * r(m,n) + gtrA_pff(m,n) / 3. 
    /*0*/  + (gAlp(m,n) * (3 * gA2(m,n) - gtrA(m,n)) * (gtrA(m,n) + gtrK(m,n))) / 3.;
}
Real BimetricEvolve::eq_fA1_t( Int m, Int n )
{
    return fRicci(m,n) + fA1_convr(m,n) + k_f * fJA1(m,n) + exp(-4 * fconf(m,n)) 
    /*0*/  * (fDers(m,n) + fAlp(m,n) * ((2 * fDconf(m,n) * fDers(m,n)) / fDAlp(m,n) - (8
    /*3*/  * pow2(fDconf(m,n))) / (3. * pow2(fA(m,n))))) + exp(-4 * fconf(m,n)) * ((-2
    /*2*/  * fDAlpr_r(m,n)) / (3. * pow2(fA(m,n))) - (4 * fDconfr_r(m,n) * fAlp(m,n))
    /*1*/  / (3. * pow2(fA(m,n)))) * r(m,n) + ftrA_pff(m,n) / 3. + (fAlp(m,n) * (3 
    /*2*/  * fA1(m,n) - ftrA(m,n)) * (ftrA(m,n) + ftrK(m,n))) / 3.;
}
Real BimetricEvolve::eq_fA2_t( Int m, Int n )
{
    return -fRicci(m,n) / 2. + fA2_convr(m,n) + k_f * fJA2(m,n) + exp(-4 
    /*1*/  * fconf(m,n)) * (-fDers(m,n) / 2. + fAlp(m,n) * (-((fDconf(m,n) * fDers(m,n))
    /*3*/  / fDAlp(m,n)) + (4 * pow2(fDconf(m,n))) / (3. * pow2(fA(m,n))))) + exp(-4 
    /*1*/  * fconf(m,n)) * (fDAlpr_r(m,n) / (3. * pow2(fA(m,n))) + (2 * fDconfr_r(m,n) 
    /*2*/  * fAlp(m,n)) / (3. * pow2(fA(m,n)))) * r(m,n) + ftrA_pff(m,n) / 3. 
    /*0*/  + (fAlp(m,n) * (3 * fA2(m,n) - ftrA(m,n)) * (ftrA(m,n) + ftrK(m,n))) / 3.;
}
Real BimetricEvolve::eq_gL_t( Int m, Int n )
{
    return gL_convr(m,n) - gBet_gL_r(m,n) + k_g * gJL(m,n) + (-(gdet_pff(m,n) 
    /*2*/  * gL(m,n)) / 3. - gdet_pff_r(m,n) / (6. * pow2(gA(m,n)))) / gdet(m,n) 
    /*0*/  + gBet_rr(m,n) / pow2(gA(m,n)) + (gdet_pff(m,n) * gdet_r(m,n)) / (6. 
    /*1*/  * pow2(gdet(m,n)) * pow2(gA(m,n))) + (2 * gBetr_r(m,n)) / pow2(gB(m,n)) + (4
    /*1*/  * gAsig(m,n) * (gDA(m,n) + gDB(m,n) + 6 * gDconf(m,n)) * gAlp(m,n) 
    /*1*/  * pow2(r(m,n))) / (3. * pow2(gA(m,n))) + (4 * gAsig(m,n) * gAlp(m,n) 
    /*1*/  * gsig(m,n) * pow3(r(m,n))) / (3. * pow2(gA(m,n))) + (gDAlp(m,n) * (-2 
    /*2*/  * gA1(m,n) + (2 * gtrA(m,n)) / 3.)) / pow2(gA(m,n)) - (4 * gAlp(m,n) 
    /*1*/  * (gtrA_r(m,n) + gtrK_r(m,n))) / (3. * pow2(gA(m,n)));
}
Real BimetricEvolve::eq_fL_t( Int m, Int n )
{
    return fL_convr(m,n) - fBet_fL_r(m,n) + k_f * fJL(m,n) + (-(fdet_pff(m,n) 
    /*2*/  * fL(m,n)) / 3. - fdet_pff_r(m,n) / (6. * pow2(fA(m,n)))) / fdet(m,n) 
    /*0*/  + fBet_rr(m,n) / pow2(fA(m,n)) + (fdet_pff(m,n) * fdet_r(m,n)) / (6. 
    /*1*/  * pow2(fdet(m,n)) * pow2(fA(m,n))) + (2 * fBetr_r(m,n)) / pow2(fB(m,n)) + (4
    /*1*/  * fAsig(m,n) * (fDA(m,n) + fDB(m,n) + 6 * fDconf(m,n)) * fAlp(m,n) 
    /*1*/  * pow2(r(m,n))) / (3. * pow2(fA(m,n))) + (4 * fAsig(m,n) * fAlp(m,n) 
    /*1*/  * fsig(m,n) * pow3(r(m,n))) / (3. * pow2(fA(m,n))) + (fDAlp(m,n) * (-2 
    /*2*/  * fA1(m,n) + (2 * ftrA(m,n)) / 3.)) / pow2(fA(m,n)) - (4 * fAlp(m,n) 
    /*1*/  * (ftrA_r(m,n) + ftrK_r(m,n))) / (3. * pow2(fA(m,n)));
}
Real BimetricEvolve::eq_gsig_t( Int m, Int n )
{
    return gsig_convr(m,n) + 2 * gBetr(m,n) * gsig(m,n) + (2 * gAsig(m,n) * gAlp(m,n)
    /*1*/  * pow2(gA(m,n))) / pow2(gB(m,n)) - (2 * gBetr_r(m,n) * pow2(gA(m,n))) 
    /*0*/  / (pow2(gB(m,n)) * r(m,n));
}
Real BimetricEvolve::eq_fsig_t( Int m, Int n )
{
    return fsig_convr(m,n) + 2 * fBetr(m,n) * fsig(m,n) + (2 * fAsig(m,n) * fAlp(m,n)
    /*1*/  * pow2(fA(m,n))) / pow2(fB(m,n)) - (2 * fBetr_r(m,n) * pow2(fA(m,n))) 
    /*0*/  / (pow2(fB(m,n)) * r(m,n));
}
Real BimetricEvolve::eq_gAsig_t( Int m, Int n )
{
    return gAsig_convr(m,n) + 2 * gAsig(m,n) * gBetr(m,n) + exp(-4 * gconf(m,n)) 
    /*0*/  * gAlp(m,n) * (((gsig_rr(m,n) / 2. + pow2(gsig(m,n))) * pow2(gB(m,n))) 
    /*1*/  / Power(gA(m,n),4) - (gsig_gL_r(m,n) * pow2(gB(m,n))) / (2. * pow2(gA(m,n))))
    /*0*/  + (k_g * (gJA1(m,n) - gJA2(m,n)) + exp(-4 * gconf(m,n)) * ((3 * gDers(m,n))
    /*2*/  / 2. + gAlp(m,n) * ((-3 * (5 * gDB(m,n) + 11 * gDconf(m,n)) * gDers(m,n)) 
    /*3*/  / gDAlp(m,n) + (40 * gDB(m,n) * gDconf(m,n) + 7 * pow2(gDB(m,n)) + 44 
    /*4*/  * pow2(gDconf(m,n))) / pow2(gA(m,n)) + (27 * pow2(gDers(m,n)) 
    /*4*/  * pow2(gA(m,n))) / (4. * pow2(gDAlp(m,n)))))) / pow2(r(m,n)) + (exp(-4 
    /*2*/  * gconf(m,n)) * (-(gDAlpr_r(m,n) / pow2(gA(m,n))) + gAlp(m,n) * (gLr_r(m,n) 
    /*3*/  + (2 * gsig_r(m,n) * pow2(gB(m,n))) / Power(gA(m,n),4) + (-2 * gDconfr_r(m,n)
    /*4*/  - 4 * gDB(m,n) * gsig(m,n) + gsig_r(m,n) - gL(m,n) * gsig(m,n) 
    /*4*/  * pow2(gB(m,n))) / pow2(gA(m,n))))) / r(m,n) + gAsig(m,n) * gAlp(m,n) 
    /*0*/  * (gtrA(m,n) + gtrK(m,n));
}
Real BimetricEvolve::eq_fAsig_t( Int m, Int n )
{
    return fAsig_convr(m,n) + 2 * fAsig(m,n) * fBetr(m,n) + exp(-4 * fconf(m,n)) 
    /*0*/  * fAlp(m,n) * (((fsig_rr(m,n) / 2. + pow2(fsig(m,n))) * pow2(fB(m,n))) 
    /*1*/  / Power(fA(m,n),4) - (fsig_fL_r(m,n) * pow2(fB(m,n))) / (2. * pow2(fA(m,n))))
    /*0*/  + (k_f * (fJA1(m,n) - fJA2(m,n)) + exp(-4 * fconf(m,n)) * ((3 * fDers(m,n))
    /*2*/  / 2. + fAlp(m,n) * ((-3 * (5 * fDB(m,n) + 11 * fDconf(m,n)) * fDers(m,n)) 
    /*3*/  / fDAlp(m,n) + (40 * fDB(m,n) * fDconf(m,n) + 7 * pow2(fDB(m,n)) + 44 
    /*4*/  * pow2(fDconf(m,n))) / pow2(fA(m,n)) + (27 * pow2(fDers(m,n)) 
    /*4*/  * pow2(fA(m,n))) / (4. * pow2(fDAlp(m,n)))))) / pow2(r(m,n)) + (exp(-4 
    /*2*/  * fconf(m,n)) * (-(fDAlpr_r(m,n) / pow2(fA(m,n))) + fAlp(m,n) * (fLr_r(m,n) 
    /*3*/  + (2 * fsig_r(m,n) * pow2(fB(m,n))) / Power(fA(m,n),4) + (-2 * fDconfr_r(m,n)
    /*4*/  - 4 * fDB(m,n) * fsig(m,n) + fsig_r(m,n) - fL(m,n) * fsig(m,n) 
    /*4*/  * pow2(fB(m,n))) / pow2(fA(m,n))))) / r(m,n) + fAsig(m,n) * fAlp(m,n) 
    /*0*/  * (ftrA(m,n) + ftrK(m,n));
}
