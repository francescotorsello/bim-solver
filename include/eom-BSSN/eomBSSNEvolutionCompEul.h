/** @file  eomBSSNEvolutionCompEul.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN Eulerian evolution equations.
 *  @version 2019-05-14T09:05:24
 *  @image html BSSNevolutionCompEul.png
 */

Real gconf_t( Int m, Int n )
{
    return gconf_r(m,n) * gBet(m,n) * pow2(r_minus(m,n)) - (gAlp(m,n) * gtrK(m,n)) / 6.;
}
Real fconf_t( Int m, Int n )
{
    return fconf_r(m,n) * fBet(m,n) * pow2(r_minus(m,n)) - (fAlp(m,n) * ftrK(m,n)) / 6.;
}
Real gtrK_t( Int m, Int n )
{
    return k_g * gJK(m,n) - (Power(r_minus(m,n),4) * exp(-4 * gconf(m,n)) * gAlp_rr(m,n)) 
    /*0*/  / pow2(gA(m,n)) + gAlp(m,n) * (-((gA1(m,n) * gAsig(m,n) * pow2(r_plus(m,n))) 
    /*2*/  / pow2(-2 + r_plus(m,n))) + 3 * pow2(gA1(m,n)) + (Power(r_plus(m,n),4) 
    /*2*/  * pow2(gAsig(m,n))) / (8. * Power(-2 + r_plus(m,n),4)) + pow2(gtrK(m,n)) / 3.) 
    /*0*/  + gAlp_r(m,n) * ((-2 * Power(r_minus(m,n),4) * exp(-4 * gconf(m,n)) * gconf_r(m,n))
    /*1*/  / pow2(gA(m,n)) - (2 * Power(r_minus(m,n),4) * exp(-4 * gconf(m,n)) * gB_r(m,n))
    /*1*/  / (gB(m,n) * pow2(gA(m,n))) - (2 * (4 + (-2 + r_minus(m,n)) * r_plus(m,n)) * exp(-4 
    /*3*/  * gconf(m,n)) * pow2(r_minus(m,n))) / (r_plus(m,n) * pow2(gA(m,n))) + (Power(r_minus(m,n),4) 
    /*2*/  * exp(-4 * gconf(m,n)) * gA_r(m,n)) / pow3(gA(m,n))) + gBet(m,n) 
    /*0*/  * pow2(r_minus(m,n)) * gtrK_r(m,n);
}
Real ftrK_t( Int m, Int n )
{
    return k_f * fJK(m,n) - (Power(r_minus(m,n),4) * exp(-4 * fconf(m,n)) * fAlp_rr(m,n)) 
    /*0*/  / pow2(fA(m,n)) + fAlp(m,n) * (-((fA1(m,n) * fAsig(m,n) * pow2(r_plus(m,n))) 
    /*2*/  / pow2(-2 + r_plus(m,n))) + 3 * pow2(fA1(m,n)) + (Power(r_plus(m,n),4) 
    /*2*/  * pow2(fAsig(m,n))) / (8. * Power(-2 + r_plus(m,n),4)) + pow2(ftrK(m,n)) / 3.) 
    /*0*/  + fAlp_r(m,n) * ((-2 * Power(r_minus(m,n),4) * exp(-4 * fconf(m,n)) * fconf_r(m,n))
    /*1*/  / pow2(fA(m,n)) - (2 * Power(r_minus(m,n),4) * exp(-4 * fconf(m,n)) * fB_r(m,n))
    /*1*/  / (fB(m,n) * pow2(fA(m,n))) - (2 * (4 + (-2 + r_minus(m,n)) * r_plus(m,n)) * exp(-4 
    /*3*/  * fconf(m,n)) * pow2(r_minus(m,n))) / (r_plus(m,n) * pow2(fA(m,n))) + (Power(r_minus(m,n),4) 
    /*2*/  * exp(-4 * fconf(m,n)) * fA_r(m,n)) / pow3(fA(m,n))) + fBet(m,n) 
    /*0*/  * pow2(r_minus(m,n)) * ftrK_r(m,n);
}
Real gA_t( Int m, Int n )
{
    return -(gA1(m,n) * gAlp(m,n) * gA(m,n)) + gBet_r(m,n) * gA(m,n) * pow2(r_minus(m,n)) 
    /*0*/  + gBet(m,n) * gA_r(m,n) * pow2(r_minus(m,n));
}
Real gB_t( Int m, Int n )
{
    return gBetr(m,n) * gB(m,n) + gBet(m,n) * gB_r(m,n) * pow2(r_minus(m,n)) + (gAlp(m,n) 
    /*1*/  * gB(m,n) * (-12 * gA1(m,n) + (3 * gAsig(m,n) * pow2(r_plus(m,n))) / pow2(-2 
    /*3*/  + r_plus(m,n)))) / 12.;
}
Real fA_t( Int m, Int n )
{
    return -(fA1(m,n) * fAlp(m,n) * fA(m,n)) + fBet_r(m,n) * fA(m,n) * pow2(r_minus(m,n)) 
    /*0*/  + fBet(m,n) * fA_r(m,n) * pow2(r_minus(m,n));
}
Real fB_t( Int m, Int n )
{
    return fBetr(m,n) * fB(m,n) + fBet(m,n) * fB_r(m,n) * pow2(r_minus(m,n)) + (fAlp(m,n) 
    /*1*/  * fB(m,n) * (-12 * fA1(m,n) + (3 * fAsig(m,n) * pow2(r_plus(m,n))) / pow2(-2 
    /*3*/  + r_plus(m,n)))) / 12.;
}
Real gA1_t( Int m, Int n )
{
    return gRicci(m,n) + k_g * gJA1(m,n) + gA1_r(m,n) * gBet(m,n) * pow2(r_minus(m,n)) 
    /*0*/  + exp(-4 * gconf(m,n)) * ((-2 * Power(r_minus(m,n),4) * gAlp_rr(m,n)) / (3. 
    /*2*/  * pow2(gA(m,n))) + gAlp(m,n) * ((-4 * Power(r_minus(m,n),4) * gconf_rr(m,n)) / (3.
    /*3*/  * pow2(gA(m,n))) + (8 * Power(r_minus(m,n),4) * pow2(gconf_r(m,n))) / (3. 
    /*3*/  * pow2(gA(m,n))) + gconf_r(m,n) * ((4 * Power(r_minus(m,n),4) * gB_r(m,n)) / (3. 
    /*4*/  * gB(m,n) * pow2(gA(m,n))) - (8 * (-2 + r_plus(m,n) + r_minus(m,n) * r_plus(m,n)) 
    /*4*/  * pow2(r_minus(m,n))) / (3. * r_plus(m,n) * pow2(gA(m,n))) + (4 * Power(r_minus(m,n),4) 
    /*4*/  * gA_r(m,n)) / (3. * pow3(gA(m,n))))) + gAlp_r(m,n) * ((8 * Power(r_minus(m,n),4) 
    /*3*/  * gconf_r(m,n)) / (3. * pow2(gA(m,n))) + (2 * Power(r_minus(m,n),4) * gB_r(m,n)) 
    /*2*/  / (3. * gB(m,n) * pow2(gA(m,n))) - (4 * (-2 + r_plus(m,n) + r_minus(m,n) * r_plus(m,n)) 
    /*3*/  * pow2(r_minus(m,n))) / (3. * r_plus(m,n) * pow2(gA(m,n))) + (2 * Power(r_minus(m,n),4) 
    /*3*/  * gA_r(m,n)) / (3. * pow3(gA(m,n))))) + gA1(m,n) * gAlp(m,n) * gtrK(m,n);
}
Real fA1_t( Int m, Int n )
{
    return fRicci(m,n) + k_f * fJA1(m,n) + fA1_r(m,n) * fBet(m,n) * pow2(r_minus(m,n)) 
    /*0*/  + exp(-4 * fconf(m,n)) * ((-2 * Power(r_minus(m,n),4) * fAlp_rr(m,n)) / (3. 
    /*2*/  * pow2(fA(m,n))) + fAlp(m,n) * ((-4 * Power(r_minus(m,n),4) * fconf_rr(m,n)) / (3.
    /*3*/  * pow2(fA(m,n))) + (8 * Power(r_minus(m,n),4) * pow2(fconf_r(m,n))) / (3. 
    /*3*/  * pow2(fA(m,n))) + fconf_r(m,n) * ((4 * Power(r_minus(m,n),4) * fB_r(m,n)) / (3. 
    /*4*/  * fB(m,n) * pow2(fA(m,n))) - (8 * (-2 + r_plus(m,n) + r_minus(m,n) * r_plus(m,n)) 
    /*4*/  * pow2(r_minus(m,n))) / (3. * r_plus(m,n) * pow2(fA(m,n))) + (4 * Power(r_minus(m,n),4) 
    /*4*/  * fA_r(m,n)) / (3. * pow3(fA(m,n))))) + fAlp_r(m,n) * ((8 * Power(r_minus(m,n),4) 
    /*3*/  * fconf_r(m,n)) / (3. * pow2(fA(m,n))) + (2 * Power(r_minus(m,n),4) * fB_r(m,n)) 
    /*2*/  / (3. * fB(m,n) * pow2(fA(m,n))) - (4 * (-2 + r_plus(m,n) + r_minus(m,n) * r_plus(m,n)) 
    /*3*/  * pow2(r_minus(m,n))) / (3. * r_plus(m,n) * pow2(fA(m,n))) + (2 * Power(r_minus(m,n),4) 
    /*3*/  * fA_r(m,n)) / (3. * pow3(fA(m,n))))) + fA1(m,n) * fAlp(m,n) * ftrK(m,n);
}
Real gL_t( Int m, Int n )
{
    return k_g * gJL(m,n) + gBet(m,n) * gL_r(m,n) * pow2(r_minus(m,n)) - (2 * gA1(m,n) 
    /*1*/  * gAlp_r(m,n) * pow2(r_minus(m,n))) / pow2(gA(m,n)) + (pow2(r_minus(m,n)) * ((6 
    /*3*/  * gBet_rr(m,n) * pow2(r_minus(m,n))) / pow2(gA(m,n)) + (12 * gBetr_r(m,n)) 
    /*2*/  / pow2(gB(m,n)))) / 6. + gBet_r(m,n) * (-(gL(m,n) * pow2(r_minus(m,n))) + (2 
    /*2*/  * pow3(r_minus(m,n))) / pow2(gA(m,n))) + gAlp(m,n) * ((2 * gAsig(m,n) 
    /*2*/  * gconf_r(m,n) * pow2(r_minus(m,n)) * pow2(r_plus(m,n))) / (pow2(-2 + r_plus(m,n)) 
    /*2*/  * pow2(gA(m,n))) - (gAsig(m,n) * gsig(m,n) * pow3(r_plus(m,n))) / (6. 
    /*2*/  * pow2(gA(m,n)) * pow3(-2 + r_plus(m,n))) + (pow2(r_minus(m,n)) * (gAsig(m,n) 
    /*3*/  * (gA_r(m,n) * gB(m,n) + gA(m,n) * gB_r(m,n)) * pow2(r_plus(m,n)) - 4 * gA(m,n) 
    /*3*/  * gB(m,n) * pow2(-2 + r_plus(m,n)) * gtrK_r(m,n))) / (3. * gB(m,n) * pow2(-2 
    /*3*/  + r_plus(m,n)) * pow3(gA(m,n))));
}
Real fL_t( Int m, Int n )
{
    return k_f * fJL(m,n) + fBet(m,n) * fL_r(m,n) * pow2(r_minus(m,n)) - (2 * fA1(m,n) 
    /*1*/  * fAlp_r(m,n) * pow2(r_minus(m,n))) / pow2(fA(m,n)) + (pow2(r_minus(m,n)) * ((6 
    /*3*/  * fBet_rr(m,n) * pow2(r_minus(m,n))) / pow2(fA(m,n)) + (12 * fBetr_r(m,n)) 
    /*2*/  / pow2(fB(m,n)))) / 6. + fBet_r(m,n) * (-(fL(m,n) * pow2(r_minus(m,n))) + (2 
    /*2*/  * pow3(r_minus(m,n))) / pow2(fA(m,n))) + fAlp(m,n) * ((2 * fAsig(m,n) 
    /*2*/  * fconf_r(m,n) * pow2(r_minus(m,n)) * pow2(r_plus(m,n))) / (pow2(-2 + r_plus(m,n)) 
    /*2*/  * pow2(fA(m,n))) - (fAsig(m,n) * fsig(m,n) * pow3(r_plus(m,n))) / (6. 
    /*2*/  * pow2(fA(m,n)) * pow3(-2 + r_plus(m,n))) + (pow2(r_minus(m,n)) * (fAsig(m,n) 
    /*3*/  * (fA_r(m,n) * fB(m,n) + fA(m,n) * fB_r(m,n)) * pow2(r_plus(m,n)) - 4 * fA(m,n) 
    /*3*/  * fB(m,n) * pow2(-2 + r_plus(m,n)) * ftrK_r(m,n))) / (3. * fB(m,n) * pow2(-2 
    /*3*/  + r_plus(m,n)) * pow3(fA(m,n))));
}
Real gsig_t( Int m, Int n )
{
    return 2 * gBetr(m,n) * gsig(m,n) + gBet(m,n) * gsig_r(m,n) * pow2(r_minus(m,n)) + (2 
    /*1*/  * gAsig(m,n) * gAlp(m,n) * pow2(gA(m,n))) / pow2(gB(m,n)) + (4 * gBetr_r(m,n)
    /*1*/  * pow2(r_minus(m,n)) * pow2(gA(m,n))) / pow2(gB(m,n)) - (8 * gBetr_r(m,n) 
    /*1*/  * pow2(r_minus(m,n)) * pow2(gA(m,n))) / (r_plus(m,n) * pow2(gB(m,n)));
}
Real fsig_t( Int m, Int n )
{
    return 2 * fBetr(m,n) * fsig(m,n) + fBet(m,n) * fsig_r(m,n) * pow2(r_minus(m,n)) + (2 
    /*1*/  * fAsig(m,n) * fAlp(m,n) * pow2(fA(m,n))) / pow2(fB(m,n)) + (4 * fBetr_r(m,n)
    /*1*/  * pow2(r_minus(m,n)) * pow2(fA(m,n))) / pow2(fB(m,n)) - (8 * fBetr_r(m,n) 
    /*1*/  * pow2(r_minus(m,n)) * pow2(fA(m,n))) / (r_plus(m,n) * pow2(fB(m,n)));
}
Real gAsig_t( Int m, Int n )
{
    return -4 * gAsig(m,n) * gBet(m,n) + k_g * (-(gJA1(m,n) * (-4 + 16 / r_plus(m,n) - 16 
    /*3*/  / pow2(r_plus(m,n)))) / 2. + gJA1(m,n) * (4 - 16 / r_plus(m,n) + 16 / pow2(r_plus(m,n)))) 
    /*0*/  + exp(-4 * gconf(m,n)) * gAlp(m,n) * ((pow2(gsig(m,n)) * pow2(gB(m,n))) 
    /*1*/  / Power(gA(m,n),4) + (2 * gL(m,n) * gsig(m,n) * pow2(gB(m,n))) 
    /*1*/  / pow2(gA(m,n))) + (8 * gAsig(m,n) * gBet(m,n) - (4 * exp(-4 * gconf(m,n)) 
    /*2*/  * gAlp(m,n) * gL(m,n) * gsig(m,n) * pow2(gB(m,n))) / pow2(gA(m,n))) / r_plus(m,n) 
    /*0*/  + pow2(r_minus(m,n)) * (gAsig_r(m,n) * gBet(m,n) + (exp(-4 * gconf(m,n)) * ((-4 
    /*4*/  * gderAlpr_r(m,n)) / pow2(gA(m,n)) + gAlp(m,n) * (4 * gLr_r(m,n) - (8 
    /*5*/  * gderconfr_r(m,n)) / pow2(gA(m,n)) - (16 * gsig(m,n) * gB_r(m,n)) / (gB(m,n)
    /*5*/  * pow2(gA(m,n))) + gsig_r(m,n) * (4 / pow2(gA(m,n)) + (8 * pow2(gB(m,n))) 
    /*5*/  / Power(gA(m,n),4))))) / r_plus(m,n) + exp(-4 * gconf(m,n)) * ((2 
    /*3*/  * gderAlpr_r(m,n)) / pow2(gA(m,n)) + gAlp(m,n) * (-2 * gLr_r(m,n) + (4 
    /*4*/  * gderconfr_r(m,n)) / pow2(gA(m,n)) + (8 * gsig(m,n) * gB_r(m,n)) / (gB(m,n)
    /*4*/  * pow2(gA(m,n))) + gsig_r(m,n) * (-2 / pow2(gA(m,n)) - ((8 + gL(m,n) 
    /*6*/  * pow2(gA(m,n))) * pow2(gB(m,n))) / (2. * Power(gA(m,n),4)))))) + (exp(-4 
    /*2*/  * gconf(m,n)) * gAlp(m,n) * gsig_r(m,n) * pow2(gB(m,n)) * pow3(r_minus(m,n))) 
    /*0*/  / Power(gA(m,n),4) + Power(r_minus(m,n),4) * ((exp(-4 * gconf(m,n)) * (gAlp_r(m,n)
    /*3*/  * ((64 * gconf_r(m,n)) / pow2(gA(m,n)) + (16 * gB_r(m,n)) / (gB(m,n) 
    /*5*/  * pow2(gA(m,n))) + (16 * gA_r(m,n)) / pow3(gA(m,n))) + gAlp(m,n) * ((64 
    /*5*/  * pow2(gconf_r(m,n))) / pow2(gA(m,n)) + (48 * pow2(gA_r(m,n))) 
    /*4*/  / Power(gA(m,n),4) + gconf_r(m,n) * ((32 * gB_r(m,n)) / (gB(m,n) 
    /*6*/  * pow2(gA(m,n))) + (32 * gA_r(m,n)) / pow3(gA(m,n))) - (64 * gA_r(m,n) 
    /*5*/  * gB_r(m,n)) / (gB(m,n) * pow3(gA(m,n)))))) / pow2(r_plus(m,n)) + exp(-4 
    /*2*/  * gconf(m,n)) * (gAlp_r(m,n) * ((16 * gconf_r(m,n)) / pow2(gA(m,n)) + (4 
    /*4*/  * gB_r(m,n)) / (gB(m,n) * pow2(gA(m,n))) + (4 * gA_r(m,n)) / pow3(gA(m,n))) 
    /*2*/  + gAlp(m,n) * ((16 * pow2(gconf_r(m,n))) / pow2(gA(m,n)) + (12 
    /*4*/  * pow2(gA_r(m,n))) / Power(gA(m,n),4) + (gsig_rr(m,n) * pow2(gB(m,n))) / (2.
    /*4*/  * Power(gA(m,n),4)) + gconf_r(m,n) * ((8 * gB_r(m,n)) / (gB(m,n) 
    /*5*/  * pow2(gA(m,n))) + (8 * gA_r(m,n)) / pow3(gA(m,n))) - (16 * gA_r(m,n) 
    /*4*/  * gB_r(m,n)) / (gB(m,n) * pow3(gA(m,n))))) + (exp(-4 * gconf(m,n)) 
    /*2*/  * (gAlp_r(m,n) * ((-64 * gconf_r(m,n)) / pow2(gA(m,n)) - (16 * gB_r(m,n)) 
    /*4*/  / (gB(m,n) * pow2(gA(m,n))) - (16 * gA_r(m,n)) / pow3(gA(m,n))) + gAlp(m,n) 
    /*3*/  * ((-64 * pow2(gconf_r(m,n))) / pow2(gA(m,n)) - (48 * pow2(gA_r(m,n))) 
    /*4*/  / Power(gA(m,n),4) + gconf_r(m,n) * ((-32 * gB_r(m,n)) / (gB(m,n) 
    /*6*/  * pow2(gA(m,n))) - (32 * gA_r(m,n)) / pow3(gA(m,n))) + (64 * gA_r(m,n) 
    /*5*/  * gB_r(m,n)) / (gB(m,n) * pow3(gA(m,n)))))) / r_plus(m,n)) + gAsig(m,n) * gAlp(m,n)
    /*0*/  * gtrK(m,n);
}
Real fAsig_t( Int m, Int n )
{
    return -4 * fAsig(m,n) * fBet(m,n) + k_f * (-(fJA1(m,n) * (-4 + 16 / r_plus(m,n) - 16 
    /*3*/  / pow2(r_plus(m,n)))) / 2. + fJA1(m,n) * (4 - 16 / r_plus(m,n) + 16 / pow2(r_plus(m,n)))) 
    /*0*/  + exp(-4 * fconf(m,n)) * fAlp(m,n) * ((pow2(fsig(m,n)) * pow2(fB(m,n))) 
    /*1*/  / Power(fA(m,n),4) + (2 * fL(m,n) * fsig(m,n) * pow2(fB(m,n))) 
    /*1*/  / pow2(fA(m,n))) + (8 * fAsig(m,n) * fBet(m,n) - (4 * exp(-4 * fconf(m,n)) 
    /*2*/  * fAlp(m,n) * fL(m,n) * fsig(m,n) * pow2(fB(m,n))) / pow2(fA(m,n))) / r_plus(m,n) 
    /*0*/  + pow2(r_minus(m,n)) * (fAsig_r(m,n) * fBet(m,n) + (exp(-4 * fconf(m,n)) * ((-4 
    /*4*/  * fderAlpr_r(m,n)) / pow2(fA(m,n)) + fAlp(m,n) * (4 * fLr_r(m,n) - (8 
    /*5*/  * fderconfr_r(m,n)) / pow2(fA(m,n)) - (16 * fsig(m,n) * fB_r(m,n)) / (fB(m,n)
    /*5*/  * pow2(fA(m,n))) + fsig_r(m,n) * (4 / pow2(fA(m,n)) + (8 * pow2(fB(m,n))) 
    /*5*/  / Power(fA(m,n),4))))) / r_plus(m,n) + exp(-4 * fconf(m,n)) * ((2 
    /*3*/  * fderAlpr_r(m,n)) / pow2(fA(m,n)) + fAlp(m,n) * (-2 * fLr_r(m,n) + (4 
    /*4*/  * fderconfr_r(m,n)) / pow2(fA(m,n)) + (8 * fsig(m,n) * fB_r(m,n)) / (fB(m,n)
    /*4*/  * pow2(fA(m,n))) + fsig_r(m,n) * (-2 / pow2(fA(m,n)) - ((8 + fL(m,n) 
    /*6*/  * pow2(fA(m,n))) * pow2(fB(m,n))) / (2. * Power(fA(m,n),4)))))) + (exp(-4 
    /*2*/  * fconf(m,n)) * fAlp(m,n) * fsig_r(m,n) * pow2(fB(m,n)) * pow3(r_minus(m,n))) 
    /*0*/  / Power(fA(m,n),4) + Power(r_minus(m,n),4) * ((exp(-4 * fconf(m,n)) * (fAlp_r(m,n)
    /*3*/  * ((64 * fconf_r(m,n)) / pow2(fA(m,n)) + (16 * fB_r(m,n)) / (fB(m,n) 
    /*5*/  * pow2(fA(m,n))) + (16 * fA_r(m,n)) / pow3(fA(m,n))) + fAlp(m,n) * ((64 
    /*5*/  * pow2(fconf_r(m,n))) / pow2(fA(m,n)) + (48 * pow2(fA_r(m,n))) 
    /*4*/  / Power(fA(m,n),4) + fconf_r(m,n) * ((32 * fB_r(m,n)) / (fB(m,n) 
    /*6*/  * pow2(fA(m,n))) + (32 * fA_r(m,n)) / pow3(fA(m,n))) - (64 * fA_r(m,n) 
    /*5*/  * fB_r(m,n)) / (fB(m,n) * pow3(fA(m,n)))))) / pow2(r_plus(m,n)) + exp(-4 
    /*2*/  * fconf(m,n)) * (fAlp_r(m,n) * ((16 * fconf_r(m,n)) / pow2(fA(m,n)) + (4 
    /*4*/  * fB_r(m,n)) / (fB(m,n) * pow2(fA(m,n))) + (4 * fA_r(m,n)) / pow3(fA(m,n))) 
    /*2*/  + fAlp(m,n) * ((16 * pow2(fconf_r(m,n))) / pow2(fA(m,n)) + (12 
    /*4*/  * pow2(fA_r(m,n))) / Power(fA(m,n),4) + (fsig_rr(m,n) * pow2(fB(m,n))) / (2.
    /*4*/  * Power(fA(m,n),4)) + fconf_r(m,n) * ((8 * fB_r(m,n)) / (fB(m,n) 
    /*5*/  * pow2(fA(m,n))) + (8 * fA_r(m,n)) / pow3(fA(m,n))) - (16 * fA_r(m,n) 
    /*4*/  * fB_r(m,n)) / (fB(m,n) * pow3(fA(m,n))))) + (exp(-4 * fconf(m,n)) 
    /*2*/  * (fAlp_r(m,n) * ((-64 * fconf_r(m,n)) / pow2(fA(m,n)) - (16 * fB_r(m,n)) 
    /*4*/  / (fB(m,n) * pow2(fA(m,n))) - (16 * fA_r(m,n)) / pow3(fA(m,n))) + fAlp(m,n) 
    /*3*/  * ((-64 * pow2(fconf_r(m,n))) / pow2(fA(m,n)) - (48 * pow2(fA_r(m,n))) 
    /*4*/  / Power(fA(m,n),4) + fconf_r(m,n) * ((-32 * fB_r(m,n)) / (fB(m,n) 
    /*6*/  * pow2(fA(m,n))) - (32 * fA_r(m,n)) / pow3(fA(m,n))) + (64 * fA_r(m,n) 
    /*5*/  * fB_r(m,n)) / (fB(m,n) * pow3(fA(m,n)))))) / r_plus(m,n)) + fAsig(m,n) * fAlp(m,n)
    /*0*/  * ftrK(m,n);
}
