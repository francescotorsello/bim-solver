/** @file  eomBSSNEvolutionCompEul.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN Eulerian evolution equations.
 *  @version 2019-05-14T16:20:45
 *  @image html BSSNevolutionCompEul.png
 */

Real gconf_t( Int m, Int n )
{
    return 2 * gconf_r(m,n) * gBet(m,n) * pow2(r_minus(m,n)) - (gAlp(m,n) * gtrK(m,n)) / 6.;
}
Real fconf_t( Int m, Int n )
{
    return 2 * fconf_r(m,n) * fBet(m,n) * pow2(r_minus(m,n)) - (fAlp(m,n) * ftrK(m,n)) / 6.;
}
Real gtrK_t( Int m, Int n )
{
    return k_g * gJK(m,n) - (4 * Power(r_minus(m,n),4) * exp(-4 * gconf(m,n)) 
    /*1*/  * gAlp_rr(m,n)) / pow2(gA(m,n)) + gAlp_r(m,n) * ((-8 * Power(r_minus(m,n),4) 
    /*2*/  * exp(-4 * gconf(m,n)) * gconf_r(m,n)) / pow2(gA(m,n)) - (8 * Power(r_minus(m,n),4)
    /*2*/  * exp(-4 * gconf(m,n)) * gB_r(m,n)) / (gB(m,n) * pow2(gA(m,n))) + (4 
    /*2*/  * Power(r_minus(m,n),4) * exp(-4 * gconf(m,n)) * gA_r(m,n)) / pow3(gA(m,n)) - (8 
    /*2*/  * Power(r_minus(m,n),4) * exp(-4 * gconf(m,n))) / (pow2(gA(m,n)) * r(m,n))) 
    /*0*/  + gAlp(m,n) * (3 * pow2(gA1(m,n)) - (gA1(m,n) * gAsig(m,n) * pow2(r(m,n))) 
    /*1*/  / pow2(r_minus(m,n)) + pow2(gtrK(m,n)) / 3. + (pow2(gAsig(m,n)) * Power(r(m,n),4))
    /*1*/  / (8. * Power(r_minus(m,n),4))) + 2 * gBet(m,n) * pow2(r_minus(m,n)) * gtrK_r(m,n);
}
Real ftrK_t( Int m, Int n )
{
    return k_f * fJK(m,n) - (4 * Power(r_minus(m,n),4) * exp(-4 * fconf(m,n)) 
    /*1*/  * fAlp_rr(m,n)) / pow2(fA(m,n)) + fAlp_r(m,n) * ((-8 * Power(r_minus(m,n),4) 
    /*2*/  * exp(-4 * fconf(m,n)) * fconf_r(m,n)) / pow2(fA(m,n)) - (8 * Power(r_minus(m,n),4)
    /*2*/  * exp(-4 * fconf(m,n)) * fB_r(m,n)) / (fB(m,n) * pow2(fA(m,n))) + (4 
    /*2*/  * Power(r_minus(m,n),4) * exp(-4 * fconf(m,n)) * fA_r(m,n)) / pow3(fA(m,n)) - (8 
    /*2*/  * Power(r_minus(m,n),4) * exp(-4 * fconf(m,n))) / (pow2(fA(m,n)) * r(m,n))) 
    /*0*/  + fAlp(m,n) * (3 * pow2(fA1(m,n)) - (fA1(m,n) * fAsig(m,n) * pow2(r(m,n))) 
    /*1*/  / pow2(r_minus(m,n)) + pow2(ftrK(m,n)) / 3. + (pow2(fAsig(m,n)) * Power(r(m,n),4))
    /*1*/  / (8. * Power(r_minus(m,n),4))) + 2 * fBet(m,n) * pow2(r_minus(m,n)) * ftrK_r(m,n);
}
Real gA_t( Int m, Int n )
{
    return -(gA1(m,n) * gAlp(m,n) * gA(m,n)) + 2 * gBet_r(m,n) * gA(m,n) 
    /*0*/  * pow2(r_minus(m,n)) + 2 * gBet(m,n) * gA_r(m,n) * pow2(r_minus(m,n));
}
Real gB_t( Int m, Int n )
{
    return gBetr(m,n) * gB(m,n) + 2 * gBet(m,n) * gB_r(m,n) * pow2(r_minus(m,n)) 
    /*0*/  + (gAlp(m,n) * gB(m,n) * (-12 * gA1(m,n) + (3 * gAsig(m,n) * pow2(1 
    /*4*/  + r_minus(m,n))) / pow2(r_minus(m,n)))) / 12.;
}
Real fA_t( Int m, Int n )
{
    return -(fA1(m,n) * fAlp(m,n) * fA(m,n)) + 2 * fBet_r(m,n) * fA(m,n) 
    /*0*/  * pow2(r_minus(m,n)) + 2 * fBet(m,n) * fA_r(m,n) * pow2(r_minus(m,n));
}
Real fB_t( Int m, Int n )
{
    return fBetr(m,n) * fB(m,n) + 2 * fBet(m,n) * fB_r(m,n) * pow2(r_minus(m,n)) 
    /*0*/  + (fAlp(m,n) * fB(m,n) * (-12 * fA1(m,n) + (3 * fAsig(m,n) * pow2(1 
    /*4*/  + r_minus(m,n))) / pow2(r_minus(m,n)))) / 12.;
}
Real gA1_t( Int m, Int n )
{
    return gRicci(m,n) + k_g * gJA1(m,n) + 2 * gA1_r(m,n) * gBet(m,n) * pow2(r_minus(m,n)) 
    /*0*/  + exp(-4 * gconf(m,n)) * ((-8 * Power(r_minus(m,n),4) * gAlp_rr(m,n)) / (3. 
    /*2*/  * pow2(gA(m,n))) + gAlp(m,n) * ((-16 * Power(r_minus(m,n),4) * gconf_rr(m,n)) / (3.
    /*3*/  * pow2(gA(m,n))) + (32 * Power(r_minus(m,n),4) * pow2(gconf_r(m,n))) / (3. 
    /*3*/  * pow2(gA(m,n))) + gconf_r(m,n) * ((16 * Power(r_minus(m,n),4) * gB_r(m,n)) / (3. 
    /*4*/  * gB(m,n) * pow2(gA(m,n))) + (16 * Power(r_minus(m,n),4) * gA_r(m,n)) / (3. 
    /*4*/  * pow3(gA(m,n))) - (16 * (3 + 2 * r_minus(m,n)) * pow3(r_minus(m,n))) / (3. 
    /*4*/  * pow2(gA(m,n)) * r(m,n)))) + gAlp_r(m,n) * ((32 * Power(r_minus(m,n),4) 
    /*3*/  * gconf_r(m,n)) / (3. * pow2(gA(m,n))) + (8 * Power(r_minus(m,n),4) * gB_r(m,n)) 
    /*2*/  / (3. * gB(m,n) * pow2(gA(m,n))) + (8 * Power(r_minus(m,n),4) * gA_r(m,n)) / (3. 
    /*3*/  * pow3(gA(m,n))) - (8 * (3 + 2 * r_minus(m,n)) * pow3(r_minus(m,n))) / (3. 
    /*3*/  * pow2(gA(m,n)) * r(m,n)))) + gA1(m,n) * gAlp(m,n) * gtrK(m,n);
}
Real fA1_t( Int m, Int n )
{
    return fRicci(m,n) + k_f * fJA1(m,n) + 2 * fA1_r(m,n) * fBet(m,n) * pow2(r_minus(m,n)) 
    /*0*/  + exp(-4 * fconf(m,n)) * ((-8 * Power(r_minus(m,n),4) * fAlp_rr(m,n)) / (3. 
    /*2*/  * pow2(fA(m,n))) + fAlp(m,n) * ((-16 * Power(r_minus(m,n),4) * fconf_rr(m,n)) / (3.
    /*3*/  * pow2(fA(m,n))) + (32 * Power(r_minus(m,n),4) * pow2(fconf_r(m,n))) / (3. 
    /*3*/  * pow2(fA(m,n))) + fconf_r(m,n) * ((16 * Power(r_minus(m,n),4) * fB_r(m,n)) / (3. 
    /*4*/  * fB(m,n) * pow2(fA(m,n))) + (16 * Power(r_minus(m,n),4) * fA_r(m,n)) / (3. 
    /*4*/  * pow3(fA(m,n))) - (16 * (3 + 2 * r_minus(m,n)) * pow3(r_minus(m,n))) / (3. 
    /*4*/  * pow2(fA(m,n)) * r(m,n)))) + fAlp_r(m,n) * ((32 * Power(r_minus(m,n),4) 
    /*3*/  * fconf_r(m,n)) / (3. * pow2(fA(m,n))) + (8 * Power(r_minus(m,n),4) * fB_r(m,n)) 
    /*2*/  / (3. * fB(m,n) * pow2(fA(m,n))) + (8 * Power(r_minus(m,n),4) * fA_r(m,n)) / (3. 
    /*3*/  * pow3(fA(m,n))) - (8 * (3 + 2 * r_minus(m,n)) * pow3(r_minus(m,n))) / (3. 
    /*3*/  * pow2(fA(m,n)) * r(m,n)))) + fA1(m,n) * fAlp(m,n) * ftrK(m,n);
}
Real gL_t( Int m, Int n )
{
    return k_g * gJL(m,n) + 2 * gBet(m,n) * gL_r(m,n) * pow2(r_minus(m,n)) - (4 * gA1(m,n) 
    /*1*/  * gAlp_r(m,n) * pow2(r_minus(m,n))) / pow2(gA(m,n)) + (pow2(r_minus(m,n)) * ((12 
    /*3*/  * gBet_rr(m,n) * pow2(r_minus(m,n))) / pow2(gA(m,n)) + (12 * gBetr_r(m,n)) 
    /*2*/  / pow2(gB(m,n)))) / 3. + gBet_r(m,n) * (-2 * gL(m,n) * pow2(r_minus(m,n)) + (8 
    /*2*/  * pow3(r_minus(m,n))) / pow2(gA(m,n))) + gAlp(m,n) * ((4 * gAsig(m,n) 
    /*2*/  * gconf_r(m,n) * pow2(1 + r_minus(m,n))) / pow2(gA(m,n)) - (gAsig(m,n) * gsig(m,n)
    /*2*/  * pow3(1 + r_minus(m,n))) / (6. * pow2(gA(m,n)) * pow3(r_minus(m,n))) + (2 
    /*2*/  * (gAsig(m,n) * (gA_r(m,n) * gB(m,n) + gA(m,n) * gB_r(m,n)) * pow2(1 
    /*4*/  + r_minus(m,n)) - 4 * gA(m,n) * gB(m,n) * pow2(r_minus(m,n)) * gtrK_r(m,n))) / (3. 
    /*2*/  * gB(m,n) * pow3(gA(m,n))));
}
Real fL_t( Int m, Int n )
{
    return k_f * fJL(m,n) + 2 * fBet(m,n) * fL_r(m,n) * pow2(r_minus(m,n)) - (4 * fA1(m,n) 
    /*1*/  * fAlp_r(m,n) * pow2(r_minus(m,n))) / pow2(fA(m,n)) + (pow2(r_minus(m,n)) * ((12 
    /*3*/  * fBet_rr(m,n) * pow2(r_minus(m,n))) / pow2(fA(m,n)) + (12 * fBetr_r(m,n)) 
    /*2*/  / pow2(fB(m,n)))) / 3. + fBet_r(m,n) * (-2 * fL(m,n) * pow2(r_minus(m,n)) + (8 
    /*2*/  * pow3(r_minus(m,n))) / pow2(fA(m,n))) + fAlp(m,n) * ((4 * fAsig(m,n) 
    /*2*/  * fconf_r(m,n) * pow2(1 + r_minus(m,n))) / pow2(fA(m,n)) - (fAsig(m,n) * fsig(m,n)
    /*2*/  * pow3(1 + r_minus(m,n))) / (6. * pow2(fA(m,n)) * pow3(r_minus(m,n))) + (2 
    /*2*/  * (fAsig(m,n) * (fA_r(m,n) * fB(m,n) + fA(m,n) * fB_r(m,n)) * pow2(1 
    /*4*/  + r_minus(m,n)) - 4 * fA(m,n) * fB(m,n) * pow2(r_minus(m,n)) * ftrK_r(m,n))) / (3. 
    /*2*/  * fB(m,n) * pow3(fA(m,n))));
}
Real gsig_t( Int m, Int n )
{
    return 2 * gBetr(m,n) * gsig(m,n) + 2 * gBet(m,n) * gsig_r(m,n) * pow2(r_minus(m,n)) 
    /*0*/  + (2 * gAsig(m,n) * gAlp(m,n) * pow2(gA(m,n))) / pow2(gB(m,n)) + (8 
    /*1*/  * gBetr_r(m,n) * pow2(gA(m,n)) * pow3(r_minus(m,n))) / (pow2(gB(m,n)) * r(m,n));
}
Real fsig_t( Int m, Int n )
{
    return 2 * fBetr(m,n) * fsig(m,n) + 2 * fBet(m,n) * fsig_r(m,n) * pow2(r_minus(m,n)) 
    /*0*/  + (2 * fAsig(m,n) * fAlp(m,n) * pow2(fA(m,n))) / pow2(fB(m,n)) + (8 
    /*1*/  * fBetr_r(m,n) * pow2(fA(m,n)) * pow3(r_minus(m,n))) / (pow2(fB(m,n)) * r(m,n));
}
Real gAsig_t( Int m, Int n )
{
    return (exp(-4 * gconf(m,n)) * gAlp(m,n) * pow2(gsig(m,n)) * pow2(gB(m,n))) 
    /*0*/  / Power(gA(m,n),4) + (6 * k_g * gJA1(m,n) * pow2(r_minus(m,n))) / pow2(r(m,n)) 
    /*0*/  + (Power(r_minus(m,n),5) * exp(-4 * gconf(m,n)) * gAlp(m,n) * ((4 * gsig_r(m,n) 
    /*3*/  * pow2(gB(m,n))) / Power(gA(m,n),4) + (4 * gsig_rr(m,n) * pow2(gB(m,n))) 
    /*2*/  / Power(gA(m,n),4)) + pow2(r_minus(m,n)) * (2 * gAsig_r(m,n) * gBet(m,n) - (exp(-4
    /*4*/  * gconf(m,n)) * gAlp(m,n) * gL(m,n) * gsig_r(m,n) * pow2(gB(m,n))) 
    /*2*/  / pow2(gA(m,n))) + Power(r_minus(m,n),4) * (2 * gAsig_r(m,n) * gBet(m,n) + exp(-4 
    /*3*/  * gconf(m,n)) * ((4 * gderAlpr_r(m,n)) / pow2(gA(m,n)) + gAlp(m,n) * (-4 
    /*4*/  * gLr_r(m,n) + (8 * gderconfr_r(m,n)) / pow2(gA(m,n)) + (2 * gsig_rr(m,n) 
    /*5*/  * pow2(gB(m,n))) / Power(gA(m,n),4) + gsig_r(m,n) * (-4 / pow2(gA(m,n)) 
    /*5*/  - (gL(m,n) * pow2(gB(m,n))) / pow2(gA(m,n)))))) + (4 * gAsig_r(m,n) 
    /*2*/  * gBet(m,n) + exp(-4 * gconf(m,n)) * ((4 * gderAlpr_r(m,n)) / pow2(gA(m,n)) 
    /*3*/  + gAlp(m,n) * (-4 * gLr_r(m,n) + (8 * gderconfr_r(m,n)) / pow2(gA(m,n)) 
    /*4*/  + gsig_r(m,n) * (-4 / pow2(gA(m,n)) - (2 * (2 + gL(m,n) * pow2(gA(m,n))) 
    /*6*/  * pow2(gB(m,n))) / Power(gA(m,n),4))))) * pow3(r_minus(m,n)) + Power(r_minus(m,n),6) 
    /*1*/  * exp(-4 * gconf(m,n)) * (gAlp_r(m,n) * ((64 * gconf_r(m,n)) / pow2(gA(m,n))
    /*3*/  + (16 * gB_r(m,n)) / (gB(m,n) * pow2(gA(m,n))) + (16 * gA_r(m,n)) 
    /*3*/  / pow3(gA(m,n))) + gAlp(m,n) * ((64 * pow2(gconf_r(m,n))) / pow2(gA(m,n)) 
    /*3*/  + (48 * pow2(gA_r(m,n))) / Power(gA(m,n),4) + (2 * gsig_rr(m,n) 
    /*4*/  * pow2(gB(m,n))) / Power(gA(m,n),4) + gconf_r(m,n) * ((32 * gB_r(m,n)) 
    /*4*/  / (gB(m,n) * pow2(gA(m,n))) + (32 * gA_r(m,n)) / pow3(gA(m,n))) - (64 
    /*4*/  * gA_r(m,n) * gB_r(m,n)) / (gB(m,n) * pow3(gA(m,n)))))) / pow2(r(m,n)) + ((16
    /*2*/  * exp(-4 * gconf(m,n)) * gAlp(m,n) * gsig(m,n) * gB_r(m,n) * pow3(r_minus(m,n)))
    /*1*/  / (gB(m,n) * pow2(gA(m,n))) + gAsig(m,n) * gAlp(m,n) * gtrK(m,n) + r_minus(m,n) 
    /*1*/  * (-4 * gAsig(m,n) * gBet(m,n) + (2 * exp(-4 * gconf(m,n)) * gAlp(m,n) 
    /*3*/  * gL(m,n) * gsig(m,n) * pow2(gB(m,n))) / pow2(gA(m,n)) + gAsig(m,n) 
    /*2*/  * gAlp(m,n) * gtrK(m,n))) / r(m,n);
}
Real fAsig_t( Int m, Int n )
{
    return (exp(-4 * fconf(m,n)) * fAlp(m,n) * pow2(fsig(m,n)) * pow2(fB(m,n))) 
    /*0*/  / Power(fA(m,n),4) + (6 * k_f * fJA1(m,n) * pow2(r_minus(m,n))) / pow2(r(m,n)) 
    /*0*/  + (Power(r_minus(m,n),5) * exp(-4 * fconf(m,n)) * fAlp(m,n) * ((4 * fsig_r(m,n) 
    /*3*/  * pow2(fB(m,n))) / Power(fA(m,n),4) + (4 * fsig_rr(m,n) * pow2(fB(m,n))) 
    /*2*/  / Power(fA(m,n),4)) + pow2(r_minus(m,n)) * (2 * fAsig_r(m,n) * fBet(m,n) - (exp(-4
    /*4*/  * fconf(m,n)) * fAlp(m,n) * fL(m,n) * fsig_r(m,n) * pow2(fB(m,n))) 
    /*2*/  / pow2(fA(m,n))) + Power(r_minus(m,n),4) * (2 * fAsig_r(m,n) * fBet(m,n) + exp(-4 
    /*3*/  * fconf(m,n)) * ((4 * fderAlpr_r(m,n)) / pow2(fA(m,n)) + fAlp(m,n) * (-4 
    /*4*/  * fLr_r(m,n) + (8 * fderconfr_r(m,n)) / pow2(fA(m,n)) + (2 * fsig_rr(m,n) 
    /*5*/  * pow2(fB(m,n))) / Power(fA(m,n),4) + fsig_r(m,n) * (-4 / pow2(fA(m,n)) 
    /*5*/  - (fL(m,n) * pow2(fB(m,n))) / pow2(fA(m,n)))))) + (4 * fAsig_r(m,n) 
    /*2*/  * fBet(m,n) + exp(-4 * fconf(m,n)) * ((4 * fderAlpr_r(m,n)) / pow2(fA(m,n)) 
    /*3*/  + fAlp(m,n) * (-4 * fLr_r(m,n) + (8 * fderconfr_r(m,n)) / pow2(fA(m,n)) 
    /*4*/  + fsig_r(m,n) * (-4 / pow2(fA(m,n)) - (2 * (2 + fL(m,n) * pow2(fA(m,n))) 
    /*6*/  * pow2(fB(m,n))) / Power(fA(m,n),4))))) * pow3(r_minus(m,n)) + Power(r_minus(m,n),6) 
    /*1*/  * exp(-4 * fconf(m,n)) * (fAlp_r(m,n) * ((64 * fconf_r(m,n)) / pow2(fA(m,n))
    /*3*/  + (16 * fB_r(m,n)) / (fB(m,n) * pow2(fA(m,n))) + (16 * fA_r(m,n)) 
    /*3*/  / pow3(fA(m,n))) + fAlp(m,n) * ((64 * pow2(fconf_r(m,n))) / pow2(fA(m,n)) 
    /*3*/  + (48 * pow2(fA_r(m,n))) / Power(fA(m,n),4) + (2 * fsig_rr(m,n) 
    /*4*/  * pow2(fB(m,n))) / Power(fA(m,n),4) + fconf_r(m,n) * ((32 * fB_r(m,n)) 
    /*4*/  / (fB(m,n) * pow2(fA(m,n))) + (32 * fA_r(m,n)) / pow3(fA(m,n))) - (64 
    /*4*/  * fA_r(m,n) * fB_r(m,n)) / (fB(m,n) * pow3(fA(m,n)))))) / pow2(r(m,n)) + ((16
    /*2*/  * exp(-4 * fconf(m,n)) * fAlp(m,n) * fsig(m,n) * fB_r(m,n) * pow3(r_minus(m,n)))
    /*1*/  / (fB(m,n) * pow2(fA(m,n))) + fAsig(m,n) * fAlp(m,n) * ftrK(m,n) + r_minus(m,n) 
    /*1*/  * (-4 * fAsig(m,n) * fBet(m,n) + (2 * exp(-4 * fconf(m,n)) * fAlp(m,n) 
    /*3*/  * fL(m,n) * fsig(m,n) * pow2(fB(m,n))) / pow2(fA(m,n)) + fAsig(m,n) 
    /*2*/  * fAlp(m,n) * ftrK(m,n))) / r(m,n);
}
