/** @file  eomBSSNEvolutionCompEul.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN Eulerian evolution equations.
 *  @version 2019-05-15T16:16:55
 *  @image html BSSNevolutionCompEul.png
 */

Real gconf_t( Int m, Int n )
{
    return (gconf_r(m,n) * gBet(m,n) * pow2(r_minus(m,n))) / r_plus(m,n) - (gAlp(m,n) * gtrK(m,n))
    /*0*/  / 6.;
}
Real fconf_t( Int m, Int n )
{
    return (fconf_r(m,n) * fBet(m,n) * pow2(r_minus(m,n))) / r_plus(m,n) - (fAlp(m,n) * ftrK(m,n))
    /*0*/  / 6.;
}
Real gtrK_t( Int m, Int n )
{
    return k_g * gJK(m,n) - (Power(r_minus(m,n),4) * exp(-4 * gconf(m,n)) * gAlp_rr(m,n)) 
    /*0*/  / (pow2(r_plus(m,n)) * pow2(gA(m,n))) + gAlp_r(m,n) * ((-2 * Power(r_minus(m,n),4) 
    /*2*/  * exp(-4 * gconf(m,n)) * gconf_r(m,n)) / (pow2(r_plus(m,n)) * pow2(gA(m,n))) - (2 
    /*2*/  * Power(r_minus(m,n),4) * exp(-4 * gconf(m,n)) * gB_r(m,n)) / (gB(m,n) 
    /*2*/  * pow2(r_plus(m,n)) * pow2(gA(m,n))) + (Power(r_minus(m,n),4) * exp(-4 * gconf(m,n)) 
    /*2*/  * gA_r(m,n)) / (pow2(r_plus(m,n)) * pow3(gA(m,n))) - (2 * Power(r_minus(m,n),4) * exp(-4
    /*3*/  * gconf(m,n))) / (pow2(gA(m,n)) * pow3(r_plus(m,n)) * r(m,n))) + gAlp(m,n) * (3 
    /*1*/  * pow2(gA1(m,n)) - (4 * gA1(m,n) * gAsig(m,n) * pow2(r(m,n))) / pow2(r_minus(m,n))
    /*1*/  + pow2(gtrK(m,n)) / 3. + (2 * pow2(gAsig(m,n)) * Power(r(m,n),4)) 
    /*1*/  / Power(r_minus(m,n),4)) + (gBet(m,n) * pow2(r_minus(m,n)) * gtrK_r(m,n)) / r_plus(m,n);
}
Real ftrK_t( Int m, Int n )
{
    return k_f * fJK(m,n) - (Power(r_minus(m,n),4) * exp(-4 * fconf(m,n)) * fAlp_rr(m,n)) 
    /*0*/  / (pow2(r_plus(m,n)) * pow2(fA(m,n))) + fAlp_r(m,n) * ((-2 * Power(r_minus(m,n),4) 
    /*2*/  * exp(-4 * fconf(m,n)) * fconf_r(m,n)) / (pow2(r_plus(m,n)) * pow2(fA(m,n))) - (2 
    /*2*/  * Power(r_minus(m,n),4) * exp(-4 * fconf(m,n)) * fB_r(m,n)) / (fB(m,n) 
    /*2*/  * pow2(r_plus(m,n)) * pow2(fA(m,n))) + (Power(r_minus(m,n),4) * exp(-4 * fconf(m,n)) 
    /*2*/  * fA_r(m,n)) / (pow2(r_plus(m,n)) * pow3(fA(m,n))) - (2 * Power(r_minus(m,n),4) * exp(-4
    /*3*/  * fconf(m,n))) / (pow2(fA(m,n)) * pow3(r_plus(m,n)) * r(m,n))) + fAlp(m,n) * (3 
    /*1*/  * pow2(fA1(m,n)) - (4 * fA1(m,n) * fAsig(m,n) * pow2(r(m,n))) / pow2(r_minus(m,n))
    /*1*/  + pow2(ftrK(m,n)) / 3. + (2 * pow2(fAsig(m,n)) * Power(r(m,n),4)) 
    /*1*/  / Power(r_minus(m,n),4)) + (fBet(m,n) * pow2(r_minus(m,n)) * ftrK_r(m,n)) / r_plus(m,n);
}
Real gA_t( Int m, Int n )
{
    return -(gA1(m,n) * gAlp(m,n) * gA(m,n)) + (gBet_r(m,n) * gA(m,n) * pow2(r_minus(m,n)))
    /*0*/  / r_plus(m,n) + (gBet(m,n) * gA_r(m,n) * pow2(r_minus(m,n))) / r_plus(m,n);
}
Real gB_t( Int m, Int n )
{
    return gBetr(m,n) * gB(m,n) + (gBet(m,n) * gB_r(m,n) * pow2(r_minus(m,n))) / r_plus(m,n) 
    /*0*/  + (gAlp(m,n) * gB(m,n) * (-3 * gA1(m,n) + (3 * gAsig(m,n) * pow2(r(m,n))) 
    /*2*/  / pow2(r_minus(m,n)))) / 3.;
}
Real fA_t( Int m, Int n )
{
    return -(fA1(m,n) * fAlp(m,n) * fA(m,n)) + (fBet_r(m,n) * fA(m,n) * pow2(r_minus(m,n)))
    /*0*/  / r_plus(m,n) + (fBet(m,n) * fA_r(m,n) * pow2(r_minus(m,n))) / r_plus(m,n);
}
Real fB_t( Int m, Int n )
{
    return fBetr(m,n) * fB(m,n) + (fBet(m,n) * fB_r(m,n) * pow2(r_minus(m,n))) / r_plus(m,n) 
    /*0*/  + (fAlp(m,n) * fB(m,n) * (-3 * fA1(m,n) + (3 * fAsig(m,n) * pow2(r(m,n))) 
    /*2*/  / pow2(r_minus(m,n)))) / 3.;
}
Real gA1_t( Int m, Int n )
{
    return gRicci(m,n) + k_g * gJA1(m,n) + (gA1_r(m,n) * gBet(m,n) * pow2(r_minus(m,n))) 
    /*0*/  / r_plus(m,n) + exp(-4 * gconf(m,n)) * ((-2 * Power(r_minus(m,n),4) * gAlp_rr(m,n)) / (3.
    /*2*/  * pow2(r_plus(m,n)) * pow2(gA(m,n))) + gAlp_r(m,n) * ((8 * Power(r_minus(m,n),4) 
    /*3*/  * gconf_r(m,n)) / (3. * pow2(r_plus(m,n)) * pow2(gA(m,n))) + (2 * Power(r_minus(m,n),4) 
    /*3*/  * gB_r(m,n)) / (3. * gB(m,n) * pow2(r_plus(m,n)) * pow2(gA(m,n))) + (2 
    /*3*/  * Power(r_minus(m,n),4) * gA_r(m,n)) / (3. * pow2(r_plus(m,n)) * pow3(gA(m,n))) - (2 
    /*3*/  * pow3(r_minus(m,n)) * (1 + 8 * pow2(r(m,n)) + 3 * Power(r(m,n),4))) / (3. 
    /*3*/  * pow2(gA(m,n)) * pow3(r_plus(m,n)) * r(m,n))) + gAlp(m,n) * ((-4 * Power(r_minus(m,n),4)
    /*3*/  * gconf_rr(m,n)) / (3. * pow2(r_plus(m,n)) * pow2(gA(m,n))) + (8 
    /*3*/  * Power(r_minus(m,n),4) * pow2(gconf_r(m,n))) / (3. * pow2(r_plus(m,n)) * pow2(gA(m,n)))
    /*2*/  + gconf_r(m,n) * ((4 * Power(r_minus(m,n),4) * gB_r(m,n)) / (3. * gB(m,n) 
    /*4*/  * pow2(r_plus(m,n)) * pow2(gA(m,n))) + (4 * Power(r_minus(m,n),4) * gA_r(m,n)) / (3. 
    /*4*/  * pow2(r_plus(m,n)) * pow3(gA(m,n))) - (4 * pow3(r_minus(m,n)) * (1 + 8 * pow2(r(m,n)) 
    /*5*/  + 3 * Power(r(m,n),4))) / (3. * pow2(gA(m,n)) * pow3(r_plus(m,n)) * r(m,n))))) 
    /*0*/  + gA1(m,n) * gAlp(m,n) * gtrK(m,n);
}
Real fA1_t( Int m, Int n )
{
    return fRicci(m,n) + k_f * fJA1(m,n) + (fA1_r(m,n) * fBet(m,n) * pow2(r_minus(m,n))) 
    /*0*/  / r_plus(m,n) + exp(-4 * fconf(m,n)) * ((-2 * Power(r_minus(m,n),4) * fAlp_rr(m,n)) / (3.
    /*2*/  * pow2(r_plus(m,n)) * pow2(fA(m,n))) + fAlp_r(m,n) * ((8 * Power(r_minus(m,n),4) 
    /*3*/  * fconf_r(m,n)) / (3. * pow2(r_plus(m,n)) * pow2(fA(m,n))) + (2 * Power(r_minus(m,n),4) 
    /*3*/  * fB_r(m,n)) / (3. * fB(m,n) * pow2(r_plus(m,n)) * pow2(fA(m,n))) + (2 
    /*3*/  * Power(r_minus(m,n),4) * fA_r(m,n)) / (3. * pow2(r_plus(m,n)) * pow3(fA(m,n))) - (2 
    /*3*/  * pow3(r_minus(m,n)) * (1 + 8 * pow2(r(m,n)) + 3 * Power(r(m,n),4))) / (3. 
    /*3*/  * pow2(fA(m,n)) * pow3(r_plus(m,n)) * r(m,n))) + fAlp(m,n) * ((-4 * Power(r_minus(m,n),4)
    /*3*/  * fconf_rr(m,n)) / (3. * pow2(r_plus(m,n)) * pow2(fA(m,n))) + (8 
    /*3*/  * Power(r_minus(m,n),4) * pow2(fconf_r(m,n))) / (3. * pow2(r_plus(m,n)) * pow2(fA(m,n)))
    /*2*/  + fconf_r(m,n) * ((4 * Power(r_minus(m,n),4) * fB_r(m,n)) / (3. * fB(m,n) 
    /*4*/  * pow2(r_plus(m,n)) * pow2(fA(m,n))) + (4 * Power(r_minus(m,n),4) * fA_r(m,n)) / (3. 
    /*4*/  * pow2(r_plus(m,n)) * pow3(fA(m,n))) - (4 * pow3(r_minus(m,n)) * (1 + 8 * pow2(r(m,n)) 
    /*5*/  + 3 * Power(r(m,n),4))) / (3. * pow2(fA(m,n)) * pow3(r_plus(m,n)) * r(m,n))))) 
    /*0*/  + fA1(m,n) * fAlp(m,n) * ftrK(m,n);
}
Real gL_t( Int m, Int n )
{
    return (6 * k_g * gJL(m,n) - (6 * gBet_r(m,n) * gL(m,n) * pow2(r_minus(m,n))) / r_plus(m,n) 
    /*1*/  + (6 * gBet(m,n) * gL_r(m,n) * pow2(r_minus(m,n))) / r_plus(m,n) + (4 * gA1(m,n) * (-3 
    /*3*/  * gAlp_r(m,n) * gA(m,n) * gB(m,n) + 2 * gAlp(m,n) * (gA_r(m,n) * gB(m,n) 
    /*4*/  + gA(m,n) * (6 * gconf_r(m,n) * gB(m,n) + gB_r(m,n)))) * pow2(r_minus(m,n))) 
    /*1*/  / (r_plus(m,n) * gB(m,n) * pow3(gA(m,n))) - (8 * gAlp(m,n) * (gA_r(m,n) * gB(m,n) 
    /*3*/  + gA(m,n) * (6 * gconf_r(m,n) * gB(m,n) + gB_r(m,n))) * (gA1(m,n) 
    /*3*/  * pow2(r_minus(m,n)) - gAsig(m,n) * pow2(r(m,n)))) / (r_plus(m,n) * gB(m,n) 
    /*2*/  * pow3(gA(m,n))) + (12 * pow2(r_minus(m,n)) * (gBet_r(m,n) * (-1 + pow2(r(m,n))) 
    /*3*/  + gBetr_r(m,n) * r(m,n))) / (pow2(gB(m,n)) * (r(m,n) + Power(r(m,n),3))) + (2
    /*2*/  * (1 - pow2(r(m,n))) * (6 * Power(r_minus(m,n),6) * gBet_r(m,n) * pow2(gA(m,n)) 
    /*3*/  + 4 * gAsig(m,n) * gAlp(m,n) * gsig(m,n) * pow2(gB(m,n)) * (1 + pow2(r(m,n)))
    /*3*/  * Power(r(m,n),4))) / (Power(r_minus(m,n),4) * pow2(gA(m,n)) * pow2(gB(m,n)) 
    /*2*/  * (r(m,n) + Power(r(m,n),3))) - (2 * pow2(r_minus(m,n)) * (-6 * gBet_r(m,n) 
    /*3*/  * r(m,n) * (-3 + 2 * pow2(r(m,n)) + Power(r(m,n),4)) + (1 + pow2(r(m,n))) 
    /*3*/  * (-3 * gBet_rr(m,n) * pow2(r_minus(m,n)) + 4 * gAlp(m,n) * (1 + pow2(r(m,n))) 
    /*4*/  * gtrK_r(m,n)))) / (pow2(gA(m,n)) * pow3(r_plus(m,n)))) / 6.;
}
Real fL_t( Int m, Int n )
{
    return (6 * k_f * fJL(m,n) - (6 * fBet_r(m,n) * fL(m,n) * pow2(r_minus(m,n))) / r_plus(m,n) 
    /*1*/  + (6 * fBet(m,n) * fL_r(m,n) * pow2(r_minus(m,n))) / r_plus(m,n) + (4 * fA1(m,n) * (-3 
    /*3*/  * fAlp_r(m,n) * fA(m,n) * fB(m,n) + 2 * fAlp(m,n) * (fA_r(m,n) * fB(m,n) 
    /*4*/  + fA(m,n) * (6 * fconf_r(m,n) * fB(m,n) + fB_r(m,n)))) * pow2(r_minus(m,n))) 
    /*1*/  / (r_plus(m,n) * fB(m,n) * pow3(fA(m,n))) - (8 * fAlp(m,n) * (fA_r(m,n) * fB(m,n) 
    /*3*/  + fA(m,n) * (6 * fconf_r(m,n) * fB(m,n) + fB_r(m,n))) * (fA1(m,n) 
    /*3*/  * pow2(r_minus(m,n)) - fAsig(m,n) * pow2(r(m,n)))) / (r_plus(m,n) * fB(m,n) 
    /*2*/  * pow3(fA(m,n))) + (12 * pow2(r_minus(m,n)) * (fBet_r(m,n) * (-1 + pow2(r(m,n))) 
    /*3*/  + fBetr_r(m,n) * r(m,n))) / (pow2(fB(m,n)) * (r(m,n) + Power(r(m,n),3))) + (2
    /*2*/  * (1 - pow2(r(m,n))) * (6 * Power(r_minus(m,n),6) * fBet_r(m,n) * pow2(fA(m,n)) 
    /*3*/  + 4 * fAsig(m,n) * fAlp(m,n) * fsig(m,n) * pow2(fB(m,n)) * (1 + pow2(r(m,n)))
    /*3*/  * Power(r(m,n),4))) / (Power(r_minus(m,n),4) * pow2(fA(m,n)) * pow2(fB(m,n)) 
    /*2*/  * (r(m,n) + Power(r(m,n),3))) - (2 * pow2(r_minus(m,n)) * (-6 * fBet_r(m,n) 
    /*3*/  * r(m,n) * (-3 + 2 * pow2(r(m,n)) + Power(r(m,n),4)) + (1 + pow2(r(m,n))) 
    /*3*/  * (-3 * fBet_rr(m,n) * pow2(r_minus(m,n)) + 4 * fAlp(m,n) * (1 + pow2(r(m,n))) 
    /*4*/  * ftrK_r(m,n)))) / (pow2(fA(m,n)) * pow3(r_plus(m,n)))) / 6.;
}
Real gsig_t( Int m, Int n )
{
    return 2 * gBetr(m,n) * gsig(m,n) + (gBet(m,n) * gsig_r(m,n) * pow2(r_minus(m,n))) 
    /*0*/  / r_plus(m,n) + (pow2(gA(m,n)) * (2 * gAsig(m,n) * gAlp(m,n) + (2 * gBetr_r(m,n) 
    /*3*/  * pow3(r_minus(m,n))) / (r(m,n) + Power(r(m,n),3)))) / pow2(gB(m,n));
}
Real fsig_t( Int m, Int n )
{
    return 2 * fBetr(m,n) * fsig(m,n) + (fBet(m,n) * fsig_r(m,n) * pow2(r_minus(m,n))) 
    /*0*/  / r_plus(m,n) + (2 * fAsig(m,n) * fAlp(m,n) * pow2(fA(m,n))) / pow2(fB(m,n)) + (2 
    /*1*/  * fBetr_r(m,n) * pow2(fA(m,n)) * pow3(r_minus(m,n))) / (pow2(fB(m,n)) * (r(m,n) 
    /*2*/  + Power(r(m,n),3)));
}
Real gAsig_t( Int m, Int n )
{
    return ((1 - pow2(r(m,n))) * (-3 * k_g * gJA1(m,n) * (-1 + pow2(r(m,n))) - (2 
    /*3*/  * exp(-4 * gconf(m,n)) * gB_r(m,n) * pow2(r_minus(m,n)) * (gAlp_r(m,n) * gA(m,n) 
    /*4*/  * pow3(r_minus(m,n)) + 2 * gAlp(m,n) * ((gconf_r(m,n) * gA(m,n) - 2 * gA_r(m,n)) 
    /*5*/  * pow3(r_minus(m,n)) + 2 * gsig(m,n) * gA(m,n) * (1 + pow2(r(m,n))) * r(m,n)))) 
    /*2*/  / (gB(m,n) * pow2(r_plus(m,n)) * pow3(gA(m,n))) + (exp(-4 * gconf(m,n)) * gAlp(m,n)
    /*3*/  * pow2(gB(m,n)) * r(m,n) * (-2 * gL(m,n) * gsig(m,n) * pow2(gA(m,n)) - (2 
    /*5*/  * pow2(gsig(m,n)) * r(m,n)) / r_minus(m,n) + (gsig_r(m,n) * (-1 + pow2(r(m,n))) 
    /*5*/  * (-4 + 4 * pow2(r(m,n)) + gL(m,n) * pow2(gA(m,n)) * r(m,n))) / r_plus(m,n) 
    /*4*/  - (pow2(r_minus(m,n)) * r(m,n) * (2 * gsig_r(m,n) * (3 + pow2(r(m,n))) * r(m,n) 
    /*6*/  + gsig_rr(m,n) * (-1 + Power(r(m,n),4)))) / pow3(r_plus(m,n)))) / Power(gA(m,n),4)
    /*2*/  + (2 * exp(-4 * gconf(m,n)) * (1 - pow2(r(m,n))) * (-((gAlp(m,n) * (-1 
    /*7*/  + pow2(r(m,n))) * (-2 * gconf_r(m,n) * gA(m,n) * gA_r(m,n) * pow3(r_minus(m,n)) - 3
    /*7*/  * pow2(gA_r(m,n)) * pow3(r_minus(m,n)) + gLr_r(m,n) * Power(gA(m,n),4) * (1 
    /*8*/  + pow2(r(m,n))) * r(m,n) + pow2(gA(m,n)) * (-4 * pow2(gconf_r(m,n)) 
    /*8*/  * pow3(r_minus(m,n)) - (2 * gderconfr_r(m,n) - gsig_r(m,n)) * (1 + pow2(r(m,n))) 
    /*8*/  * r(m,n)))) / pow2(r_plus(m,n))) + (gA(m,n) * (Power(r_minus(m,n),4) * gAlp_r(m,n) 
    /*6*/  * gA_r(m,n) + exp(4 * gconf(m,n)) * gAsig_r(m,n) * gBet(m,n) * pow2(r(m,n)) 
    /*6*/  * (1 + pow2(r(m,n))) * pow3(gA(m,n)) + gA(m,n) * (-1 + pow2(r(m,n))) * (4 
    /*7*/  * gconf_r(m,n) * gAlp_r(m,n) * pow3(r_minus(m,n)) + gderAlpr_r(m,n) * (pow3(r(m,n))
    /*8*/  + r(m,n))))) / pow2(r_plus(m,n)) - (exp(4 * gconf(m,n)) * gAsig(m,n) 
    /*5*/  * Power(gA(m,n),4) * r(m,n) * (2 * gBet(m,n) * (-1 + pow2(r(m,n))) 
    /*6*/  - gAlp(m,n) * r(m,n) * gtrK(m,n))) / pow2(r_minus(m,n)))) / Power(gA(m,n),4))) 
    /*0*/  / (2. * pow2(r(m,n)));
}
Real fAsig_t( Int m, Int n )
{
    return ((1 - pow2(r(m,n))) * (-3 * k_f * fJA1(m,n) * (-1 + pow2(r(m,n))) - (2 
    /*3*/  * exp(-4 * fconf(m,n)) * fB_r(m,n) * pow2(r_minus(m,n)) * (fAlp_r(m,n) * fA(m,n) 
    /*4*/  * pow3(r_minus(m,n)) + 2 * fAlp(m,n) * ((fconf_r(m,n) * fA(m,n) - 2 * fA_r(m,n)) 
    /*5*/  * pow3(r_minus(m,n)) + 2 * fsig(m,n) * fA(m,n) * (1 + pow2(r(m,n))) * r(m,n)))) 
    /*2*/  / (fB(m,n) * pow2(r_plus(m,n)) * pow3(fA(m,n))) + (exp(-4 * fconf(m,n)) * fAlp(m,n)
    /*3*/  * pow2(fB(m,n)) * r(m,n) * (-2 * fL(m,n) * fsig(m,n) * pow2(fA(m,n)) - (2 
    /*5*/  * pow2(fsig(m,n)) * r(m,n)) / r_minus(m,n) + (fsig_r(m,n) * (-1 + pow2(r(m,n))) 
    /*5*/  * (-4 + 4 * pow2(r(m,n)) + fL(m,n) * pow2(fA(m,n)) * r(m,n))) / r_plus(m,n) 
    /*4*/  - (pow2(r_minus(m,n)) * r(m,n) * (2 * fsig_r(m,n) * (3 + pow2(r(m,n))) * r(m,n) 
    /*6*/  + fsig_rr(m,n) * (-1 + Power(r(m,n),4)))) / pow3(r_plus(m,n)))) / Power(fA(m,n),4)
    /*2*/  + (2 * exp(-4 * fconf(m,n)) * (1 - pow2(r(m,n))) * (-((fAlp(m,n) * (-1 
    /*7*/  + pow2(r(m,n))) * (-2 * fconf_r(m,n) * fA(m,n) * fA_r(m,n) * pow3(r_minus(m,n)) - 3
    /*7*/  * pow2(fA_r(m,n)) * pow3(r_minus(m,n)) + fLr_r(m,n) * Power(fA(m,n),4) * (1 
    /*8*/  + pow2(r(m,n))) * r(m,n) + pow2(fA(m,n)) * (-4 * pow2(fconf_r(m,n)) 
    /*8*/  * pow3(r_minus(m,n)) - (2 * fderconfr_r(m,n) - fsig_r(m,n)) * (1 + pow2(r(m,n))) 
    /*8*/  * r(m,n)))) / pow2(r_plus(m,n))) + (fA(m,n) * (Power(r_minus(m,n),4) * fAlp_r(m,n) 
    /*6*/  * fA_r(m,n) + exp(4 * fconf(m,n)) * fAsig_r(m,n) * fBet(m,n) * pow2(r(m,n)) 
    /*6*/  * (1 + pow2(r(m,n))) * pow3(fA(m,n)) + fA(m,n) * (-1 + pow2(r(m,n))) * (4 
    /*7*/  * fconf_r(m,n) * fAlp_r(m,n) * pow3(r_minus(m,n)) + fderAlpr_r(m,n) * (pow3(r(m,n))
    /*8*/  + r(m,n))))) / pow2(r_plus(m,n)) - (exp(4 * fconf(m,n)) * fAsig(m,n) 
    /*5*/  * Power(fA(m,n),4) * r(m,n) * (2 * fBet(m,n) * (-1 + pow2(r(m,n))) 
    /*6*/  - fAlp(m,n) * r(m,n) * ftrK(m,n))) / pow2(r_minus(m,n)))) / Power(fA(m,n),4))) 
    /*0*/  / (2. * pow2(r(m,n)));
}
