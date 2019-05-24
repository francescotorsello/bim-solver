/** @file  eomBSSNEvolutionCompEul.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN Eulerian evolution equations.
 *  @version 2019-05-24T15:31:53
 *  @image html BSSNevolutionCompEul.png
 */

Real gconf_t( Int m, Int n )
{
    return (gconf_r(m,n) * gBet(m,n) * pow2(r_minus(m,n))) / r_plus(m,n)
    /*0*/  - (gAlp(m,n) * gtrK(m,n)) / 6.;
}
Real fconf_t( Int m, Int n )
{
    return (fconf_r(m,n) * fBet(m,n) * pow2(r_minus(m,n))) / r_plus(m,n)
    /*0*/  - (fAlp(m,n) * ftrK(m,n)) / 6.;
}
Real gtrK_t( Int m, Int n )
{
    return -((Power(r_minus(m,n),4) * gAlp_rr(m,n)) / (TINY_Real 
    /*2*/  + Power(r_plus(m,n),2) * exp(4 * gconf(m,n)) 
    /*2*/  * Power(gA(m,n),2))) + k_g * gJK(m,n) + gAlp(m,n) 
    /*0*/  * (3 * pow2(gA1(m,n)) - (4 * gA1(m,n) * gAsig(m,n) 
    /*2*/  * pow2(r(m,n))) / (TINY_Real + Power(r_minus(m,n),2)) 
    /*1*/  + pow2(gtrK(m,n)) / (3 + TINY_Real) + (2 
    /*2*/  * pow2(gAsig(m,n)) * Power(r(m,n),4)) / (TINY_Real 
    /*2*/  + Power(r_minus(m,n),4))) + gAlp_r(m,n) * ((-2 
    /*2*/  * Power(r_minus(m,n),4) * gconf_r(m,n)) / (TINY_Real 
    /*2*/  + Power(r_plus(m,n),2) * exp(4 * gconf(m,n)) 
    /*2*/  * Power(gA(m,n),2)) + (Power(r_minus(m,n),4) * gA_r(m,n))
    /*1*/  / (TINY_Real + Power(r_plus(m,n),2) * exp(4 
    /*3*/  * gconf(m,n)) * Power(gA(m,n),3)) - (2 
    /*2*/  * Power(r_minus(m,n),4) * gB_r(m,n)) / (TINY_Real 
    /*2*/  + Power(r_plus(m,n),2) * exp(4 * gconf(m,n)) 
    /*2*/  * Power(gA(m,n),2) * gB(m,n)) - (2 
    /*2*/  * Power(r_minus(m,n),4)) / (TINY_Real + Power(r_plus(m,n),3) 
    /*2*/  * exp(4 * gconf(m,n)) * Power(gA(m,n),2) * r(m,n)))
    /*0*/  + (gBet(m,n) * pow2(r_minus(m,n)) * gtrK_r(m,n)) / (1 
    /*1*/  + TINY_Real + Power(r(m,n),2));
}
Real ftrK_t( Int m, Int n )
{
    return -((Power(r_minus(m,n),4) * fAlp_rr(m,n)) / (TINY_Real 
    /*2*/  + Power(r_plus(m,n),2) * exp(4 * fconf(m,n)) 
    /*2*/  * Power(fA(m,n),2))) + k_f * fJK(m,n) + fAlp(m,n) 
    /*0*/  * (3 * pow2(fA1(m,n)) - (4 * fA1(m,n) * fAsig(m,n) 
    /*2*/  * pow2(r(m,n))) / (TINY_Real + Power(r_minus(m,n),2)) 
    /*1*/  + pow2(ftrK(m,n)) / (3 + TINY_Real) + (2 
    /*2*/  * pow2(fAsig(m,n)) * Power(r(m,n),4)) / (TINY_Real 
    /*2*/  + Power(r_minus(m,n),4))) + fAlp_r(m,n) * ((-2 
    /*2*/  * Power(r_minus(m,n),4) * fconf_r(m,n)) / (TINY_Real 
    /*2*/  + Power(r_plus(m,n),2) * exp(4 * fconf(m,n)) 
    /*2*/  * Power(fA(m,n),2)) + (Power(r_minus(m,n),4) * fA_r(m,n))
    /*1*/  / (TINY_Real + Power(r_plus(m,n),2) * exp(4 
    /*3*/  * fconf(m,n)) * Power(fA(m,n),3)) - (2 
    /*2*/  * Power(r_minus(m,n),4) * fB_r(m,n)) / (TINY_Real 
    /*2*/  + Power(r_plus(m,n),2) * exp(4 * fconf(m,n)) 
    /*2*/  * Power(fA(m,n),2) * fB(m,n)) - (2 
    /*2*/  * Power(r_minus(m,n),4)) / (TINY_Real + Power(r_plus(m,n),3) 
    /*2*/  * exp(4 * fconf(m,n)) * Power(fA(m,n),2) * r(m,n)))
    /*0*/  + (fBet(m,n) * pow2(r_minus(m,n)) * ftrK_r(m,n)) / (1 
    /*1*/  + TINY_Real + Power(r(m,n),2));
}
Real gA_t( Int m, Int n )
{
    return -(gA1(m,n) * gAlp(m,n) * gA(m,n)) + (gBet_r(m,n) 
    /*1*/  * gA(m,n) * pow2(r_minus(m,n))) / r_plus(m,n) + (gBet(m,n) 
    /*1*/  * gA_r(m,n) * pow2(r_minus(m,n))) / r_plus(m,n);
}
Real gB_t( Int m, Int n )
{
    return gBetr(m,n) * gB(m,n) + gAlp(m,n) * (-(gA1(m,n) 
    /*2*/  * gB(m,n)) + (gAsig(m,n) * gB(m,n) * pow2(r(m,n))) 
    /*1*/  / (TINY_Real + Power(r_minus(m,n),2))) + (gBet(m,n) 
    /*1*/  * gB_r(m,n) * pow2(r_minus(m,n))) / (1 + TINY_Real 
    /*1*/  + Power(r(m,n),2));
}
Real fA_t( Int m, Int n )
{
    return -(fA1(m,n) * fAlp(m,n) * fA(m,n)) + (fBet_r(m,n) 
    /*1*/  * fA(m,n) * pow2(r_minus(m,n))) / r_plus(m,n) + (fBet(m,n) 
    /*1*/  * fA_r(m,n) * pow2(r_minus(m,n))) / r_plus(m,n);
}
Real fB_t( Int m, Int n )
{
    return fBetr(m,n) * fB(m,n) + fAlp(m,n) * (-(fA1(m,n) 
    /*2*/  * fB(m,n)) + (fAsig(m,n) * fB(m,n) * pow2(r(m,n))) 
    /*1*/  / (TINY_Real + Power(r_minus(m,n),2))) + (fBet(m,n) 
    /*1*/  * fB_r(m,n) * pow2(r_minus(m,n))) / (1 + TINY_Real 
    /*1*/  + Power(r(m,n),2));
}
Real gA1_t( Int m, Int n )
{
    return (3 * gRicci(m,n) + (2 * Power(r_minus(m,n),4) * (2 
    /*3*/  * gconf_r(m,n) * gAlp(m,n) + gAlp_r(m,n)) 
    /*2*/  * gA_r(m,n)) / (TINY_Real + Power(r_plus(m,n),2) * exp(4 
    /*3*/  * gconf(m,n)) * Power(gA(m,n),3)) + 3 * k_g 
    /*1*/  * gJA1(m,n) + (3 * gA1_r(m,n) * gBet(m,n) 
    /*2*/  * pow2(r_minus(m,n))) / (1 + TINY_Real + Power(r(m,n),2))
    /*1*/  - (2 * pow3(r_minus(m,n)) * (-(gAlp_r(m,n) * gB_r(m,n) 
    /*4*/  * r(m,n) * (-1 + Power(r(m,n),4))) + gB(m,n) 
    /*3*/  * (gAlp_rr(m,n) * r(m,n) * (-1 + Power(r(m,n),4)) 
    /*4*/  + gAlp_r(m,n) * (1 + 8 * pow2(r(m,n)) + 3 
    /*5*/  * Power(r(m,n),4) - 4 * gconf_r(m,n) * r(m,n) * (-1
    /*6*/  + Power(r(m,n),4)))) + 2 * gAlp(m,n) 
    /*3*/  * (-(gconf_r(m,n) * gB_r(m,n) * r(m,n) * (-1 
    /*6*/  + Power(r(m,n),4))) + gB(m,n) * (gconf_rr(m,n) 
    /*5*/  * r(m,n) * (-1 + Power(r(m,n),4)) - 2 
    /*5*/  * pow2(gconf_r(m,n)) * r(m,n) * (-1 
    /*6*/  + Power(r(m,n),4)) + gconf_r(m,n) * (1 + 8 
    /*6*/  * pow2(r(m,n)) + 3 * Power(r(m,n),4)))))) 
    /*1*/  / (TINY_Real + Power(r_plus(m,n),3) * exp(4 * gconf(m,n))
    /*2*/  * Power(gA(m,n),2) * gB(m,n) * r(m,n)) + 3 
    /*1*/  * gA1(m,n) * gAlp(m,n) * gtrK(m,n)) / (3 
    /*1*/  + TINY_Real);
}
Real fA1_t( Int m, Int n )
{
    return (3 * fRicci(m,n) + (2 * Power(r_minus(m,n),4) * (2 
    /*3*/  * fconf_r(m,n) * fAlp(m,n) + fAlp_r(m,n)) 
    /*2*/  * fA_r(m,n)) / (TINY_Real + Power(r_plus(m,n),2) * exp(4 
    /*3*/  * fconf(m,n)) * Power(fA(m,n),3)) + 3 * k_f 
    /*1*/  * fJA1(m,n) + (3 * fA1_r(m,n) * fBet(m,n) 
    /*2*/  * pow2(r_minus(m,n))) / (1 + TINY_Real + Power(r(m,n),2))
    /*1*/  - (2 * (2 * fconf_r(m,n) * fAlp(m,n) 
    /*3*/  + fAlp_r(m,n)) * pow3(r_minus(m,n))) / (TINY_Real + exp(4
    /*3*/  * fconf(m,n)) * Power(fA(m,n),2) * (r(m,n) 
    /*3*/  + Power(r(m,n),3))) + (2 * pow3(r_minus(m,n)) 
    /*2*/  * (fAlp_r(m,n) * fB_r(m,n) * (-1 + Power(r(m,n),4))
    /*3*/  + fB(m,n) * (-(fAlp_rr(m,n) * (-1 
    /*6*/  + Power(r(m,n),4))) + fAlp_r(m,n) * (-2 * (3 
    /*6*/  + pow2(r(m,n))) * r(m,n) + 4 * fconf_r(m,n) * (-1 
    /*6*/  + Power(r(m,n),4)))) + 2 * fAlp(m,n) * (fconf_r(m,n)
    /*4*/  * fB_r(m,n) * (-1 + Power(r(m,n),4)) + fB(m,n) 
    /*4*/  * (-2 * fconf_r(m,n) * (3 + pow2(r(m,n))) * r(m,n) 
    /*5*/  - fconf_rr(m,n) * (-1 + Power(r(m,n),4)) + 2 
    /*5*/  * pow2(fconf_r(m,n)) * (-1 + Power(r(m,n),4)))))) 
    /*1*/  / (TINY_Real + Power(r_plus(m,n),3) * exp(4 * fconf(m,n))
    /*2*/  * Power(fA(m,n),2) * fB(m,n)) + 3 * fA1(m,n) 
    /*1*/  * fAlp(m,n) * ftrK(m,n)) / (3 + TINY_Real);
}
Real gL_t( Int m, Int n )
{
    return (6 * k_g * gJL(m,n) - (6 * gBet_r(m,n) * gL(m,n) 
    /*2*/  * pow2(r_minus(m,n))) / (1 + TINY_Real + Power(r(m,n),2))
    /*1*/  + (6 * gBet(m,n) * gL_r(m,n) * pow2(r_minus(m,n))) / (1
    /*2*/  + TINY_Real + Power(r(m,n),2)) + (4 * gA1(m,n) 
    /*2*/  * (-3 * gAlp_r(m,n) * gA(m,n) * gB(m,n) + 2 
    /*3*/  * gAlp(m,n) * (gA_r(m,n) * gB(m,n) + gA(m,n) * (6 
    /*5*/  * gconf_r(m,n) * gB(m,n) + gB_r(m,n)))) 
    /*2*/  * pow2(r_minus(m,n))) / (TINY_Real + Power(gA(m,n),3) 
    /*2*/  * gB(m,n) * (1 + Power(r(m,n),2))) - (8 * gAlp(m,n)
    /*2*/  * (gA_r(m,n) * gB(m,n) + gA(m,n) * (6 
    /*4*/  * gconf_r(m,n) * gB(m,n) + gB_r(m,n))) * (gA1(m,n) 
    /*3*/  * pow2(r_minus(m,n)) - gAsig(m,n) * pow2(r(m,n)))) 
    /*1*/  / (TINY_Real + Power(gA(m,n),3) * gB(m,n) * (1 
    /*3*/  + Power(r(m,n),2))) + (12 * pow2(r_minus(m,n)) 
    /*2*/  * (gBet_r(m,n) * (-1 + pow2(r(m,n))) + gBetr_r(m,n)
    /*3*/  * r(m,n))) / (TINY_Real + Power(gB(m,n),2) 
    /*2*/  * (r(m,n) + Power(r(m,n),3))) + (2 * (1 
    /*3*/  - pow2(r(m,n))) * (6 * Power(r_minus(m,n),6) * gBet_r(m,n)
    /*3*/  * pow2(gA(m,n)) + 4 * gAsig(m,n) * gAlp(m,n) 
    /*3*/  * gsig(m,n) * pow2(gB(m,n)) * (1 + pow2(r(m,n))) 
    /*3*/  * Power(r(m,n),4))) / (TINY_Real + Power(r_minus(m,n),4) 
    /*2*/  * Power(gA(m,n),2) * Power(gB(m,n),2) * (r(m,n) 
    /*3*/  + Power(r(m,n),3))) - (2 * pow2(r_minus(m,n)) * (-6 
    /*3*/  * gBet_r(m,n) * r(m,n) * (-3 + 2 * pow2(r(m,n)) 
    /*4*/  + Power(r(m,n),4)) + (1 + pow2(r(m,n))) * (-3 
    /*4*/  * gBet_rr(m,n) * pow2(r_minus(m,n)) + 4 * gAlp(m,n) * (1 
    /*5*/  + pow2(r(m,n))) * gtrK_r(m,n)))) / (TINY_Real 
    /*2*/  + Power(r_plus(m,n),3) * Power(gA(m,n),2))) / (6 
    /*1*/  + TINY_Real);
}
Real fL_t( Int m, Int n )
{
    return (6 * k_f * fJL(m,n) - (6 * fBet_r(m,n) * fL(m,n) 
    /*2*/  * pow2(r_minus(m,n))) / (1 + TINY_Real + Power(r(m,n),2))
    /*1*/  + (6 * fBet(m,n) * fL_r(m,n) * pow2(r_minus(m,n))) / (1
    /*2*/  + TINY_Real + Power(r(m,n),2)) + (4 * fA1(m,n) 
    /*2*/  * (-3 * fAlp_r(m,n) * fA(m,n) * fB(m,n) + 2 
    /*3*/  * fAlp(m,n) * (fA_r(m,n) * fB(m,n) + fA(m,n) * (6 
    /*5*/  * fconf_r(m,n) * fB(m,n) + fB_r(m,n)))) 
    /*2*/  * pow2(r_minus(m,n))) / (TINY_Real + Power(fA(m,n),3) 
    /*2*/  * fB(m,n) * (1 + Power(r(m,n),2))) - (8 * fAlp(m,n)
    /*2*/  * (fA_r(m,n) * fB(m,n) + fA(m,n) * (6 
    /*4*/  * fconf_r(m,n) * fB(m,n) + fB_r(m,n))) * (fA1(m,n) 
    /*3*/  * pow2(r_minus(m,n)) - fAsig(m,n) * pow2(r(m,n)))) 
    /*1*/  / (TINY_Real + Power(fA(m,n),3) * fB(m,n) * (1 
    /*3*/  + Power(r(m,n),2))) + (12 * pow2(r_minus(m,n)) 
    /*2*/  * (fBet_r(m,n) * (-1 + pow2(r(m,n))) + fBetr_r(m,n)
    /*3*/  * r(m,n))) / (TINY_Real + Power(fB(m,n),2) 
    /*2*/  * (r(m,n) + Power(r(m,n),3))) + (2 * (1 
    /*3*/  - pow2(r(m,n))) * (6 * Power(r_minus(m,n),6) * fBet_r(m,n)
    /*3*/  * pow2(fA(m,n)) + 4 * fAsig(m,n) * fAlp(m,n) 
    /*3*/  * fsig(m,n) * pow2(fB(m,n)) * (1 + pow2(r(m,n))) 
    /*3*/  * Power(r(m,n),4))) / (TINY_Real + Power(r_minus(m,n),4) 
    /*2*/  * Power(fA(m,n),2) * Power(fB(m,n),2) * (r(m,n) 
    /*3*/  + Power(r(m,n),3))) - (2 * pow2(r_minus(m,n)) * (-6 
    /*3*/  * fBet_r(m,n) * r(m,n) * (-3 + 2 * pow2(r(m,n)) 
    /*4*/  + Power(r(m,n),4)) + (1 + pow2(r(m,n))) * (-3 
    /*4*/  * fBet_rr(m,n) * pow2(r_minus(m,n)) + 4 * fAlp(m,n) * (1 
    /*5*/  + pow2(r(m,n))) * ftrK_r(m,n)))) / (TINY_Real 
    /*2*/  + Power(r_plus(m,n),3) * Power(fA(m,n),2))) / (6 
    /*1*/  + TINY_Real);
}
Real gsig_t( Int m, Int n )
{
    return 2 * gBetr(m,n) * gsig(m,n) + (gBet(m,n) 
    /*1*/  * gsig_r(m,n) * pow2(r_minus(m,n))) / (1 + TINY_Real 
    /*1*/  + Power(r(m,n),2)) + (pow2(gA(m,n)) * (2 
    /*2*/  * gAsig(m,n) * gAlp(m,n) + (2 * gBetr_r(m,n) 
    /*3*/  * pow3(r_minus(m,n))) / (r(m,n) + Power(r(m,n),3)))) 
    /*0*/  / (TINY_Real + Power(gB(m,n),2));
}
Real fsig_t( Int m, Int n )
{
    return 2 * fBetr(m,n) * fsig(m,n) + (2 * fAsig(m,n) 
    /*1*/  * fAlp(m,n) * pow2(fA(m,n))) / (TINY_Real 
    /*1*/  + Power(fB(m,n),2)) + (fBet(m,n) * fsig_r(m,n) 
    /*1*/  * pow2(r_minus(m,n))) / (1 + TINY_Real + Power(r(m,n),2))
    /*0*/  + (2 * fBetr_r(m,n) * pow2(fA(m,n)) 
    /*1*/  * pow3(r_minus(m,n))) / ((TINY_Real + Power(fB(m,n),2)) 
    /*1*/  * (r(m,n) + Power(r(m,n),3)));
}
Real gAsig_t( Int m, Int n )
{
    return ((1 - pow2(r(m,n))) * (-3 * k_g * gJA1(m,n) * (-1
    /*3*/  + pow2(r(m,n))) - (2 * gB_r(m,n) * pow2(r_minus(m,n)) 
    /*3*/  * (gAlp_r(m,n) * gA(m,n) * pow3(r_minus(m,n)) + 2 
    /*4*/  * gAlp(m,n) * ((gconf_r(m,n) * gA(m,n) - 2 
    /*6*/  * gA_r(m,n)) * pow3(r_minus(m,n)) + 2 * gsig(m,n) 
    /*5*/  * gA(m,n) * (1 + pow2(r(m,n))) * r(m,n)))) 
    /*2*/  / (TINY_Real + Power(r_plus(m,n),2) * exp(4 * gconf(m,n))
    /*3*/  * Power(gA(m,n),3) * gB(m,n)) + (gAlp(m,n) 
    /*3*/  * pow2(gB(m,n)) * r(m,n) * (-2 * gL(m,n) * gsig(m,n)
    /*4*/  * pow2(gA(m,n)) - (2 * pow2(gsig(m,n)) * r(m,n)) 
    /*4*/  / (TINY_Real + r_minus(m,n)) + (gsig_r(m,n) * (-1 
    /*6*/  + pow2(r(m,n))) * (-4 + 4 * pow2(r(m,n)) + gL(m,n) 
    /*6*/  * pow2(gA(m,n)) * r(m,n))) / r_plus(m,n) - (pow2(r_minus(m,n)) 
    /*5*/  * r(m,n) * (2 * gsig_r(m,n) * (3 + pow2(r(m,n))) 
    /*6*/  * r(m,n) + gsig_rr(m,n) * (-1 + Power(r(m,n),4)))) 
    /*4*/  / pow3(r_plus(m,n)))) / (TINY_Real + exp(4 * gconf(m,n)) 
    /*3*/  * Power(gA(m,n),4)) + (2 * (1 - pow2(r(m,n))) 
    /*3*/  * (-((gAlp(m,n) * (-1 + pow2(r(m,n))) * (-2 
    /*7*/  * gconf_r(m,n) * gA(m,n) * gA_r(m,n) * pow3(r_minus(m,n))
    /*7*/  - 3 * pow2(gA_r(m,n)) * pow3(r_minus(m,n)) + gLr_r(m,n)
    /*7*/  * Power(gA(m,n),4) * (1 + pow2(r(m,n))) * r(m,n) 
    /*7*/  + pow2(gA(m,n)) * (-4 * pow2(gconf_r(m,n)) 
    /*8*/  * pow3(r_minus(m,n)) - (2 * gderconfr_r(m,n) 
    /*9*/  - gsig_r(m,n)) * (1 + pow2(r(m,n))) * r(m,n)))) 
    /*5*/  / (TINY_Real + Power(r_plus(m,n),2))) + (gA(m,n) 
    /*5*/  * (Power(r_minus(m,n),4) * gAlp_r(m,n) * gA_r(m,n) + exp(4
    /*7*/  * gconf(m,n)) * gAsig_r(m,n) * gBet(m,n) 
    /*6*/  * pow2(r(m,n)) * (1 + pow2(r(m,n))) * pow3(gA(m,n))
    /*6*/  + gA(m,n) * (-1 + pow2(r(m,n))) * (4 
    /*7*/  * gconf_r(m,n) * gAlp_r(m,n) * pow3(r_minus(m,n)) 
    /*7*/  + gderAlpr_r(m,n) * (pow3(r(m,n)) + r(m,n))))) 
    /*4*/  / (TINY_Real + Power(r_plus(m,n),2)) - (exp(4 
    /*6*/  * gconf(m,n)) * gAsig(m,n) * Power(gA(m,n),4) 
    /*5*/  * r(m,n) * (2 * gBet(m,n) * (-1 + pow2(r(m,n))) 
    /*6*/  - gAlp(m,n) * r(m,n) * gtrK(m,n))) / (TINY_Real 
    /*5*/  + Power(r_minus(m,n),2)))) / (TINY_Real + exp(4 
    /*4*/  * gconf(m,n)) * Power(gA(m,n),4)))) / (2. 
    /*1*/  * pow2(r(m,n)));
}
Real fAsig_t( Int m, Int n )
{
    return ((1 - pow2(r(m,n))) * (-3 * k_f * fJA1(m,n) * (-1
    /*3*/  + pow2(r(m,n))) - (2 * fB_r(m,n) * pow2(r_minus(m,n)) 
    /*3*/  * (fAlp_r(m,n) * fA(m,n) * pow3(r_minus(m,n)) + 2 
    /*4*/  * fAlp(m,n) * ((fconf_r(m,n) * fA(m,n) - 2 
    /*6*/  * fA_r(m,n)) * pow3(r_minus(m,n)) + 2 * fsig(m,n) 
    /*5*/  * fA(m,n) * (1 + pow2(r(m,n))) * r(m,n)))) 
    /*2*/  / (TINY_Real + Power(r_plus(m,n),2) * exp(4 * fconf(m,n))
    /*3*/  * Power(fA(m,n),3) * fB(m,n)) + (fAlp(m,n) 
    /*3*/  * pow2(fB(m,n)) * r(m,n) * (-2 * fL(m,n) * fsig(m,n)
    /*4*/  * pow2(fA(m,n)) - (2 * pow2(fsig(m,n)) * r(m,n)) 
    /*4*/  / (TINY_Real + r_minus(m,n)) + (fsig_r(m,n) * (-1 
    /*6*/  + pow2(r(m,n))) * (-4 + 4 * pow2(r(m,n)) + fL(m,n) 
    /*6*/  * pow2(fA(m,n)) * r(m,n))) / r_plus(m,n) - (pow2(r_minus(m,n)) 
    /*5*/  * r(m,n) * (2 * fsig_r(m,n) * (3 + pow2(r(m,n))) 
    /*6*/  * r(m,n) + fsig_rr(m,n) * (-1 + Power(r(m,n),4)))) 
    /*4*/  / pow3(r_plus(m,n)))) / (TINY_Real + exp(4 * fconf(m,n)) 
    /*3*/  * Power(fA(m,n),4)) + (2 * (1 - pow2(r(m,n))) 
    /*3*/  * (-((fAlp(m,n) * (-1 + pow2(r(m,n))) * (-2 
    /*7*/  * fconf_r(m,n) * fA(m,n) * fA_r(m,n) * pow3(r_minus(m,n))
    /*7*/  - 3 * pow2(fA_r(m,n)) * pow3(r_minus(m,n)) + fLr_r(m,n)
    /*7*/  * Power(fA(m,n),4) * (1 + pow2(r(m,n))) * r(m,n) 
    /*7*/  + pow2(fA(m,n)) * (-4 * pow2(fconf_r(m,n)) 
    /*8*/  * pow3(r_minus(m,n)) - (2 * fderconfr_r(m,n) 
    /*9*/  - fsig_r(m,n)) * (1 + pow2(r(m,n))) * r(m,n)))) 
    /*5*/  / (TINY_Real + Power(r_plus(m,n),2))) + (fA(m,n) 
    /*5*/  * (Power(r_minus(m,n),4) * fAlp_r(m,n) * fA_r(m,n) + exp(4
    /*7*/  * fconf(m,n)) * fAsig_r(m,n) * fBet(m,n) 
    /*6*/  * pow2(r(m,n)) * (1 + pow2(r(m,n))) * pow3(fA(m,n))
    /*6*/  + fA(m,n) * (-1 + pow2(r(m,n))) * (4 
    /*7*/  * fconf_r(m,n) * fAlp_r(m,n) * pow3(r_minus(m,n)) 
    /*7*/  + fderAlpr_r(m,n) * (pow3(r(m,n)) + r(m,n))))) 
    /*4*/  / (TINY_Real + Power(r_plus(m,n),2)) - (exp(4 
    /*6*/  * fconf(m,n)) * fAsig(m,n) * Power(fA(m,n),4) 
    /*5*/  * r(m,n) * (2 * fBet(m,n) * (-1 + pow2(r(m,n))) 
    /*6*/  - fAlp(m,n) * r(m,n) * ftrK(m,n))) / (TINY_Real 
    /*5*/  + Power(r_minus(m,n),2)))) / (TINY_Real + exp(4 
    /*4*/  * fconf(m,n)) * Power(fA(m,n),4)))) / (2. 
    /*1*/  * pow2(r(m,n)));
}
