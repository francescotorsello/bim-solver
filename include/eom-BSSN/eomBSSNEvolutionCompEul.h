/** @file  eomBSSNEvolutionCompEul.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN Eulerian evolution equations.
 *  @version 2019-05-06T15:30:57
 *  @image html BSSNevolutionCompEul.png
 */

Real BimetricEvolve::eq_gconf_t( Int m, Int n )
{
    return gconf_r(m,n) * gBet(m,n) * pow2(-1 + r(m,n)) - (gAlp(m,n) * gtrK(m,n)) / (6
    /*1*/  + TINY_Real);
}
Real BimetricEvolve::eq_fconf_t( Int m, Int n )
{
    return fconf_r(m,n) * fBet(m,n) * pow2(-1 + r(m,n)) - (fAlp(m,n) * ftrK(m,n)) / (6
    /*1*/  + TINY_Real);
}
Real BimetricEvolve::eq_gtrK_t( Int m, Int n )
{
    return k_g * gJK(m,n) + gAlp(m,n) * ((3 * pow2(gA1(m,n))) / 2. + pow2(gtrK(m,n)) 
    /*1*/  / 3.) + ((gAlp_r(m,n) * gA_r(m,n) * Power(-1 + r(m,n),4)) / (TINY_Real 
    /*2*/  + Power(gA(m,n),3)) - (Power(-1 + r(m,n),4) * (2 * gAlp_r(m,n) * gB_r(m,n) 
    /*3*/  * (1 + r(m,n)) + gB(m,n) * (gAlp_rr(m,n) * (1 + r(m,n)) + 2 * gAlp_r(m,n) 
    /*4*/  * (1 + gconf_r(m,n) * (1 + r(m,n)))))) / (TINY_Real + Power(gA(m,n),2) 
    /*2*/  * gB(m,n) * (1 + r(m,n)))) / (TINY_Real + exp(4 * gconf(m,n))) + gBet(m,n) 
    /*0*/  * pow2(-1 + r(m,n)) * gtrK_r(m,n);
}
Real BimetricEvolve::eq_ftrK_t( Int m, Int n )
{
    return k_f * fJK(m,n) + fAlp(m,n) * ((3 * pow2(fA1(m,n))) / 2. + pow2(ftrK(m,n)) 
    /*1*/  / 3.) + ((fAlp_r(m,n) * fA_r(m,n) * Power(-1 + r(m,n),4)) / (TINY_Real 
    /*2*/  + Power(fA(m,n),3)) - (Power(-1 + r(m,n),4) * (2 * fAlp_r(m,n) * fB_r(m,n) 
    /*3*/  * (1 + r(m,n)) + fB(m,n) * (fAlp_rr(m,n) * (1 + r(m,n)) + 2 * fAlp_r(m,n) 
    /*4*/  * (1 + fconf_r(m,n) * (1 + r(m,n)))))) / (TINY_Real + Power(fA(m,n),2) 
    /*2*/  * fB(m,n) * (1 + r(m,n)))) / (TINY_Real + exp(4 * fconf(m,n))) + fBet(m,n) 
    /*0*/  * pow2(-1 + r(m,n)) * ftrK_r(m,n);
}
Real BimetricEvolve::eq_gA_t( Int m, Int n )
{
    return -(gA1(m,n) * gAlp(m,n) * gA(m,n)) + gBet_r(m,n) * gA(m,n) * pow2(-1 
    /*1*/  + r(m,n)) + gBet(m,n) * gA_r(m,n) * pow2(-1 + r(m,n));
}
Real BimetricEvolve::eq_gB_t( Int m, Int n )
{
    return (gA1(m,n) * gAlp(m,n) * gB(m,n)) / (2 + TINY_Real) + (gBet(m,n) * (-2 
    /*2*/  * gB(m,n) + gB_r(m,n) * (-1 + pow2(r(m,n)))) * (-1 + r(m,n))) / (1 
    /*1*/  + TINY_Real + r(m,n));
}
Real BimetricEvolve::eq_fA_t( Int m, Int n )
{
    return -(fA1(m,n) * fAlp(m,n) * fA(m,n)) + fBet_r(m,n) * fA(m,n) * pow2(-1 
    /*1*/  + r(m,n)) + fBet(m,n) * fA_r(m,n) * pow2(-1 + r(m,n));
}
Real BimetricEvolve::eq_fB_t( Int m, Int n )
{
    return (fA1(m,n) * fAlp(m,n) * fB(m,n)) / (2 + TINY_Real) + (fBet(m,n) * (-2 
    /*2*/  * fB(m,n) + fB_r(m,n) * (-1 + pow2(r(m,n)))) * (-1 + r(m,n))) / (1 
    /*1*/  + TINY_Real + r(m,n));
}
Real BimetricEvolve::eq_gA1_t( Int m, Int n )
{
    return gRicci(m,n) + k_g * gJA1(m,n) + gA1_r(m,n) * gBet(m,n) * pow2(-1 + r(m,n))
    /*0*/  + ((2 * gAlp_r(m,n) * gA_r(m,n) * Power(-1 + r(m,n),4)) / (TINY_Real + 3 
    /*2*/  * Power(gA(m,n),3)) + (2 * pow3(-1 + r(m,n)) * (gAlp_r(m,n) * gB_r(m,n) * (-1
    /*4*/  + pow2(r(m,n))) + gB(m,n) * (-(gAlp_rr(m,n) * (-1 + pow2(r(m,n)))) 
    /*4*/  + gAlp_r(m,n) * (4 * gconf_r(m,n) * (-1 + pow2(r(m,n))) - 2 * (2 
    /*6*/  + r(m,n)))))) / (TINY_Real + 3 * Power(gA(m,n),2) * gB(m,n) * (1 + r(m,n))))
    /*0*/  / (TINY_Real + exp(4 * gconf(m,n))) + gAlp(m,n) * (exp(-4 * gconf(m,n)) 
    /*1*/  * ((4 * gconf_r(m,n) * gA_r(m,n) * Power(-1 + r(m,n),4)) / (3. 
    /*3*/  * pow3(gA(m,n))) + (4 * pow3(-1 + r(m,n)) * (gconf_r(m,n) * gB_r(m,n) * (-1 
    /*5*/  + pow2(r(m,n))) + gB(m,n) * (-(gconf_rr(m,n) * (-1 + pow2(r(m,n)))) + 2 
    /*5*/  * pow2(gconf_r(m,n)) * (-1 + pow2(r(m,n))) - 2 * gconf_r(m,n) * (2 
    /*6*/  + r(m,n))))) / (3. * gB(m,n) * pow2(gA(m,n)) * (1 + r(m,n)))) + gA1(m,n) 
    /*1*/  * gtrK(m,n));
}
Real BimetricEvolve::eq_fA1_t( Int m, Int n )
{
    return fRicci(m,n) + k_f * fJA1(m,n) + fA1_r(m,n) * fBet(m,n) * pow2(-1 + r(m,n))
    /*0*/  + ((2 * fAlp_r(m,n) * fA_r(m,n) * Power(-1 + r(m,n),4)) / (TINY_Real + 3 
    /*2*/  * Power(fA(m,n),3)) + (2 * pow3(-1 + r(m,n)) * (fAlp_r(m,n) * fB_r(m,n) * (-1
    /*4*/  + pow2(r(m,n))) + fB(m,n) * (-(fAlp_rr(m,n) * (-1 + pow2(r(m,n)))) 
    /*4*/  + fAlp_r(m,n) * (4 * fconf_r(m,n) * (-1 + pow2(r(m,n))) - 2 * (2 
    /*6*/  + r(m,n)))))) / (TINY_Real + 3 * Power(fA(m,n),2) * fB(m,n) * (1 + r(m,n))))
    /*0*/  / (TINY_Real + exp(4 * fconf(m,n))) + fAlp(m,n) * (exp(-4 * fconf(m,n)) 
    /*1*/  * ((4 * fconf_r(m,n) * fA_r(m,n) * Power(-1 + r(m,n),4)) / (3. 
    /*3*/  * pow3(fA(m,n))) + (4 * pow3(-1 + r(m,n)) * (fconf_r(m,n) * fB_r(m,n) * (-1 
    /*5*/  + pow2(r(m,n))) + fB(m,n) * (-(fconf_rr(m,n) * (-1 + pow2(r(m,n)))) + 2 
    /*5*/  * pow2(fconf_r(m,n)) * (-1 + pow2(r(m,n))) - 2 * fconf_r(m,n) * (2 
    /*6*/  + r(m,n))))) / (3. * fB(m,n) * pow2(fA(m,n)) * (1 + r(m,n)))) + fA1(m,n) 
    /*1*/  * ftrK(m,n));
}
Real BimetricEvolve::eq_gL_t( Int m, Int n )
{
    return k_g * gJL(m,n) - (2 * gA1(m,n) * gAlp_r(m,n) * pow2(-1 + r(m,n))) 
    /*0*/  / (TINY_Real + Power(gA(m,n),2)) + (pow3(-1 + r(m,n)) * (2 * gBet_r(m,n) 
    /*2*/  + gBet_rr(m,n) * (-1 + r(m,n)))) / (TINY_Real + Power(gA(m,n),2)) + (pow2(-1
    /*2*/  + r(m,n)) * (gBet(m,n) * (-8 + gL_r(m,n) * pow2(gB(m,n)) * pow2(1 
    /*4*/  + r(m,n))) - gBet_r(m,n) * (1 + r(m,n)) * (4 * (-1 + r(m,n)) + gL(m,n) 
    /*3*/  * pow2(gB(m,n)) * (1 + r(m,n))))) / (TINY_Real + Power(gB(m,n),2) * Power(1 
    /*2*/  + r(m,n),2)) + gAlp(m,n) * (gA1(m,n) * ((2 * gA_r(m,n) * pow2(-1 + r(m,n))) 
    /*2*/  / pow3(gA(m,n)) + (4 * (-1 + r(m,n))) / (pow2(gB(m,n)) * (1 + r(m,n))) + (2 
    /*3*/  * (gB(m,n) * (-2 + 6 * gconf_r(m,n) * (-1 + pow2(r(m,n)))) + gB_r(m,n) * (-1
    /*5*/  + pow2(r(m,n)))) * (-1 + r(m,n))) / (gB(m,n) * pow2(gA(m,n)) * (1 
    /*4*/  + r(m,n)))) - (4 * pow2(-1 + r(m,n)) * gtrK_r(m,n)) / (3. * pow2(gA(m,n))));
}
Real BimetricEvolve::eq_fL_t( Int m, Int n )
{
    return k_f * fJL(m,n) - (2 * fA1(m,n) * fAlp_r(m,n) * pow2(-1 + r(m,n))) 
    /*0*/  / (TINY_Real + Power(fA(m,n),2)) + (pow3(-1 + r(m,n)) * (2 * fBet_r(m,n) 
    /*2*/  + fBet_rr(m,n) * (-1 + r(m,n)))) / (TINY_Real + Power(fA(m,n),2)) + (pow2(-1
    /*2*/  + r(m,n)) * (fBet(m,n) * (-8 + fL_r(m,n) * pow2(fB(m,n)) * pow2(1 
    /*4*/  + r(m,n))) - fBet_r(m,n) * (1 + r(m,n)) * (4 * (-1 + r(m,n)) + fL(m,n) 
    /*3*/  * pow2(fB(m,n)) * (1 + r(m,n))))) / (TINY_Real + Power(fB(m,n),2) * Power(1 
    /*2*/  + r(m,n),2)) + fAlp(m,n) * (fA1(m,n) * ((2 * fA_r(m,n) * pow2(-1 + r(m,n))) 
    /*2*/  / pow3(fA(m,n)) + (4 * (-1 + r(m,n))) / (pow2(fB(m,n)) * (1 + r(m,n))) + (2 
    /*3*/  * (fB(m,n) * (-2 + 6 * fconf_r(m,n) * (-1 + pow2(r(m,n)))) + fB_r(m,n) * (-1
    /*5*/  + pow2(r(m,n)))) * (-1 + r(m,n))) / (fB(m,n) * pow2(fA(m,n)) * (1 
    /*4*/  + r(m,n)))) - (4 * pow2(-1 + r(m,n)) * ftrK_r(m,n)) / (3. * pow2(fA(m,n))));
}
