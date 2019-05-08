/** @file  eomBSSNEvolutionCompEul.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN Eulerian evolution equations.
 *  @version 2019-05-08T16:05:33
 *  @image html BSSNevolutionCompEul.png
 */

Real gconf_t( Int m, Int n )
{
    return gconf_r(m,n) * gBet(m,n) * pow2(-1 + r(m,n)) - (gAlp(m,n) * gtrK(m,n)) / 6.;
}
Real fconf_t( Int m, Int n )
{
    return fconf_r(m,n) * fBet(m,n) * pow2(-1 + r(m,n)) - (fAlp(m,n) * ftrK(m,n)) / 6.;
}
Real gtrK_t( Int m, Int n )
{
    return k_g * gJK(m,n) + gAlp(m,n) * ((3 * pow2(gA1(m,n))) / 2. + pow2(gtrK(m,n)) 
    /*1*/  / 3.) + exp(-4 * gconf(m,n)) * ((gAlp_r(m,n) * gA_r(m,n) * Power(-1 
    /*3*/  + r(m,n),4)) / pow3(gA(m,n)) - (Power(-1 + r(m,n),4) * (2 * gAlp_r(m,n) 
    /*3*/  * gB_r(m,n) * (1 + TINY_Real + r(m,n)) + gB(m,n) * (gAlp_rr(m,n) * (1 
    /*5*/  + TINY_Real + r(m,n)) + 2 * gAlp_r(m,n) * (1 + gconf_r(m,n) * (1 + TINY_Real
    /*6*/  + r(m,n)))))) / (gB(m,n) * pow2(gA(m,n)) * (1 + TINY_Real + r(m,n)))) 
    /*0*/  + gBet(m,n) * pow2(-1 + r(m,n)) * gtrK_r(m,n);
}
Real ftrK_t( Int m, Int n )
{
    return k_f * fJK(m,n) + fAlp(m,n) * ((3 * pow2(fA1(m,n))) / 2. + pow2(ftrK(m,n)) 
    /*1*/  / 3.) + exp(-4 * fconf(m,n)) * ((fAlp_r(m,n) * fA_r(m,n) * Power(-1 
    /*3*/  + r(m,n),4)) / pow3(fA(m,n)) - (Power(-1 + r(m,n),4) * (2 * fAlp_r(m,n) 
    /*3*/  * fB_r(m,n) * (1 + TINY_Real + r(m,n)) + fB(m,n) * (fAlp_rr(m,n) * (1 
    /*5*/  + TINY_Real + r(m,n)) + 2 * fAlp_r(m,n) * (1 + fconf_r(m,n) * (1 + TINY_Real
    /*6*/  + r(m,n)))))) / (fB(m,n) * pow2(fA(m,n)) * (1 + TINY_Real + r(m,n)))) 
    /*0*/  + fBet(m,n) * pow2(-1 + r(m,n)) * ftrK_r(m,n);
}
Real gA_t( Int m, Int n )
{
    return -(gA1(m,n) * gAlp(m,n) * gA(m,n)) + gBet_r(m,n) * gA(m,n) * pow2(-1 
    /*1*/  + r(m,n)) + gBet(m,n) * gA_r(m,n) * pow2(-1 + r(m,n));
}
Real gB_t( Int m, Int n )
{
    return gBet(m,n) * gB_r(m,n) * pow2(-1 + r(m,n)) + gB(m,n) * ((gA1(m,n) 
    /*2*/  * gAlp(m,n)) / 2. + (gBet(m,n) * (2 - 2 * r(m,n))) / (1 + TINY_Real 
    /*2*/  + r(m,n)));
}
Real fA_t( Int m, Int n )
{
    return -(fA1(m,n) * fAlp(m,n) * fA(m,n)) + fBet_r(m,n) * fA(m,n) * pow2(-1 
    /*1*/  + r(m,n)) + fBet(m,n) * fA_r(m,n) * pow2(-1 + r(m,n));
}
Real fB_t( Int m, Int n )
{
    return fBet(m,n) * fB_r(m,n) * pow2(-1 + r(m,n)) + fB(m,n) * ((fA1(m,n) 
    /*2*/  * fAlp(m,n)) / 2. + (fBet(m,n) * (2 - 2 * r(m,n))) / (1 + TINY_Real 
    /*2*/  + r(m,n)));
}
Real gA1_t( Int m, Int n )
{
    return gRicci(m,n) + k_g * gJA1(m,n) + gA1_r(m,n) * gBet(m,n) * pow2(-1 + r(m,n))
    /*0*/  + gB_r(m,n) * ((4 * exp(-4 * gconf(m,n)) * gconf_r(m,n) * gAlp(m,n) 
    /*2*/  * Power(-1 + r(m,n),4)) / (3. * gB(m,n) * pow2(gA(m,n))) + (2 * exp(-4 
    /*3*/  * gconf(m,n)) * gAlp_r(m,n) * Power(-1 + r(m,n),4)) / (3. * gB(m,n) 
    /*2*/  * pow2(gA(m,n)))) + gA_r(m,n) * ((4 * exp(-4 * gconf(m,n)) * gconf_r(m,n) 
    /*2*/  * gAlp(m,n) * Power(-1 + r(m,n),4)) / (3. * pow3(gA(m,n))) + (2 * exp(-4 
    /*3*/  * gconf(m,n)) * gAlp_r(m,n) * Power(-1 + r(m,n),4)) / (3. * pow3(gA(m,n)))) 
    /*0*/  + (8 * exp(-4 * gconf(m,n)) * gconf_r(m,n) * gAlp_r(m,n) * Power(-1 
    /*2*/  + r(m,n),4)) / (3. * pow2(gA(m,n))) - (2 * exp(-4 * gconf(m,n)) 
    /*1*/  * gAlp_rr(m,n) * Power(-1 + r(m,n),4)) / (3. * pow2(gA(m,n))) - (4 * exp(-4 
    /*2*/  * gconf(m,n)) * gAlp_r(m,n) * pow3(-1 + r(m,n)) * (2 + r(m,n))) / (3. 
    /*1*/  * pow2(gA(m,n)) * (1 + TINY_Real + r(m,n))) + gAlp(m,n) * ((-4 * exp(-4 
    /*3*/  * gconf(m,n)) * gconf_rr(m,n) * Power(-1 + r(m,n),4)) / (3. * pow2(gA(m,n)))
    /*1*/  + (8 * exp(-4 * gconf(m,n)) * pow2(gconf_r(m,n)) * Power(-1 + r(m,n),4)) 
    /*1*/  / (3. * pow2(gA(m,n))) - (8 * exp(-4 * gconf(m,n)) * gconf_r(m,n) * pow3(-1 
    /*3*/  + r(m,n)) * (2 + r(m,n))) / (3. * pow2(gA(m,n)) * (1 + TINY_Real + r(m,n))) 
    /*1*/  + gA1(m,n) * gtrK(m,n));
}
Real fA1_t( Int m, Int n )
{
    return fRicci(m,n) + k_f * fJA1(m,n) + fA1_r(m,n) * fBet(m,n) * pow2(-1 + r(m,n))
    /*0*/  + fB_r(m,n) * ((4 * exp(-4 * fconf(m,n)) * fconf_r(m,n) * fAlp(m,n) 
    /*2*/  * Power(-1 + r(m,n),4)) / (3. * fB(m,n) * pow2(fA(m,n))) + (2 * exp(-4 
    /*3*/  * fconf(m,n)) * fAlp_r(m,n) * Power(-1 + r(m,n),4)) / (3. * fB(m,n) 
    /*2*/  * pow2(fA(m,n)))) + fA_r(m,n) * ((4 * exp(-4 * fconf(m,n)) * fconf_r(m,n) 
    /*2*/  * fAlp(m,n) * Power(-1 + r(m,n),4)) / (3. * pow3(fA(m,n))) + (2 * exp(-4 
    /*3*/  * fconf(m,n)) * fAlp_r(m,n) * Power(-1 + r(m,n),4)) / (3. * pow3(fA(m,n)))) 
    /*0*/  + (8 * exp(-4 * fconf(m,n)) * fconf_r(m,n) * fAlp_r(m,n) * Power(-1 
    /*2*/  + r(m,n),4)) / (3. * pow2(fA(m,n))) - (2 * exp(-4 * fconf(m,n)) 
    /*1*/  * fAlp_rr(m,n) * Power(-1 + r(m,n),4)) / (3. * pow2(fA(m,n))) - (4 * exp(-4 
    /*2*/  * fconf(m,n)) * fAlp_r(m,n) * pow3(-1 + r(m,n)) * (2 + r(m,n))) / (3. 
    /*1*/  * pow2(fA(m,n)) * (1 + TINY_Real + r(m,n))) + fAlp(m,n) * ((-4 * exp(-4 
    /*3*/  * fconf(m,n)) * fconf_rr(m,n) * Power(-1 + r(m,n),4)) / (3. * pow2(fA(m,n)))
    /*1*/  + (8 * exp(-4 * fconf(m,n)) * pow2(fconf_r(m,n)) * Power(-1 + r(m,n),4)) 
    /*1*/  / (3. * pow2(fA(m,n))) - (8 * exp(-4 * fconf(m,n)) * fconf_r(m,n) * pow3(-1 
    /*3*/  + r(m,n)) * (2 + r(m,n))) / (3. * pow2(fA(m,n)) * (1 + TINY_Real + r(m,n))) 
    /*1*/  + fA1(m,n) * ftrK(m,n));
}
Real gL_t( Int m, Int n )
{
    return k_g * gJL(m,n) - gBet_r(m,n) * gL(m,n) * pow2(-1 + r(m,n)) + gBet(m,n) 
    /*0*/  * gL_r(m,n) * pow2(-1 + r(m,n)) - (2 * gA1(m,n) * gAlp_r(m,n) * pow2(-1 
    /*2*/  + r(m,n))) / pow2(gA(m,n)) + (2 * gA1(m,n) * gAlp(m,n) * gB_r(m,n) * pow2(-1
    /*2*/  + r(m,n))) / (gB(m,n) * pow2(gA(m,n))) - (8 * gBet(m,n) * pow2(-1 
    /*2*/  + r(m,n))) / (pow2(gB(m,n)) * pow2(1 + TINY_Real + r(m,n))) + (2 * gA1(m,n) 
    /*1*/  * gAlp(m,n) * gA_r(m,n) * pow2(-1 + r(m,n))) / pow3(gA(m,n)) + (gBet_rr(m,n)
    /*1*/  * Power(-1 + r(m,n),4)) / pow2(gA(m,n)) + gBet_r(m,n) * ((2 * pow3(-1 
    /*3*/  + r(m,n))) / pow2(gA(m,n)) - (4 * pow3(-1 + r(m,n))) / (pow2(gB(m,n)) * (1 
    /*3*/  + TINY_Real + r(m,n)))) + gAlp(m,n) * ((12 * gA1(m,n) * gconf_r(m,n) 
    /*2*/  * pow2(-1 + r(m,n))) / pow2(gA(m,n)) - (4 * gA1(m,n) * (-1 + r(m,n))) 
    /*1*/  / (pow2(gA(m,n)) * (1 + TINY_Real + r(m,n))) + (4 * gA1(m,n) * (-1 + r(m,n)))
    /*1*/  / (pow2(gB(m,n)) * (1 + TINY_Real + r(m,n))) - (4 * pow2(-1 + r(m,n)) 
    /*2*/  * gtrK_r(m,n)) / (3. * pow2(gA(m,n))));
}
Real fL_t( Int m, Int n )
{
    return k_f * fJL(m,n) - fBet_r(m,n) * fL(m,n) * pow2(-1 + r(m,n)) + fBet(m,n) 
    /*0*/  * fL_r(m,n) * pow2(-1 + r(m,n)) - (2 * fA1(m,n) * fAlp_r(m,n) * pow2(-1 
    /*2*/  + r(m,n))) / pow2(fA(m,n)) + (2 * fA1(m,n) * fAlp(m,n) * fB_r(m,n) * pow2(-1
    /*2*/  + r(m,n))) / (fB(m,n) * pow2(fA(m,n))) - (8 * fBet(m,n) * pow2(-1 
    /*2*/  + r(m,n))) / (pow2(fB(m,n)) * pow2(1 + TINY_Real + r(m,n))) + (2 * fA1(m,n) 
    /*1*/  * fAlp(m,n) * fA_r(m,n) * pow2(-1 + r(m,n))) / pow3(fA(m,n)) + (fBet_rr(m,n)
    /*1*/  * Power(-1 + r(m,n),4)) / pow2(fA(m,n)) + fBet_r(m,n) * ((2 * pow3(-1 
    /*3*/  + r(m,n))) / pow2(fA(m,n)) - (4 * pow3(-1 + r(m,n))) / (pow2(fB(m,n)) * (1 
    /*3*/  + TINY_Real + r(m,n)))) + fAlp(m,n) * ((12 * fA1(m,n) * fconf_r(m,n) 
    /*2*/  * pow2(-1 + r(m,n))) / pow2(fA(m,n)) - (4 * fA1(m,n) * (-1 + r(m,n))) 
    /*1*/  / (pow2(fA(m,n)) * (1 + TINY_Real + r(m,n))) + (4 * fA1(m,n) * (-1 + r(m,n)))
    /*1*/  / (pow2(fB(m,n)) * (1 + TINY_Real + r(m,n))) - (4 * pow2(-1 + r(m,n)) 
    /*2*/  * ftrK_r(m,n)) / (3. * pow2(fA(m,n))));
}
