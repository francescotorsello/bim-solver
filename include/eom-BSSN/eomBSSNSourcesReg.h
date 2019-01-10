/** @file  eomBSSNSourcesReg.h
 *  @author Francesco Torsello
 *  @brief The regularized BSSN sources contained in the evolution equations.
 *  @version 2019-01-09T13:39:39
 *  @image html BSSNsourcesReg.png
 */

Real BimetricEvolve::eq_gJK( Int m, Int n )
{
    return fAlp(m,n) * ((exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) 
    /*2*/  * P_1_2(R(m,n))) / gA(m,n) + (2 * P_1_1(R(m,n)) + P_2_1(R(m,n))) / (2. 
    /*2*/  * Lt(m,n))) + gAlp(m,n) * ((gJ11(m,n) + 2 * gJ22(m,n) + grho(m,n) + 2 
    /*2*/  * P_1_0(R(m,n))) / 2. + (exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * (2
    /*3*/  * P_1_1(R(m,n)) + (1 - 2 * Lt2(m,n)) * P_2_1(R(m,n)))) / (2. * gA(m,n) 
    /*2*/  * Lt(m,n)));
}
Real BimetricEvolve::eq_fJK( Int m, Int n )
{
    return fAlp(m,n) * (-(exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * (2 
    /*3*/  * P_1_1(R(m,n)) + (-3 + 2 * Lt2(m,n)) * P_2_1(R(m,n)))) / (2. * fA(m,n) 
    /*2*/  * Lt(m,n) * pow2(R(m,n))) + ((fJ11(m,n) + 2 * fJ22(m,n) + frho(m,n)) 
    /*2*/  * pow2(R(m,n)) - 2 * P_1_2(R(m,n)) + 2 * P_2_2(R(m,n))) / (2. 
    /*2*/  * pow2(R(m,n)))) + gAlp(m,n) * ((exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*2*/  * gA(m,n) * P_1_1(R(m,n))) / (fA(m,n) * R(m,n)) + (P_2_1(R(m,n)) + 2 
    /*2*/  * P_1_2(R(m,n)) * R(m,n)) / (2. * Lt(m,n) * pow2(R(m,n))));
}
Real BimetricEvolve::eq_gJA1( Int m, Int n )
{
    return fAlp(m,n) * ((2 * exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) 
    /*2*/  * P_1_2(R(m,n))) / (3. * gA(m,n)) + (2 * (P_1_1(R(m,n)) - P_2_1(R(m,n)))) 
    /*1*/  / (3. * Lt(m,n))) + gAlp(m,n) * ((-2 * (gJ11(m,n) - gJ22(m,n) - P_1_0(R(m,n))
    /*3*/  + P_2_0(R(m,n)))) / 3. + (2 * exp(2 * fconf(m,n) - 2 * gconf(m,n)) 
    /*2*/  * fA(m,n) * (P_1_1(R(m,n)) + pow2(p(m,n)) * P_2_1(R(m,n)))) / (3. * gA(m,n) 
    /*2*/  * Lt(m,n)));
}
Real BimetricEvolve::eq_gJA2( Int m, Int n )
{
    return fAlp(m,n) * (-(exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) 
    /*2*/  * P_1_2(R(m,n))) / (3. * gA(m,n)) + (-P_1_1(R(m,n)) + P_2_1(R(m,n))) / (3. 
    /*2*/  * Lt(m,n))) + gAlp(m,n) * ((gJ11(m,n) - gJ22(m,n) - P_1_0(R(m,n)) 
    /*2*/  + P_2_0(R(m,n))) / 3. - (exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) 
    /*2*/  * (P_1_1(R(m,n)) + pow2(p(m,n)) * P_2_1(R(m,n)))) / (3. * gA(m,n) * Lt(m,n)));
}
Real BimetricEvolve::eq_fJA1( Int m, Int n )
{
    return fAlp(m,n) * ((2 * (-fJ11(m,n) + fJ22(m,n) - P_1_2(R(m,n)) / pow2(R(m,n))))
    /*1*/  / 3. - (2 * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) 
    /*2*/  * (P_1_1(R(m,n)) - Lt2(m,n) * P_2_1(R(m,n)))) / (3. * fA(m,n) * Lt(m,n) 
    /*2*/  * pow2(R(m,n)))) + gAlp(m,n) * ((2 * exp(-2 * fconf(m,n) + 2 * gconf(m,n)) 
    /*2*/  * gA(m,n) * P_1_1(R(m,n))) / (3. * fA(m,n) * R(m,n)) + (2 * (-P_2_1(R(m,n)) 
    /*3*/  + P_1_2(R(m,n)) * R(m,n))) / (3. * Lt(m,n) * pow2(R(m,n))));
}
Real BimetricEvolve::eq_fJA2( Int m, Int n )
{
    return fAlp(m,n) * ((fJ11(m,n) - fJ22(m,n) + P_1_2(R(m,n)) / pow2(R(m,n))) / 3. 
    /*1*/  + (exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * (P_1_1(R(m,n)) 
    /*3*/  - Lt2(m,n) * P_2_1(R(m,n)))) / (3. * fA(m,n) * Lt(m,n) * pow2(R(m,n)))) 
    /*0*/  + gAlp(m,n) * (-(exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) 
    /*2*/  * P_1_1(R(m,n))) / (3. * fA(m,n) * R(m,n)) + (P_2_1(R(m,n)) - P_1_2(R(m,n)) 
    /*2*/  * R(m,n)) / (3. * Lt(m,n) * pow2(R(m,n))));
}
Real BimetricEvolve::eq_gJL( Int m, Int n )
{
    return gAlp(m,n) * (-2 * exp(4 * gconf(m,n)) * gj(m,n) + (2 * exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * p(m,n) * P_2_1(R(m,n))) / pow2(gA(m,n)));
}
Real BimetricEvolve::eq_fJL( Int m, Int n )
{
    return fAlp(m,n) * (-2 * exp(4 * fconf(m,n)) * fj(m,n) - (2 * exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * p(m,n) * P_2_1(R(m,n))) / (pow2(fA(m,n)) * pow2(R(m,n))));
}
