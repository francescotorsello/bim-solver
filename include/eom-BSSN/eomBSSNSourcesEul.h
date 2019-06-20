/** @file  eomBSSNSourcesEul.h
 *  @author Francesco Torsello
 *  @brief The regularized and simplified BSSN sources contained in the Eulerian evolution equations.
 *  @version 2019-04-04T11:00:52
 *  @image html BSSNsourcesEul.png
 */

Real BimetricEvolve::eq_gJK( Int m, Int n )
{
    return (gAlp(m,n) * (gJ11(m,n) + 2 * gJ22(m,n) + grho(m,n) + 2 * P_1_0(R(m,n)))) 
    /*0*/  / (2 + TINY_Real) + (exp(2 * fconf(m,n)) * fAlp(m,n) * fA(m,n) 
    /*1*/  * P_1_2(R(m,n))) / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n)) + (-((exp(2 
    /*4*/  * fconf(m,n)) * fA(m,n) * gAlp(m,n) * Lt2(m,n) * P_2_1(R(m,n))) / (TINY_Real
    /*3*/  + exp(2 * gconf(m,n)) * gA(m,n))) + (fAlp(m,n) * (2 * P_1_1(R(m,n)) 
    /*3*/  + P_2_1(R(m,n)))) / (2 + TINY_Real) + (exp(2 * fconf(m,n)) * fA(m,n) 
    /*2*/  * gAlp(m,n) * (2 * P_1_1(R(m,n)) + P_2_1(R(m,n)))) / (TINY_Real + 2 * exp(2 
    /*3*/  * gconf(m,n)) * gA(m,n))) / (TINY_Real + Lt(m,n));
}
Real BimetricEvolve::eq_fJK( Int m, Int n )
{
    return (exp(2 * gconf(m,n)) * gAlp(m,n) * gA(m,n) * P_1_1(R(m,n))) / (TINY_Real 
    /*1*/  + exp(2 * fconf(m,n)) * fA(m,n) * R(m,n)) + (fAlp(m,n) * ((fJ11(m,n) + 2 
    /*3*/  * fJ22(m,n) + frho(m,n)) * pow2(R(m,n)) - 2 * P_1_2(R(m,n)) + 2 
    /*2*/  * P_2_2(R(m,n)))) / (TINY_Real + 2 * Power(R(m,n),2)) + ((gAlp(m,n) 
    /*2*/  * (P_2_1(R(m,n)) + 2 * P_1_2(R(m,n)) * R(m,n))) / (TINY_Real + 2 
    /*2*/  * Power(R(m,n),2)) - (exp(2 * gconf(m,n)) * fAlp(m,n) * gA(m,n) * Lt2(m,n) 
    /*2*/  * P_2_1(R(m,n))) / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n) 
    /*2*/  * Power(R(m,n),2)) - (exp(2 * gconf(m,n)) * fAlp(m,n) * gA(m,n) * (2 
    /*3*/  * P_1_1(R(m,n)) - 3 * P_2_1(R(m,n)))) / (TINY_Real + 2 * exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Power(R(m,n),2))) / (TINY_Real + Lt(m,n));
}
Real BimetricEvolve::eq_gJA1( Int m, Int n )
{
    return (2 * exp(2 * fconf(m,n)) * fAlp(m,n) * fA(m,n) * P_1_2(R(m,n))) 
    /*0*/  / (TINY_Real + 3 * exp(2 * gconf(m,n)) * gA(m,n)) - (2 * gAlp(m,n) 
    /*1*/  * (gJ11(m,n) - gJ22(m,n) - P_1_0(R(m,n)) + P_2_0(R(m,n)))) / (3 + TINY_Real)
    /*0*/  + ((2 * exp(2 * fconf(m,n)) * fA(m,n) * gAlp(m,n) * P_1_1(R(m,n))) 
    /*1*/  / (TINY_Real + 3 * exp(2 * gconf(m,n)) * gA(m,n)) + (2 * fAlp(m,n) 
    /*2*/  * (P_1_1(R(m,n)) - P_2_1(R(m,n)))) / (3 + TINY_Real)) / (TINY_Real + Lt(m,n))
    /*0*/  + (2 * exp(2 * fconf(m,n)) * fA(m,n) * gAlp(m,n) * pow2(p(m,n)) 
    /*1*/  * P_2_1(R(m,n))) / (TINY_Real + 3 * exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n));
}
Real BimetricEvolve::eq_fJA1( Int m, Int n )
{
    return (2 * exp(2 * gconf(m,n)) * gAlp(m,n) * gA(m,n) * P_1_1(R(m,n))) 
    /*0*/  / (TINY_Real + 3 * exp(2 * fconf(m,n)) * fA(m,n) * R(m,n)) - (2 * fAlp(m,n) 
    /*1*/  * ((fJ11(m,n) - fJ22(m,n)) * pow2(R(m,n)) + P_1_2(R(m,n)))) / (TINY_Real + 3
    /*1*/  * Power(R(m,n),2)) + ((2 * gAlp(m,n) * (-P_2_1(R(m,n)) + P_1_2(R(m,n)) 
    /*3*/  * R(m,n))) / (TINY_Real + 3 * Power(R(m,n),2)) - (2 * exp(2 * gconf(m,n)) 
    /*2*/  * fAlp(m,n) * gA(m,n) * P_1_1(R(m,n))) / (TINY_Real + 3 * exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Power(R(m,n),2)) + (2 * exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*2*/  * gA(m,n) * Lt2(m,n) * P_2_1(R(m,n))) / (TINY_Real + 3 * exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * Power(R(m,n),2))) / (TINY_Real + Lt(m,n));
}
Real BimetricEvolve::eq_gJL( Int m, Int n )
{
    return -2 * exp(4 * gconf(m,n)) * gj(m,n) * gAlp(m,n) + (2 * exp(2 * fconf(m,n)) 
    /*1*/  * fA(m,n) * gAlp(m,n) * p(m,n) * P_2_1(R(m,n))) / (TINY_Real 
    /*1*/  + Power(gA(m,n),2));
}
Real BimetricEvolve::eq_fJL( Int m, Int n )
{
    return -2 * exp(4 * fconf(m,n)) * fj(m,n) * fAlp(m,n) - (2 * exp(2 * gconf(m,n)) 
    /*1*/  * fAlp(m,n) * gA(m,n) * p(m,n) * P_2_1(R(m,n))) / (TINY_Real 
    /*1*/  + Power(fA(m,n),2) * Power(R(m,n),2));
}
