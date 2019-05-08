/** @file  eomBSSNSourcesCompEul.h
 *  @author Francesco Torsello
 *  @brief The compactified BSSN sources contained in the Eulerian evolution equations.
 *  @version 2019-05-08T11:14:49
 *  @image html BSSNsourcesCompEul.png
 */

Real gJK( Int m, Int n )
{
    return (exp(-2 * gconf(m,n)) * (exp(2 * gconf(m,n)) * gA(m,n) * (gAlp(m,n) 
    /*3*/  * Lt(m,n) * (gJ11(m,n) + 2 * gJ22(m,n) + grho(m,n) + 2 * P_1_0(R(m,n))) 
    /*3*/  + fAlp(m,n) * (2 * P_1_1(R(m,n)) + P_2_1(R(m,n)))) + exp(2 * fconf(m,n)) 
    /*2*/  * fA(m,n) * (2 * fAlp(m,n) * Lt(m,n) * P_1_2(R(m,n)) + gAlp(m,n) * (2 
    /*4*/  * P_1_1(R(m,n)) - (Lt2(m,n) + pow2(p(m,n))) * P_2_1(R(m,n)))))) / (2. 
    /*1*/  * gA(m,n) * Lt(m,n));
}
Real fJK( Int m, Int n )
{
    return (exp(-2 * fconf(m,n)) * (fAlp(m,n) * (-(exp(2 * gconf(m,n)) * gA(m,n) * (2
    /*5*/  * P_1_1(R(m,n)) + (-2 + Lt2(m,n) + pow2(p(m,n))) * P_2_1(R(m,n)))) - exp(2
    /*4*/  * fconf(m,n)) * fA(m,n) * Lt(m,n) * (-((fJ11(m,n) + 2 * fJ22(m,n) 
    /*6*/  + frho(m,n)) * pow2(R(m,n))) + 2 * P_1_2(R(m,n)) - 2 * P_2_2(R(m,n)))) 
    /*2*/  + gAlp(m,n) * (2 * exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n) * P_1_1(R(m,n)) 
    /*3*/  * R(m,n) + exp(2 * fconf(m,n)) * fA(m,n) * (P_2_1(R(m,n)) + 2 * P_1_2(R(m,n))
    /*4*/  * R(m,n))))) / (2. * fA(m,n) * Lt(m,n) * pow2(R(m,n)));
}
Real gJA1( Int m, Int n )
{
    return (2 * (-(gJ11(m,n) * gAlp(m,n)) + (exp(-2 * gconf(m,n)) * (exp(2 
    /*5*/  * gconf(m,n)) * gJ22(m,n) * gAlp(m,n) * gA(m,n) * Lt(m,n) + fAlp(m,n) 
    /*4*/  * (exp(2 * fconf(m,n)) * fA(m,n) * Lt(m,n) * P_1_2(R(m,n)) + exp(2 
    /*6*/  * gconf(m,n)) * gA(m,n) * (P_1_1(R(m,n)) - P_2_1(R(m,n)))) + gAlp(m,n) 
    /*4*/  * (exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n) * (P_1_0(R(m,n)) - P_2_0(R(m,n)))
    /*5*/  + exp(2 * fconf(m,n)) * fA(m,n) * (P_1_1(R(m,n)) + pow2(p(m,n)) 
    /*6*/  * P_2_1(R(m,n)))))) / (gA(m,n) * Lt(m,n)))) / 3.;
}
Real fJA1( Int m, Int n )
{
    return (-2 * exp(-2 * fconf(m,n)) * (fAlp(m,n) * (exp(2 * fconf(m,n)) * fA(m,n) 
    /*3*/  * Lt(m,n) * ((fJ11(m,n) - fJ22(m,n)) * pow2(R(m,n)) + P_1_2(R(m,n))) + exp(2
    /*4*/  * gconf(m,n)) * gA(m,n) * (P_1_1(R(m,n)) - Lt2(m,n) * P_2_1(R(m,n)))) 
    /*2*/  + gAlp(m,n) * (-(exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n) * P_1_1(R(m,n)) 
    /*4*/  * R(m,n)) + exp(2 * fconf(m,n)) * fA(m,n) * (P_2_1(R(m,n)) - P_1_2(R(m,n)) 
    /*4*/  * R(m,n))))) / (3. * fA(m,n) * Lt(m,n) * pow2(R(m,n)));
}
Real gJL( Int m, Int n )
{
    return 2 * gAlp(m,n) * (-(exp(4 * gconf(m,n)) * gj(m,n)) + (exp(2 * fconf(m,n)) 
    /*2*/  * fA(m,n) * p(m,n) * P_2_1(R(m,n))) / pow2(gA(m,n)));
}
Real fJL( Int m, Int n )
{
    return 2 * fAlp(m,n) * (-(exp(4 * fconf(m,n)) * fj(m,n)) - (exp(2 * gconf(m,n)) 
    /*2*/  * gA(m,n) * p(m,n) * P_2_1(R(m,n))) / (pow2(fA(m,n)) * pow2(R(m,n))));
}
