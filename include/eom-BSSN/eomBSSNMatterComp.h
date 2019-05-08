/** @file  eomBSSNMatterComp.h
 *  @author Francesco Torsello
 *  @brief The compactified Valencia formulation of the hydrodynamical equations for a perfect fluid.
 *  @version 2019-05-08T11:58:34
 *  @image html BSSNevolutionComp.png
 */

Real pfv( Int m, Int n )
{
    return pfS(m,n) / (pfD(m,n) + pftau(m,n));
}
Real pfv_r( Int m, Int n )
{
    return ((pfS_r(m,n) * (pfD(m,n) + pftau(m,n)) - pfS(m,n) * (pfD_r(m,n) 
    /*3*/  + pftau_r(m,n))) * pow2(-1 + r(m,n))) / pow2(pfD(m,n) + pftau(m,n));
}
Real gD_t( Int m, Int n )
{
    return ((pfD(m,n) * (-4 * gBet(m,n) + gAlp(m,n) * (4 * pfv(m,n) - pfv_r(m,n) * (-1
    /*5*/  + pow2(r(m,n)))) - (-gBet_r(m,n) + gAlp_r(m,n) * pfv(m,n)) * (-1 
    /*4*/  + pow2(r(m,n)))) + pfD_r(m,n) * (gBet(m,n) - gAlp(m,n) * pfv(m,n)) * (-1 
    /*3*/  + pow2(r(m,n)))) * (-1 + r(m,n))) / (1 + r(m,n));
}
Real gS_t( Int m, Int n )
{
    return ((-(gBet(m,n) * gA(m,n) * (4 * pfS(m,n) - pfS_r(m,n) * (-1 
    /*5*/  + pow2(r(m,n))))) + gAlp(m,n) * (pfS(m,n) * (pfv(m,n) * (2 * gA(m,n) * (2 
    /*6*/  + gconf_r(m,n) * (-1 + pow2(r(m,n)))) + gA_r(m,n) * (-1 + pow2(r(m,n)))) 
    /*4*/  - pfv_r(m,n) * gA(m,n) * (-1 + pow2(r(m,n)))) - pfS_r(m,n) * pfv(m,n) 
    /*3*/  * gA(m,n) * (-1 + pow2(r(m,n)))) - (pfD(m,n) * gAlp_r(m,n) + pfS(m,n) * (-2 
    /*4*/  * gBet_r(m,n) + gAlp_r(m,n) * pfv(m,n)) + gAlp_r(m,n) * pftau(m,n)) * gA(m,n)
    /*2*/  * (-1 + pow2(r(m,n)))) * (-1 + r(m,n))) / (gA(m,n) * (1 + r(m,n)));
}
Real gtau_t( Int m, Int n )
{
    return (exp(-4 * gconf(m,n)) * ((exp(4 * gconf(m,n)) * gAlp(m,n) * pow2(gA(m,n)) 
    /*3*/  * (pfv(m,n) * (4 * pftau(m,n) - pftau_r(m,n) * (-1 + pow2(r(m,n)))) 
    /*4*/  - pfv_r(m,n) * pftau(m,n) * (-1 + pow2(r(m,n)))) - exp(4 * gconf(m,n)) 
    /*3*/  * gBet(m,n) * pow2(gA(m,n)) * (4 * pftau(m,n) - pftau_r(m,n) * (-1 
    /*5*/  + pow2(r(m,n)))) - (gAlp_r(m,n) * pfS(m,n) + exp(4 * gconf(m,n)) 
    /*4*/  * (-gBet_r(m,n) + gAlp_r(m,n) * pfv(m,n)) * pftau(m,n) * pow2(gA(m,n))) * (-1
    /*4*/  + pow2(r(m,n)))) * (-1 + r(m,n)) + exp(4 * gconf(m,n)) * gK1(m,n) 
    /*2*/  * gAlp(m,n) * pfS(m,n) * pfv(m,n) * pow2(gA(m,n)) * (1 + r(m,n)))) 
    /*0*/  / (pow2(gA(m,n)) * (1 + r(m,n)));
}
Real gW( Int m, Int n )
{
    return (pfD(m,n) + pftau(m,n)) / sqrt(Power(pfD(m,n) + pftau(m,n),2) - (exp(-4 
    /*3*/  * gconf(m,n)) * Power(pfS(m,n),2)) / Power(gA(m,n),2));
}
Real grhobar( Int m, Int n )
{
    return (exp(-6 * gconf(m,n)) * pfD(m,n) * sqrt(pow2(pfD(m,n) + pftau(m,n)) 
    /*2*/  - (exp(-4 * gconf(m,n)) * pow2(pfS(m,n))) / pow2(gA(m,n)))) / ((pfD(m,n) 
    /*2*/  + pftau(m,n)) * gA(m,n) * pow2(gB(m,n)));
}
Real grho( Int m, Int n )
{
    return (exp(-6 * gconf(m,n)) * (pfD(m,n) + pftau(m,n))) / (gA(m,n) 
    /*1*/  * pow2(gB(m,n)));
}
Real gj( Int m, Int n )
{
    return (exp(-6 * gconf(m,n)) * pfS(m,n)) / (gA(m,n) * pow2(gB(m,n)));
}
Real gJ11( Int m, Int n )
{
    return (exp(-2 * gconf(m,n)) * pfS(m,n) * pfv(m,n) * gA(m,n)) / pow2(gB(m,n));
}
Real gJ22( Int m, Int n )
{
    return 0;
}
Real frho( Int m, Int n )
{
    return 0;
}
Real fj( Int m, Int n )
{
    return 0;
}
Real fJ11( Int m, Int n )
{
    return 0;
}
Real fJ22( Int m, Int n )
{
    return 0;
}
