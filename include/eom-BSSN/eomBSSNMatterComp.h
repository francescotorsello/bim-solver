/** @file  eomBSSNMatterComp.h
 *  @author Francesco Torsello
 *  @brief The compactified Valencia formulation of the hydrodynamical equations for a perfect fluid.
 *  @version 2019-05-14T16:25:43
 *  @image html BSSNevolutionComp.png
 */

Real pfv( Int m, Int n )
{
    return pfS(m,n) / (pfD(m,n) + pftau(m,n));
}
Real pfv_r( Int m, Int n )
{
    return (2 * (pfS_r(m,n) * (pfD(m,n) + pftau(m,n)) - pfS(m,n) * (pfD_r(m,n) 
    /*3*/  + pftau_r(m,n))) * pow2(r_minus(m,n))) / pow2(pfD(m,n) + pftau(m,n));
}
Real gD_t( Int m, Int n )
{
    return (-2 * r_minus(m,n) * (r_minus(m,n) * pfD_r(m,n) * (-gBet(m,n) + gAlp(m,n) * pfv(m,n)) 
    /*2*/  * r(m,n) + pfD(m,n) * (2 * gBet(m,n) + r_minus(m,n) * (-gBet_r(m,n) + gAlp_r(m,n) 
    /*4*/  * pfv(m,n)) * r(m,n) + gAlp(m,n) * (-2 * pfv(m,n) + r_minus(m,n) * pfv_r(m,n) 
    /*4*/  * r(m,n))))) / r(m,n);
}
Real gS_t( Int m, Int n )
{
    return 2 * r_minus(m,n) * (r_minus(m,n) * pfS_r(m,n) * (gBet(m,n) - gAlp(m,n) * pfv(m,n)) 
    /*1*/  - r_minus(m,n) * gAlp_r(m,n) * (pfD(m,n) + pftau(m,n)) + pfS(m,n) * (2 * r_minus(m,n) 
    /*2*/  * gBet_r(m,n) - r_minus(m,n) * gAlp_r(m,n) * pfv(m,n) + gAlp(m,n) * (-(r_minus(m,n) 
    /*4*/  * pfv_r(m,n)) + pfv(m,n) * (2 * r_minus(m,n) * gconf_r(m,n) + (r_minus(m,n) * gA_r(m,n))
    /*4*/  / gA(m,n) + 2 / r(m,n))) - (2 * gBet(m,n)) / r(m,n)));
}
Real gtau_t( Int m, Int n )
{
    return gK1(m,n) * gAlp(m,n) * pfS(m,n) * pfv(m,n) + 2 * (gBet(m,n) - gAlp(m,n) 
    /*1*/  * pfv(m,n)) * pftau_r(m,n) * pow2(r_minus(m,n)) - (2 * exp(-4 * gconf(m,n)) 
    /*1*/  * gAlp_r(m,n) * pfS(m,n) * pow2(r_minus(m,n))) / pow2(gA(m,n)) - (2 * r_minus(m,n) 
    /*1*/  * pftau(m,n) * (2 * gBet(m,n) + r_minus(m,n) * (-gBet_r(m,n) + gAlp_r(m,n) 
    /*3*/  * pfv(m,n)) * r(m,n) + gAlp(m,n) * (-2 * pfv(m,n) + r_minus(m,n) * pfv_r(m,n) 
    /*3*/  * r(m,n)))) / r(m,n);
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
