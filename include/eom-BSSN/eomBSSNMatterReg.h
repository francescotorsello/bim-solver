/** @file  eomBSSNMatterReg.h
 *  @author Francesco Torsello
 *  @brief The regularized Valencia formulation of the hydrodynamical equations for a perfect fluid.
 *  @version 2018-10-01T16:45:45
 *  @image html BSSNevolutionReg.png
 */

Real BimetricEvolve::eq_pf_gD_t( Int m, Int n )
{
    return (pfD_r(m,n) * (gBet(m,n) - gAlp(m,n) * pfv(m,n)) * (TINY_Real + r(m,n)) 
    /*1*/  + pfD(m,n) * (2 * gBet(m,n) + gBet_r(m,n) * r(m,n) - gAlp(m,n) * (pfv_r(m,n)
    /*3*/  * r(m,n) + pfv(m,n) * (2 + gDAlp(m,n) * r(m,n))))) / (TINY_Real + r(m,n));
}
Real BimetricEvolve::eq_pf_gS_t( Int m, Int n )
{
    return gBet_r(m,n) * pfS(m,n) + gBet(m,n) * pfS_r(m,n) + gDA(m,n) * gAlp(m,n) 
    /*0*/  * pfS(m,n) * pfv(m,n) + 2 * gconf(m,n) * gDconf(m,n) * gAlp(m,n) * pfS(m,n) 
    /*0*/  * pfv(m,n) - gAlp(m,n) * pfS_r(m,n) * pfv(m,n) - gAlp(m,n) * pfS(m,n) 
    /*0*/  * pfv_r(m,n) - gDAlp(m,n) * gAlp(m,n) * (pfD(m,n) + pfS(m,n) * pfv(m,n) 
    /*1*/  + pftau(m,n)) + (2 * pfS(m,n) * (gBet(m,n) - gAlp(m,n) * pfv(m,n))) / r(m,n);
}
Real BimetricEvolve::eq_pf_gtau_t( Int m, Int n )
{
    return gDAlp(m,n) * pfD(m,n) * gAlp(m,n) * pfv(m,n) + gBet(m,n) * pftau_r(m,n) 
    /*0*/  - (2 * gDAlp(m,n) * gAlp(m,n) * pfS(m,n)) / (TINY_Real + exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),2)) + (pftau(m,n) * (2 * gBet(m,n) + gBet_r(m,n) * r(m,n)))
    /*0*/  / (TINY_Real + r(m,n)) + gAlp(m,n) * (pfD_r(m,n) * pfv(m,n) + pfD(m,n) 
    /*1*/  * pfv_r(m,n) + pfv(m,n) * (gK1(m,n) * pfS(m,n) + (2 * pfD(m,n)) / r(m,n)) 
    /*1*/  + (exp(-4 * gconf(m,n)) * (-(pfS_r(m,n) * r(m,n)) + 2 * pfS(m,n) * (-1 
    /*4*/  + gDA(m,n) * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)))) 
    /*1*/  / (pow2(gA(m,n)) * r(m,n)));
}
Real BimetricEvolve::eq_pf_gv( Int m, Int n )
{
    return pfS(m,n) / (TINY_Real + pfD(m,n) + pftau(m,n));
}
Real BimetricEvolve::eq_pf_gv_r( Int m, Int n )
{
    return (pfS_r(m,n) * (pfD(m,n) + pftau(m,n)) - pfS(m,n) * (pfD_r(m,n) 
    /*2*/  + pftau_r(m,n))) / (TINY_Real + Power(pfD(m,n),2) + 2 * pfD(m,n) * pftau(m,n)
    /*1*/  + Power(pftau(m,n),2));
}
Real BimetricEvolve::eq_pf_gW( Int m, Int n )
{
    return (pfD(m,n) + pftau(m,n)) / (TINY_Real + Sqrt(Power(pfD(m,n) + pftau(m,n),2)
    /*2*/  - (exp(-4 * gconf(m,n)) * Power(pfS(m,n),2)) / Power(gA(m,n),2)));
}
Real BimetricEvolve::eq_pf_grhobar( Int m, Int n )
{
    return (pfD(m,n) * sqrt(pow2(pfD(m,n) + pftau(m,n)) - (exp(-4 * gconf(m,n)) 
    /*3*/  * pow2(pfS(m,n))) / pow2(gA(m,n)))) / (TINY_Real + exp(6 * gconf(m,n)) 
    /*1*/  * (pfD(m,n) + pftau(m,n)) * gA(m,n) * Power(gB(m,n),2));
}
Real BimetricEvolve::eq_pf_grho( Int m, Int n )
{
    return (pfD(m,n) + pftau(m,n)) / (TINY_Real + exp(6 * gconf(m,n)) * gA(m,n) 
    /*1*/  * Power(gB(m,n),2));
}
Real BimetricEvolve::eq_pf_gj( Int m, Int n )
{
    return pfS(m,n) / (TINY_Real + exp(6 * gconf(m,n)) * gA(m,n) * Power(gB(m,n),2));
}
Real BimetricEvolve::eq_pf_gJ11( Int m, Int n )
{
    return (pfS(m,n) * pfv(m,n) * gA(m,n)) / (TINY_Real + exp(2 * gconf(m,n)) 
    /*1*/  * Power(gB(m,n),2));
}
Real BimetricEvolve::eq_pf_gJ22( Int m, Int n )
{
    return 0;
}
Real BimetricEvolve::eq_pf_frho( Int m, Int n )
{
    return 0;
}
Real BimetricEvolve::eq_pf_fj( Int m, Int n )
{
    return 0;
}
Real BimetricEvolve::eq_pf_fJ11( Int m, Int n )
{
    return 0;
}
Real BimetricEvolve::eq_pf_fJ22( Int m, Int n )
{
    return 0;
}
