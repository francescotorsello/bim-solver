/** @file  eomBSSNMatterReg.h
 *  @author Francesco Torsello
 *  @brief The regularized Valencia formulation of the hydrodynamical equations for a perfect fluid.
 *  @version 2018-11-07T10:13:17
 *  @image html BSSNevolutionReg.png
 */

Real BimetricEvolve::eq_pf_gD_t( Int m, Int n )
{
    return pfD_convr(m,n) - pfD_r(m,n) * gAlp(m,n) * pfv(m,n) + pfD(m,n) 
    /*0*/  * (gBet_r(m,n) + 2 * gBetr(m,n) + gAlp(m,n) * (-pfv_r(m,n) + pfv(m,n) 
    /*2*/  * (-gDAlp(m,n) - 2 / (TINY_Real + r(m,n)))));
}
Real BimetricEvolve::eq_pf_gS_t( Int m, Int n )
{
    return pfS_convr(m,n) - gAlp(m,n) * pfS_r(m,n) * pfv(m,n) + pfS(m,n) * (2 
    /*1*/  * gBet_r(m,n) + 2 * gBetr(m,n) - gDAlp(m,n) * gAlp(m,n) * pfv(m,n) 
    /*1*/  + gAlp(m,n) * ((gDA(m,n) + 2 * gconf(m,n) * gDconf(m,n)) * pfv(m,n) 
    /*2*/  - pfv_r(m,n))) - gDAlp(m,n) * gAlp(m,n) * (pfD(m,n) + pftau(m,n)) - (2 
    /*1*/  * gAlp(m,n) * pfS(m,n) * pfv(m,n)) / (TINY_Real + r(m,n));
}
Real BimetricEvolve::eq_pf_gtau_t( Int m, Int n )
{
    return pftau_convr(m,n) + gK1(m,n) * gAlp(m,n) * pfS(m,n) * pfv(m,n) 
    /*0*/  + (gBet_r(m,n) + 2 * gBetr(m,n) - gDAlp(m,n) * gAlp(m,n) * pfv(m,n) 
    /*1*/  - gAlp(m,n) * pfv_r(m,n)) * pftau(m,n) - gAlp(m,n) * pfv(m,n) * pftau_r(m,n)
    /*0*/  - (gDAlp(m,n) * gAlp(m,n) * pfS(m,n)) / (TINY_Real + exp(4 * gconf(m,n)) 
    /*1*/  * Power(gA(m,n),2)) - (2 * gAlp(m,n) * pfv(m,n) * pftau(m,n)) / (TINY_Real 
    /*1*/  + r(m,n));
}
Real BimetricEvolve::eq_pf_gv( Int m, Int n )
{
    return pfS(m,n) / (TINY_Real + pfD(m,n) + pftau(m,n));
}
Real BimetricEvolve::eq_pf_gv_r( Int m, Int n )
{
    return (pfS_r(m,n) * (pfD(m,n) + pftau(m,n)) - pfS(m,n) * (pfD_r(m,n) 
    /*2*/  + pftau_r(m,n))) / (TINY_Real + Power(pfD(m,n) + pftau(m,n),2));
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
