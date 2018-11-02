/** @file  eomMatter.h
 *  @brief The equations of motion for a perfect fluid.
 *  @image html matter.png
 */

Real BimetricEvolve::eq_pf_D_t( Int m, Int n )
{
    return gBet_r(m,n) * pfD(m,n) + gBet(m,n) * pfD_r(m,n) + gAlp(m,n) * (-(gDAlp(m,n)
    /*2*/  * pfD(m,n) * pfv(m,n)) - pfD_r(m,n) * pfv(m,n) - pfD(m,n) * pfv_r(m,n)) 
    /*0*/  + (2 * gBet(m,n) * pfD(m,n) - 2 * gAlp(m,n) * pfD(m,n) * pfv(m,n)) / r(m,n);
}

Real BimetricEvolve::eq_pf_S_t( Int m, Int n )
{
    return 2 * gBet_r(m,n) * pfS(m,n) + gBet(m,n) * pfS_r(m,n) + gAlp(m,n) * (gDA(m,n)
    /*1*/  * pfS(m,n) * pfv(m,n) - pfS_r(m,n) * pfv(m,n) - pfS(m,n) * pfv_r(m,n) 
    /*1*/  - gDAlp(m,n) * (pfD(m,n) + pfS(m,n) * pfv(m,n) + pftau(m,n))) + (2 
    /*1*/  * gBet(m,n) * pfS(m,n) - 2 * gAlp(m,n) * pfS(m,n) * pfv(m,n)) / r(m,n);
}

Real BimetricEvolve::eq_pf_tau_t( Int m, Int n )
{
    return gBet_r(m,n) * pftau(m,n) + gBet(m,n) * pftau_r(m,n) + gAlp(m,n) * (gK1(m,n)
    /*1*/  * pfS(m,n) * pfv(m,n) - gDAlp(m,n) * pfv(m,n) * pftau(m,n) - pfv_r(m,n) 
    /*1*/  * pftau(m,n) - pfv(m,n) * pftau_r(m,n) - (gDAlp(m,n) * pfS(m,n)) 
    /*1*/  / pow2(gA(m,n))) + (2 * gBet(m,n) * pftau(m,n) - 2 * gAlp(m,n) * pfv(m,n) 
    /*1*/  * pftau(m,n)) / r(m,n);
}

Real BimetricEvolve::eq_pf_v( Int m, Int n )
{
    return pfS(m,n) / ((TINY_Real + pfD(m,n) + pftau(m,n)) * pow2(gA(m,n)));
}

Real BimetricEvolve::eq_pf_v_r( Int m, Int n )
{
    return (-2 * gDA(m,n) * pfS(m,n) * (pfD(m,n) + pftau(m,n)) + pfS_r(m,n) 
    /*1*/  * (pfD(m,n) + pftau(m,n)) - pfS(m,n) * (pfD_r(m,n) + pftau_r(m,n))) 
    /*0*/  / (pow2(gA(m,n)) * pow2(TINY_Real + pfD(m,n) + pftau(m,n)));
}

Real BimetricEvolve::eq_pf_rho( Int m, Int n )
{
    return (pfD(m,n) + pftau(m,n)) / (gA(m,n) * pow2(gB(m,n)));
}

Real BimetricEvolve::eq_pf_j( Int m, Int n )
{
    return pfS(m,n) / (gA(m,n) * pow2(gB(m,n)));
}

Real BimetricEvolve::eq_pf_J1( Int m, Int n )
{
    return (pfS(m,n) * pfv(m,n)) / (gA(m,n) * pow2(gB(m,n)));
}

Real BimetricEvolve::eq_pf_J2( Int m, Int n )
{
    return 0;
}
