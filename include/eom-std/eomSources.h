/** @file  eomSources.h
 *  @brief The stress-energy tensor components.
 *  @image html sources.png
 */

Real BimetricEvolve::eq_g_rho( Int m, Int n )
{
    return (pfD(m,n) + pftau(m,n)) / (gA(m,n) * pow2(gB(m,n))) - P_2_0(R(m,n)) 
    /*0*/  - (fA(m,n) * Lt(m,n) * P_2_1(R(m,n))) / gA(m,n);
}

Real BimetricEvolve::eq_g_j( Int m, Int n )
{
    return pfS(m,n) / (gA(m,n) * pow2(gB(m,n))) - fA(m,n) * p(m,n) * P_2_1(R(m,n));
}

Real BimetricEvolve::eq_g_JK1( Int m, Int n )
{
    return -P_1_0(R(m,n)) - (fA(m,n) * pow2(p(m,n)) * P_2_1(R(m,n))) / (gA(m,n) 
    /*1*/  * Lt(m,n)) + (fAlp(m,n) * (-((fA(m,n) * P_1_2(R(m,n))) / gA(m,n)) + (-2 
    /*3*/  * P_1_1(R(m,n)) + P_2_1(R(m,n))) / (2. * Lt(m,n)))) / gAlp(m,n) + ((pfD(m,n)
    /*2*/  + pfS(m,n) * pfv(m,n) + pftau(m,n)) / (2. * pow2(gB(m,n))) - (fA(m,n) * (2
    /*3*/  * P_1_1(R(m,n)) + P_2_1(R(m,n)))) / (2. * Lt(m,n))) / gA(m,n);
}

Real BimetricEvolve::eq_g_JK2( Int m, Int n )
{
    return -((pfD(m,n) - pfS(m,n) * pfv(m,n) + pftau(m,n)) * P_2_0(R(m,n))) / (2. 
    /*1*/  * gA(m,n) * pow2(gB(m,n))) - (fAlp(m,n) * (pfD(m,n) - pfS(m,n) * pfv(m,n) 
    /*2*/  + pftau(m,n)) * P_2_1(R(m,n))) / (4. * gA(m,n) * gAlp(m,n) * Lt(m,n) 
    /*1*/  * pow2(gB(m,n))) - (fA(m,n) * (pfD(m,n) - pfS(m,n) * pfv(m,n) + pftau(m,n)) 
    /*1*/  * P_2_1(R(m,n))) / (4. * Lt(m,n) * pow2(gA(m,n)) * pow2(gB(m,n)));
}

Real BimetricEvolve::eq_f_rho( Int m, Int n )
{
    return -((gA(m,n) * Lt(m,n) * P_2_1(R(m,n))) / (fA(m,n) * pow2(R(m,n)))) 
    /*0*/  - P_2_2(R(m,n)) / pow2(R(m,n));
}

Real BimetricEvolve::eq_f_j( Int m, Int n )
{
    return (gA(m,n) * p(m,n) * pow2(gB(m,n)) * P_2_1(R(m,n))) / pow2(fB(m,n));
}

Real BimetricEvolve::eq_f_JK1( Int m, Int n )
{
    return -((gA(m,n) * pow2(p(m,n)) * P_2_1(R(m,n))) / (fA(m,n) * Lt(m,n) 
    /*2*/  * pow2(R(m,n)))) - P_1_3(R(m,n)) / R(m,n) - (gA(m,n) * (P_2_1(R(m,n)) + 2 
    /*2*/  * P_1_2(R(m,n)) * R(m,n))) / (2. * fA(m,n) * Lt(m,n) * pow2(R(m,n))) 
    /*0*/  + (gAlp(m,n) * (-((gA(m,n) * P_1_1(R(m,n))) / (fA(m,n) * R(m,n))) 
    /*2*/  + (P_2_1(R(m,n)) - 2 * P_1_2(R(m,n)) * R(m,n)) / (2. * Lt(m,n) 
    /*3*/  * pow2(R(m,n))))) / fAlp(m,n);
}

Real BimetricEvolve::eq_f_JK2( Int m, Int n )
{
    return -(gA(m,n) * P_2_1(R(m,n))) / (2. * fA(m,n) * Lt(m,n) * pow2(R(m,n))) 
    /*0*/  - (gAlp(m,n) * P_2_1(R(m,n))) / (2. * fAlp(m,n) * Lt(m,n) * pow2(R(m,n))) 
    /*0*/  - P_2_2(R(m,n)) / pow2(R(m,n));
}

Real BimetricEvolve::eq_g_JK( Int m, Int n )
{
    return (gAlp(m,n) * (3 * pfD(m,n) - pfS(m,n) * pfv(m,n) + 3 * pftau(m,n))) / (2. 
    /*1*/  * gA(m,n) * pow2(gB(m,n))) - (fA(m,n) * gAlp(m,n) * pow2(p(m,n)) 
    /*1*/  * P_2_1(R(m,n))) / (gA(m,n) * Lt(m,n)) + fAlp(m,n) * (-((fA(m,n) 
    /*3*/  * P_1_2(R(m,n))) / gA(m,n)) - (2 * P_1_1(R(m,n)) + P_2_1(R(m,n))) / (2. 
    /*2*/  * Lt(m,n))) + gAlp(m,n) * (-P_1_0(R(m,n)) - 2 * P_2_0(R(m,n)) - (fA(m,n) * (2
    /*3*/  * P_1_1(R(m,n)) + 3 * P_2_1(R(m,n)))) / (2. * gA(m,n) * Lt(m,n)));
}

Real BimetricEvolve::eq_g_JKD( Int m, Int n )
{
    return fAlp(m,n) * P_1_2(R(m,n)) * (-(fA(m,n) / gA(m,n)) + R(m,n) / Lt(m,n)) 
    /*0*/  + gAlp(m,n) * ((pfS(m,n) * pfv(m,n)) / (gA(m,n) * pow2(gB(m,n))) - (fA(m,n) 
    /*2*/  * pow2(p(m,n)) * P_2_1(R(m,n))) / (gA(m,n) * Lt(m,n)) + P_1_1(R(m,n)) 
    /*1*/  * (-(fA(m,n) / (gA(m,n) * Lt(m,n))) + R(m,n)));
}

Real BimetricEvolve::eq_f_JK( Int m, Int n )
{
    return -((fAlp(m,n) * gA(m,n) * pow2(p(m,n)) * P_2_1(R(m,n))) / (fA(m,n) * Lt(m,n)
    /*2*/  * pow2(R(m,n)))) + gAlp(m,n) * (-((gA(m,n) * P_1_1(R(m,n))) / (fA(m,n) 
    /*3*/  * R(m,n))) - (P_2_1(R(m,n)) + 2 * P_1_2(R(m,n)) * R(m,n)) / (2. * Lt(m,n) 
    /*2*/  * pow2(R(m,n)))) + fAlp(m,n) * (-(gA(m,n) * (3 * P_2_1(R(m,n)) + 2 
    /*3*/  * P_1_2(R(m,n)) * R(m,n))) / (2. * fA(m,n) * Lt(m,n) * pow2(R(m,n))) - (2 
    /*2*/  * P_2_2(R(m,n)) + P_1_3(R(m,n)) * R(m,n)) / pow2(R(m,n)));
}

Real BimetricEvolve::eq_f_JKD( Int m, Int n )
{
    return (gAlp(m,n) * P_1_1(R(m,n)) * (1 / Lt(m,n) - (gA(m,n) * R(m,n)) / fA(m,n)))
    /*0*/  / pow2(R(m,n)) + fAlp(m,n) * (-((gA(m,n) * pow2(p(m,n)) * P_2_1(R(m,n))) 
    /*2*/  / (fA(m,n) * Lt(m,n) * pow2(R(m,n)))) + (P_1_2(R(m,n)) * (1 - (gA(m,n) 
    /*4*/  * R(m,n)) / (fA(m,n) * Lt(m,n)))) / pow2(R(m,n)));
}
