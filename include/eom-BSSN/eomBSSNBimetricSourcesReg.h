/** @file  eomBSSNBimetricSourcesReg.h
 *  @author Francesco Torsello
 *  @brief The regularized bimetric sources in the cBSSN formalism.
 *  @version 2019-01-09T13:42:13
 *  @image html cBSSNBimetricSourcesReg.png
 */

Real BimetricEvolve::eq_grhob( Int m, Int n )
{
    return -P_2_0(R(m,n)) - (exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * Lt(m,n)
    /*1*/  * P_2_1(R(m,n))) / gA(m,n);
}
Real BimetricEvolve::eq_frhob( Int m, Int n )
{
    return -(((exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * Lt(m,n) 
    /*3*/  * P_2_1(R(m,n))) / fA(m,n) + P_2_2(R(m,n))) / pow2(R(m,n)));
}
Real BimetricEvolve::eq_gjb_u( Int m, Int n )
{
    return -((exp(2 * fconf(m,n) - 4 * gconf(m,n)) * fA(m,n) * p(m,n) * P_2_1(R(m,n)))
    /*1*/  / pow2(gA(m,n)));
}
Real BimetricEvolve::eq_fjb_u( Int m, Int n )
{
    return (exp(-4 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n) * p(m,n) * P_2_1(R(m,n)))
    /*0*/  / (pow2(fA(m,n)) * pow2(R(m,n)));
}
Real BimetricEvolve::eq_gjb_d( Int m, Int n )
{
    return -(exp(2 * fconf(m,n)) * fA(m,n) * p(m,n) * P_2_1(R(m,n)));
}
Real BimetricEvolve::eq_fjb_d( Int m, Int n )
{
    return (exp(2 * gconf(m,n)) * gA(m,n) * p(m,n) * P_2_1(R(m,n))) / pow2(R(m,n));
}
Real BimetricEvolve::eq_gJb1_dd( Int m, Int n )
{
    return exp(4 * gconf(m,n)) * pow2(gA(m,n)) * (P_2_0(R(m,n)) + ((fAlp(m,n) 
    /*3*/  / gAlp(m,n) - (exp(2 * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * pow2(p(m,n)))
    /*3*/  / gA(m,n)) * P_2_1(R(m,n))) / Lt(m,n));
}
Real BimetricEvolve::eq_gJb2_dd( Int m, Int n )
{
    return (exp(2 * gconf(m,n)) * pow2(gB(m,n)) * pow2(r(m,n)) * (exp(2 * gconf(m,n))
    /*2*/  * gA(m,n) * (gAlp(m,n) * Lt(m,n) * P_1_0(R(m,n)) + fAlp(m,n) 
    /*3*/  * P_1_1(R(m,n))) + exp(2 * fconf(m,n)) * fA(m,n) * (gAlp(m,n) * P_1_1(R(m,n))
    /*3*/  + fAlp(m,n) * Lt(m,n) * P_1_2(R(m,n))))) / (gAlp(m,n) * gA(m,n) * Lt(m,n));
}
Real BimetricEvolve::eq_fJb1_dd( Int m, Int n )
{
    return (exp(2 * fconf(m,n)) * fA(m,n) * (-(exp(2 * gconf(m,n)) * fAlp(m,n) 
    /*3*/  * gA(m,n) * pow2(p(m,n)) * P_2_1(R(m,n))) + exp(2 * fconf(m,n)) * fA(m,n) 
    /*2*/  * (gAlp(m,n) * P_2_1(R(m,n)) + fAlp(m,n) * Lt(m,n) * P_2_2(R(m,n))))) 
    /*0*/  / (fAlp(m,n) * Lt(m,n) * pow2(R(m,n)));
}
Real BimetricEvolve::eq_fJb2_dd( Int m, Int n )
{
    return (exp(-2 * fconf(m,n) + 4 * gconf(m,n)) * pow2(gB(m,n)) * pow2(r(m,n)) 
    /*1*/  * (fAlp(m,n) * (-(exp(2 * gconf(m,n)) * gA(m,n) * (P_1_1(R(m,n)) 
    /*5*/  - P_2_1(R(m,n)))) - exp(2 * fconf(m,n)) * fA(m,n) * Lt(m,n) * (P_1_2(R(m,n))
    /*4*/  - P_2_2(R(m,n)))) + gAlp(m,n) * (exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n) 
    /*3*/  * P_1_1(R(m,n)) + exp(2 * fconf(m,n)) * fA(m,n) * P_1_2(R(m,n))) * R(m,n))) 
    /*0*/  / (fAlp(m,n) * fA(m,n) * Lt(m,n));
}
Real BimetricEvolve::eq_gJb1_ud( Int m, Int n )
{
    return P_2_0(R(m,n)) + ((fAlp(m,n) / gAlp(m,n) - (exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * fA(m,n) * pow2(p(m,n))) / gA(m,n)) * P_2_1(R(m,n))) / Lt(m,n);
}
Real BimetricEvolve::eq_gJb2_ud( Int m, Int n )
{
    return (exp(-2 * gconf(m,n)) * (exp(2 * gconf(m,n)) * gA(m,n) * (gAlp(m,n) 
    /*3*/  * Lt(m,n) * P_1_0(R(m,n)) + fAlp(m,n) * P_1_1(R(m,n))) + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * (gAlp(m,n) * P_1_1(R(m,n)) + fAlp(m,n) * Lt(m,n) 
    /*3*/  * P_1_2(R(m,n))))) / (gAlp(m,n) * gA(m,n) * Lt(m,n));
}
Real BimetricEvolve::eq_fJb1_ud( Int m, Int n )
{
    return (((gAlp(m,n) / fAlp(m,n) - (exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n)
    /*4*/  * pow2(p(m,n))) / fA(m,n)) * P_2_1(R(m,n))) / Lt(m,n) + P_2_2(R(m,n))) 
    /*0*/  / pow2(R(m,n));
}
Real BimetricEvolve::eq_fJb2_ud( Int m, Int n )
{
    return (exp(-2 * fconf(m,n)) * (fAlp(m,n) * (-(exp(2 * gconf(m,n)) * gA(m,n) 
    /*4*/  * (P_1_1(R(m,n)) - P_2_1(R(m,n)))) - exp(2 * fconf(m,n)) * fA(m,n) * Lt(m,n)
    /*3*/  * (P_1_2(R(m,n)) - P_2_2(R(m,n)))) + gAlp(m,n) * (exp(2 * gconf(m,n)) 
    /*3*/  * gA(m,n) * Lt(m,n) * P_1_1(R(m,n)) + exp(2 * fconf(m,n)) * fA(m,n) 
    /*3*/  * P_1_2(R(m,n))) * R(m,n))) / (fAlp(m,n) * fA(m,n) * Lt(m,n) * pow2(R(m,n)));
}
Real BimetricEvolve::eq_gJb1_du( Int m, Int n )
{
    return P_2_0(R(m,n)) + ((fAlp(m,n) / gAlp(m,n) - (exp(2 * fconf(m,n) - 2 
    /*4*/  * gconf(m,n)) * fA(m,n) * pow2(p(m,n))) / gA(m,n)) * P_2_1(R(m,n))) / Lt(m,n);
}
Real BimetricEvolve::eq_gJb2_du( Int m, Int n )
{
    return (exp(-2 * gconf(m,n)) * (exp(2 * gconf(m,n)) * gA(m,n) * (gAlp(m,n) 
    /*3*/  * Lt(m,n) * P_1_0(R(m,n)) + fAlp(m,n) * P_1_1(R(m,n))) + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * (gAlp(m,n) * P_1_1(R(m,n)) + fAlp(m,n) * Lt(m,n) 
    /*3*/  * P_1_2(R(m,n))))) / (gAlp(m,n) * gA(m,n) * Lt(m,n));
}
Real BimetricEvolve::eq_fJb1_du( Int m, Int n )
{
    return (((gAlp(m,n) / fAlp(m,n) - (exp(-2 * fconf(m,n) + 2 * gconf(m,n)) * gA(m,n)
    /*4*/  * pow2(p(m,n))) / fA(m,n)) * P_2_1(R(m,n))) / Lt(m,n) + P_2_2(R(m,n))) 
    /*0*/  / pow2(R(m,n));
}
Real BimetricEvolve::eq_fJb2_du( Int m, Int n )
{
    return (exp(-2 * fconf(m,n)) * (fAlp(m,n) * (-(exp(2 * gconf(m,n)) * gA(m,n) 
    /*4*/  * (P_1_1(R(m,n)) - P_2_1(R(m,n)))) - exp(2 * fconf(m,n)) * fA(m,n) * Lt(m,n)
    /*3*/  * (P_1_2(R(m,n)) - P_2_2(R(m,n)))) + gAlp(m,n) * (exp(2 * gconf(m,n)) 
    /*3*/  * gA(m,n) * Lt(m,n) * P_1_1(R(m,n)) + exp(2 * fconf(m,n)) * fA(m,n) 
    /*3*/  * P_1_2(R(m,n))) * R(m,n))) / (fAlp(m,n) * fA(m,n) * Lt(m,n) * pow2(R(m,n)));
}
Real BimetricEvolve::eq_gJb1_uu( Int m, Int n )
{
    return (exp(-4 * gconf(m,n)) * (P_2_0(R(m,n)) + ((fAlp(m,n) / gAlp(m,n) - (exp(2 
    /*6*/  * fconf(m,n) - 2 * gconf(m,n)) * fA(m,n) * pow2(p(m,n))) / gA(m,n)) 
    /*3*/  * P_2_1(R(m,n))) / Lt(m,n))) / pow2(gA(m,n));
}
Real BimetricEvolve::eq_gJb2_uu( Int m, Int n )
{
    return (exp(-6 * gconf(m,n)) * (exp(2 * gconf(m,n)) * gA(m,n) * (gAlp(m,n) 
    /*3*/  * Lt(m,n) * P_1_0(R(m,n)) + fAlp(m,n) * P_1_1(R(m,n))) + exp(2 * fconf(m,n))
    /*2*/  * fA(m,n) * (gAlp(m,n) * P_1_1(R(m,n)) + fAlp(m,n) * Lt(m,n) 
    /*3*/  * P_1_2(R(m,n))))) / (gAlp(m,n) * gA(m,n) * Lt(m,n) * pow2(gB(m,n)) 
    /*1*/  * pow2(r(m,n)));
}
Real BimetricEvolve::eq_fJb1_uu( Int m, Int n )
{
    return (exp(-6 * fconf(m,n)) * (-(exp(2 * gconf(m,n)) * fAlp(m,n) * gA(m,n) 
    /*3*/  * pow2(p(m,n)) * P_2_1(R(m,n))) + exp(2 * fconf(m,n)) * fA(m,n) * (gAlp(m,n)
    /*3*/  * P_2_1(R(m,n)) + fAlp(m,n) * Lt(m,n) * P_2_2(R(m,n))))) / (fAlp(m,n) 
    /*1*/  * Lt(m,n) * pow2(R(m,n)) * pow3(fA(m,n)));
}
Real BimetricEvolve::eq_fJb2_uu( Int m, Int n )
{
    return (exp(-2 * (fconf(m,n) + 2 * gconf(m,n))) * (fAlp(m,n) * (-(exp(2 
    /*5*/  * gconf(m,n)) * gA(m,n) * (P_1_1(R(m,n)) - P_2_1(R(m,n)))) - exp(2 
    /*4*/  * fconf(m,n)) * fA(m,n) * Lt(m,n) * (P_1_2(R(m,n)) - P_2_2(R(m,n)))) 
    /*2*/  + gAlp(m,n) * (exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n) * P_1_1(R(m,n)) 
    /*3*/  + exp(2 * fconf(m,n)) * fA(m,n) * P_1_2(R(m,n))) * R(m,n))) / (fAlp(m,n) 
    /*1*/  * fA(m,n) * Lt(m,n) * pow2(gB(m,n)) * pow2(r(m,n)) * Power(R(m,n),4));
}
