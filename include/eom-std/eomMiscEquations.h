/** @file  eomMiscEquations.h
 *  @brief Contains p, gHorz, fHorz, gBet, fBet and sigma eqs.
 */

Real BimetricEvolve::eq_invr_p_g( Int m, Int n )
{
    return (-2 * gKD(m,n)) / (k_g * fA(m,n) * P_2_1(R(m,n)));
}

Real BimetricEvolve::eq_base_p_g( Int m, Int n )
{
    return (-2 * (3 * gDB(m,n) * gKD(m,n) + gKD_r(m,n) - gK_r(m,n))) / (3. * k_g 
    /*1*/  * fA(m,n) * P_2_1(R(m,n))) + pfS(m,n) / (fA(m,n) * gA(m,n) * pow2(gB(m,n)) 
    /*1*/  * P_2_1(R(m,n)));
}

Real BimetricEvolve::eq_invr_p_f( Int m, Int n )
{
    return (2 * fKD(m,n) * pow2(fB(m,n))) / (k_f * gA(m,n) * pow2(gB(m,n)) 
    /*1*/  * P_2_1(R(m,n)));
}

Real BimetricEvolve::eq_base_p_f( Int m, Int n )
{
    return (2 * (3 * fDB(m,n) * fKD(m,n) + fKD_r(m,n) - fK_r(m,n)) * pow2(fB(m,n))) 
    /*0*/  / (3. * k_f * gA(m,n) * pow2(gB(m,n)) * P_2_1(R(m,n)));
}

Real BimetricEvolve::eq_p_t( Int m, Int n )
{
    return eq_pr(m,n) * p(m,n) * ((-2 * fB(m,n) * gAlp(m,n) * P_1_1(R(m,n))) 
    /*1*/  / (fA(m,n) * gA(m,n) * (1 + Lt(m,n)) * P_2_1(R(m,n))) + (2 * fAlp(m,n) 
    /*2*/  * fB(m,n) * P_1_2(R(m,n))) / (fA(m,n) * gA(m,n) * (1 + Lt(m,n)) 
    /*2*/  * P_2_1(R(m,n)))) + p(m,n) * ((gAlp(m,n) * (gK(m,n) + 2 * gKD(m,n) + (2 
    /*4*/  * (gK(m,n) - gKD(m,n)) * P_1_1(R(m,n))) / P_2_1(R(m,n)))) / 3. + (fAlp(m,n) 
    /*2*/  * (2 * fKD(m,n) * P_1_1(R(m,n)) + fK(m,n) * (-2 * P_1_1(R(m,n)) + 3 
    /*4*/  * P_2_1(R(m,n))))) / (3. * P_2_1(R(m,n)))) + p_r(m,n) * q(m,n) + gAlp(m,n) 
    /*0*/  * ((Lt(m,n) * (-gDAlp(m,n) - (2 * gDB(m,n) * P_1_1(R(m,n))) / P_2_1(R(m,n)))
    /*2*/  + (2 * fB(m,n) * (gSig(m,n) - fSig(m,n) * Lt(m,n)) * P_1_1(R(m,n))) 
    /*2*/  / (fA(m,n) * P_2_1(R(m,n)))) / gA(m,n) + (2 * fDB(m,n) * P_1_1(R(m,n)) 
    /*2*/  * R(m,n)) / (fA(m,n) * P_2_1(R(m,n)))) + fAlp(m,n) * (((-2 * gDB(m,n) 
    /*3*/  * P_1_2(R(m,n))) / P_2_1(R(m,n)) - (2 * fB(m,n) * (fSig(m,n) - gSig(m,n) 
    /*4*/  * Lt(m,n)) * P_1_2(R(m,n))) / (fA(m,n) * P_2_1(R(m,n)))) / gA(m,n) + (Lt(m,n)
    /*2*/  * (fDAlp(m,n) + (2 * fDB(m,n) * P_1_2(R(m,n)) * R(m,n)) / P_2_1(R(m,n)))) 
    /*1*/  / fA(m,n));
}

Real BimetricEvolve::eq_gBet( Int m, Int n )
{
    return (gAlp(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n)) + q(m,n);
}

Real BimetricEvolve::eq_fBet( Int m, Int n )
{
    return -((fAlp(m,n) * p(m,n)) / (fA(m,n) * Lt(m,n))) + q(m,n);
}

Real BimetricEvolve::eq_gBet_r( Int m, Int n )
{
    return (gAlp(m,n) * (-gDA(m,n) + gDAlp(m,n)) * p(m,n)) / (gA(m,n) * Lt(m,n)) 
    /*0*/  + (gAlp(m,n) * p_r(m,n)) / (gA(m,n) * Lt(m,n) * (1 + pow2(p(m,n)))) 
    /*0*/  + q_r(m,n);
}

Real BimetricEvolve::eq_fBet_r( Int m, Int n )
{
    return (fAlp(m,n) * (fDA(m,n) - fDAlp(m,n)) * p(m,n)) / (fA(m,n) * Lt(m,n)) 
    /*0*/  - (fAlp(m,n) * p_r(m,n)) / (fA(m,n) * Lt(m,n) * (1 + pow2(p(m,n)))) 
    /*0*/  + q_r(m,n);
}

Real BimetricEvolve::eq_gBetr_r( Int m, Int n )
{
    return eq_qr_r(m,n) + (eq_pr_r(m,n) * gAlp(m,n)) / (gA(m,n) * Lt(m,n) * (1 
    /*2*/  + pow2(p(m,n)))) - (gAlp(m,n) * pow3(p(m,n))) / (gA(m,n) * Lt(m,n) * (1 
    /*2*/  + pow2(p(m,n))) * pow2(r(m,n))) + (gAlp(m,n) * (-gDA(m,n) + gDAlp(m,n)) 
    /*1*/  * p(m,n)) / (gA(m,n) * Lt(m,n) * r(m,n));
}

Real BimetricEvolve::eq_fBetr_r( Int m, Int n )
{
    return eq_qr_r(m,n) - (eq_pr_r(m,n) * fAlp(m,n)) / (fA(m,n) * Lt(m,n) * (1 
    /*2*/  + pow2(p(m,n)))) + (fAlp(m,n) * pow3(p(m,n))) / (fA(m,n) * Lt(m,n) * (1 
    /*2*/  + pow2(p(m,n))) * pow2(r(m,n))) + (fAlp(m,n) * (fDA(m,n) - fDAlp(m,n)) 
    /*1*/  * p(m,n)) / (fA(m,n) * Lt(m,n) * r(m,n));
}

Real BimetricEvolve::eq_gHorz( Int m, Int n )
{
    return gDB(m,n) + (gA(m,n) * (-gK(m,n) + gKD(m,n))) / 3. + 1 / r(m,n);
}

Real BimetricEvolve::eq_fHorz( Int m, Int n )
{
    return fDB(m,n) + (fA(m,n) * (-fK(m,n) + fKD(m,n))) / 3. + 1 / r(m,n);
}
