/** @file  eomBSSNMiscEquationsReg.h
 *  @author Francesco Torsello
 *  @brief Contains the regularized p, gHorz, fHorz, gBet, fBet.
 *  @version 2019-01-09T13:43:05
 *  @image html BSSNevolutionReg.png
 */

Real BimetricEvolve::eq_gBet( Int m, Int n )
{
    return (gAlp(m,n) * p(m,n)) / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n) 
    /*1*/  * Lt(m,n)) + q(m,n);
}
Real BimetricEvolve::eq_fBet( Int m, Int n )
{
    return -((fAlp(m,n) * p(m,n)) / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n) 
    /*2*/  * Lt(m,n))) + q(m,n);
}
Real BimetricEvolve::eq_gBet_r( Int m, Int n )
{
    return (((gDAlp(m,n) / gA(m,n) - ((gDA(m,n) + 2 * gDconf(m,n)) * gAlp(m,n)) 
    /*3*/  / gA(m,n)) * p(m,n)) / (TINY_Real + exp(2 * gconf(m,n))) + (gAlp(m,n) 
    /*2*/  * p_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n) * Lt2(m,n))) 
    /*0*/  / (TINY_Real + Lt(m,n)) + q_r(m,n);
}
Real BimetricEvolve::eq_fBet_r( Int m, Int n )
{
    return (((-(fDAlp(m,n) / fA(m,n)) + ((fDA(m,n) + 2 * fDconf(m,n)) * fAlp(m,n)) 
    /*3*/  / fA(m,n)) * p(m,n)) / (TINY_Real + exp(2 * fconf(m,n))) - (fAlp(m,n) 
    /*2*/  * p_r(m,n)) / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n) * Lt2(m,n))) 
    /*0*/  / (TINY_Real + Lt(m,n)) + q_r(m,n);
}
Real BimetricEvolve::eq_gBet_rr( Int m, Int n )
{
    return (((((2 * gDAlp(m,n)) / gA(m,n) - (2 * (gDA(m,n) + 2 * gDconf(m,n)) 
    /*5*/  * gAlp(m,n)) / gA(m,n)) * p_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n))) 
    /*2*/  + (gAlp(m,n) * p_rr(m,n)) / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n))) 
    /*1*/  / (TINY_Real + Lt2(m,n)) + (p(m,n) * ((-2 * gDA(m,n) * gDAlp(m,n) - 4 
    /*4*/  * gDconf(m,n) * gDAlp(m,n) + gDAlp_r(m,n)) / gA(m,n) + (gAlp(m,n) 
    /*4*/  * (-gDA_r(m,n) + 4 * gDA(m,n) * gDconf(m,n) - 2 * gDconf_r(m,n) 
    /*5*/  + pow2(gDA(m,n)) + 4 * pow2(gDconf(m,n)))) / gA(m,n))) / (TINY_Real + exp(2 
    /*3*/  * gconf(m,n))) - (3 * gAlp(m,n) * p(m,n) * pow2(p_r(m,n))) / (TINY_Real 
    /*2*/  + exp(2 * gconf(m,n)) * gA(m,n) * Power(Lt2(m,n),2))) / (TINY_Real + Lt(m,n))
    /*0*/  + q_rr(m,n);
}
Real BimetricEvolve::eq_fBet_rr( Int m, Int n )
{
    return (((((-2 * fDAlp(m,n)) / fA(m,n) + (2 * (fDA(m,n) + 2 * fDconf(m,n)) 
    /*5*/  * fAlp(m,n)) / fA(m,n)) * p_r(m,n)) / (TINY_Real + exp(2 * fconf(m,n))) 
    /*2*/  - (fAlp(m,n) * p_rr(m,n)) / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n))) 
    /*1*/  / (TINY_Real + Lt2(m,n)) + (p(m,n) * ((2 * fDA(m,n) * fDAlp(m,n) + 4 
    /*4*/  * fDconf(m,n) * fDAlp(m,n) - fDAlp_r(m,n)) / fA(m,n) + (fAlp(m,n) 
    /*4*/  * (fDA_r(m,n) - 4 * fDA(m,n) * fDconf(m,n) + 2 * fDconf_r(m,n) 
    /*5*/  - pow2(fDA(m,n)) - 4 * pow2(fDconf(m,n)))) / fA(m,n))) / (TINY_Real + exp(2 
    /*3*/  * fconf(m,n))) + (3 * fAlp(m,n) * p(m,n) * pow2(p_r(m,n))) / (TINY_Real 
    /*2*/  + exp(2 * fconf(m,n)) * fA(m,n) * Power(Lt2(m,n),2))) / (TINY_Real + Lt(m,n))
    /*0*/  + q_rr(m,n);
}
Real BimetricEvolve::eq_gBetr( Int m, Int n )
{
    return (gAlp(m,n) * pr(m,n)) / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n) 
    /*1*/  * Lt(m,n)) + qr(m,n);
}
Real BimetricEvolve::eq_fBetr( Int m, Int n )
{
    return -((fAlp(m,n) * pr(m,n)) / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n) 
    /*2*/  * Lt(m,n))) + qr(m,n);
}
Real BimetricEvolve::eq_gBetr_r( Int m, Int n )
{
    return ((exp(-2 * gconf(m,n)) * (gDAlp(m,n) / gA(m,n) - ((gDA(m,n) + 2 
    /*5*/  * gDconf(m,n)) * gAlp(m,n)) / gA(m,n)) - (exp(-2 * gconf(m,n)) * gAlp(m,n) 
    /*3*/  * p(m,n) * p_r(m,n)) / (gA(m,n) * Lt2(m,n))) * pr(m,n)) / (TINY_Real 
    /*1*/  + Lt(m,n)) + (gAlp(m,n) * pr_r(m,n)) / (TINY_Real + exp(2 * gconf(m,n)) 
    /*1*/  * gA(m,n) * Lt(m,n)) + qr_r(m,n);
}
Real BimetricEvolve::eq_fBetr_r( Int m, Int n )
{
    return ((exp(-2 * fconf(m,n)) * (-(fDAlp(m,n) / fA(m,n)) + ((fDA(m,n) + 2 
    /*5*/  * fDconf(m,n)) * fAlp(m,n)) / fA(m,n)) + (exp(-2 * fconf(m,n)) * fAlp(m,n) 
    /*3*/  * p(m,n) * p_r(m,n)) / (fA(m,n) * Lt2(m,n))) * pr(m,n)) / (TINY_Real 
    /*1*/  + Lt(m,n)) - (fAlp(m,n) * pr_r(m,n)) / (TINY_Real + exp(2 * fconf(m,n)) 
    /*1*/  * fA(m,n) * Lt(m,n)) + qr_r(m,n);
}
Real BimetricEvolve::eq_gBetr_rr( Int m, Int n )
{
    return ((exp(-2 * gconf(m,n)) * ((-2 * gDA(m,n) * gDAlp(m,n) - 4 * gDconf(m,n) 
    /*4*/  * gDAlp(m,n) + gDAlp_r(m,n)) / gA(m,n) + (gAlp(m,n) * (-gDA_r(m,n) + 4 
    /*5*/  * gDA(m,n) * gDconf(m,n) - 2 * gDconf_r(m,n) + pow2(gDA(m,n)) + 4 
    /*5*/  * pow2(gDconf(m,n)))) / gA(m,n)) - (3 * exp(-2 * gconf(m,n)) * gAlp(m,n) 
    /*3*/  * pow2(p_r(m,n))) / (gA(m,n) * pow2(Lt2(m,n))) + (exp(-2 * gconf(m,n)) * ((-2
    /*5*/  * gDAlp(m,n) * p(m,n) * p_r(m,n)) / gA(m,n) + (gAlp(m,n) * (2 * gDA(m,n) 
    /*6*/  * p(m,n) * p_r(m,n) + 4 * gDconf(m,n) * p(m,n) * p_r(m,n) - p(m,n) 
    /*6*/  * p_rr(m,n) + 2 * pow2(p_r(m,n)))) / gA(m,n))) / Lt2(m,n)) * pr(m,n)) 
    /*0*/  / (TINY_Real + Lt(m,n)) + ((exp(-2 * gconf(m,n)) * ((2 * gDAlp(m,n)) 
    /*3*/  / gA(m,n) - (2 * (gDA(m,n) + 2 * gDconf(m,n)) * gAlp(m,n)) / gA(m,n)) - (2 
    /*3*/  * exp(-2 * gconf(m,n)) * gAlp(m,n) * p(m,n) * p_r(m,n)) / (gA(m,n) 
    /*3*/  * Lt2(m,n))) * pr_r(m,n)) / (TINY_Real + Lt(m,n)) + (gAlp(m,n) * pr_rr(m,n))
    /*0*/  / (TINY_Real + exp(2 * gconf(m,n)) * gA(m,n) * Lt(m,n)) + qr_rr(m,n);
}
Real BimetricEvolve::eq_fBetr_rr( Int m, Int n )
{
    return ((exp(-2 * fconf(m,n)) * ((2 * fDA(m,n) * fDAlp(m,n) + 4 * fDconf(m,n) 
    /*4*/  * fDAlp(m,n) - fDAlp_r(m,n)) / fA(m,n) + (fAlp(m,n) * (fDA_r(m,n) - 4 
    /*5*/  * fDA(m,n) * fDconf(m,n) + 2 * fDconf_r(m,n) - pow2(fDA(m,n)) - 4 
    /*5*/  * pow2(fDconf(m,n)))) / fA(m,n)) + (exp(-2 * fconf(m,n)) * ((2 * fDAlp(m,n) 
    /*5*/  * p(m,n) * p_r(m,n)) / fA(m,n) + (fAlp(m,n) * (-2 * fDA(m,n) * p(m,n) 
    /*6*/  * p_r(m,n) - 4 * fDconf(m,n) * p(m,n) * p_r(m,n) + p(m,n) * p_rr(m,n) - 2 
    /*6*/  * pow2(p_r(m,n)))) / fA(m,n))) / Lt2(m,n) + (3 * exp(-2 * fconf(m,n)) 
    /*3*/  * fAlp(m,n) * pow2(p_r(m,n))) / (fA(m,n) * pow2(Lt2(m,n)))) * pr(m,n)) 
    /*0*/  / (TINY_Real + Lt(m,n)) + ((exp(-2 * fconf(m,n)) * ((-2 * fDAlp(m,n)) 
    /*3*/  / fA(m,n) + (2 * (fDA(m,n) + 2 * fDconf(m,n)) * fAlp(m,n)) / fA(m,n)) + (2 
    /*3*/  * exp(-2 * fconf(m,n)) * fAlp(m,n) * p(m,n) * p_r(m,n)) / (fA(m,n) 
    /*3*/  * Lt2(m,n))) * pr_r(m,n)) / (TINY_Real + Lt(m,n)) - (fAlp(m,n) * pr_rr(m,n))
    /*0*/  / (TINY_Real + exp(2 * fconf(m,n)) * fA(m,n) * Lt(m,n)) + qr_rr(m,n);
}
Real BimetricEvolve::eq_base_p_g( Int m, Int n )
{
    return (exp(-2 * fconf(m,n)) * (6 * gA2(m,n) * (1 + gDB(m,n) * r(m,n)) - 6 
    /*2*/  * gA1(m,n) * (1 + gDB(m,n) * r(m,n) + 3 * gDconf(m,n) * r(m,n)) + r(m,n) 
    /*2*/  * (-3 * gA1_r(m,n) + 3 * k_g * exp(4 * gconf(m,n)) * gj(m,n) * pow2(gA(m,n))
    /*3*/  + 6 * gDconf(m,n) * gtrA(m,n) + 3 * gtrA_r(m,n) + 2 * gtrK_r(m,n)))) / (3.
    /*1*/  * k_g * fA(m,n) * P_2_1(R(m,n)) * r(m,n));
}
Real BimetricEvolve::eq_base_p_f( Int m, Int n )
{
    return -(exp(4 * fconf(m,n) - 6 * gconf(m,n)) * pow2(fB(m,n)) * (6 * fA2(m,n) * (1
    /*3*/  + fDB(m,n) * r(m,n)) - 6 * fA1(m,n) * (1 + fDB(m,n) * r(m,n) + 3 
    /*3*/  * fDconf(m,n) * r(m,n)) + r(m,n) * (-3 * fA1_r(m,n) + 3 * k_f * exp(4 
    /*4*/  * fconf(m,n)) * fj(m,n) * pow2(fA(m,n)) + 6 * fDconf(m,n) * ftrA(m,n) + 3 
    /*3*/  * ftrA_r(m,n) + 2 * ftrK_r(m,n)))) / (3. * k_f * gA(m,n) * pow2(gB(m,n)) 
    /*1*/  * P_2_1(R(m,n)) * r(m,n));
}
Real BimetricEvolve::eq_invr_p_g( Int m, Int n )
{
    return -((exp(-2 * fconf(m,n) + 4 * gconf(m,n)) * gjb_u(m,n) * pow2(gA(m,n))) 
    /*1*/  / (fA(m,n) * P_2_1(R(m,n))));
}
Real BimetricEvolve::eq_invr_p_f( Int m, Int n )
{
    return (exp(4 * fconf(m,n) - 2 * gconf(m,n)) * fjb_u(m,n) * pow2(fA(m,n)) 
    /*1*/  * pow2(R(m,n))) / (gA(m,n) * P_2_1(R(m,n)));
}
Real BimetricEvolve::eq_p_t( Int m, Int n )
{
    return p_r(m,n) * q(m,n) + (2 * exp(-2 * gconf(m,n)) * (fB(m,n) * gAlp(m,n) 
    /*2*/  * gA(m,n) * P_1_1(R(m,n)) * (1 + fDB(m,n) * r(m,n) + 2 * fDconf(m,n) 
    /*3*/  * r(m,n)) - fAlp(m,n) * fA(m,n) * gB(m,n) * P_1_2(R(m,n)) * (1 + gDB(m,n) 
    /*3*/  * r(m,n) + 2 * gDconf(m,n) * r(m,n)))) / (fA(m,n) * gA(m,n) * gB(m,n) 
    /*1*/  * P_2_1(R(m,n)) * r(m,n)) + (exp(-2 * (fconf(m,n) + gconf(m,n))) * Lt(m,n) 
    /*1*/  * (gA(m,n) * (exp(2 * gconf(m,n)) * fDAlp(m,n) * gB(m,n) * P_2_1(R(m,n)) 
    /*3*/  * r(m,n) + 2 * exp(2 * fconf(m,n)) * fAlp(m,n) * fB(m,n) * P_1_2(R(m,n)) * (1
    /*4*/  + fDB(m,n) * r(m,n) + 2 * fDconf(m,n) * r(m,n))) - exp(2 * fconf(m,n)) 
    /*2*/  * fA(m,n) * gB(m,n) * (gDAlp(m,n) * P_2_1(R(m,n)) * r(m,n) + 2 * gAlp(m,n) 
    /*3*/  * P_1_1(R(m,n)) * (1 + gDB(m,n) * r(m,n) + 2 * gDconf(m,n) * r(m,n))))) 
    /*0*/  / (fA(m,n) * gA(m,n) * gB(m,n) * P_2_1(R(m,n)) * r(m,n)) + (p(m,n) 
    /*1*/  * (exp(-fconf(m,n)) * fAlp(m,n) * (3 * fA1(m,n) - ftrA(m,n) + exp(fconf(m,n))
    /*3*/  * ftrK(m,n) + 2 * (1 - P_1_1(R(m,n)) / P_2_1(R(m,n))) * (3 * fA2(m,n) 
    /*4*/  - ftrA(m,n) + exp(fconf(m,n)) * ftrK(m,n))) + (exp(-gconf(m,n)) * gAlp(m,n) 
    /*3*/  * (6 * gA2(m,n) * P_1_1(R(m,n)) + 3 * gA1(m,n) * P_2_1(R(m,n)) - (2 
    /*5*/  * P_1_1(R(m,n)) + P_2_1(R(m,n))) * (gtrA(m,n) - exp(gconf(m,n)) 
    /*5*/  * gtrK(m,n)))) / P_2_1(R(m,n)))) / 3.;
}
Real BimetricEvolve::eq_p_r( Int m, Int n )
{
    return (exp(-2 * fconf(m,n)) * (-3 * exp(4 * gconf(m,n)) * pow2(gA(m,n)) 
    /*2*/  * (fDA(m,n) * gj(m,n) * P_2_1(R(m,n)) + 2 * fDconf(m,n) * gj(m,n) 
    /*3*/  * P_2_1(R(m,n)) - 2 * gDA(m,n) * gj(m,n) * P_2_1(R(m,n)) - 4 * gDconf(m,n) 
    /*3*/  * gj(m,n) * P_2_1(R(m,n)) - gj_r(m,n) * P_2_1(R(m,n)) + 2 * gj(m,n) 
    /*3*/  * P_1_2(R(m,n)) * R_r(m,n)) + (-6 * gA2(m,n) * (P_2_1(R(m,n)) * (1 
    /*5*/  - gDB_r(m,n) * pow2(r(m,n)) + fDA(m,n) * r(m,n) * (1 + gDB(m,n) * r(m,n)) + 2
    /*5*/  * fDconf(m,n) * r(m,n) * (1 + gDB(m,n) * r(m,n))) + 2 * P_1_2(R(m,n)) 
    /*4*/  * r(m,n) * (1 + gDB(m,n) * r(m,n)) * R_r(m,n)) + 6 * gA1(m,n) 
    /*3*/  * (P_2_1(R(m,n)) * (1 - gDB_r(m,n) * pow2(r(m,n)) - 3 * gDconf_r(m,n) 
    /*5*/  * pow2(r(m,n)) + fDA(m,n) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 3 * gDconf(m,n)
    /*6*/  * r(m,n)) + 2 * fDconf(m,n) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 3 
    /*6*/  * gDconf(m,n) * r(m,n))) + 2 * P_1_2(R(m,n)) * r(m,n) * (1 + gDB(m,n) 
    /*5*/  * r(m,n) + 3 * gDconf(m,n) * r(m,n)) * R_r(m,n)) - r(m,n) * (6 * gA1_r(m,n) 
    /*4*/  * P_2_1(R(m,n)) - 6 * gA2_r(m,n) * P_2_1(R(m,n)) + 3 * gA1_rr(m,n) 
    /*4*/  * P_2_1(R(m,n)) * r(m,n) + 6 * gA1_r(m,n) * gDB(m,n) * P_2_1(R(m,n)) * r(m,n)
    /*4*/  - 6 * gA2_r(m,n) * gDB(m,n) * P_2_1(R(m,n)) * r(m,n) + 18 * gA1_r(m,n) 
    /*4*/  * gDconf(m,n) * P_2_1(R(m,n)) * r(m,n) - 6 * gA1_r(m,n) * P_1_2(R(m,n)) 
    /*4*/  * r(m,n) * R_r(m,n) - 6 * gDconf_r(m,n) * P_2_1(R(m,n)) * r(m,n) * gtrA(m,n)
    /*4*/  + 12 * gDconf(m,n) * P_1_2(R(m,n)) * r(m,n) * R_r(m,n) * gtrA(m,n) - 6 
    /*4*/  * gDconf(m,n) * P_2_1(R(m,n)) * r(m,n) * gtrA_r(m,n) + 6 * P_1_2(R(m,n)) 
    /*4*/  * r(m,n) * R_r(m,n) * gtrA_r(m,n) - 3 * P_2_1(R(m,n)) * r(m,n) * gtrA_rr(m,n)
    /*4*/  + 4 * P_1_2(R(m,n)) * r(m,n) * R_r(m,n) * gtrK_r(m,n) + fDA(m,n) 
    /*4*/  * P_2_1(R(m,n)) * r(m,n) * (-3 * gA1_r(m,n) + 6 * gDconf(m,n) * gtrA(m,n) + 3
    /*5*/  * gtrA_r(m,n) + 2 * gtrK_r(m,n)) + 2 * fDconf(m,n) * P_2_1(R(m,n)) * r(m,n)
    /*4*/  * (-3 * gA1_r(m,n) + 6 * gDconf(m,n) * gtrA(m,n) + 3 * gtrA_r(m,n) + 2 
    /*5*/  * gtrK_r(m,n)) - 2 * P_2_1(R(m,n)) * r(m,n) * gtrK_rr(m,n))) / (k_g 
    /*3*/  * pow2(r(m,n))))) / (3. * fA(m,n) * pow2(P_2_1(R(m,n))));
}
Real BimetricEvolve::eq_gHorz( Int m, Int n )
{
    return 2 * (gDB(m,n) + 2 * gDconf(m,n) - exp(2 * gconf(m,n)) * gK2(m,n) * gA(m,n)
    /*1*/  + 1 / r(m,n));
}
Real BimetricEvolve::eq_fHorz( Int m, Int n )
{
    return 2 * (fDB(m,n) + 2 * fDconf(m,n) - exp(2 * fconf(m,n)) * fK2(m,n) * fA(m,n)
    /*1*/  + 1 / r(m,n));
}
