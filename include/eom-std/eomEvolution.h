/** @file  eomEvolution.h
 *  @brief The evolution equations for metrics.
 *  @image html evolution.png
 */

Real BimetricEvolve::eq_gA_t( Int m, Int n )
{
    return (gAlp(m,n) * gDAlp(m,n) * p(m,n)) / Lt(m,n) - (gAlp(m,n) * pow2(p(m,n)) 
    /*1*/  * p_r(m,n)) / (Lt(m,n) * (1 + pow2(p(m,n)))) + gAlp(m,n) * ((gA(m,n) 
    /*2*/  * (-gK(m,n) - 2 * gKD(m,n))) / 3. + p_r(m,n) / Lt(m,n)) + gA(m,n) * gDA(m,n)
    /*0*/  * q(m,n) + gA(m,n) * q_r(m,n);
}

Real BimetricEvolve::eq_gB_t( Int m, Int n )
{
    return eq_qr(m,n) * gB(m,n) + (gAlp(m,n) * gB(m,n) * (-gK(m,n) + gKD(m,n))) / 3. 
    /*0*/  + (eq_pr(m,n) * gAlp(m,n) * gB(m,n)) / (gA(m,n) * Lt(m,n)) + (gAlp(m,n) 
    /*1*/  * gB(m,n) * gDB(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n)) + gB(m,n) * gDB(m,n) 
    /*0*/  * q(m,n);
}

Real BimetricEvolve::eq_invr_gK1_t( Int m, Int n )
{
    return (2 * gAlp(m,n) * (gDA(m,n) - 2 * gDB(m,n))) / pow2(gA(m,n));
}

Real BimetricEvolve::eq_base_gK1_t( Int m, Int n )
{
    return k_g * (-g_JK(m,n) / 3. - (2 * g_JKD(m,n)) / 3.) + (gAlp(m,n) * (2 
    /*2*/  * gKD_r(m,n) + gK_r(m,n)) * p(m,n)) / (3. * gA(m,n) * Lt(m,n)) + gAlp(m,n) 
    /*0*/  * ((gK(m,n) * (gK(m,n) + 2 * gKD(m,n))) / 3. + (-gDAlp_r(m,n) + gDA(m,n) 
    /*2*/  * (gDAlp(m,n) + 2 * gDB(m,n)) - 2 * gDB_r(m,n) - pow2(gDAlp(m,n)) - 2 
    /*2*/  * pow2(gDB(m,n))) / pow2(gA(m,n))) + ((2 * gKD_r(m,n) + gK_r(m,n)) * q(m,n))
    /*0*/  / 3.;
}

Real BimetricEvolve::eq_invr_gK2_t( Int m, Int n )
{
    return gAlp(m,n) * (gSig(m,n) / (gA(m,n) * gB(m,n)) + (gDA(m,n) - gDAlp(m,n) - 4 
    /*2*/  * gDB(m,n) + gSig(m,n)) / pow2(gA(m,n)));
}

Real BimetricEvolve::eq_base_gK2_t( Int m, Int n )
{
    return k_g * (-g_JK(m,n) / 3. + g_JKD(m,n) / 3.) + (gAlp(m,n) * (-gKD_r(m,n) 
    /*2*/  + gK_r(m,n)) * p(m,n)) / (3. * gA(m,n) * Lt(m,n)) + gAlp(m,n) * ((gK(m,n) 
    /*2*/  * (gK(m,n) - gKD(m,n))) / 3. + (gDA(m,n) * gDB(m,n) - gDAlp(m,n) * gDB(m,n) 
    /*2*/  - gDB_r(m,n) - 2 * pow2(gDB(m,n))) / pow2(gA(m,n))) + ((-gKD_r(m,n) 
    /*2*/  + gK_r(m,n)) * q(m,n)) / 3.;
}

Real BimetricEvolve::eq_gDA_t( Int m, Int n )
{
    return (gAlp(m,n) * (-(gDAlp(m,n) * (gK(m,n) + 2 * gKD(m,n))) - 2 * gKD_r(m,n) 
    /*2*/  - gK_r(m,n))) / 3. + (gAlp(m,n) * p(m,n) * ((-(gDA(m,n) * gDAlp(m,n)) 
    /*3*/  + gDAlp_r(m,n) + pow2(gDAlp(m,n))) / Lt(m,n) - (3 * pow2(p_r(m,n))) 
    /*2*/  / (Lt(m,n) * pow2(1 + pow2(p(m,n)))))) / gA(m,n) + (gAlp(m,n) * (-((gDA(m,n)
    /*4*/  - 2 * gDAlp(m,n)) * p_r(m,n)) + p_rr(m,n))) / (gA(m,n) * Lt(m,n) * (1 
    /*2*/  + pow2(p(m,n)))) + gDA_r(m,n) * q(m,n) + gDA(m,n) * q_r(m,n) + q_rr(m,n);
}

Real BimetricEvolve::eq_gDB_t( Int m, Int n )
{
    return eq_qr_r(m,n) + (eq_pr(m,n) * gAlp(m,n) * (-gDA(m,n) + gDAlp(m,n))) 
    /*0*/  / (gA(m,n) * Lt(m,n)) + (gAlp(m,n) * (-(gDA(m,n) * gDB(m,n)) + gDAlp(m,n) 
    /*2*/  * gDB(m,n) + gDB_r(m,n)) * p(m,n)) / (gA(m,n) * Lt(m,n)) - (gAlp(m,n) 
    /*1*/  * p(m,n) * pow2(eq_pr(m,n))) / (gA(m,n) * Lt(m,n) * (1 + pow2(p(m,n)))) 
    /*0*/  + gAlp(m,n) * ((gDAlp(m,n) * (-gK(m,n) + gKD(m,n)) + gKD_r(m,n) - gK_r(m,n))
    /*1*/  / 3. + (eq_pr_r(m,n) + gDB(m,n) * p_r(m,n)) / (gA(m,n) * Lt(m,n) * (1 
    /*3*/  + pow2(p(m,n))))) + gDB_r(m,n) * q(m,n) + gDB(m,n) * q_r(m,n);
}

Real BimetricEvolve::eq_invr_gK_t( Int m, Int n )
{
    return gAlp(m,n) * ((2 * gSig(m,n)) / (gA(m,n) * gB(m,n)) + (2 * (2 * gDA(m,n) 
    /*3*/  - gDAlp(m,n) - 6 * gDB(m,n) + gSig(m,n))) / pow2(gA(m,n)));
}

Real BimetricEvolve::eq_base_gK_t( Int m, Int n )
{
    return -(k_g * g_JK(m,n)) + (gAlp(m,n) * gK_r(m,n) * p(m,n)) / (gA(m,n) * Lt(m,n))
    /*0*/  + gAlp(m,n) * ((-gDAlp_r(m,n) - 2 * gDAlp(m,n) * gDB(m,n) + gDA(m,n) 
    /*2*/  * (gDAlp(m,n) + 4 * gDB(m,n)) - 4 * gDB_r(m,n) - pow2(gDAlp(m,n)) - 6 
    /*2*/  * pow2(gDB(m,n))) / pow2(gA(m,n)) + pow2(gK(m,n))) + gK_r(m,n) * q(m,n);
}

Real BimetricEvolve::eq_invr_gKD_t( Int m, Int n )
{
    return gAlp(m,n) * (-(gSig(m,n) / (gA(m,n) * gB(m,n))) + (gDA(m,n) + gDAlp(m,n) 
    /*2*/  - gSig(m,n)) / pow2(gA(m,n)));
}

Real BimetricEvolve::eq_base_gKD_t( Int m, Int n )
{
    return -(k_g * g_JKD(m,n)) + (gAlp(m,n) * gKD_r(m,n) * p(m,n)) / (gA(m,n) 
    /*1*/  * Lt(m,n)) + gAlp(m,n) * (gK(m,n) * gKD(m,n) + (-gDAlp_r(m,n) + gDAlp(m,n) 
    /*2*/  * gDB(m,n) + gDA(m,n) * (gDAlp(m,n) + gDB(m,n)) - gDB_r(m,n) 
    /*2*/  - pow2(gDAlp(m,n))) / pow2(gA(m,n))) + gKD_r(m,n) * q(m,n);
}

Real BimetricEvolve::eq_gSig_t( Int m, Int n )
{
    return (eq_qr_r(m,n) * gA(m,n)) / gB(m,n) + (eq_qr(m,n) * gA(m,n) * (gDA(m,n) 
    /*2*/  - gDB(m,n))) / gB(m,n) + (eq_pr(m,n) * gAlp(m,n) * (gDAlp(m,n) - gDB(m,n))) 
    /*0*/  / (gB(m,n) * Lt(m,n)) - (gAlp(m,n) * p(m,n) * pow2(eq_pr(m,n))) / (gB(m,n) 
    /*1*/  * Lt(m,n) * (1 + pow2(p(m,n)))) + gAlp(m,n) * ((gA(m,n) * (3 * gDB(m,n) 
    /*3*/  * gKD(m,n) + gKD_r(m,n) - gK_r(m,n))) / (3. * gB(m,n)) + eq_pr_r(m,n) 
    /*1*/  / (gB(m,n) * Lt(m,n) * (1 + pow2(p(m,n))))) - (k_g * gAlp(m,n) * pfS(m,n)) 
    /*0*/  / (2. * pow3(gB(m,n))) + (k_g * fA(m,n) * gA(m,n) * gAlp(m,n) * p(m,n) 
    /*1*/  * P_2_1(R(m,n))) / (2. * gB(m,n));
}

Real BimetricEvolve::eq_fA_t( Int m, Int n )
{
    return -((fAlp(m,n) * fDAlp(m,n) * p(m,n)) / Lt(m,n)) + (fAlp(m,n) * pow2(p(m,n))
    /*1*/  * p_r(m,n)) / (Lt(m,n) * (1 + pow2(p(m,n)))) + fAlp(m,n) * ((fA(m,n) 
    /*2*/  * (-fK(m,n) - 2 * fKD(m,n))) / 3. - p_r(m,n) / Lt(m,n)) + fA(m,n) * fDA(m,n)
    /*0*/  * q(m,n) + fA(m,n) * q_r(m,n);
}

Real BimetricEvolve::eq_fB_t( Int m, Int n )
{
    return eq_qr(m,n) * fB(m,n) + (fAlp(m,n) * fB(m,n) * (-fK(m,n) + fKD(m,n))) / 3. 
    /*0*/  - (eq_pr(m,n) * fAlp(m,n) * fB(m,n)) / (fA(m,n) * Lt(m,n)) - (fAlp(m,n) 
    /*1*/  * fB(m,n) * fDB(m,n) * p(m,n)) / (fA(m,n) * Lt(m,n)) + fB(m,n) * fDB(m,n) 
    /*0*/  * q(m,n);
}

Real BimetricEvolve::eq_invr_fK1_t( Int m, Int n )
{
    return (2 * fAlp(m,n) * (fDA(m,n) - 2 * fDB(m,n))) / pow2(fA(m,n));
}

Real BimetricEvolve::eq_base_fK1_t( Int m, Int n )
{
    return k_f * (-f_JK(m,n) / 3. - (2 * f_JKD(m,n)) / 3.) - (fAlp(m,n) * (2 
    /*2*/  * fKD_r(m,n) + fK_r(m,n)) * p(m,n)) / (3. * fA(m,n) * Lt(m,n)) + fAlp(m,n) 
    /*0*/  * ((fK(m,n) * (fK(m,n) + 2 * fKD(m,n))) / 3. + (-fDAlp_r(m,n) + fDA(m,n) 
    /*2*/  * (fDAlp(m,n) + 2 * fDB(m,n)) - 2 * fDB_r(m,n) - pow2(fDAlp(m,n)) - 2 
    /*2*/  * pow2(fDB(m,n))) / pow2(fA(m,n))) + ((2 * fKD_r(m,n) + fK_r(m,n)) * q(m,n))
    /*0*/  / 3.;
}

Real BimetricEvolve::eq_invr_fK2_t( Int m, Int n )
{
    return fAlp(m,n) * (fSig(m,n) / (fA(m,n) * fB(m,n)) + (fDA(m,n) - fDAlp(m,n) - 4 
    /*2*/  * fDB(m,n) + fSig(m,n)) / pow2(fA(m,n)));
}

Real BimetricEvolve::eq_base_fK2_t( Int m, Int n )
{
    return k_f * (-f_JK(m,n) / 3. + f_JKD(m,n) / 3.) + (fAlp(m,n) * (fKD_r(m,n) 
    /*2*/  - fK_r(m,n)) * p(m,n)) / (3. * fA(m,n) * Lt(m,n)) + fAlp(m,n) * ((fK(m,n) 
    /*2*/  * (fK(m,n) - fKD(m,n))) / 3. + (fDA(m,n) * fDB(m,n) - fDAlp(m,n) * fDB(m,n) 
    /*2*/  - fDB_r(m,n) - 2 * pow2(fDB(m,n))) / pow2(fA(m,n))) + ((-fKD_r(m,n) 
    /*2*/  + fK_r(m,n)) * q(m,n)) / 3.;
}

Real BimetricEvolve::eq_fDA_t( Int m, Int n )
{
    return -(fAlp(m,n) * (fDAlp(m,n) * (fK(m,n) + 2 * fKD(m,n)) + 2 * fKD_r(m,n) 
    /*2*/  + fK_r(m,n))) / 3. + p(m,n) * ((fAlp(m,n) * ((fDA(m,n) - fDAlp(m,n)) 
    /*3*/  * fDAlp(m,n) - fDAlp_r(m,n))) / (fA(m,n) * Lt(m,n)) + (3 * fAlp(m,n) 
    /*2*/  * pow2(p_r(m,n))) / (fA(m,n) * Lt(m,n) * pow2(1 + pow2(p(m,n))))) 
    /*0*/  + (fAlp(m,n) * ((fDA(m,n) - 2 * fDAlp(m,n)) * p_r(m,n) - p_rr(m,n))) 
    /*0*/  / (fA(m,n) * Lt(m,n) * (1 + pow2(p(m,n)))) + fDA_r(m,n) * q(m,n) + fDA(m,n) 
    /*0*/  * q_r(m,n) + q_rr(m,n);
}

Real BimetricEvolve::eq_fDB_t( Int m, Int n )
{
    return eq_qr_r(m,n) + (eq_pr(m,n) * fAlp(m,n) * (fDA(m,n) - fDAlp(m,n))) 
    /*0*/  / (fA(m,n) * Lt(m,n)) - (fAlp(m,n) * (-(fDA(m,n) * fDB(m,n)) + fDAlp(m,n) 
    /*2*/  * fDB(m,n) + fDB_r(m,n)) * p(m,n)) / (fA(m,n) * Lt(m,n)) + (fAlp(m,n) 
    /*1*/  * p(m,n) * pow2(eq_pr(m,n))) / (fA(m,n) * Lt(m,n) * (1 + pow2(p(m,n)))) 
    /*0*/  + fAlp(m,n) * ((fDAlp(m,n) * (-fK(m,n) + fKD(m,n)) + fKD_r(m,n) - fK_r(m,n))
    /*1*/  / 3. - (eq_pr_r(m,n) + fDB(m,n) * p_r(m,n)) / (fA(m,n) * Lt(m,n) * (1 
    /*3*/  + pow2(p(m,n))))) + fDB_r(m,n) * q(m,n) + fDB(m,n) * q_r(m,n);
}

Real BimetricEvolve::eq_invr_fK_t( Int m, Int n )
{
    return fAlp(m,n) * ((2 * fSig(m,n)) / (fA(m,n) * fB(m,n)) + (2 * (2 * fDA(m,n) 
    /*3*/  - fDAlp(m,n) - 6 * fDB(m,n) + fSig(m,n))) / pow2(fA(m,n)));
}

Real BimetricEvolve::eq_base_fK_t( Int m, Int n )
{
    return -(k_f * f_JK(m,n)) - (fAlp(m,n) * fK_r(m,n) * p(m,n)) / (fA(m,n) * Lt(m,n))
    /*0*/  + fAlp(m,n) * ((-fDAlp_r(m,n) - 2 * fDAlp(m,n) * fDB(m,n) + fDA(m,n) 
    /*2*/  * (fDAlp(m,n) + 4 * fDB(m,n)) - 4 * fDB_r(m,n) - pow2(fDAlp(m,n)) - 6 
    /*2*/  * pow2(fDB(m,n))) / pow2(fA(m,n)) + pow2(fK(m,n))) + fK_r(m,n) * q(m,n);
}

Real BimetricEvolve::eq_invr_fKD_t( Int m, Int n )
{
    return fAlp(m,n) * (-(fSig(m,n) / (fA(m,n) * fB(m,n))) + (fDA(m,n) + fDAlp(m,n) 
    /*2*/  - fSig(m,n)) / pow2(fA(m,n)));
}

Real BimetricEvolve::eq_base_fKD_t( Int m, Int n )
{
    return -(k_f * f_JKD(m,n)) - (fAlp(m,n) * fKD_r(m,n) * p(m,n)) / (fA(m,n) 
    /*1*/  * Lt(m,n)) + fAlp(m,n) * (fK(m,n) * fKD(m,n) + (-fDAlp_r(m,n) + fDAlp(m,n) 
    /*2*/  * fDB(m,n) + fDA(m,n) * (fDAlp(m,n) + fDB(m,n)) - fDB_r(m,n) 
    /*2*/  - pow2(fDAlp(m,n))) / pow2(fA(m,n))) + fKD_r(m,n) * q(m,n);
}

Real BimetricEvolve::eq_fSig_t( Int m, Int n )
{
    return (eq_qr_r(m,n) * fA(m,n)) / fB(m,n) + (eq_qr(m,n) * fA(m,n) * (fDA(m,n) 
    /*2*/  - fDB(m,n))) / fB(m,n) + (eq_pr(m,n) * fAlp(m,n) * (-fDAlp(m,n) + fDB(m,n)))
    /*0*/  / (fB(m,n) * Lt(m,n)) + (fAlp(m,n) * p(m,n) * pow2(eq_pr(m,n))) / (fB(m,n)
    /*1*/  * Lt(m,n) * (1 + pow2(p(m,n)))) + fAlp(m,n) * ((fA(m,n) * (3 * fDB(m,n) 
    /*3*/  * fKD(m,n) + fKD_r(m,n) - fK_r(m,n))) / (3. * fB(m,n)) - eq_pr_r(m,n) 
    /*1*/  / (fB(m,n) * Lt(m,n) * (1 + pow2(p(m,n))))) - (k_f * fA(m,n) * fAlp(m,n) 
    /*1*/  * gA(m,n) * p(m,n) * P_2_1(R(m,n))) / (2. * fB(m,n) * pow2(R(m,n)));
}

