/** @file  eomConstraints.h
 *  @brief The constraint equations (without sources).
 *  @image html constraints.png
 */

Real BimetricEvolve::eq_gC_rho( Int m, Int n )
{
    return (2 * gDA(m,n) * gDB(m,n) - 2 * gDB_r(m,n) - 3 * pow2(gDB(m,n))) 
    /*0*/  / pow2(gA(m,n)) + (pow2(gK(m,n)) - pow2(gKD(m,n))) / 3. + (gSig(m,n) 
    /*1*/  / (gA(m,n) * gB(m,n)) + (2 * gDA(m,n) - 6 * gDB(m,n) + gSig(m,n)) 
    /*1*/  / pow2(gA(m,n))) / r(m,n);
}

Real BimetricEvolve::eq_fC_rho( Int m, Int n )
{
    return (2 * fDA(m,n) * fDB(m,n) - 2 * fDB_r(m,n) - 3 * pow2(fDB(m,n))) 
    /*0*/  / pow2(fA(m,n)) + (pow2(fK(m,n)) - pow2(fKD(m,n))) / 3. + (fSig(m,n) 
    /*1*/  / (fA(m,n) * fB(m,n)) + (2 * fDA(m,n) - 6 * fDB(m,n) + fSig(m,n)) 
    /*1*/  / pow2(fA(m,n))) / r(m,n);
}

Real BimetricEvolve::eq_gC_j( Int m, Int n )
{
    return (2 * (3 * gDB(m,n) * gKD(m,n) + gKD_r(m,n) - gK_r(m,n))) / 3. + (2 
    /*1*/  * gKD(m,n)) / r(m,n);
}

Real BimetricEvolve::eq_fC_j( Int m, Int n )
{
    return (2 * (3 * fDB(m,n) * fKD(m,n) + fKD_r(m,n) - fK_r(m,n))) / 3. + (2 
    /*1*/  * fKD(m,n)) / r(m,n);
}

Real BimetricEvolve::eq_C_bimConsLaw( Int m, Int n )
{
    return (fA(m,n) * (2 * (fKD(m,n) * P_1_1(R(m,n)) + (-gK(m,n) + gKD(m,n)) * Lt(m,n)
    /*3*/  * P_1_2(R(m,n))) + fK(m,n) * (-2 * P_1_1(R(m,n)) + 3 * P_2_1(R(m,n))))) 
    /*0*/  / 3. + (P_2_1(R(m,n)) * p_r(m,n)) / Lt(m,n) + (p(m,n) * ((2 * fA(m,n) 
    /*3*/  * P_1_2(R(m,n))) / gA(m,n) + (2 * gA(m,n) * P_1_1(R(m,n)) * R(m,n)) 
    /*2*/  / fA(m,n))) / r(m,n) + p(m,n) * ((2 * fA(m,n) * gDB(m,n) * P_1_2(R(m,n))) 
    /*1*/  / gA(m,n) + (2 * fDB(m,n) * gA(m,n) * P_1_1(R(m,n)) * R(m,n)) / fA(m,n)) 
    /*0*/  + gA(m,n) * (-((gK(m,n) + 2 * gKD(m,n)) * P_2_1(R(m,n))) / 3. - (2 
    /*2*/  * P_1_1(R(m,n)) * (gK(m,n) - gKD(m,n) + (-fK(m,n) + fKD(m,n)) * Lt(m,n) 
    /*3*/  * R(m,n))) / 3.);
}

Real BimetricEvolve::eq_C_1( Int m, Int n )
{
    return -(k_g * g_rho(m,n)) + (2 * gDA(m,n) * gDB(m,n) - 2 * gDB_r(m,n) - 3 
    /*1*/  * pow2(gDB(m,n))) / pow2(gA(m,n)) + (pow2(gK(m,n)) - pow2(gKD(m,n))) / 3. 
    /*0*/  + (gSig(m,n) / (gA(m,n) * gB(m,n)) + (2 * gDA(m,n) - 6 * gDB(m,n) 
    /*2*/  + gSig(m,n)) / pow2(gA(m,n))) / r(m,n);
}

Real BimetricEvolve::eq_C_2( Int m, Int n )
{
    return -(k_f * f_rho(m,n)) + (2 * fDA(m,n) * fDB(m,n) - 2 * fDB_r(m,n) - 3 
    /*1*/  * pow2(fDB(m,n))) / pow2(fA(m,n)) + (pow2(fK(m,n)) - pow2(fKD(m,n))) / 3. 
    /*0*/  + (fSig(m,n) / (fA(m,n) * fB(m,n)) + (2 * fDA(m,n) - 6 * fDB(m,n) 
    /*2*/  + fSig(m,n)) / pow2(fA(m,n))) / r(m,n);
}

Real BimetricEvolve::eq_C_3( Int m, Int n )
{
    return k_f * ((gA(m,n) * (3 * gDB(m,n) * gKD(m,n) + gKD_r(m,n) - gK_r(m,n))) / 3.
    /*1*/  - (k_g * pfS(m,n)) / (2. * pow2(gB(m,n)))) + (k_g * fA(m,n) * (3 * fDB(m,n)
    /*2*/  * fKD(m,n) + fKD_r(m,n) - fK_r(m,n)) * pow2(R(m,n))) / 3. + (k_f * gA(m,n)
    /*1*/  * gKD(m,n) + k_g * fA(m,n) * fKD(m,n) * pow2(R(m,n))) / r(m,n);
}

Real BimetricEvolve::eq_C_4( Int m, Int n )
{
    return (fA(m,n) * (2 * (fKD(m,n) * P_1_1(R(m,n)) + (-gK(m,n) + gKD(m,n)) * Lt(m,n)
    /*3*/  * P_1_2(R(m,n))) + fK(m,n) * (-2 * P_1_1(R(m,n)) + 3 * P_2_1(R(m,n))))) 
    /*0*/  / 3. + (P_2_1(R(m,n)) * p_r(m,n)) / Lt(m,n) + (p(m,n) * ((2 * fA(m,n) 
    /*3*/  * P_1_2(R(m,n))) / gA(m,n) + (2 * gA(m,n) * P_1_1(R(m,n)) * R(m,n)) 
    /*2*/  / fA(m,n))) / r(m,n) + p(m,n) * ((2 * fA(m,n) * gDB(m,n) * P_1_2(R(m,n))) 
    /*1*/  / gA(m,n) + (2 * fDB(m,n) * gA(m,n) * P_1_1(R(m,n)) * R(m,n)) / fA(m,n)) 
    /*0*/  + gA(m,n) * (-((gK(m,n) + 2 * gKD(m,n)) * P_2_1(R(m,n))) / 3. - (2 
    /*2*/  * P_1_1(R(m,n)) * (gK(m,n) - gKD(m,n) + (-fK(m,n) + fKD(m,n)) * Lt(m,n) 
    /*3*/  * R(m,n))) / 3.);
}
