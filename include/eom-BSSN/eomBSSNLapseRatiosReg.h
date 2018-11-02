/** @file  eomBSSNLapseRatiosReg.h
 *  @author Francesco Torsello
 *  @brief The regularized ratio of the lapses in the BSSN formulation.
 *  @version 2018-10-15T11:08:53
 *  @image html BSSNLapseRatiosReg.png
 */

Real BimetricEvolve::eq_gW( Int m, Int n )
{
    return -(k_g * (Lt2(m,n) * (-12 * exp(6 * fconf(m,n) + 4 * gconf(m,n)) 
    /*3*/  * Power(fA(m,n),4) * gA(m,n) * pow2(gB(m,n)) * pow2(p(m,n)) 
    /*3*/  * pow2(P_2_1(R(m,n))) * pow2(r(m,n)) * pow2(R(m,n)) * P_1_2(R(m,n)) - 6 
    /*3*/  * exp(2 * fconf(m,n) + 8 * gconf(m,n)) * pow2(fA(m,n)) * pow2(gB(m,n)) 
    /*3*/  * pow2(r(m,n)) * pow2(R(m,n)) * pow3(gA(m,n)) * P_2_1(R(m,n)) * (2 
    /*4*/  * P_1_1(R(m,n)) * (eq_pf_gJ22(m,n) - P_2_0(R(m,n))) + (eq_pf_gJ11(m,n) - eq_pf_grho(m,n) 
    /*5*/  - P_1_0(R(m,n)) + P_2_0(R(m,n))) * P_2_1(R(m,n)))) + Lt(m,n) * (-12 * exp(4 
    /*4*/  * fconf(m,n) + 6 * gconf(m,n)) * eq_pf_gJ22(m,n) * Lt2(m,n) * pow2(gA(m,n)) 
    /*3*/  * pow2(gB(m,n)) * pow2(r(m,n)) * pow2(R(m,n)) * pow3(fA(m,n)) * P_1_2(R(m,n))
    /*3*/  * P_2_1(R(m,n)) + exp(4 * fconf(m,n) + 6 * gconf(m,n)) * pow2(gA(m,n)) 
    /*3*/  * pow2(gB(m,n)) * pow3(fA(m,n)) * (12 * eq_pf_grho(m,n) * pow2(p(m,n)) 
    /*4*/  * pow2(r(m,n)) * pow2(R(m,n)) * P_1_2(R(m,n)) * P_2_1(R(m,n)) + 3 
    /*4*/  * pow2(r(m,n)) * pow2(R(m,n)) * P_2_1(R(m,n)) * (2 * eq_pf_grho(m,n) 
    /*5*/  * P_1_2(R(m,n)) + 2 * P_1_2(R(m,n)) * P_2_0(R(m,n)) + (4 * P_1_1(R(m,n)) 
    /*6*/  - P_2_1(R(m,n))) * P_2_1(R(m,n))))))) - 12 * exp(2 * fconf(m,n) + 4 
    /*1*/  * gconf(m,n)) * gA(m,n) * p(m,n) * p_r(m,n) * pow2(fA(m,n)) * pow2(gB(m,n)) 
    /*0*/  * pow3(R(m,n)) * P_1_2(R(m,n)) * P_2_1(R(m,n)) * r(m,n) * (1 + gDB(m,n) 
    /*1*/  * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) - k_f * (-18 * exp(2 
    /*2*/  * fconf(m,n) + 8 * gconf(m,n)) * Lt2(m,n) * pow2(fA(m,n)) * pow2(gB(m,n)) 
    /*1*/  * pow2(P_2_1(R(m,n))) * pow2(r(m,n)) * pow3(gA(m,n)) * P_1_1(R(m,n)) * R(m,n)
    /*1*/  + exp(4 * fconf(m,n) + 6 * gconf(m,n)) * Lt(m,n) * pow2(gA(m,n)) 
    /*1*/  * pow2(gB(m,n)) * pow3(fA(m,n)) * (3 * pow2(P_2_1(R(m,n))) * pow2(r(m,n)) 
    /*2*/  * (2 * P_1_1(R(m,n)) - P_2_1(R(m,n))) + 6 * eq_pf_frho(m,n) * pow2(r(m,n)) 
    /*2*/  * pow3(R(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) - 6 * pow2(r(m,n)) 
    /*2*/  * P_2_1(R(m,n)) * (P_1_2(R(m,n)) * P_2_1(R(m,n)) + P_1_1(R(m,n)) 
    /*3*/  * P_2_2(R(m,n))) * R(m,n))) - exp(2 * fconf(m,n) + 6 * gconf(m,n)) * p_r(m,n)
    /*0*/  * pow2(fA(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n)) * pow3(R(m,n)) * (-12 
    /*1*/  * fA1(m,n) * pow2(r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) - 4 * pow2(r(m,n))
    /*1*/  * P_1_1(R(m,n)) * P_2_1(R(m,n)) * ftrK(m,n)) - Lt(m,n) * (exp(4 
    /*2*/  * fconf(m,n) + 2 * gconf(m,n)) * pow2(gB(m,n)) * pow3(fA(m,n)) * (-6 * pow2(1
    /*3*/  + gDB(m,n) * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) * pow2(R(m,n))
    /*2*/  * P_1_2(R(m,n)) * P_2_1(R(m,n)) + pow2(p(m,n)) * (12 * pow2(1 + gDB(m,n) 
    /*4*/  * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) * pow3(R(m,n)) 
    /*3*/  * P_0_3(R(m,n)) * P_2_1(R(m,n)) + 12 * pow2(1 + gDB(m,n) * r(m,n) + 2 
    /*4*/  * gconf(m,n) * gDconf(m,n) * r(m,n)) * pow2(R(m,n)) * P_1_2(R(m,n)) 
    /*3*/  * P_2_1(R(m,n)))) + exp(4 * fconf(m,n) + 4 * gconf(m,n)) * gA(m,n) * p(m,n) 
    /*1*/  * pow2(gB(m,n)) * pow3(fA(m,n)) * pow3(R(m,n)) * (-12 * fA1(m,n) 
    /*2*/  * (P_0_2(R(m,n)) - 3 * P_1_2(R(m,n))) * P_2_1(R(m,n)) * r(m,n) * (1 
    /*3*/  + gDB(m,n) * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) - 4 
    /*2*/  * (P_0_2(R(m,n)) - 3 * P_1_2(R(m,n))) * P_2_1(R(m,n)) * r(m,n) * (1 
    /*3*/  + gDB(m,n) * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) * ftrK(m,n)) 
    /*1*/  + pow3(gA(m,n)) * (exp(2 * fconf(m,n) + 6 * gconf(m,n)) * fA(m,n) * fB(m,n) 
    /*2*/  * gB(m,n) * p(m,n) * pow3(R(m,n)) * (-24 * fA1(m,n) * pow2(P_1_1(R(m,n))) 
    /*3*/  * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) 
    /*3*/  - 8 * pow2(P_1_1(R(m,n))) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n)
    /*4*/  * fDconf(m,n) * r(m,n)) * ftrK(m,n)) + exp(2 * fconf(m,n) + 8 * gconf(m,n))
    /*2*/  * pow2(fA(m,n)) * pow2(gB(m,n)) * pow2(p(m,n)) * pow3(R(m,n)) * (-12 
    /*3*/  * fA1(m,n) * pow2(r(m,n)) * P_1_1(R(m,n)) * (2 * gA2(m,n) * P_1_1(R(m,n)) 
    /*4*/  + gA1(m,n) * P_2_1(R(m,n))) - 4 * pow2(r(m,n)) * P_1_1(R(m,n)) * (2 
    /*4*/  * gA2(m,n) * P_1_1(R(m,n)) + gA1(m,n) * P_2_1(R(m,n))) * ftrK(m,n) + (-4 
    /*4*/  * fA1(m,n) * pow2(r(m,n)) * P_1_1(R(m,n)) * (2 * P_1_1(R(m,n)) 
    /*5*/  + P_2_1(R(m,n))) - (4 * pow2(r(m,n)) * P_1_1(R(m,n)) * (2 * P_1_1(R(m,n)) 
    /*6*/  + P_2_1(R(m,n))) * ftrK(m,n)) / 3.) * gtrK(m,n))) + Lt2(m,n) * (24 * exp(4 
    /*3*/  * fconf(m,n) + 2 * gconf(m,n)) * pow2(gB(m,n)) * pow2(1 + gDB(m,n) * r(m,n) 
    /*3*/  + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) * pow2(R(m,n)) * pow3(fA(m,n)) 
    /*2*/  * P_1_1(R(m,n)) * P_1_2(R(m,n)) + pow2(gA(m,n)) * (exp(2 * fconf(m,n) + 4 
    /*4*/  * gconf(m,n)) * fA(m,n) * fB(m,n) * gB(m,n) * (12 * pow3(R(m,n)) 
    /*4*/  * P_0_2(R(m,n)) * P_2_1(R(m,n)) * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) 
    /*5*/  * fDconf(m,n) * r(m,n)) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) 
    /*5*/  * gDconf(m,n) * r(m,n)) + 12 * pow2(R(m,n)) * P_1_1(R(m,n)) * (2 
    /*5*/  * P_1_1(R(m,n)) + P_2_1(R(m,n))) * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) 
    /*5*/  * fDconf(m,n) * r(m,n)) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) 
    /*5*/  * gDconf(m,n) * r(m,n))) + exp(4 * fconf(m,n) + 6 * gconf(m,n)) 
    /*3*/  * pow2(gB(m,n)) * pow3(fA(m,n)) * pow3(R(m,n)) * (12 * pow2(gA2(m,n)) 
    /*4*/  * pow2(r(m,n)) * P_0_3(R(m,n)) * P_2_1(R(m,n)) + (4 * pow2(r(m,n)) 
    /*5*/  * pow2(gtrK(m,n)) * P_0_3(R(m,n)) * P_2_1(R(m,n))) / 3. + 8 * gA2(m,n) 
    /*4*/  * pow2(r(m,n)) * P_0_3(R(m,n)) * P_2_1(R(m,n)) * gtrK(m,n))) + exp(2 
    /*3*/  * fconf(m,n) + 8 * gconf(m,n)) * pow2(fA(m,n)) * pow2(gB(m,n)) 
    /*2*/  * pow3(gA(m,n)) * (pow3(R(m,n)) * (12 * fA1(m,n) * (gA1(m,n) - gA2(m,n)) 
    /*4*/  * pow2(r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) + 4 * (gA1(m,n) - gA2(m,n)) 
    /*4*/  * pow2(r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) * ftrK(m,n)) + Power(R(m,n),4)
    /*3*/  * (-12 * fA1(m,n) * gA2(m,n) * pow2(r(m,n)) * P_0_2(R(m,n)) * P_2_1(R(m,n))
    /*4*/  - 4 * gA2(m,n) * pow2(r(m,n)) * P_0_2(R(m,n)) * P_2_1(R(m,n)) * ftrK(m,n) 
    /*4*/  + (-4 * fA1(m,n) * pow2(r(m,n)) * P_0_2(R(m,n)) * P_2_1(R(m,n)) - (4 
    /*6*/  * pow2(r(m,n)) * P_0_2(R(m,n)) * P_2_1(R(m,n)) * ftrK(m,n)) / 3.) 
    /*4*/  * gtrK(m,n)))) + pow2(gA(m,n)) * (fA(m,n) * (-24 * exp(2 * fconf(m,n) + 4 
    /*4*/  * gconf(m,n)) * fB(m,n) * gB(m,n) * pow3(R(m,n)) * P_1_1(R(m,n)) 
    /*3*/  * P_1_2(R(m,n)) * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) 
    /*4*/  * r(m,n)) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) 
    /*3*/  + exp(4 * fconf(m,n) + 2 * gconf(m,n)) * pow2(fB(m,n)) * (12 * pow2(1 
    /*5*/  + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) * pow2(R(m,n)) 
    /*4*/  * (2 * P_1_1(R(m,n)) * P_1_2(R(m,n)) - P_0_2(R(m,n)) * P_2_1(R(m,n))) + 6 
    /*4*/  * pow2(1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) 
    /*4*/  * P_1_1(R(m,n)) * P_2_1(R(m,n)) * R(m,n))) + pow3(fA(m,n)) * (exp(4 
    /*4*/  * fconf(m,n) + 2 * gconf(m,n)) * (6 * pow2(R(m,n)) * P_1_2(R(m,n)) 
    /*4*/  * P_2_1(R(m,n)) - 6 * P_1_1(R(m,n)) * P_2_1(R(m,n)) * R(m,n)) + exp(4 
    /*4*/  * fconf(m,n) + 6 * gconf(m,n)) * pow2(gB(m,n)) * (pow3(R(m,n)) * (-18 
    /*5*/  * pow2(fA1(m,n)) * pow2(r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) - 2 
    /*5*/  * pow2(r(m,n)) * pow2(ftrK(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) - 12 
    /*5*/  * fA1(m,n) * pow2(r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) * ftrK(m,n)) 
    /*4*/  + pow2(R(m,n)) * (18 * pow2(gA2(m,n)) * pow2(r(m,n)) * P_1_2(R(m,n)) 
    /*5*/  * P_2_1(R(m,n)) + 2 * pow2(r(m,n)) * pow2(gtrK(m,n)) * P_1_2(R(m,n)) 
    /*5*/  * P_2_1(R(m,n)) + 12 * gA2(m,n) * pow2(r(m,n)) * P_1_2(R(m,n)) 
    /*5*/  * P_2_1(R(m,n)) * gtrK(m,n)) + pow2(p(m,n)) * pow2(R(m,n)) * (12 
    /*5*/  * pow2(gA2(m,n)) * pow2(r(m,n)) * P_1_2(R(m,n)) * (2 * P_1_1(R(m,n)) 
    /*6*/  + P_2_1(R(m,n))) + (4 * pow2(r(m,n)) * pow2(gtrK(m,n)) * P_1_2(R(m,n)) * (2 
    /*7*/  * P_1_1(R(m,n)) + P_2_1(R(m,n)))) / 3. + 8 * gA2(m,n) * pow2(r(m,n)) 
    /*5*/  * P_1_2(R(m,n)) * (2 * P_1_1(R(m,n)) + P_2_1(R(m,n))) * gtrK(m,n)))) 
    /*2*/  + pow2(fA(m,n)) * (exp(4 * fconf(m,n) + 4 * gconf(m,n)) * fB(m,n) * gB(m,n) 
    /*3*/  * p(m,n) * pow2(R(m,n)) * (12 * gA2(m,n) * (4 * P_1_1(R(m,n)) * P_1_2(R(m,n))
    /*5*/  - P_0_2(R(m,n)) * P_2_1(R(m,n))) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 
    /*5*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) + 4 * (4 * P_1_1(R(m,n)) * P_1_2(R(m,n))
    /*5*/  - P_0_2(R(m,n)) * P_2_1(R(m,n))) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 
    /*5*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) * gtrK(m,n)) + exp(2 * fconf(m,n) + 6 
    /*4*/  * gconf(m,n)) * pow2(gB(m,n)) * (p_r(m,n) * pow3(R(m,n)) * (-12 * gA2(m,n) 
    /*5*/  * pow2(r(m,n)) * P_1_2(R(m,n)) * P_2_1(R(m,n)) - 4 * pow2(r(m,n)) 
    /*5*/  * P_1_2(R(m,n)) * P_2_1(R(m,n)) * gtrK(m,n)) + p(m,n) * pow3(R(m,n)) * (-12 
    /*5*/  * P_1_2(R(m,n)) * (2 * gA2(m,n) * P_1_1(R(m,n)) + gA1(m,n) * P_2_1(R(m,n))) 
    /*5*/  * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) 
    /*5*/  - 4 * P_1_2(R(m,n)) * (2 * P_1_1(R(m,n)) + P_2_1(R(m,n))) * r(m,n) * (1 
    /*6*/  + gDB(m,n) * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) * gtrK(m,n))))))
    /*0*/  - Lt2(m,n) * (pow2(gA(m,n)) * (exp(2 * fconf(m,n) + 6 * gconf(m,n)) 
    /*2*/  * p(m,n) * pow2(fA(m,n)) * pow2(gB(m,n)) * (Power(R(m,n),4) * (12 * fA1(m,n)
    /*4*/  * P_0_2(R(m,n)) * P_2_1(R(m,n)) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 
    /*5*/  * gconf(m,n) * gDconf(m,n) * r(m,n)) + 4 * P_0_2(R(m,n)) * P_2_1(R(m,n)) 
    /*4*/  * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) 
    /*4*/  * ftrK(m,n)) + pow3(R(m,n)) * (12 * fA1(m,n) * P_1_1(R(m,n)) * (2 
    /*5*/  * P_1_1(R(m,n)) + P_2_1(R(m,n))) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 
    /*5*/  * gconf(m,n) * gDconf(m,n) * r(m,n)) + 4 * P_1_1(R(m,n)) * (2 * P_1_1(R(m,n))
    /*5*/  + P_2_1(R(m,n))) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) 
    /*5*/  * gDconf(m,n) * r(m,n)) * ftrK(m,n))) + exp(4 * fconf(m,n) + 6 * gconf(m,n))
    /*2*/  * pow2(gB(m,n)) * pow3(fA(m,n)) * pow3(R(m,n)) * (12 * fA1(m,n) * gA2(m,n)
    /*3*/  * pow2(r(m,n)) * (P_0_2(R(m,n)) - 3 * P_1_2(R(m,n))) * P_2_1(R(m,n)) + 4 
    /*3*/  * gA2(m,n) * pow2(r(m,n)) * (P_0_2(R(m,n)) - 3 * P_1_2(R(m,n))) 
    /*3*/  * P_2_1(R(m,n)) * ftrK(m,n) + (4 * fA1(m,n) * pow2(r(m,n)) * (P_0_2(R(m,n)) 
    /*5*/  - 3 * P_1_2(R(m,n))) * P_2_1(R(m,n)) + (4 * pow2(r(m,n)) * (P_0_2(R(m,n)) - 3
    /*6*/  * P_1_2(R(m,n))) * P_2_1(R(m,n)) * ftrK(m,n)) / 3.) * gtrK(m,n))) 
    /*1*/  + pow3(gA(m,n)) * (-24 * exp(2 * fconf(m,n) + 4 * gconf(m,n)) * pow2(fB(m,n))
    /*2*/  * pow2(P_1_1(R(m,n))) * pow2(1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) 
    /*3*/  * fDconf(m,n) * r(m,n)) * pow2(R(m,n)) + exp(2 * fconf(m,n) + 6 * gconf(m,n))
    /*2*/  * fA(m,n) * fB(m,n) * gB(m,n) * p(m,n) * (pow3(R(m,n)) * (-12 * gA2(m,n) 
    /*4*/  * P_0_2(R(m,n)) * P_2_1(R(m,n)) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 
    /*5*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) - 4 * P_0_2(R(m,n)) * P_2_1(R(m,n)) 
    /*4*/  * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) 
    /*4*/  * gtrK(m,n)) + pow2(R(m,n)) * (-12 * gA2(m,n) * P_1_1(R(m,n)) * (2 
    /*5*/  * P_1_1(R(m,n)) + P_2_1(R(m,n))) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 
    /*5*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) - 4 * P_1_1(R(m,n)) * (2 * P_1_1(R(m,n))
    /*5*/  + P_2_1(R(m,n))) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) 
    /*5*/  * fDconf(m,n) * r(m,n)) * gtrK(m,n))) + pow2(fA(m,n)) * (6 * exp(2 
    /*4*/  * fconf(m,n) + 4 * gconf(m,n)) * pow2(R(m,n)) * (2 * P_1_1(R(m,n)) 
    /*4*/  - P_2_1(R(m,n))) * P_2_1(R(m,n)) + exp(2 * fconf(m,n) + 8 * gconf(m,n)) 
    /*3*/  * pow2(gB(m,n)) * (pow3(R(m,n)) * ((4 * pow2(r(m,n)) * pow2(gtrK(m,n)) 
    /*6*/  * (P_0_2(R(m,n)) + P_1_2(R(m,n))) * P_2_1(R(m,n))) / 3. + 12 * gA2(m,n) 
    /*5*/  * pow2(r(m,n)) * (gA2(m,n) * P_0_2(R(m,n)) + gA1(m,n) * P_1_2(R(m,n))) 
    /*5*/  * P_2_1(R(m,n)) + 4 * pow2(r(m,n)) * (gA1(m,n) * P_1_2(R(m,n)) + gA2(m,n) 
    /*6*/  * (2 * P_0_2(R(m,n)) + P_1_2(R(m,n)))) * P_2_1(R(m,n)) * gtrK(m,n)) 
    /*4*/  + pow2(R(m,n)) * (6 * pow2(gA2(m,n)) * pow2(r(m,n)) * (4 * P_1_1(R(m,n)) 
    /*6*/  - P_2_1(R(m,n))) * P_2_1(R(m,n)) + (2 * pow2(r(m,n)) * pow2(gtrK(m,n)) * (4 
    /*7*/  * P_1_1(R(m,n)) - P_2_1(R(m,n))) * P_2_1(R(m,n))) / 3. + 4 * gA2(m,n) 
    /*5*/  * pow2(r(m,n)) * (4 * P_1_1(R(m,n)) - P_2_1(R(m,n))) * P_2_1(R(m,n)) 
    /*5*/  * gtrK(m,n))))) + gA(m,n) * (pow2(fA(m,n)) * (exp(2 * fconf(m,n) + 4 
    /*4*/  * gconf(m,n)) * pow2(gB(m,n)) * (-6 * pow2(1 + gDB(m,n) * r(m,n) + 2 
    /*5*/  * gconf(m,n) * gDconf(m,n) * r(m,n)) * pow2(R(m,n)) * (4 * P_1_1(R(m,n)) 
    /*5*/  - P_2_1(R(m,n))) * P_2_1(R(m,n)) + 12 * pow2(1 + gDB(m,n) * r(m,n) + 2 
    /*5*/  * gconf(m,n) * gDconf(m,n) * r(m,n)) * pow3(R(m,n)) * (2 * P_1_1(R(m,n)) 
    /*5*/  * P_1_2(R(m,n)) - P_0_2(R(m,n)) * P_2_1(R(m,n)))) - 12 * exp(4 * fconf(m,n) 
    /*4*/  + 2 * gconf(m,n)) * fB(m,n) * gB(m,n) * pow2(R(m,n)) * (4 * P_1_1(R(m,n)) 
    /*4*/  * P_1_2(R(m,n)) - P_0_2(R(m,n)) * P_2_1(R(m,n))) * (1 + fDB(m,n) * r(m,n) + 2
    /*4*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) * (1 + gDB(m,n) * r(m,n) + 2 
    /*4*/  * gconf(m,n) * gDconf(m,n) * r(m,n))) + exp(4 * fconf(m,n) + 4 * gconf(m,n))
    /*2*/  * p(m,n) * pow2(gB(m,n)) * pow3(fA(m,n)) * (pow3(R(m,n)) * (-24 * gA2(m,n)
    /*4*/  * P_0_3(R(m,n)) * P_2_1(R(m,n)) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 
    /*5*/  * gconf(m,n) * gDconf(m,n) * r(m,n)) - 8 * P_0_3(R(m,n)) * P_2_1(R(m,n)) 
    /*4*/  * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) 
    /*4*/  * gtrK(m,n)) + pow2(R(m,n)) * (-8 * P_1_2(R(m,n)) * (2 * P_1_1(R(m,n)) 
    /*5*/  + P_2_1(R(m,n))) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) 
    /*5*/  * gDconf(m,n) * r(m,n)) * gtrK(m,n) - 8 * P_1_2(R(m,n)) * r(m,n) * (6 
    /*5*/  * gA2(m,n) * P_1_1(R(m,n)) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) 
    /*6*/  * gDconf(m,n) * r(m,n)) + P_2_1(R(m,n)) * (3 * gA1(m,n) * (1 + gDB(m,n) 
    /*7*/  * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) - r(m,n) * (3 * gA2_r(m,n)
    /*7*/  + gtrK_r(m,n))))))));
}
Real BimetricEvolve::eq_fW( Int m, Int n )
{
    return k_g * (-18 * exp(4 * fconf(m,n) + 6 * gconf(m,n)) * Lt2(m,n) 
    /*1*/  * pow2(gA(m,n)) * pow2(gB(m,n)) * pow2(P_2_1(R(m,n))) * pow2(r(m,n)) 
    /*1*/  * pow2(R(m,n)) * pow3(fA(m,n)) * P_1_2(R(m,n)) + 3 * exp(2 * fconf(m,n) + 8 
    /*2*/  * gconf(m,n)) * Lt(m,n) * pow2(fA(m,n)) * pow2(gB(m,n)) * pow2(r(m,n)) 
    /*1*/  * pow2(R(m,n)) * pow3(gA(m,n)) * P_2_1(R(m,n)) * (pow2(P_2_1(R(m,n))) + 2 
    /*2*/  * P_1_2(R(m,n)) * (eq_pf_grho(m,n) - P_2_0(R(m,n))) - 4 * P_1_1(R(m,n)) 
    /*2*/  * P_2_1(R(m,n)))) + k_f * (Lt2(m,n) * (-12 * exp(10 * gconf(m,n)) * fA(m,n) 
    /*2*/  * Power(gA(m,n),4) * pow2(gB(m,n)) * pow2(p(m,n)) * pow2(P_2_1(R(m,n))) 
    /*2*/  * pow2(r(m,n)) * P_1_1(R(m,n)) * R(m,n) + exp(4 * fconf(m,n) + 6 
    /*3*/  * gconf(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n)) * pow3(fA(m,n)) * (-6 
    /*3*/  * eq_pf_frho(m,n) * pow2(r(m,n)) * pow3(R(m,n)) * P_1_2(R(m,n)) * P_2_1(R(m,n)) - 6
    /*3*/  * pow2(r(m,n)) * pow2(R(m,n)) * P_2_1(R(m,n)) * (eq_pf_frho(m,n) * (P_1_1(R(m,n))
    /*5*/  - 2 * P_2_1(R(m,n))) - 2 * eq_pf_fJ22(m,n) * (P_1_1(R(m,n)) - P_2_1(R(m,n))) 
    /*4*/  + eq_pf_fJ11(m,n) * P_2_1(R(m,n))) - 6 * pow2(r(m,n)) * P_1_1(R(m,n)) 
    /*3*/  * P_2_1(R(m,n)) * P_2_2(R(m,n)) + 6 * pow2(r(m,n)) * P_2_1(R(m,n)) 
    /*3*/  * (P_1_3(R(m,n)) * P_2_1(R(m,n)) + P_1_2(R(m,n)) * P_2_2(R(m,n))) * R(m,n)))
    /*1*/  + Lt(m,n) * (-12 * exp(2 * fconf(m,n) + 8 * gconf(m,n)) * eq_pf_fJ22(m,n) 
    /*2*/  * Lt2(m,n) * pow2(fA(m,n)) * pow2(gB(m,n)) * pow2(r(m,n)) * pow3(gA(m,n)) 
    /*2*/  * pow3(R(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) + exp(2 * fconf(m,n) + 8 
    /*3*/  * gconf(m,n)) * pow2(fA(m,n)) * pow2(gB(m,n)) * pow3(gA(m,n)) * (-3 
    /*3*/  * pow2(r(m,n)) * pow3(P_2_1(R(m,n))) + 6 * eq_pf_frho(m,n) * pow2(r(m,n)) 
    /*3*/  * pow3(R(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) + 6 * pow2(r(m,n)) 
    /*3*/  * P_2_1(R(m,n)) * (2 * P_1_2(R(m,n)) * P_2_1(R(m,n)) + P_1_1(R(m,n)) 
    /*4*/  * P_2_2(R(m,n))) * R(m,n) + pow2(p(m,n)) * (6 * pow2(P_2_1(R(m,n))) 
    /*4*/  * pow2(r(m,n)) * (P_1_1(R(m,n)) - P_2_1(R(m,n))) + 12 * eq_pf_frho(m,n) 
    /*4*/  * pow2(r(m,n)) * pow3(R(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) + 6 
    /*4*/  * pow2(P_2_1(R(m,n))) * pow2(r(m,n)) * P_1_2(R(m,n)) * R(m,n))))) 
    /*0*/  + pow2(gA(m,n)) * (-12 * exp(2 * fconf(m,n) + 4 * gconf(m,n)) * fA(m,n) 
    /*1*/  * fB(m,n) * gB(m,n) * p(m,n) * p_r(m,n) * pow2(R(m,n)) * P_1_2(R(m,n)) 
    /*1*/  * P_2_1(R(m,n)) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) 
    /*2*/  * fDconf(m,n) * r(m,n)) + exp(2 * fconf(m,n) + 6 * gconf(m,n)) * p_r(m,n) 
    /*1*/  * pow2(fA(m,n)) * pow2(gB(m,n)) * pow2(R(m,n)) * (12 * gA2(m,n) 
    /*2*/  * pow2(r(m,n)) * P_1_2(R(m,n)) * P_2_1(R(m,n)) + 4 * pow2(r(m,n)) 
    /*2*/  * P_1_2(R(m,n)) * P_2_1(R(m,n)) * gtrK(m,n))) + Lt2(m,n) * (-24 * exp(4 
    /*2*/  * fconf(m,n) + 2 * gconf(m,n)) * pow2(gB(m,n)) * pow2(P_1_2(R(m,n))) * pow2(1
    /*2*/  + gDB(m,n) * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) * pow2(R(m,n))
    /*1*/  * pow3(fA(m,n)) + exp(4 * fconf(m,n) + 4 * gconf(m,n)) * gA(m,n) * p(m,n) 
    /*1*/  * pow2(gB(m,n)) * pow3(fA(m,n)) * (pow2(R(m,n)) * (-24 * fA1(m,n) 
    /*3*/  * P_1_2(R(m,n)) * (P_1_1(R(m,n)) - P_2_1(R(m,n))) * r(m,n) * (1 + gDB(m,n) 
    /*4*/  * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) - 8 * P_1_2(R(m,n)) 
    /*3*/  * (P_1_1(R(m,n)) - P_2_1(R(m,n))) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 
    /*4*/  * gconf(m,n) * gDconf(m,n) * r(m,n)) * ftrK(m,n)) + pow3(R(m,n)) * (-12 
    /*3*/  * fA1(m,n) * P_0_3(R(m,n)) * P_2_1(R(m,n)) * r(m,n) * (1 + gDB(m,n) * r(m,n)
    /*4*/  + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) - 4 * P_0_3(R(m,n)) 
    /*3*/  * P_2_1(R(m,n)) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) 
    /*4*/  * gDconf(m,n) * r(m,n)) * ftrK(m,n))) + pow2(gA(m,n)) * (fA(m,n) * (exp(2 
    /*4*/  * fconf(m,n) + 4 * gconf(m,n)) * fB(m,n) * gB(m,n) * (-12 * pow2(R(m,n)) 
    /*4*/  * P_1_2(R(m,n)) * (2 * P_1_1(R(m,n)) + P_2_1(R(m,n))) * (1 + fDB(m,n) 
    /*5*/  * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) * (1 + gDB(m,n) * r(m,n) 
    /*5*/  + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) + 12 * pow3(R(m,n)) * (2 
    /*5*/  * pow2(P_1_2(R(m,n))) - P_0_3(R(m,n)) * P_2_1(R(m,n))) * (1 + fDB(m,n) 
    /*5*/  * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) * (1 + gDB(m,n) * r(m,n) 
    /*5*/  + 2 * gconf(m,n) * gDconf(m,n) * r(m,n))) + exp(4 * fconf(m,n) + 2 
    /*4*/  * gconf(m,n)) * pow2(fB(m,n)) * (6 * pow2(1 + fDB(m,n) * r(m,n) + 2 
    /*5*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) - 12 
    /*4*/  * pow2(1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) 
    /*4*/  * pow2(R(m,n)) * (2 * pow2(P_1_2(R(m,n))) - P_0_3(R(m,n)) * P_2_1(R(m,n))) 
    /*4*/  - 6 * pow2(1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) 
    /*4*/  * P_1_2(R(m,n)) * P_2_1(R(m,n)) * R(m,n))) + pow3(fA(m,n)) * (exp(4 
    /*4*/  * fconf(m,n) + 2 * gconf(m,n)) * (-6 * P_1_1(R(m,n)) * P_2_1(R(m,n)) + 6 
    /*4*/  * P_1_2(R(m,n)) * P_2_1(R(m,n)) * R(m,n)) + exp(4 * fconf(m,n) + 6 
    /*4*/  * gconf(m,n)) * pow2(gB(m,n)) * (pow2(R(m,n)) * (-6 * pow2(fA1(m,n)) 
    /*5*/  * pow2(r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) - (2 * pow2(r(m,n)) 
    /*6*/  * pow2(ftrK(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n))) / 3. - 4 * fA1(m,n) 
    /*5*/  * pow2(r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) * ftrK(m,n)) + pow3(R(m,n)) 
    /*4*/  * (6 * pow2(fA1(m,n)) * pow2(r(m,n)) * (2 * P_0_2(R(m,n)) - 3 
    /*6*/  * P_1_2(R(m,n))) * P_2_1(R(m,n)) + (2 * pow2(r(m,n)) * pow2(ftrK(m,n)) * (2 
    /*7*/  * P_0_2(R(m,n)) - 3 * P_1_2(R(m,n))) * P_2_1(R(m,n))) / 3. + 4 * fA1(m,n) 
    /*5*/  * pow2(r(m,n)) * (2 * P_0_2(R(m,n)) - 3 * P_1_2(R(m,n))) * P_2_1(R(m,n)) 
    /*5*/  * ftrK(m,n)))) + exp(4 * fconf(m,n) + 4 * gconf(m,n)) * fB(m,n) * gB(m,n) 
    /*2*/  * p(m,n) * pow2(fA(m,n)) * pow2(R(m,n)) * (-12 * gA2(m,n) * (2 
    /*4*/  * pow2(P_1_2(R(m,n))) - P_0_3(R(m,n)) * P_2_1(R(m,n))) * r(m,n) * (1 
    /*4*/  + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) - 4 * (2 
    /*4*/  * pow2(P_1_2(R(m,n))) - P_0_3(R(m,n)) * P_2_1(R(m,n))) * r(m,n) * (1 
    /*4*/  + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) * gtrK(m,n))) 
    /*1*/  + pow3(gA(m,n)) * (fA(m,n) * (exp(2 * fconf(m,n) + 6 * gconf(m,n)) * fB(m,n)
    /*3*/  * gB(m,n) * p(m,n) * (pow2(R(m,n)) * (-24 * fA1(m,n) * P_1_1(R(m,n)) 
    /*5*/  * (P_1_1(R(m,n)) - P_2_1(R(m,n))) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 
    /*6*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) - 8 * P_1_1(R(m,n)) * (P_1_1(R(m,n)) 
    /*6*/  - P_2_1(R(m,n))) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) 
    /*6*/  * fDconf(m,n) * r(m,n)) * ftrK(m,n)) + pow3(R(m,n)) * (24 * fA1(m,n) 
    /*5*/  * (P_1_1(R(m,n)) * P_1_2(R(m,n)) - P_0_2(R(m,n)) * P_2_1(R(m,n))) * r(m,n) 
    /*5*/  * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) + 8 
    /*5*/  * (P_1_1(R(m,n)) * P_1_2(R(m,n)) - P_0_2(R(m,n)) * P_2_1(R(m,n))) * r(m,n) 
    /*5*/  * (1 + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) 
    /*5*/  * ftrK(m,n))) - 8 * exp(8 * gconf(m,n)) * p(m,n) * pow2(gB(m,n)) 
    /*3*/  * pow2(r(m,n)) * pow3(R(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) * (3 
    /*4*/  * fA1_r(m,n) + ftrK_r(m,n))) + exp(2 * fconf(m,n) + 8 * gconf(m,n)) 
    /*2*/  * pow2(fA(m,n)) * pow2(gB(m,n)) * pow3(R(m,n)) * (12 * fA1(m,n) 
    /*3*/  * pow2(r(m,n)) * (gA2(m,n) * P_0_2(R(m,n)) + gA1(m,n) * P_1_2(R(m,n))) 
    /*3*/  * P_2_1(R(m,n)) + 4 * pow2(r(m,n)) * (gA2(m,n) * P_0_2(R(m,n)) + gA1(m,n) 
    /*4*/  * P_1_2(R(m,n))) * P_2_1(R(m,n)) * ftrK(m,n) + (4 * fA1(m,n) * pow2(r(m,n)) 
    /*4*/  * (P_0_2(R(m,n)) + P_1_2(R(m,n))) * P_2_1(R(m,n)) + (4 * pow2(r(m,n)) 
    /*5*/  * (P_0_2(R(m,n)) + P_1_2(R(m,n))) * P_2_1(R(m,n)) * ftrK(m,n)) / 3.) 
    /*3*/  * gtrK(m,n)))) + Lt(m,n) * (gA(m,n) * (pow2(fA(m,n)) * (exp(2 * fconf(m,n) 
    /*4*/  + 4 * gconf(m,n)) * pow2(gB(m,n)) * (18 * pow2(1 + gDB(m,n) * r(m,n) + 2 
    /*5*/  * gconf(m,n) * gDconf(m,n) * r(m,n)) * pow2(R(m,n)) * P_1_2(R(m,n)) 
    /*4*/  * P_2_1(R(m,n)) - 12 * pow2(1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) 
    /*5*/  * gDconf(m,n) * r(m,n)) * pow3(R(m,n)) * (2 * pow2(P_1_2(R(m,n))) 
    /*5*/  - P_0_3(R(m,n)) * P_2_1(R(m,n)))) + exp(4 * fconf(m,n) + 2 * gconf(m,n)) 
    /*3*/  * fB(m,n) * gB(m,n) * (48 * pow2(P_1_2(R(m,n))) * pow2(R(m,n)) * (1 
    /*5*/  + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) * (1 + gDB(m,n)
    /*5*/  * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) + 24 * pow2(p(m,n)) 
    /*4*/  * pow2(P_1_2(R(m,n))) * pow2(R(m,n)) * (1 + fDB(m,n) * r(m,n) + 2 
    /*5*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) * (1 + gDB(m,n) * r(m,n) + 2 
    /*5*/  * gconf(m,n) * gDconf(m,n) * r(m,n)))) + exp(4 * fconf(m,n) + 4 * gconf(m,n))
    /*2*/  * p(m,n) * pow2(gB(m,n)) * pow2(R(m,n)) * pow3(fA(m,n)) * (24 * gA2(m,n) 
    /*3*/  * pow2(P_1_2(R(m,n))) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) 
    /*4*/  * gDconf(m,n) * r(m,n)) + 8 * pow2(P_1_2(R(m,n))) * r(m,n) * (1 + gDB(m,n) 
    /*4*/  * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) * gtrK(m,n))) 
    /*1*/  + pow2(gA(m,n)) * (pow2(fA(m,n)) * (exp(4 * fconf(m,n) + 4 * gconf(m,n)) 
    /*3*/  * fB(m,n) * gB(m,n) * p(m,n) * pow2(R(m,n)) * (12 * fA1(m,n) * P_1_2(R(m,n))
    /*4*/  * (2 * P_1_1(R(m,n)) - 3 * P_2_1(R(m,n))) * r(m,n) * (1 + fDB(m,n) * r(m,n)
    /*5*/  + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) + 4 * P_1_2(R(m,n)) * (2 
    /*5*/  * P_1_1(R(m,n)) - 3 * P_2_1(R(m,n))) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 
    /*5*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) * ftrK(m,n)) + exp(2 * fconf(m,n) + 6 
    /*4*/  * gconf(m,n)) * pow2(gB(m,n)) * (p_r(m,n) * pow3(R(m,n)) * (-12 * fA1(m,n) 
    /*5*/  * pow2(r(m,n)) * P_1_2(R(m,n)) * P_2_1(R(m,n)) - 4 * pow2(r(m,n)) 
    /*5*/  * P_1_2(R(m,n)) * P_2_1(R(m,n)) * ftrK(m,n)) + p(m,n) * pow3(R(m,n)) * (-12 
    /*5*/  * fA1(m,n) * (4 * P_1_1(R(m,n)) * P_1_2(R(m,n)) - P_0_2(R(m,n)) 
    /*6*/  * P_2_1(R(m,n))) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 * gconf(m,n) 
    /*6*/  * gDconf(m,n) * r(m,n)) - 4 * (4 * P_1_1(R(m,n)) * P_1_2(R(m,n)) 
    /*6*/  - P_0_2(R(m,n)) * P_2_1(R(m,n))) * r(m,n) * (1 + gDB(m,n) * r(m,n) + 2 
    /*6*/  * gconf(m,n) * gDconf(m,n) * r(m,n)) * ftrK(m,n)))) + exp(4 * fconf(m,n) + 6
    /*3*/  * gconf(m,n)) * pow2(gB(m,n)) * pow2(p(m,n)) * pow2(R(m,n)) * pow3(fA(m,n))
    /*2*/  * (12 * fA1(m,n) * gA2(m,n) * pow2(r(m,n)) * P_1_2(R(m,n)) * (2 
    /*4*/  * P_1_1(R(m,n)) - 3 * P_2_1(R(m,n))) + 4 * gA2(m,n) * pow2(r(m,n)) 
    /*3*/  * P_1_2(R(m,n)) * (2 * P_1_1(R(m,n)) - 3 * P_2_1(R(m,n))) * ftrK(m,n) + (4 
    /*4*/  * fA1(m,n) * pow2(r(m,n)) * P_1_2(R(m,n)) * (2 * P_1_1(R(m,n)) - 3 
    /*5*/  * P_2_1(R(m,n))) + (4 * pow2(r(m,n)) * P_1_2(R(m,n)) * (2 * P_1_1(R(m,n)) - 3
    /*6*/  * P_2_1(R(m,n))) * ftrK(m,n)) / 3.) * gtrK(m,n))) + Lt2(m,n) * (-12 * exp(4
    /*3*/  * fconf(m,n) + 2 * gconf(m,n)) * fB(m,n) * gA(m,n) * gB(m,n) 
    /*2*/  * pow2(fA(m,n)) * pow2(R(m,n)) * P_0_3(R(m,n)) * P_2_1(R(m,n)) * (1 
    /*3*/  + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) * (1 + gDB(m,n)
    /*3*/  * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n)) + pow3(gA(m,n)) * (24 
    /*3*/  * exp(2 * fconf(m,n) + 4 * gconf(m,n)) * pow2(fB(m,n)) * pow2(1 + fDB(m,n) 
    /*4*/  * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) * pow2(R(m,n)) 
    /*3*/  * P_1_1(R(m,n)) * P_1_2(R(m,n)) + exp(2 * fconf(m,n) + 8 * gconf(m,n)) 
    /*3*/  * pow2(fA(m,n)) * pow2(gB(m,n)) * Power(R(m,n),4) * (-12 * pow2(fA1(m,n)) 
    /*4*/  * pow2(r(m,n)) * P_0_2(R(m,n)) * P_2_1(R(m,n)) - (4 * pow2(r(m,n)) 
    /*5*/  * pow2(ftrK(m,n)) * P_0_2(R(m,n)) * P_2_1(R(m,n))) / 3. - 8 * fA1(m,n) 
    /*4*/  * pow2(r(m,n)) * P_0_2(R(m,n)) * P_2_1(R(m,n)) * ftrK(m,n))) + exp(4 
    /*3*/  * fconf(m,n) + 6 * gconf(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n)) 
    /*2*/  * pow3(fA(m,n)) * (pow3(R(m,n)) * (12 * fA1(m,n) * gA2(m,n) * pow2(r(m,n)) 
    /*4*/  * P_0_3(R(m,n)) * P_2_1(R(m,n)) + 4 * gA2(m,n) * pow2(r(m,n)) * P_0_3(R(m,n))
    /*4*/  * P_2_1(R(m,n)) * ftrK(m,n) + (4 * fA1(m,n) * pow2(r(m,n)) * P_0_3(R(m,n))
    /*5*/  * P_2_1(R(m,n)) + (4 * pow2(r(m,n)) * P_0_3(R(m,n)) * P_2_1(R(m,n)) 
    /*6*/  * ftrK(m,n)) / 3.) * gtrK(m,n)) + pow2(R(m,n)) * (12 * fA1(m,n) * gA2(m,n) 
    /*4*/  * pow2(r(m,n)) * P_1_2(R(m,n)) * P_2_1(R(m,n)) + 4 * gA2(m,n) * pow2(r(m,n))
    /*4*/  * P_1_2(R(m,n)) * P_2_1(R(m,n)) * ftrK(m,n) + (4 * fA1(m,n) * pow2(r(m,n))
    /*5*/  * P_1_2(R(m,n)) * P_2_1(R(m,n)) + (4 * pow2(r(m,n)) * P_1_2(R(m,n)) 
    /*6*/  * P_2_1(R(m,n)) * ftrK(m,n)) / 3.) * gtrK(m,n)))) + pow3(gA(m,n)) * (exp(2 
    /*3*/  * fconf(m,n) + 4 * gconf(m,n)) * pow2(fB(m,n)) * (-12 * pow2(p(m,n)) * pow2(1
    /*4*/  + fDB(m,n) * r(m,n) + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) * pow2(R(m,n))
    /*3*/  * P_0_2(R(m,n)) * P_2_1(R(m,n)) - 6 * pow2(1 + fDB(m,n) * r(m,n) + 2 
    /*4*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) 
    /*3*/  * R(m,n)) + exp(2 * fconf(m,n) + 6 * gconf(m,n)) * fA(m,n) * fB(m,n) 
    /*2*/  * gB(m,n) * p(m,n) * pow2(R(m,n)) * (12 * (gA2(m,n) * P_0_2(R(m,n)) 
    /*4*/  + gA1(m,n) * P_1_2(R(m,n))) * P_2_1(R(m,n)) * r(m,n) * (1 + fDB(m,n) * r(m,n)
    /*4*/  + 2 * fconf(m,n) * fDconf(m,n) * r(m,n)) + 4 * (P_0_2(R(m,n)) 
    /*4*/  + P_1_2(R(m,n))) * P_2_1(R(m,n)) * r(m,n) * (1 + fDB(m,n) * r(m,n) + 2 
    /*4*/  * fconf(m,n) * fDconf(m,n) * r(m,n)) * gtrK(m,n)) + pow2(fA(m,n)) * (exp(2 
    /*4*/  * fconf(m,n) + 4 * gconf(m,n)) * (-6 * pow2(R(m,n)) * P_1_2(R(m,n)) 
    /*4*/  * P_2_1(R(m,n)) + 6 * P_1_1(R(m,n)) * P_2_1(R(m,n)) * R(m,n)) + exp(2 
    /*4*/  * fconf(m,n) + 8 * gconf(m,n)) * pow2(gB(m,n)) * (pow3(R(m,n)) * (6 
    /*5*/  * pow2(fA1(m,n)) * pow2(r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) + (2 
    /*6*/  * pow2(r(m,n)) * pow2(ftrK(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n))) / 3. + 4 
    /*5*/  * fA1(m,n) * pow2(r(m,n)) * P_1_1(R(m,n)) * P_2_1(R(m,n)) * ftrK(m,n)) 
    /*4*/  + pow2(p(m,n)) * pow3(R(m,n)) * (24 * pow2(fA1(m,n)) * pow2(r(m,n)) 
    /*5*/  * P_1_1(R(m,n)) * (-P_1_1(R(m,n)) + P_2_1(R(m,n))) + (8 * pow2(r(m,n)) 
    /*6*/  * pow2(ftrK(m,n)) * P_1_1(R(m,n)) * (-P_1_1(R(m,n)) + P_2_1(R(m,n)))) / 3. 
    /*5*/  + 16 * fA1(m,n) * pow2(r(m,n)) * P_1_1(R(m,n)) * (-P_1_1(R(m,n)) 
    /*6*/  + P_2_1(R(m,n))) * ftrK(m,n)) + pow2(R(m,n)) * (-6 * gA2(m,n) * (2 * gA1(m,n)
    /*6*/  + gA2(m,n)) * pow2(r(m,n)) * P_1_2(R(m,n)) * P_2_1(R(m,n)) - 2 
    /*5*/  * pow2(r(m,n)) * pow2(gtrK(m,n)) * P_1_2(R(m,n)) * P_2_1(R(m,n)) - 4 
    /*5*/  * (gA1(m,n) + 2 * gA2(m,n)) * pow2(r(m,n)) * P_1_2(R(m,n)) * P_2_1(R(m,n)) 
    /*5*/  * gtrK(m,n))))));
}
