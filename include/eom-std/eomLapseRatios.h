/** @file  eomLapseRatios.h
 *  @brief The lapse factors BimetricEvolve::eq_gW and BimetricEvolve::eq_fW.
 */

Real BimetricEvolve::eq_gW( Int m, Int n )
{
    return (2 * fA(m,n) * gA(m,n) * Lt(m,n) * p(m,n) * pow2(gB(m,n)) * pow2(R(m,n)) 
    /*1*/  * r(m,n) * (2 * fDB(m,n) * (gK(m,n) - gKD(m,n)) * pow2(gA(m,n)) 
    /*2*/  * pow2(p(m,n)) * pow2(b_1 + b_2 * R(m,n)) * R(m,n) + fDB(m,n) * (gK(m,n) 
    /*3*/  - gKD(m,n)) * pow2(gA(m,n)) * R(m,n) * (2 * pow2(b_1 + b_2 * R(m,n)) + (1 
    /*4*/  + pow2(p(m,n))) * (pow2(b_1) + 2 * pow2(b_2) * pow2(R(m,n)) + b_1 * R(m,n) 
    /*4*/  * (2 * b_2 - b_3 * R(m,n)))) + fA(m,n) * gA(m,n) * Lt(m,n) * R(m,n) 
    /*2*/  * (-((fK(m,n) - fKD(m,n)) * gDB(m,n) * Lt(m,n) * (3 * pow2(b_1) + 4 
    /*5*/  * pow2(b_2) * pow2(R(m,n)) + b_1 * R(m,n) * (6 * b_2 - b_3 * R(m,n)))) 
    /*3*/  + fDB(m,n) * (gK(m,n) - gKD(m,n)) * (-(b_2 * R(m,n) * (2 * b_2 + 3 * b_3 
    /*6*/  * R(m,n))) - b_1 * (3 * b_2 + 4 * b_3 * R(m,n)))) + gDB(m,n) * Lt(m,n) 
    /*2*/  * pow2(fA(m,n)) * ((fK(m,n) - fKD(m,n)) * R(m,n) * (-(b_2 * R(m,n) * (2 * b_2
    /*6*/  + 3 * b_3 * R(m,n))) - b_1 * (3 * b_2 + 4 * b_3 * R(m,n))) + 2 * (gK(m,n) 
    /*4*/  - gKD(m,n)) * Lt(m,n) * (b_2 * R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n)) + b_1 
    /*4*/  * (3 * b_2 + 4 * b_3 * R(m,n)))))) / 3. + (1 + pow2(p(m,n))) * R(m,n) 
    /*0*/  * (pow2(fDB(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n)) * pow2(R(m,n)) * (4 
    /*2*/  * gA(m,n) * Lt(m,n) * pow2(b_1 + b_2 * R(m,n)) * R(m,n) - fA(m,n) 
    /*2*/  * (pow2(b_1) + 5 * b_1 * R(m,n) * (b_2 + b_3 * R(m,n)) + b_2 * pow2(R(m,n)) 
    /*3*/  * (2 * b_2 + 3 * b_3 * R(m,n)))) - 2 * fA(m,n) * fDB(m,n) * gA(m,n) 
    /*1*/  * gDB(m,n) * Lt(m,n) * pow2(gB(m,n)) * pow2(R(m,n)) * (gA(m,n) * Lt(m,n) * (3
    /*3*/  * pow2(b_1) + 4 * pow2(b_2) * pow2(R(m,n)) + b_1 * R(m,n) * (6 * b_2 - b_3
    /*4*/  * R(m,n))) + fA(m,n) * (-(b_2 * R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n))) 
    /*3*/  - b_1 * (3 * b_2 + 4 * b_3 * R(m,n)))) + pow2(fA(m,n)) * (-(gA(m,n) * Lt(m,n)
    /*3*/  * pow2(gB(m,n)) * pow2(gDB(m,n)) * R(m,n) * (-3 * pow2(b_1) - 4 * pow2(b_2)
    /*4*/  * pow2(R(m,n)) + pow2(b_3) * pow4(R(m,n)) - 6 * b_1 * b_2 * R(m,n) + 2 
    /*4*/  * b_3 * pow2(R(m,n)) * (b_1 + b_2 * R(m,n)))) + fA(m,n) * pow2(gA(m,n)) 
    /*2*/  * (b_1 - b_3 * pow2(R(m,n))) * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) 
    /*2*/  + Lt(m,n) * (-b_1 + b_3 * pow2(R(m,n))) * pow3(gA(m,n)) * R(m,n) * (b_1 
    /*3*/  + R(m,n) * (2 * b_2 + b_3 * R(m,n))) + fA(m,n) * pow2(gB(m,n)) 
    /*2*/  * pow2(gDB(m,n)) * R(m,n) * (pow2(p(m,n)) * (-(pow2(b_3) * pow3(R(m,n))) - 5
    /*4*/  * b_3 * R(m,n) * (b_1 + b_2 * R(m,n)) - b_2 * (3 * b_1 + 2 * b_2 * R(m,n)))
    /*3*/  + (1 + pow2(p(m,n))) * (b_2 + b_3 * R(m,n)) * (-3 * b_1 + R(m,n) * (-2 
    /*5*/  * b_2 + b_3 * R(m,n)))))) + pow2(r(m,n)) * ((k_g * pow2(fA(m,n)) * (1 
    /*3*/  + pow2(p(m,n))) * pow2(R(m,n)) * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) 
    /*2*/  * (2 * Lt(m,n) * pow2(gB(m,n)) * pow3(gA(m,n)) * (2 * b_0 - b_3 
    /*4*/  * pow3(R(m,n)) + 3 * b_1 * R(m,n)) * (b_1 + b_2 * R(m,n)) + 4 * fA(m,n) 
    /*3*/  * Lt(m,n) * p(m,n) * pfS(m,n) * (b_2 + b_3 * R(m,n)) - 2 * fA(m,n) * gA(m,n)
    /*3*/  * ((pfD(m,n) + pftau(m,n)) * pow2(p(m,n)) + pfS(m,n) * pfv(m,n) * (1 
    /*5*/  + pow2(p(m,n)))) * (b_2 + b_3 * R(m,n)) + pow2(gA(m,n)) * (Lt(m,n) 
    /*4*/  * (pfD(m,n) - pfS(m,n) * pfv(m,n) + pftau(m,n)) * (b_1 - b_3 * pow2(R(m,n)))
    /*4*/  + fA(m,n) * pow2(gB(m,n)) * (3 * pow2(b_1) + 2 * b_2 * (b_0 + 3 * b_2 
    /*6*/  * pow2(R(m,n))) - pow2(b_3) * pow4(R(m,n)) + 2 * b_3 * (b_0 + b_2 
    /*6*/  * pow2(R(m,n))) * R(m,n) + 6 * b_1 * R(m,n) * (2 * b_2 + b_3 * R(m,n)))))) 
    /*1*/  / 2. + (pow2(fA(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n)) * (1 + pow2(p(m,n))) 
    /*2*/  * pow2(R(m,n)) * (-(fA(m,n) * pow2(gK(m,n) - gKD(m,n)) * pow2(p(m,n)) * (b_2
    /*5*/  + b_3 * R(m,n)) * (3 * b_1 + R(m,n) * (2 * b_2 - b_3 * R(m,n)))) + gA(m,n)
    /*3*/  * Lt(m,n) * (-4 * pow2(fK(m,n) - fKD(m,n)) * pow2(R(m,n)) * pow2(b_1 + b_2
    /*5*/  * R(m,n)) + pow2(gK(m,n) - gKD(m,n)) * (-3 * pow2(b_1) - 4 * pow2(b_2) 
    /*5*/  * pow2(R(m,n)) + pow2(b_3) * pow4(R(m,n)) - 6 * b_1 * b_2 * R(m,n) + 2 * b_3
    /*5*/  * pow2(R(m,n)) * (b_1 + b_2 * R(m,n))) + 2 * (fK(m,n) - fKD(m,n)) 
    /*4*/  * (gK(m,n) - gKD(m,n)) * Lt(m,n) * R(m,n) * (3 * pow2(b_1) + 4 * pow2(b_2) 
    /*5*/  * pow2(R(m,n)) + b_1 * R(m,n) * (6 * b_2 - b_3 * R(m,n)))) + fA(m,n) 
    /*3*/  * (pow2(gK(m,n) - gKD(m,n)) * (1 + pow2(p(m,n))) * (-(pow2(b_3) 
    /*6*/  * pow3(R(m,n))) - 5 * b_3 * R(m,n) * (b_1 + b_2 * R(m,n)) - b_2 * (3 * b_1 
    /*6*/  + 2 * b_2 * R(m,n))) + pow2(fK(m,n) - fKD(m,n)) * R(m,n) * (b_1 + b_2 
    /*5*/  * R(m,n)) * (b_1 - R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n))) + 2 * (fK(m,n) 
    /*5*/  - fKD(m,n)) * (gK(m,n) - gKD(m,n)) * Lt(m,n) * R(m,n) * (b_2 * R(m,n) * (2 
    /*6*/  * b_2 + 3 * b_3 * R(m,n)) + b_1 * (3 * b_2 + 4 * b_3 * R(m,n)))))) / 9. 
    /*1*/  - (k_f * pow2(fA(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n)) * (1 + pow2(p(m,n))) 
    /*2*/  * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) * (6 * gA(m,n) * Lt(m,n) * R(m,n)
    /*3*/  * (b_1 + b_2 * R(m,n)) * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) 
    /*3*/  + fA(m,n) * (-pow2(b_1) + 2 * b_1 * R(m,n) * (b_2 + R(m,n) * (3 * b_3 + b_4 
    /*6*/  * R(m,n))) + pow2(R(m,n)) * (6 * pow2(b_2) + 3 * pow2(b_3) * pow2(R(m,n)) + 2
    /*5*/  * b_2 * R(m,n) * (6 * b_3 + b_4 * R(m,n)))))) / 2.);
}

Real BimetricEvolve::eq_fW( Int m, Int n )
{
    return (2 * fA(m,n) * gA(m,n) * Lt(m,n) * p(m,n) * pow2(gB(m,n)) * pow2(R(m,n)) 
    /*1*/  * r(m,n) * (2 * (-fK(m,n) + fKD(m,n)) * gDB(m,n) * pow2(fA(m,n)) 
    /*2*/  * pow2(p(m,n)) * (b_2 + b_3 * R(m,n)) * (b_1 + R(m,n) * (2 * b_2 + b_3 
    /*4*/  * R(m,n))) + (fK(m,n) - fKD(m,n)) * gDB(m,n) * pow2(fA(m,n)) * ((1 
    /*4*/  + pow2(p(m,n))) * (5 * pow2(b_3) * pow3(R(m,n)) + 2 * b_2 * (b_1 + 4 * b_2 
    /*5*/  * R(m,n)) + b_3 * R(m,n) * (b_1 + 12 * b_2 * R(m,n))) - 2 * (b_2 + b_3 
    /*4*/  * R(m,n)) * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n)))) + fA(m,n) * gA(m,n) 
    /*2*/  * Lt(m,n) * R(m,n) * (fDB(m,n) * (gK(m,n) - gKD(m,n)) * Lt(m,n) * (-4 
    /*4*/  * pow2(b_2) + b_3 * (b_1 - 3 * b_3 * pow2(R(m,n))) - 6 * b_2 * b_3 * R(m,n))
    /*3*/  + (fK(m,n) - fKD(m,n)) * gDB(m,n) * (-(b_2 * R(m,n) * (2 * b_2 + 3 * b_3 
    /*6*/  * R(m,n))) - b_1 * (3 * b_2 + 4 * b_3 * R(m,n)))) + fDB(m,n) * Lt(m,n) 
    /*2*/  * pow2(gA(m,n)) * R(m,n) * ((gK(m,n) - gKD(m,n)) * (-(b_2 * R(m,n) * (2 * b_2
    /*6*/  + 3 * b_3 * R(m,n))) - b_1 * (3 * b_2 + 4 * b_3 * R(m,n))) + 2 * (fK(m,n) 
    /*4*/  - fKD(m,n)) * Lt(m,n) * R(m,n) * (b_2 * R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n))
    /*4*/  + b_1 * (3 * b_2 + 4 * b_3 * R(m,n)))))) / 3. + (1 + pow2(p(m,n))) 
    /*0*/  * (Lt(m,n) * pow3(fA(m,n)) * (-4 * pow2(gB(m,n)) * pow2(gDB(m,n)) 
    /*2*/  * pow2(R(m,n)) * pow2(b_2 + b_3 * R(m,n)) + pow2(gA(m,n)) * (-b_1 + b_3 
    /*3*/  * pow2(R(m,n))) * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n)))) - fA(m,n) 
    /*1*/  * fDB(m,n) * Lt(m,n) * pow2(gA(m,n)) * pow2(gB(m,n)) * pow2(R(m,n)) * (4 
    /*2*/  * R(m,n) * (b_2 + b_3 * R(m,n)) * (gDB(m,n) * (b_1 + b_2 * R(m,n)) 
    /*3*/  + (fDB(m,n) - gDB(m,n)) * R(m,n) * (b_2 + b_3 * R(m,n))) + (b_1 + R(m,n) * (2
    /*4*/  * b_2 + b_3 * R(m,n))) * (-(fDB(m,n) * (b_1 + b_3 * pow2(R(m,n)))) + 2 
    /*3*/  * gDB(m,n) * R(m,n) * (b_2 + 2 * b_3 * R(m,n)))) + pow2(fDB(m,n)) 
    /*1*/  * pow2(gB(m,n)) * pow3(gA(m,n)) * pow3(R(m,n)) * (pow2(p(m,n)) * (pow2(b_1) 
    /*3*/  + 5 * b_1 * R(m,n) * (b_2 + b_3 * R(m,n)) + b_2 * pow2(R(m,n)) * (2 * b_2 + 3
    /*4*/  * b_3 * R(m,n))) + (1 + pow2(p(m,n))) * (b_1 + b_2 * R(m,n)) * (-b_1 
    /*3*/  + R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n)))) + gA(m,n) * pow2(fA(m,n)) * R(m,n)
    /*1*/  * (pow2(gA(m,n)) * (b_1 + b_2 * R(m,n) - R(m,n) * (b_2 + b_3 * R(m,n))) 
    /*2*/  * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) + gDB(m,n) * pow2(gB(m,n)) 
    /*2*/  * R(m,n) * (-4 * (gDB(m,n) - 2 * fDB(m,n) * (1 + pow2(p(m,n)))) * pow2(b_2 
    /*4*/  + b_3 * R(m,n)) * R(m,n) + (2 * b_3 * (gDB(m,n) - fDB(m,n) * (1 
    /*6*/  + pow2(p(m,n)))) * R(m,n) + 3 * gDB(m,n) * (b_2 + b_3 * R(m,n))) * (b_1 
    /*4*/  + R(m,n) * (2 * b_2 + b_3 * R(m,n)))))) + pow2(r(m,n)) * ((k_g 
    /*2*/  * pow2(fA(m,n)) * pow2(gA(m,n)) * (1 + pow2(p(m,n))) * pow2(R(m,n)) * (b_1 
    /*3*/  + R(m,n) * (2 * b_2 + b_3 * R(m,n))) * (gA(m,n) * pow2(gB(m,n)) * (3 
    /*4*/  * pow2(b_1) + 2 * b_2 * (b_0 + 3 * b_2 * pow2(R(m,n))) - pow2(b_3) 
    /*4*/  * pow4(R(m,n)) + 2 * b_3 * (b_0 + b_2 * pow2(R(m,n))) * R(m,n) + 6 * b_1 
    /*4*/  * R(m,n) * (2 * b_2 + b_3 * R(m,n))) + 2 * (b_2 + b_3 * R(m,n)) * (pfD(m,n) 
    /*4*/  + pftau(m,n) + 3 * fA(m,n) * Lt(m,n) * pow2(gB(m,n)) * (b_1 + R(m,n) * (2 
    /*6*/  * b_2 + b_3 * R(m,n)))))) / 2. + (pow2(fA(m,n)) * pow2(gA(m,n)) 
    /*2*/  * pow2(gB(m,n)) * (1 + pow2(p(m,n))) * pow2(R(m,n)) * (-(gA(m,n) 
    /*4*/  * pow2(fK(m,n) - fKD(m,n)) * pow2(p(m,n)) * R(m,n) * (b_1 + b_2 * R(m,n)) 
    /*4*/  * (b_1 - R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n)))) + gA(m,n) * (pow2(gK(m,n) 
    /*5*/  - gKD(m,n)) * (b_2 + b_3 * R(m,n)) * (3 * b_1 + R(m,n) * (2 * b_2 - b_3 
    /*6*/  * R(m,n))) + pow2(fK(m,n) - fKD(m,n)) * (1 + pow2(p(m,n))) * R(m,n) 
    /*4*/  * (pow2(b_1) + 5 * b_1 * R(m,n) * (b_2 + b_3 * R(m,n)) + b_2 * pow2(R(m,n)) 
    /*5*/  * (2 * b_2 + 3 * b_3 * R(m,n))) + 2 * (fK(m,n) - fKD(m,n)) * (gK(m,n) 
    /*5*/  - gKD(m,n)) * Lt(m,n) * R(m,n) * (-(b_2 * R(m,n) * (2 * b_2 + 3 * b_3 
    /*7*/  * R(m,n))) - b_1 * (3 * b_2 + 4 * b_3 * R(m,n)))) + fA(m,n) * Lt(m,n) * (4 
    /*4*/  * pow2(gK(m,n) - gKD(m,n)) * pow2(b_2 + b_3 * R(m,n)) + 2 * (fK(m,n) 
    /*5*/  - fKD(m,n)) * (gK(m,n) - gKD(m,n)) * Lt(m,n) * R(m,n) * (-4 * pow2(b_2) + b_3
    /*5*/  * (b_1 - 3 * b_3 * pow2(R(m,n))) - 6 * b_2 * b_3 * R(m,n)) + pow2(fK(m,n) 
    /*5*/  - fKD(m,n)) * (-pow2(b_1) - 2 * b_1 * R(m,n) * (b_2 + b_3 * R(m,n)) 
    /*5*/  + pow2(R(m,n)) * (4 * pow2(b_2) + 3 * pow2(b_3) * pow2(R(m,n)) + 6 * b_2 
    /*6*/  * b_3 * R(m,n)))))) / 9. - (k_f * pow2(fA(m,n)) * pow2(gA(m,n)) 
    /*2*/  * pow2(gB(m,n)) * (1 + pow2(p(m,n))) * (b_1 + R(m,n) * (2 * b_2 + b_3 
    /*4*/  * R(m,n))) * (2 * fA(m,n) * Lt(m,n) * (b_2 + b_3 * R(m,n)) * (-b_1 
    /*4*/  + pow2(R(m,n)) * (3 * b_3 + 2 * b_4 * R(m,n))) + gA(m,n) * (-((b_1 + R(m,n) 
    /*6*/  * (2 * b_2 + b_3 * R(m,n))) * (b_1 - R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n))))
    /*4*/  - 2 * pow2(p(m,n)) * R(m,n) * (b_1 + b_2 * R(m,n)) * (b_2 + R(m,n) * (2 
    /*6*/  * b_3 + b_4 * R(m,n))) + 2 * (1 + pow2(p(m,n))) * R(m,n) * (b_1 + b_2 
    /*5*/  * R(m,n)) * (b_2 + R(m,n) * (2 * b_3 + b_4 * R(m,n)))))) / 2.);
}

Real BimetricEvolve::eq_cW( Int m, Int n )
{
    return -(Lt(m,n) * pow2(fA(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n)) * (1 
    /*2*/  + pow2(p(m,n))) * pow2(r(m,n)) * pow2(R(m,n)) * (b_1 + R(m,n) * (2 * b_2 
    /*3*/  + b_3 * R(m,n))));
}

Real BimetricEvolve::eq_gW_0( Int m, Int n )
{
    return (1 + pow2(p(m,n))) * R(m,n) * (pow2(fDB(m,n)) * pow2(gA(m,n)) 
    /*1*/  * pow2(gB(m,n)) * pow2(R(m,n)) * (4 * gA(m,n) * Lt(m,n) * pow2(b_1 + b_2 
    /*3*/  * R(m,n)) * R(m,n) - fA(m,n) * (pow2(b_1) + 5 * b_1 * R(m,n) * (b_2 + b_3 
    /*4*/  * R(m,n)) + b_2 * pow2(R(m,n)) * (2 * b_2 + 3 * b_3 * R(m,n)))) - 2 * fA(m,n)
    /*1*/  * fDB(m,n) * gA(m,n) * gDB(m,n) * Lt(m,n) * pow2(gB(m,n)) * pow2(R(m,n)) 
    /*1*/  * (gA(m,n) * Lt(m,n) * (3 * pow2(b_1) + 4 * pow2(b_2) * pow2(R(m,n)) + b_1 
    /*3*/  * R(m,n) * (6 * b_2 - b_3 * R(m,n))) + fA(m,n) * (-(b_2 * R(m,n) * (2 * b_2 
    /*5*/  + 3 * b_3 * R(m,n))) - b_1 * (3 * b_2 + 4 * b_3 * R(m,n)))) + pow2(fA(m,n)) 
    /*1*/  * (-(gA(m,n) * Lt(m,n) * pow2(gB(m,n)) * pow2(gDB(m,n)) * R(m,n) * (-3 
    /*4*/  * pow2(b_1) - 4 * pow2(b_2) * pow2(R(m,n)) + pow2(b_3) * pow4(R(m,n)) - 6 
    /*4*/  * b_1 * b_2 * R(m,n) + 2 * b_3 * pow2(R(m,n)) * (b_1 + b_2 * R(m,n)))) 
    /*2*/  + fA(m,n) * pow2(gA(m,n)) * (b_1 - b_3 * pow2(R(m,n))) * (b_1 + R(m,n) * (2 
    /*4*/  * b_2 + b_3 * R(m,n))) + Lt(m,n) * (-b_1 + b_3 * pow2(R(m,n))) 
    /*2*/  * pow3(gA(m,n)) * R(m,n) * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) 
    /*2*/  + fA(m,n) * pow2(gB(m,n)) * pow2(gDB(m,n)) * R(m,n) * (pow2(p(m,n)) 
    /*3*/  * (-(pow2(b_3) * pow3(R(m,n))) - 5 * b_3 * R(m,n) * (b_1 + b_2 * R(m,n)) 
    /*4*/  - b_2 * (3 * b_1 + 2 * b_2 * R(m,n))) + (1 + pow2(p(m,n))) * (b_2 + b_3 
    /*4*/  * R(m,n)) * (-3 * b_1 + R(m,n) * (-2 * b_2 + b_3 * R(m,n))))));
}

Real BimetricEvolve::eq_gW_1( Int m, Int n )
{
    return (2 * fA(m,n) * gA(m,n) * Lt(m,n) * p(m,n) * pow2(gB(m,n)) * pow2(R(m,n)) 
    /*1*/  * (2 * fDB(m,n) * (gK(m,n) - gKD(m,n)) * pow2(gA(m,n)) * pow2(p(m,n)) 
    /*2*/  * pow2(b_1 + b_2 * R(m,n)) * R(m,n) + fDB(m,n) * (gK(m,n) - gKD(m,n)) 
    /*2*/  * pow2(gA(m,n)) * R(m,n) * (2 * pow2(b_1 + b_2 * R(m,n)) + (1 + pow2(p(m,n)))
    /*3*/  * (pow2(b_1) + 2 * pow2(b_2) * pow2(R(m,n)) + b_1 * R(m,n) * (2 * b_2 - b_3
    /*5*/  * R(m,n)))) + fA(m,n) * gA(m,n) * Lt(m,n) * R(m,n) * (-((fK(m,n) 
    /*5*/  - fKD(m,n)) * gDB(m,n) * Lt(m,n) * (3 * pow2(b_1) + 4 * pow2(b_2) 
    /*5*/  * pow2(R(m,n)) + b_1 * R(m,n) * (6 * b_2 - b_3 * R(m,n)))) + fDB(m,n) 
    /*3*/  * (gK(m,n) - gKD(m,n)) * (-(b_2 * R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n))) 
    /*4*/  - b_1 * (3 * b_2 + 4 * b_3 * R(m,n)))) + gDB(m,n) * Lt(m,n) * pow2(fA(m,n)) 
    /*2*/  * ((fK(m,n) - fKD(m,n)) * R(m,n) * (-(b_2 * R(m,n) * (2 * b_2 + 3 * b_3 
    /*6*/  * R(m,n))) - b_1 * (3 * b_2 + 4 * b_3 * R(m,n))) + 2 * (gK(m,n) - gKD(m,n)) 
    /*3*/  * Lt(m,n) * (b_2 * R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n)) + b_1 * (3 * b_2 + 4
    /*5*/  * b_3 * R(m,n)))))) / 3.;
}

Real BimetricEvolve::eq_gW_2( Int m, Int n )
{
    return (k_g * pow2(fA(m,n)) * (1 + pow2(p(m,n))) * pow2(R(m,n)) * (b_1 + R(m,n) 
    /*2*/  * (2 * b_2 + b_3 * R(m,n))) * (2 * Lt(m,n) * pow2(gB(m,n)) * pow3(gA(m,n)) 
    /*2*/  * (2 * b_0 - b_3 * pow3(R(m,n)) + 3 * b_1 * R(m,n)) * (b_1 + b_2 * R(m,n)) 
    /*2*/  + 4 * fA(m,n) * Lt(m,n) * p(m,n) * pfS(m,n) * (b_2 + b_3 * R(m,n)) - 2 
    /*2*/  * fA(m,n) * gA(m,n) * ((pfD(m,n) + pftau(m,n)) * pow2(p(m,n)) + pfS(m,n) 
    /*3*/  * pfv(m,n) * (1 + pow2(p(m,n)))) * (b_2 + b_3 * R(m,n)) + pow2(gA(m,n)) 
    /*2*/  * (Lt(m,n) * (pfD(m,n) - pfS(m,n) * pfv(m,n) + pftau(m,n)) * (b_1 - b_3 
    /*4*/  * pow2(R(m,n))) + fA(m,n) * pow2(gB(m,n)) * (3 * pow2(b_1) + 2 * b_2 * (b_0 
    /*5*/  + 3 * b_2 * pow2(R(m,n))) - pow2(b_3) * pow4(R(m,n)) + 2 * b_3 * (b_0 + b_2 
    /*5*/  * pow2(R(m,n))) * R(m,n) + 6 * b_1 * R(m,n) * (2 * b_2 + b_3 * R(m,n)))))) 
    /*0*/  / 2. + (pow2(fA(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n)) * (1 + pow2(p(m,n))) 
    /*1*/  * pow2(R(m,n)) * (-(fA(m,n) * pow2(gK(m,n) - gKD(m,n)) * pow2(p(m,n)) * (b_2
    /*4*/  + b_3 * R(m,n)) * (3 * b_1 + R(m,n) * (2 * b_2 - b_3 * R(m,n)))) + gA(m,n)
    /*2*/  * Lt(m,n) * (-4 * pow2(fK(m,n) - fKD(m,n)) * pow2(R(m,n)) * pow2(b_1 + b_2
    /*4*/  * R(m,n)) + pow2(gK(m,n) - gKD(m,n)) * (-3 * pow2(b_1) - 4 * pow2(b_2) 
    /*4*/  * pow2(R(m,n)) + pow2(b_3) * pow4(R(m,n)) - 6 * b_1 * b_2 * R(m,n) + 2 * b_3
    /*4*/  * pow2(R(m,n)) * (b_1 + b_2 * R(m,n))) + 2 * (fK(m,n) - fKD(m,n)) 
    /*3*/  * (gK(m,n) - gKD(m,n)) * Lt(m,n) * R(m,n) * (3 * pow2(b_1) + 4 * pow2(b_2) 
    /*4*/  * pow2(R(m,n)) + b_1 * R(m,n) * (6 * b_2 - b_3 * R(m,n)))) + fA(m,n) 
    /*2*/  * (pow2(gK(m,n) - gKD(m,n)) * (1 + pow2(p(m,n))) * (-(pow2(b_3) 
    /*5*/  * pow3(R(m,n))) - 5 * b_3 * R(m,n) * (b_1 + b_2 * R(m,n)) - b_2 * (3 * b_1 
    /*5*/  + 2 * b_2 * R(m,n))) + pow2(fK(m,n) - fKD(m,n)) * R(m,n) * (b_1 + b_2 
    /*4*/  * R(m,n)) * (b_1 - R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n))) + 2 * (fK(m,n) 
    /*4*/  - fKD(m,n)) * (gK(m,n) - gKD(m,n)) * Lt(m,n) * R(m,n) * (b_2 * R(m,n) * (2 
    /*5*/  * b_2 + 3 * b_3 * R(m,n)) + b_1 * (3 * b_2 + 4 * b_3 * R(m,n)))))) / 9. 
    /*0*/  - (k_f * pow2(fA(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n)) * (1 + pow2(p(m,n))) 
    /*1*/  * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) * (6 * gA(m,n) * Lt(m,n) * R(m,n)
    /*2*/  * (b_1 + b_2 * R(m,n)) * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) 
    /*2*/  + fA(m,n) * (-pow2(b_1) + 2 * b_1 * R(m,n) * (b_2 + R(m,n) * (3 * b_3 + b_4 
    /*5*/  * R(m,n))) + pow2(R(m,n)) * (6 * pow2(b_2) + 3 * pow2(b_3) * pow2(R(m,n)) + 2
    /*4*/  * b_2 * R(m,n) * (6 * b_3 + b_4 * R(m,n)))))) / 2.;
}

Real BimetricEvolve::eq_fW_0( Int m, Int n )
{
    return (1 + pow2(p(m,n))) * (Lt(m,n) * pow3(fA(m,n)) * (-4 * pow2(gB(m,n)) 
    /*2*/  * pow2(gDB(m,n)) * pow2(R(m,n)) * pow2(b_2 + b_3 * R(m,n)) + pow2(gA(m,n)) 
    /*2*/  * (-b_1 + b_3 * pow2(R(m,n))) * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n)))) 
    /*1*/  - fA(m,n) * fDB(m,n) * Lt(m,n) * pow2(gA(m,n)) * pow2(gB(m,n)) * pow2(R(m,n))
    /*1*/  * (4 * R(m,n) * (b_2 + b_3 * R(m,n)) * (gDB(m,n) * (b_1 + b_2 * R(m,n)) 
    /*3*/  + (fDB(m,n) - gDB(m,n)) * R(m,n) * (b_2 + b_3 * R(m,n))) + (b_1 + R(m,n) * (2
    /*4*/  * b_2 + b_3 * R(m,n))) * (-(fDB(m,n) * (b_1 + b_3 * pow2(R(m,n)))) + 2 
    /*3*/  * gDB(m,n) * R(m,n) * (b_2 + 2 * b_3 * R(m,n)))) + pow2(fDB(m,n)) 
    /*1*/  * pow2(gB(m,n)) * pow3(gA(m,n)) * pow3(R(m,n)) * (pow2(p(m,n)) * (pow2(b_1) 
    /*3*/  + 5 * b_1 * R(m,n) * (b_2 + b_3 * R(m,n)) + b_2 * pow2(R(m,n)) * (2 * b_2 + 3
    /*4*/  * b_3 * R(m,n))) + (1 + pow2(p(m,n))) * (b_1 + b_2 * R(m,n)) * (-b_1 
    /*3*/  + R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n)))) + gA(m,n) * pow2(fA(m,n)) * R(m,n)
    /*1*/  * (pow2(gA(m,n)) * (b_1 + b_2 * R(m,n) - R(m,n) * (b_2 + b_3 * R(m,n))) 
    /*2*/  * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) + gDB(m,n) * pow2(gB(m,n)) 
    /*2*/  * R(m,n) * (-4 * (gDB(m,n) - 2 * fDB(m,n) * (1 + pow2(p(m,n)))) * pow2(b_2 
    /*4*/  + b_3 * R(m,n)) * R(m,n) + (2 * b_3 * (gDB(m,n) - fDB(m,n) * (1 
    /*6*/  + pow2(p(m,n)))) * R(m,n) + 3 * gDB(m,n) * (b_2 + b_3 * R(m,n))) * (b_1 
    /*4*/  + R(m,n) * (2 * b_2 + b_3 * R(m,n))))));
}

Real BimetricEvolve::eq_fW_1( Int m, Int n )
{
    return (2 * fA(m,n) * gA(m,n) * Lt(m,n) * p(m,n) * pow2(gB(m,n)) * pow2(R(m,n)) 
    /*1*/  * (2 * (-fK(m,n) + fKD(m,n)) * gDB(m,n) * pow2(fA(m,n)) * pow2(p(m,n)) * (b_2
    /*3*/  + b_3 * R(m,n)) * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) + (fK(m,n) 
    /*3*/  - fKD(m,n)) * gDB(m,n) * pow2(fA(m,n)) * ((1 + pow2(p(m,n))) * (5 * pow2(b_3)
    /*4*/  * pow3(R(m,n)) + 2 * b_2 * (b_1 + 4 * b_2 * R(m,n)) + b_3 * R(m,n) * (b_1 
    /*5*/  + 12 * b_2 * R(m,n))) - 2 * (b_2 + b_3 * R(m,n)) * (b_1 + R(m,n) * (2 * b_2 
    /*5*/  + b_3 * R(m,n)))) + fA(m,n) * gA(m,n) * Lt(m,n) * R(m,n) * (fDB(m,n) 
    /*3*/  * (gK(m,n) - gKD(m,n)) * Lt(m,n) * (-4 * pow2(b_2) + b_3 * (b_1 - 3 * b_3 
    /*5*/  * pow2(R(m,n))) - 6 * b_2 * b_3 * R(m,n)) + (fK(m,n) - fKD(m,n)) * gDB(m,n) 
    /*3*/  * (-(b_2 * R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n))) - b_1 * (3 * b_2 + 4 * b_3
    /*5*/  * R(m,n)))) + fDB(m,n) * Lt(m,n) * pow2(gA(m,n)) * R(m,n) * ((gK(m,n) 
    /*4*/  - gKD(m,n)) * (-(b_2 * R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n))) - b_1 * (3 
    /*5*/  * b_2 + 4 * b_3 * R(m,n))) + 2 * (fK(m,n) - fKD(m,n)) * Lt(m,n) * R(m,n) 
    /*3*/  * (b_2 * R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n)) + b_1 * (3 * b_2 + 4 * b_3 
    /*5*/  * R(m,n)))))) / 3.;
}

Real BimetricEvolve::eq_fW_2( Int m, Int n )
{
    return (k_g * pow2(fA(m,n)) * pow2(gA(m,n)) * (1 + pow2(p(m,n))) * pow2(R(m,n)) 
    /*1*/  * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))) * (gA(m,n) * pow2(gB(m,n)) * (3 
    /*3*/  * pow2(b_1) + 2 * b_2 * (b_0 + 3 * b_2 * pow2(R(m,n))) - pow2(b_3) 
    /*3*/  * pow4(R(m,n)) + 2 * b_3 * (b_0 + b_2 * pow2(R(m,n))) * R(m,n) + 6 * b_1 
    /*3*/  * R(m,n) * (2 * b_2 + b_3 * R(m,n))) + 2 * (b_2 + b_3 * R(m,n)) * (pfD(m,n) 
    /*3*/  + pftau(m,n) + 3 * fA(m,n) * Lt(m,n) * pow2(gB(m,n)) * (b_1 + R(m,n) * (2 
    /*5*/  * b_2 + b_3 * R(m,n)))))) / 2. + (pow2(fA(m,n)) * pow2(gA(m,n)) 
    /*1*/  * pow2(gB(m,n)) * (1 + pow2(p(m,n))) * pow2(R(m,n)) * (-(gA(m,n) 
    /*3*/  * pow2(fK(m,n) - fKD(m,n)) * pow2(p(m,n)) * R(m,n) * (b_1 + b_2 * R(m,n)) 
    /*3*/  * (b_1 - R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n)))) + gA(m,n) * (pow2(gK(m,n) 
    /*4*/  - gKD(m,n)) * (b_2 + b_3 * R(m,n)) * (3 * b_1 + R(m,n) * (2 * b_2 - b_3 
    /*5*/  * R(m,n))) + pow2(fK(m,n) - fKD(m,n)) * (1 + pow2(p(m,n))) * R(m,n) 
    /*3*/  * (pow2(b_1) + 5 * b_1 * R(m,n) * (b_2 + b_3 * R(m,n)) + b_2 * pow2(R(m,n)) 
    /*4*/  * (2 * b_2 + 3 * b_3 * R(m,n))) + 2 * (fK(m,n) - fKD(m,n)) * (gK(m,n) 
    /*4*/  - gKD(m,n)) * Lt(m,n) * R(m,n) * (-(b_2 * R(m,n) * (2 * b_2 + 3 * b_3 
    /*6*/  * R(m,n))) - b_1 * (3 * b_2 + 4 * b_3 * R(m,n)))) + fA(m,n) * Lt(m,n) * (4 
    /*3*/  * pow2(gK(m,n) - gKD(m,n)) * pow2(b_2 + b_3 * R(m,n)) + 2 * (fK(m,n) 
    /*4*/  - fKD(m,n)) * (gK(m,n) - gKD(m,n)) * Lt(m,n) * R(m,n) * (-4 * pow2(b_2) + b_3
    /*4*/  * (b_1 - 3 * b_3 * pow2(R(m,n))) - 6 * b_2 * b_3 * R(m,n)) + pow2(fK(m,n) 
    /*4*/  - fKD(m,n)) * (-pow2(b_1) - 2 * b_1 * R(m,n) * (b_2 + b_3 * R(m,n)) 
    /*4*/  + pow2(R(m,n)) * (4 * pow2(b_2) + 3 * pow2(b_3) * pow2(R(m,n)) + 6 * b_2 
    /*5*/  * b_3 * R(m,n)))))) / 9. - (k_f * pow2(fA(m,n)) * pow2(gA(m,n)) 
    /*1*/  * pow2(gB(m,n)) * (1 + pow2(p(m,n))) * (b_1 + R(m,n) * (2 * b_2 + b_3 
    /*3*/  * R(m,n))) * (2 * fA(m,n) * Lt(m,n) * (b_2 + b_3 * R(m,n)) * (-b_1 
    /*3*/  + pow2(R(m,n)) * (3 * b_3 + 2 * b_4 * R(m,n))) + gA(m,n) * (-((b_1 + R(m,n) 
    /*5*/  * (2 * b_2 + b_3 * R(m,n))) * (b_1 - R(m,n) * (2 * b_2 + 3 * b_3 * R(m,n))))
    /*3*/  - 2 * pow2(p(m,n)) * R(m,n) * (b_1 + b_2 * R(m,n)) * (b_2 + R(m,n) * (2 
    /*5*/  * b_3 + b_4 * R(m,n))) + 2 * (1 + pow2(p(m,n))) * R(m,n) * (b_1 + b_2 
    /*4*/  * R(m,n)) * (b_2 + R(m,n) * (2 * b_3 + b_4 * R(m,n)))))) / 2.;
}

Real BimetricEvolve::eq_cW_0( Int m, Int n )
{
    return 0;
}

Real BimetricEvolve::eq_cW_1( Int m, Int n )
{
    return 0;
}

Real BimetricEvolve::eq_cW_2( Int m, Int n )
{
    return -(Lt(m,n) * pow2(fA(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n)) * (1 
    /*2*/  + pow2(p(m,n))) * pow2(R(m,n)) * (b_1 + R(m,n) * (2 * b_2 + b_3 * R(m,n))));
}
