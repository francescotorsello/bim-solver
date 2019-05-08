/** @file  eomBSSNRicciComp.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN terms involving the Ricci tensor in the evolution equations.
 *  @version 2019-05-08T11:23:56
 *  @image html BSSNRicciComp.png
 */

Real gRicci( Int m, Int n )
{
    return (2 * exp(-4 * gconf(m,n)) * gAlp(m,n) * (-1 + r(m,n)) * (3 
    /*2*/  * pow2(gA_r(m,n)) * pow2(1 + r(m,n)) * pow3(gB(m,n)) * pow3(-1 + r(m,n)) 
    /*2*/  + Power(gA(m,n),4) * (12 * gB(m,n) * (-1 + r(m,n)) - 8 * gB_r(m,n) * pow2(-1
    /*4*/  + r(m,n)) * (1 + r(m,n)) + (2 * gL(m,n) + gL_r(m,n) * (-1 + pow2(r(m,n))))
    /*3*/  * pow3(gB(m,n)) * (1 + r(m,n))) - gA(m,n) * pow2(gB(m,n)) * pow3(-1 
    /*3*/  + r(m,n)) * (1 + r(m,n)) * (3 * gA_r(m,n) * gB_r(m,n) * (1 + r(m,n)) 
    /*3*/  + gB(m,n) * (2 * gA_r(m,n) + gA_rr(m,n) * (1 + r(m,n)))) + gB(m,n) 
    /*2*/  * pow2(gA(m,n)) * (-1 + r(m,n)) * (-12 * pow2(gB(m,n)) - pow2(gB_r(m,n)) 
    /*3*/  * pow2(-1 + pow2(r(m,n))) + gB(m,n) * (-1 + pow2(r(m,n))) * (gB_rr(m,n) * (-1
    /*5*/  + pow2(r(m,n))) + 2 * gB_r(m,n) * (3 + r(m,n)))))) / (3. * Power(gA(m,n),4)
    /*1*/  * pow2(1 + r(m,n)) * pow3(gB(m,n)));
}
Real fRicci( Int m, Int n )
{
    return (2 * exp(-4 * fconf(m,n)) * fAlp(m,n) * (-1 + r(m,n)) * (3 
    /*2*/  * pow2(fA_r(m,n)) * pow2(1 + r(m,n)) * pow3(fB(m,n)) * pow3(-1 + r(m,n)) 
    /*2*/  + Power(fA(m,n),4) * (12 * fB(m,n) * (-1 + r(m,n)) - 8 * fB_r(m,n) * pow2(-1
    /*4*/  + r(m,n)) * (1 + r(m,n)) + (2 * fL(m,n) + fL_r(m,n) * (-1 + pow2(r(m,n))))
    /*3*/  * pow3(fB(m,n)) * (1 + r(m,n))) - fA(m,n) * pow2(fB(m,n)) * pow3(-1 
    /*3*/  + r(m,n)) * (1 + r(m,n)) * (3 * fA_r(m,n) * fB_r(m,n) * (1 + r(m,n)) 
    /*3*/  + fB(m,n) * (2 * fA_r(m,n) + fA_rr(m,n) * (1 + r(m,n)))) + fB(m,n) 
    /*2*/  * pow2(fA(m,n)) * (-1 + r(m,n)) * (-12 * pow2(fB(m,n)) - pow2(fB_r(m,n)) 
    /*3*/  * pow2(-1 + pow2(r(m,n))) + fB(m,n) * (-1 + pow2(r(m,n))) * (fB_rr(m,n) * (-1
    /*5*/  + pow2(r(m,n))) + 2 * fB_r(m,n) * (3 + r(m,n)))))) / (3. * Power(fA(m,n),4)
    /*1*/  * pow2(1 + r(m,n)) * pow3(fB(m,n)));
}
