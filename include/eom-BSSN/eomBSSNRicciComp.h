/** @file  eomBSSNRicciComp.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN terms involving the Ricci tensor in the evolution equations.
 *  @version 2019-05-14T16:24:00
 *  @image html BSSNRicciComp.png
 */

Real gRicci( Int m, Int n )
{
    return exp(-4 * gconf(m,n)) * gAlp(m,n) * ((4 * gL_r(m,n) * pow2(r_minus(m,n))) / 3. 
    /*1*/  - (2 * gsig(m,n)) / pow2(gA(m,n)) + (8 * Power(r_minus(m,n),4) * gB_rr(m,n)) / (3.
    /*2*/  * gB(m,n) * pow2(gA(m,n))) + (8 * Power(r_minus(m,n),4) * pow2(gA_r(m,n))) 
    /*1*/  / Power(gA(m,n),4) - (8 * Power(r_minus(m,n),4) * pow2(gB_r(m,n))) / (3. 
    /*2*/  * pow2(gA(m,n)) * pow2(gB(m,n))) - (8 * Power(r_minus(m,n),4) * gA_rr(m,n)) / (3. 
    /*2*/  * pow3(gA(m,n))) + (4 * r_minus(m,n) * gL(m,n)) / (3. * r(m,n)) - (16 
    /*2*/  * Power(r_minus(m,n),4) * gA_r(m,n)) / (3. * pow3(gA(m,n)) * r(m,n)) + gB_r(m,n) 
    /*1*/  * ((-8 * Power(r_minus(m,n),4) * gA_r(m,n)) / (gB(m,n) * pow3(gA(m,n))) + (4 
    /*3*/  * r_minus(m,n) * (-4 * pow2(r_minus(m,n)) * pow2(gA(m,n)) + pow2(gB(m,n)) * r(m,n) * (4 
    /*5*/  * pow2(r_minus(m,n)) + gsig(m,n) * r(m,n)))) / (3. * pow2(gA(m,n)) * pow3(gB(m,n))
    /*3*/  * r(m,n))));
}
Real fRicci( Int m, Int n )
{
    return exp(-4 * fconf(m,n)) * fAlp(m,n) * ((4 * fL_r(m,n) * pow2(r_minus(m,n))) / 3. 
    /*1*/  - (2 * fsig(m,n)) / pow2(fA(m,n)) + (8 * Power(r_minus(m,n),4) * fB_rr(m,n)) / (3.
    /*2*/  * fB(m,n) * pow2(fA(m,n))) + (8 * Power(r_minus(m,n),4) * pow2(fA_r(m,n))) 
    /*1*/  / Power(fA(m,n),4) - (8 * Power(r_minus(m,n),4) * pow2(fB_r(m,n))) / (3. 
    /*2*/  * pow2(fA(m,n)) * pow2(fB(m,n))) - (8 * Power(r_minus(m,n),4) * fA_rr(m,n)) / (3. 
    /*2*/  * pow3(fA(m,n))) + (4 * r_minus(m,n) * fL(m,n)) / (3. * r(m,n)) - (16 
    /*2*/  * Power(r_minus(m,n),4) * fA_r(m,n)) / (3. * pow3(fA(m,n)) * r(m,n)) + fB_r(m,n) 
    /*1*/  * ((-8 * Power(r_minus(m,n),4) * fA_r(m,n)) / (fB(m,n) * pow3(fA(m,n))) + (4 
    /*3*/  * r_minus(m,n) * (-4 * pow2(r_minus(m,n)) * pow2(fA(m,n)) + pow2(fB(m,n)) * r(m,n) * (4 
    /*5*/  * pow2(r_minus(m,n)) + fsig(m,n) * r(m,n)))) / (3. * pow2(fA(m,n)) * pow3(fB(m,n))
    /*3*/  * r(m,n))));
}
