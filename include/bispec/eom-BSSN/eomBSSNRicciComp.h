/** @file  eomBSSNRicciComp.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN terms involving the Ricci tensor in the evolution equations.
 *  @version 2019-05-24T12:57:16
 *  @image html BSSNRicciComp.png
 */

Real gRicci( Int m, Int n )
{
    return exp(-4 * gconf(m,n)) * gAlp(m,n) * ((2 * gL_r(m,n)
    /*2*/  * pow2(r_minus(m,n))) / (3. * r_plus(m,n)) - (2 * gsig(m,n)) 
    /*1*/  / pow2(gA(m,n)) + (2 * Power(r_minus(m,n),4) * gB_rr(m,n))
    /*1*/  / (3. * gB(m,n) * pow2(r_plus(m,n)) * pow2(gA(m,n))) 
    /*1*/  + (2 * Power(r_minus(m,n),4) * pow2(gA_r(m,n))) 
    /*1*/  / (Power(gA(m,n),4) * pow2(r_plus(m,n))) - (2 
    /*2*/  * Power(r_minus(m,n),4) * pow2(gB_r(m,n))) / (3. 
    /*2*/  * pow2(r_plus(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n))) - (2
    /*2*/  * Power(r_minus(m,n),4) * gA_rr(m,n)) / (3. 
    /*2*/  * pow2(r_plus(m,n)) * pow3(gA(m,n))) + gB_r(m,n) * ((-2 
    /*3*/  * Power(r_minus(m,n),4) * gA_r(m,n)) / (gB(m,n) 
    /*3*/  * pow2(r_plus(m,n)) * pow3(gA(m,n))) - (4 * (-1 
    /*4*/  + pow2(r(m,n))) * (-(pow2(gB(m,n)) * pow2(r(m,n)) 
    /*5*/  * (gsig(m,n) * pow2(r_plus(m,n)) + pow2(r_minus(m,n)) * (3 
    /*7*/  + pow2(r(m,n))))) + pow2(gA(m,n)) * pow2(-1 
    /*5*/  + Power(r(m,n),4)))) / (3. * pow2(gA(m,n)) 
    /*3*/  * pow3(r_plus(m,n)) * pow3(gB(m,n)) * r(m,n))) - (4 
    /*2*/  * Power(r_minus(m,n),4) * gA_r(m,n)) / (3. * pow3(r_plus(m,n)) 
    /*2*/  * pow3(gA(m,n)) * r(m,n)) + (2 * gL(m,n) * (-(1 
    /*4*/  / r(m,n)) + r(m,n))) / 3.);
}
Real fRicci( Int m, Int n )
{
    return exp(-4 * fconf(m,n)) * fAlp(m,n) * ((2 * fL_r(m,n)
    /*2*/  * pow2(r_minus(m,n))) / (3. * r_plus(m,n)) - (2 * fsig(m,n)) 
    /*1*/  / pow2(fA(m,n)) + (2 * Power(r_minus(m,n),4) * fB_rr(m,n))
    /*1*/  / (3. * fB(m,n) * pow2(r_plus(m,n)) * pow2(fA(m,n))) 
    /*1*/  + (2 * Power(r_minus(m,n),4) * pow2(fA_r(m,n))) 
    /*1*/  / (Power(fA(m,n),4) * pow2(r_plus(m,n))) - (2 
    /*2*/  * Power(r_minus(m,n),4) * pow2(fB_r(m,n))) / (3. 
    /*2*/  * pow2(r_plus(m,n)) * pow2(fA(m,n)) * pow2(fB(m,n))) - (2
    /*2*/  * Power(r_minus(m,n),4) * fA_rr(m,n)) / (3. 
    /*2*/  * pow2(r_plus(m,n)) * pow3(fA(m,n))) + fB_r(m,n) * ((-2 
    /*3*/  * Power(r_minus(m,n),4) * fA_r(m,n)) / (fB(m,n) 
    /*3*/  * pow2(r_plus(m,n)) * pow3(fA(m,n))) - (4 * (-1 
    /*4*/  + pow2(r(m,n))) * (-(pow2(fB(m,n)) * pow2(r(m,n)) 
    /*5*/  * (fsig(m,n) * pow2(r_plus(m,n)) + pow2(r_minus(m,n)) * (3 
    /*7*/  + pow2(r(m,n))))) + pow2(fA(m,n)) * pow2(-1 
    /*5*/  + Power(r(m,n),4)))) / (3. * pow2(fA(m,n)) 
    /*3*/  * pow3(r_plus(m,n)) * pow3(fB(m,n)) * r(m,n))) - (4 
    /*2*/  * Power(r_minus(m,n),4) * fA_r(m,n)) / (3. * pow3(r_plus(m,n)) 
    /*2*/  * pow3(fA(m,n)) * r(m,n)) + (2 * fL(m,n) * (-(1 
    /*4*/  / r(m,n)) + r(m,n))) / 3.);
}
