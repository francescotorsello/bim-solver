/** @file  eomBSSNRicciComp.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN terms involving the Ricci tensor in the evolution equations.
 *  @version 2019-05-15T15:50:45
 *  @image html BSSNRicciComp.png
 */

Real gRicci( Int m, Int n )
{
    return exp(-4 * gconf(m,n)) * gAlp(m,n) * ((2 * gL_r(m,n) * pow2(r_minus(m,n))) / (3. 
    /*2*/  * r_plus(m,n)) - (2 * gsig(m,n)) / pow2(gA(m,n)) + (2 * Power(r_minus(m,n),4) 
    /*2*/  * pow2(gA_r(m,n))) / (Power(gA(m,n),4) * pow2(r_plus(m,n))) - (2 * Power(r_minus(m,n),4)
    /*2*/  * pow2(gB_r(m,n))) / (3. * pow2(r_plus(m,n)) * pow2(gA(m,n)) * pow2(gB(m,n))) 
    /*1*/  + (2 * gL(m,n) * (-(1 / r(m,n)) + r(m,n))) / 3. + (2 * gB_rr(m,n) 
    /*2*/  * pow3(r_minus(m,n)) * (-1 + Power(r(m,n),4))) / (3. * gB(m,n) * pow2(gA(m,n)) 
    /*2*/  * pow3(r_plus(m,n))) - (2 * gA_rr(m,n) * pow3(r_minus(m,n)) * (-1 + Power(r(m,n),4))) 
    /*1*/  / (3. * pow3(r_plus(m,n)) * pow3(gA(m,n))) - (4 * gA_r(m,n) * pow3(r_minus(m,n)) * (3 
    /*3*/  * pow2(r(m,n)) - pow3(r_plus(m,n)) + 4 * Power(r(m,n),4) + Power(r(m,n),6))) / (3.
    /*2*/  * pow3(r_plus(m,n)) * pow3(gA(m,n)) * (r(m,n) + Power(r(m,n),3))) + gB_r(m,n) 
    /*1*/  * ((-2 * gA_r(m,n) * (-1 + pow2(r(m,n))) * pow3(r_minus(m,n))) / (gB(m,n) 
    /*3*/  * pow2(r_plus(m,n)) * pow3(gA(m,n))) + (4 * (-(r_plus(m,n) * pow2(r_minus(m,n)) * pow2(gA(m,n))
    /*5*/  * (-1 + Power(r(m,n),4))) + pow2(gB(m,n)) * pow2(r(m,n)) * ((3 
    /*6*/  + pow2(r(m,n))) * pow3(r_minus(m,n)) + r_plus(m,n) * gsig(m,n) * (-1 
    /*6*/  + Power(r(m,n),4))))) / (3. * pow2(gA(m,n)) * pow3(r_plus(m,n)) * pow3(gB(m,n)) 
    /*3*/  * r(m,n))));
}
Real fRicci( Int m, Int n )
{
    return exp(-4 * fconf(m,n)) * fAlp(m,n) * ((2 * fL_r(m,n) * pow2(r_minus(m,n))) / (3. 
    /*2*/  * r_plus(m,n)) - (2 * fsig(m,n)) / pow2(fA(m,n)) + (2 * Power(r_minus(m,n),4) 
    /*2*/  * pow2(fA_r(m,n))) / (Power(fA(m,n),4) * pow2(r_plus(m,n))) - (2 * Power(r_minus(m,n),4)
    /*2*/  * pow2(fB_r(m,n))) / (3. * pow2(r_plus(m,n)) * pow2(fA(m,n)) * pow2(fB(m,n))) 
    /*1*/  + (2 * fL(m,n) * (-(1 / r(m,n)) + r(m,n))) / 3. + (2 * fB_rr(m,n) 
    /*2*/  * pow3(r_minus(m,n)) * (-1 + Power(r(m,n),4))) / (3. * fB(m,n) * pow2(fA(m,n)) 
    /*2*/  * pow3(r_plus(m,n))) - (2 * fA_rr(m,n) * pow3(r_minus(m,n)) * (-1 + Power(r(m,n),4))) 
    /*1*/  / (3. * pow3(r_plus(m,n)) * pow3(fA(m,n))) - (4 * fA_r(m,n) * pow3(r_minus(m,n)) * (3 
    /*3*/  * pow2(r(m,n)) - pow3(r_plus(m,n)) + 4 * Power(r(m,n),4) + Power(r(m,n),6))) / (3.
    /*2*/  * pow3(r_plus(m,n)) * pow3(fA(m,n)) * (r(m,n) + Power(r(m,n),3))) + fB_r(m,n) 
    /*1*/  * ((-2 * fA_r(m,n) * (-1 + pow2(r(m,n))) * pow3(r_minus(m,n))) / (fB(m,n) 
    /*3*/  * pow2(r_plus(m,n)) * pow3(fA(m,n))) + (4 * (-(r_plus(m,n) * pow2(r_minus(m,n)) * pow2(fA(m,n))
    /*5*/  * (-1 + Power(r(m,n),4))) + pow2(fB(m,n)) * pow2(r(m,n)) * ((3 
    /*6*/  + pow2(r(m,n))) * pow3(r_minus(m,n)) + r_plus(m,n) * fsig(m,n) * (-1 
    /*6*/  + Power(r(m,n),4))))) / (3. * pow2(fA(m,n)) * pow3(r_plus(m,n)) * pow3(fB(m,n)) 
    /*3*/  * r(m,n))));
}
