/** @file  eomBSSNRicciComp.h
 *  @author Francesco Torsello
 *  @brief The compactified cBSSN terms involving the Ricci tensor in the evolution equations.
 *  @version 2019-05-24T13:15:04
 *  @image html BSSNRicciComp.png
 */

Real gRicci( Int m, Int n )
{
    return (-2 * gAlp(m,n) * gsig(m,n)) / (TINY_Real + exp(4
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2)) - (2 
    /*1*/  * Power(r_minus(m,n),4) * gAlp(m,n) * gA_rr(m,n)) 
    /*0*/  / (TINY_Real + 3 * Power(r_plus(m,n),2) * exp(4 
    /*2*/  * gconf(m,n)) * Power(gA(m,n),3)) - (2 
    /*1*/  * Power(r_minus(m,n),4) * gAlp(m,n) * gA_r(m,n) 
    /*1*/  * gB_r(m,n)) / (TINY_Real + Power(r_plus(m,n),2) * exp(4 
    /*2*/  * gconf(m,n)) * Power(gA(m,n),3) * gB(m,n)) + (2 
    /*1*/  * Power(r_minus(m,n),4) * gAlp(m,n) * gB_rr(m,n)) 
    /*0*/  / (TINY_Real + 3 * Power(r_plus(m,n),2) * exp(4 
    /*2*/  * gconf(m,n)) * Power(gA(m,n),2) * gB(m,n)) + (2 
    /*1*/  * gAlp(m,n) * gL_r(m,n) * pow2(r_minus(m,n))) / (TINY_Real
    /*1*/  + 3 * r_plus(m,n) * exp(4 * gconf(m,n))) + (2 
    /*1*/  * Power(r_minus(m,n),4) * gAlp(m,n) * pow2(gA_r(m,n))) 
    /*0*/  / (TINY_Real + Power(r_plus(m,n),2) * exp(4 * gconf(m,n))
    /*1*/  * Power(gA(m,n),4)) - (2 * Power(r_minus(m,n),4) 
    /*1*/  * gAlp(m,n) * pow2(gB_r(m,n))) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),2) * exp(4 * gconf(m,n)) 
    /*1*/  * Power(gA(m,n),2) * Power(gB(m,n),2)) + (4 
    /*1*/  * gAlp(m,n) * gsig(m,n) * gB_r(m,n) * pow3(r(m,n)))
    /*0*/  / (TINY_Real + 3 * r_plus(m,n) * exp(4 * gconf(m,n)) 
    /*1*/  * Power(gA(m,n),2) * gB(m,n)) - (8 * gAlp(m,n) 
    /*1*/  * gB_r(m,n) * pow3(r(m,n))) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * gconf(m,n)) 
    /*1*/  * Power(gB(m,n),3)) + (8 * gAlp(m,n) * gB_r(m,n) 
    /*1*/  * pow2(r_minus(m,n)) * pow3(r(m,n))) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * gconf(m,n)) 
    /*1*/  * Power(gA(m,n),2) * gB(m,n)) + (2 * gAlp(m,n) 
    /*1*/  * gL(m,n) * r(m,n)) / (TINY_Real + 3 * exp(4 
    /*2*/  * gconf(m,n))) - (4 * gAlp(m,n) * gsig(m,n) 
    /*1*/  * gB_r(m,n) * r(m,n)) / (TINY_Real + 3 * r_plus(m,n) 
    /*1*/  * exp(4 * gconf(m,n)) * Power(gA(m,n),2) * gB(m,n))
    /*0*/  - (4 * gAlp(m,n) * gB_r(m,n) * r(m,n)) 
    /*0*/  / (TINY_Real + 3 * Power(r_plus(m,n),3) * exp(4 
    /*2*/  * gconf(m,n)) * Power(gB(m,n),3)) - (4 * gAlp(m,n) 
    /*1*/  * gB_r(m,n) * pow2(r_minus(m,n)) * r(m,n)) / (TINY_Real 
    /*1*/  + Power(r_plus(m,n),3) * exp(4 * gconf(m,n)) 
    /*1*/  * Power(gA(m,n),2) * gB(m,n)) + (8 * gAlp(m,n) 
    /*1*/  * gB_r(m,n) * Power(r(m,n),5)) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * gconf(m,n)) 
    /*1*/  * Power(gB(m,n),3)) + (4 * gAlp(m,n) * gB_r(m,n) 
    /*1*/  * pow2(r_minus(m,n)) * Power(r(m,n),5)) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * gconf(m,n)) 
    /*1*/  * Power(gA(m,n),2) * gB(m,n)) + (4 * gAlp(m,n) 
    /*1*/  * gB_r(m,n) * Power(r(m,n),7)) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * gconf(m,n)) 
    /*1*/  * Power(gB(m,n),3)) - (4 * gAlp(m,n) * gB_r(m,n) 
    /*1*/  * Power(r(m,n),9)) / (TINY_Real + 3 * Power(r_plus(m,n),3)
    /*1*/  * exp(4 * gconf(m,n)) * Power(gB(m,n),3)) - (2 
    /*1*/  * gAlp(m,n) * gL(m,n)) / (TINY_Real + 3 * exp(4 
    /*2*/  * gconf(m,n)) * r(m,n)) - (4 * Power(r_minus(m,n),4) 
    /*1*/  * gAlp(m,n) * gA_r(m,n)) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * gconf(m,n)) 
    /*1*/  * Power(gA(m,n),3) * r(m,n)) + (4 * gAlp(m,n) 
    /*1*/  * gB_r(m,n)) / (TINY_Real + 3 * Power(r_plus(m,n),3) 
    /*1*/  * exp(4 * gconf(m,n)) * Power(gB(m,n),3) * r(m,n));
}
Real fRicci( Int m, Int n )
{
    return (-2 * fAlp(m,n) * fsig(m,n)) / (TINY_Real + exp(4
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2)) - (2 
    /*1*/  * Power(r_minus(m,n),4) * fAlp(m,n) * fA_rr(m,n)) 
    /*0*/  / (TINY_Real + 3 * Power(r_plus(m,n),2) * exp(4 
    /*2*/  * fconf(m,n)) * Power(fA(m,n),3)) - (2 
    /*1*/  * Power(r_minus(m,n),4) * fAlp(m,n) * fA_r(m,n) 
    /*1*/  * fB_r(m,n)) / (TINY_Real + Power(r_plus(m,n),2) * exp(4 
    /*2*/  * fconf(m,n)) * Power(fA(m,n),3) * fB(m,n)) + (2 
    /*1*/  * Power(r_minus(m,n),4) * fAlp(m,n) * fB_rr(m,n)) 
    /*0*/  / (TINY_Real + 3 * Power(r_plus(m,n),2) * exp(4 
    /*2*/  * fconf(m,n)) * Power(fA(m,n),2) * fB(m,n)) + (2 
    /*1*/  * fAlp(m,n) * fL_r(m,n) * pow2(r_minus(m,n))) / (TINY_Real
    /*1*/  + 3 * r_plus(m,n) * exp(4 * fconf(m,n))) + (2 
    /*1*/  * Power(r_minus(m,n),4) * fAlp(m,n) * pow2(fA_r(m,n))) 
    /*0*/  / (TINY_Real + Power(r_plus(m,n),2) * exp(4 * fconf(m,n))
    /*1*/  * Power(fA(m,n),4)) - (2 * Power(r_minus(m,n),4) 
    /*1*/  * fAlp(m,n) * pow2(fB_r(m,n))) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),2) * exp(4 * fconf(m,n)) 
    /*1*/  * Power(fA(m,n),2) * Power(fB(m,n),2)) + (4 
    /*1*/  * fAlp(m,n) * fsig(m,n) * fB_r(m,n) * pow3(r(m,n)))
    /*0*/  / (TINY_Real + 3 * r_plus(m,n) * exp(4 * fconf(m,n)) 
    /*1*/  * Power(fA(m,n),2) * fB(m,n)) - (8 * fAlp(m,n) 
    /*1*/  * fB_r(m,n) * pow3(r(m,n))) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * fconf(m,n)) 
    /*1*/  * Power(fB(m,n),3)) + (8 * fAlp(m,n) * fB_r(m,n) 
    /*1*/  * pow2(r_minus(m,n)) * pow3(r(m,n))) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * fconf(m,n)) 
    /*1*/  * Power(fA(m,n),2) * fB(m,n)) + (2 * fAlp(m,n) 
    /*1*/  * fL(m,n) * r(m,n)) / (TINY_Real + 3 * exp(4 
    /*2*/  * fconf(m,n))) - (4 * fAlp(m,n) * fsig(m,n) 
    /*1*/  * fB_r(m,n) * r(m,n)) / (TINY_Real + 3 * r_plus(m,n) 
    /*1*/  * exp(4 * fconf(m,n)) * Power(fA(m,n),2) * fB(m,n))
    /*0*/  - (4 * fAlp(m,n) * fB_r(m,n) * r(m,n)) 
    /*0*/  / (TINY_Real + 3 * Power(r_plus(m,n),3) * exp(4 
    /*2*/  * fconf(m,n)) * Power(fB(m,n),3)) - (4 * fAlp(m,n) 
    /*1*/  * fB_r(m,n) * pow2(r_minus(m,n)) * r(m,n)) / (TINY_Real 
    /*1*/  + Power(r_plus(m,n),3) * exp(4 * fconf(m,n)) 
    /*1*/  * Power(fA(m,n),2) * fB(m,n)) + (8 * fAlp(m,n) 
    /*1*/  * fB_r(m,n) * Power(r(m,n),5)) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * fconf(m,n)) 
    /*1*/  * Power(fB(m,n),3)) + (4 * fAlp(m,n) * fB_r(m,n) 
    /*1*/  * pow2(r_minus(m,n)) * Power(r(m,n),5)) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * fconf(m,n)) 
    /*1*/  * Power(fA(m,n),2) * fB(m,n)) + (4 * fAlp(m,n) 
    /*1*/  * fB_r(m,n) * Power(r(m,n),7)) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * fconf(m,n)) 
    /*1*/  * Power(fB(m,n),3)) - (4 * fAlp(m,n) * fB_r(m,n) 
    /*1*/  * Power(r(m,n),9)) / (TINY_Real + 3 * Power(r_plus(m,n),3)
    /*1*/  * exp(4 * fconf(m,n)) * Power(fB(m,n),3)) - (2 
    /*1*/  * fAlp(m,n) * fL(m,n)) / (TINY_Real + 3 * exp(4 
    /*2*/  * fconf(m,n)) * r(m,n)) - (4 * Power(r_minus(m,n),4) 
    /*1*/  * fAlp(m,n) * fA_r(m,n)) / (TINY_Real + 3 
    /*1*/  * Power(r_plus(m,n),3) * exp(4 * fconf(m,n)) 
    /*1*/  * Power(fA(m,n),3) * r(m,n)) + (4 * fAlp(m,n) 
    /*1*/  * fB_r(m,n)) / (TINY_Real + 3 * Power(r_plus(m,n),3) 
    /*1*/  * exp(4 * fconf(m,n)) * Power(fB(m,n),3) * r(m,n));
}
