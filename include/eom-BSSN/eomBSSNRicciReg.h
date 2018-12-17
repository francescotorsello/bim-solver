/** @file  eomBSSNRicciReg.h
 *  @author Francesco Torsello
 *  @brief The regularized cBSSN terms involving the Ricci tensor in the evolution equations.
 *  @version 2018-12-17T10:09:27
 *  @image html BSSNRicciReg.png
 */

Real BimetricEvolve::eq_gRicci( Int m, Int n )
{
    return exp(-4 * gconf(m,n)) * gAlp(m,n) * ((2 * gL_r(m,n)) / 3. + ((-2 
    /*3*/  * gDA_r(m,n)) / 3. + (2 * gDB_r(m,n)) / 3. + (4 * pow2(gDA(m,n))) / 3. 
    /*2*/  + gDA(m,n) * (-2 * gDB(m,n) - 4 / (3. * r(m,n)))) / pow2(gA(m,n)) - (2 
    /*2*/  * gL(m,n)) / (3. * r(m,n)) + (4 * gDB(m,n)) / (3. * pow2(gB(m,n)) * r(m,n)) 
    /*1*/  + (gsig(m,n) * (-2 - (4 * gDB(m,n) * r(m,n)) / 3.)) / pow2(gA(m,n)));
}
Real BimetricEvolve::eq_fRicci( Int m, Int n )
{
    return exp(-4 * fconf(m,n)) * fAlp(m,n) * ((2 * fL_r(m,n)) / 3. + ((-2 
    /*3*/  * fDA_r(m,n)) / 3. + (2 * fDB_r(m,n)) / 3. + (4 * pow2(fDA(m,n))) / 3. 
    /*2*/  + fDA(m,n) * (-2 * fDB(m,n) - 4 / (3. * r(m,n)))) / pow2(fA(m,n)) - (2 
    /*2*/  * fL(m,n)) / (3. * r(m,n)) + (4 * fDB(m,n)) / (3. * pow2(fB(m,n)) * r(m,n)) 
    /*1*/  + (fsig(m,n) * (-2 - (4 * fDB(m,n) * r(m,n)) / 3.)) / pow2(fA(m,n)));
}
Real BimetricEvolve::eq_gDers( Int m, Int n )
{
    return (2 * gDA(m,n) * gDAlp(m,n)) / (3. * pow2(gA(m,n))) + (2 * gDB(m,n) 
    /*1*/  * gDAlp(m,n)) / (3. * pow2(gA(m,n))) + (8 * gDconf(m,n) * gDAlp(m,n)) / (3. 
    /*1*/  * pow2(gA(m,n)));
}
Real BimetricEvolve::eq_fDers( Int m, Int n )
{
    return (2 * fDA(m,n) * fDAlp(m,n)) / (3. * pow2(fA(m,n))) + (2 * fDB(m,n) 
    /*1*/  * fDAlp(m,n)) / (3. * pow2(fA(m,n))) + (8 * fDconf(m,n) * fDAlp(m,n)) / (3. 
    /*1*/  * pow2(fA(m,n)));
}
