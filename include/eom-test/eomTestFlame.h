/** @file  eomTest.h
 *  @author Francesco Torsello
 *  @brief The evolution equations for the test mode.
 *  @version 2019-06-14T09:32:41
 *  @image html eomTest.png
 */

Real BimetricEvolve::eq_u_t( Int m, Int n )
{
    return pow2(u(m,n)) - pow3(u(m,n));
}
Real BimetricEvolve::eq_v_t( Int m, Int n )
{
    return 0;
}
