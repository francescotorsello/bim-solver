/** @file  eomTest.h
 *  @author Francesco Torsello
 *  @brief The evolution equations for the test mode.
 *  @version 2019-06-14T09:32:41
 *  @image html eomTest.png
 */

Real BimetricEvolve::eq_u_t( Int m, Int n )
{
    return -80.6 * u(m,n) + 119.4 * v(m,n);
}
Real BimetricEvolve::eq_v_t( Int m, Int n )
{
    return 79.6 * u(m,n) -120.4 * v(m,n);
}
