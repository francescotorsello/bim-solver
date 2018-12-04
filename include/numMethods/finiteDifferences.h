/**
 *  @file      finiteDifferences.h
 *  @brief     Macros to emit various finite differences (e.g., for approximating
 *             derivatives, or for the Keisser-Oliger term).
 *  @authors   Mikica Kocic, Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _FINITE_DIFFERENCES_H_INCLUDED
#define _FINITE_DIFFERENCES_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g11a Finite differences (FD)
 *
 *  Here are the macros to emit the spatial derivatives (centered, up-, or downwind)
 *
 *  Convention: If f() denotes a function, then f_t() denotes a time derivative,
 *  f_r() denotes a spatial derivative, f_rr() denotes the 2nd derivative in r, etc.
 *  The centered FD scheme is used by default.
 *  Upwind derivatives are denoted as f_up_r() and downwind as f_down_r().
 *
 *  @warning The code in this group is generated in the Mathematica notebook
 *          "Finite differences.nb" which is in the directory "algorithms".
 */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
#ifndef CFDS_ORDER
    /** The order of the centered final difference scheme (CFDS).
     *  @todo A run-time CFDS_ORDER?
     */
    #define CFDS_ORDER 4
#endif

#if CFDS_ORDER == 2

    #define GF_r(f,m,n)  ( ( \
        1./2. * f(m,n+1) - 1./2. * f(m,n-1) \
    ) * inv_delta_r )

    #define GF_rr(f,m,n)  ( ( \
        f(m,n+1) - 2 * f(m,n) + f(m,n-1) \
    ) * inv_delta_rr )

    #define GF_up_r(f,m,n)  ( ( \
        1./2. * f(m,n+1) - 1./2. * f(m,n-1) \
    ) * inv_delta_r )

    #define GF_up_rr(f,m,n)  ( ( \
        f(m,n+1) - 2 * f(m,n) + f(m,n-1) \
    ) * inv_delta_rr )

    #define GF_down_r(f,m,n)  ( ( \
        1./2. * f(m,n+1) - 1./2. * f(m,n-1) \
    ) * inv_delta_r )

    #define GF_down_rr(f,m,n)  ( ( \
        f(m,n+1) - 2 * f(m,n) + f(m,n-1) \
    ) * inv_delta_rr )

#elif CFDS_ORDER == 4

    #define GF_r(f,m,n)  ( ( \
        - 1./12. * f(m,n+2) + 2./3. * f(m,n+1) - 2./3. * f(m,n-1) + 1./12.  \
        * f(m,n-2) \
    ) * inv_delta_r )

    #define GF_rr(f,m,n)  ( ( \
        - 1./12. * f(m,n+2) + 4./3. * f(m,n+1) - 5./2. * f(m,n) + 4./3. * f(m,n-1) \
        - 1./12. * f(m,n-2) \
    ) * inv_delta_rr )

    #define GF_up_r(f,m,n)  ( ( \
        1./12. * f(m,n+3) - 1./2. * f(m,n+2) + 3./2. * f(m,n+1) - 5./6. * f(m,n) \
        - 1./4. * f(m,n-1) \
    ) * inv_delta_r )

    #define GF_up_rr(f,m,n)  ( ( \
        - 1./12. * f(m,n+3) + 1./3. * f(m,n+2) + 1./2. * f(m,n+1) - 5./3. * f(m,n) \
        + 11./12. * f(m,n-1) \
    ) * inv_delta_rr )

    #define GF_down_r(f,m,n)  ( ( \
        1./4. * f(m,n+1) + 5./6. * f(m,n) - 3./2. * f(m,n-1) + 1./2. * f(m,n-2) \
        - 1./12. * f(m,n-3) \
    ) * inv_delta_r )

    #define GF_down_rr(f,m,n)  ( ( \
        11./12. * f(m,n+1) - 5./3. * f(m,n) + 1./2. * f(m,n-1) + 1./3.  \
        * f(m,n-2) - 1./12. * f(m,n-3) \
    ) * inv_delta_rr )

#elif CFDS_ORDER == 6

    #define GF_r(f,m,n)  ( ( \
        1./60. * f(m,n+3) - 3./20. * f(m,n+2) + 3./4. * f(m,n+1) - 3./4.  \
        * f(m,n-1) + 3./20. * f(m,n-2) - 1./60. * f(m,n-3) \
    ) * inv_delta_r )

    #define GF_rr(f,m,n)  ( ( \
        1./90. * f(m,n+3) - 3./20. * f(m,n+2) + 3./2. * f(m,n+1) - 49./18.  \
        * f(m,n) + 3./2. * f(m,n-1) - 3./20. * f(m,n-2) + 1./90. * f(m,n-3) \
    ) * inv_delta_rr )

    #define GF_up_r(f,m,n)  ( ( \
        1./30. * f(m,n+5) - 1./4. * f(m,n+4) + 5./6. * f(m,n+3) - 5./3.  \
        * f(m,n+2) + 5./2. * f(m,n+1) - 77./60. * f(m,n) - 1./6. * f(m,n-1) \
    ) * inv_delta_r )

    #define GF_up_rr(f,m,n)  ( ( \
        - 13./180. * f(m,n+5) + 31./60. * f(m,n+4) - 19./12. * f(m,n+3) + 47./18. \
        * f(m,n+2) - 17./12. * f(m,n+1) - 49./60. * f(m,n) + 137./180. * f(m,n-1) \
    ) * inv_delta_rr )

    #define GF_down_r(f,m,n)  ( ( \
        1./6. * f(m,n+1) + 77./60. * f(m,n) - 5./2. * f(m,n-1) + 5./3.  \
        * f(m,n-2) - 5./6. * f(m,n-3) + 1./4. * f(m,n-4) - 1./30. * f(m,n-5) \
    ) * inv_delta_r )

    #define GF_down_rr(f,m,n)  ( ( \
        137./180. * f(m,n+1) - 49./60. * f(m,n) - 17./12. * f(m,n-1) + 47./18.  \
        * f(m,n-2) - 19./12. * f(m,n-3) + 31./60. * f(m,n-4) - 13./180. * f(m,n-5) \
    ) * inv_delta_rr )

#elif CFDS_ORDER == 8

    #define GF_r(f,m,n)  ( ( \
        - 1./280. * f(m,n+4) + 4./105. * f(m,n+3) - 1./5. * f(m,n+2) + 4./5.  \
        * f(m,n+1) - 4./5. * f(m,n-1) + 1./5. * f(m,n-2) - 4./105. * f(m,n-3)  \
        + 1./280. * f(m,n-4) \
    ) * inv_delta_r )

    #define GF_rr(f,m,n)  ( ( \
        - 1./560. * f(m,n+4) + 8./315. * f(m,n+3) - 1./5. * f(m,n+2) + 8./5.  \
        * f(m,n+1) - 205./72. * f(m,n) + 8./5. * f(m,n-1) - 1./5. * f(m,n-2)  \
        + 8./315. * f(m,n-3) - 1./560. * f(m,n-4) \
    ) * inv_delta_rr )

    #define GF_up_r(f,m,n)  ( ( \
        1./56. * f(m,n+7) - 1./6. * f(m,n+6) + 7./10. * f(m,n+5) - 7./4.  \
        * f(m,n+4) + 35./12. * f(m,n+3) - 7./2. * f(m,n+2) + 7./2. * f(m,n+1)  \
        - 223./140. * f(m,n) - 1./8. * f(m,n-1) \
    ) * inv_delta_r )

    #define GF_up_rr(f,m,n)  ( ( \
        - 29./560. * f(m,n+7) + 599./1260. * f(m,n+6) - 39./20. * f(m,n+5)  \
        + 47./10. * f(m,n+4) - 529./72. * f(m,n+3) + 153./20. * f(m,n+2) - 83./20.  \
        * f(m,n+1) + 8./315. * f(m,n) + 363./560. * f(m,n-1) \
    ) * inv_delta_rr )

    #define GF_down_r(f,m,n)  ( ( \
        1./8. * f(m,n+1) + 223./140. * f(m,n) - 7./2. * f(m,n-1) + 7./2.  \
        * f(m,n-2) - 35./12. * f(m,n-3) + 7./4. * f(m,n-4) - 7./10. * f(m,n-5)  \
        + 1./6. * f(m,n-6) - 1./56. * f(m,n-7) \
    ) * inv_delta_r )

    #define GF_down_rr(f,m,n)  ( ( \
        363./560. * f(m,n+1) + 8./315. * f(m,n) - 83./20. * f(m,n-1) + 153./20. \
        * f(m,n-2) - 529./72. * f(m,n-3) + 47./10. * f(m,n-4) - 39./20. * f(m,n-5) \
        + 599./1260. * f(m,n-6) - 29./560. * f(m,n-7) \
    ) * inv_delta_rr )

#else
    #error "CFDS_ORDER must be either 2, 4, 6, or 8"
#endif
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g11b Extrapolation macros                                                */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
/** extrapolate_TO4 is an optimized version of 4th order in accuracy extrapolation
 *  using the 4th order Taylor expansion.
 */
#define extrapolate_TO4(f,m,n)  GF(f,m,n) = \
    -   7./48. * GF(f,m,n-8) +  209./144. * GF(f,m,n-7) - 103./16. * GF(f,m,n-6) \
    + 793./48. * GF(f,m,n-5) - 3835./144. * GF(f,m,n-4) + 439./16. * GF(f,m,n-3) \
    - 281./16. * GF(f,m,n-2) +  917./144. * GF(f,m,n-1)

/** extrapolate_LIN4 is a linear extrapolation of the 4th order in accuracy
 *  (used for time derivatives just before the right ghost zone).
 */
#define extrapolate_LIN4(f,m,n)  GF(f,m,n) = \
    1./4. * GF(f,m,n-5) -   4./3. * GF(f,m,n-4) + 3. * GF(f,m,n-3) \
     - 4. * GF(f,m,n-2) + 37./12. * GF(f,m,n-1)

// The following are supplied by FT

#define extrapolate_TO2(f,m,n)  GF(f,m,n) = \
    - 1./2. * GF(f,m,n-4) + 5./2. * GF(f,m,n-3) - 9./2. * GF(f,m,n-2) \
    + 7./2. * GF(f,m,n-1)

#define extrapolate_TO6(f,m,n)  GF(f,m,n) = \
    -     781./34560. * GF(f,m,n-12) +   58793./172800. * GF(f,m,n-11) \
    -   83609./34560. * GF(f,m,n-10) +     13749./1280. * GF(f,m,n-9) \
    -   190463./5760. * GF(f,m,n-8)  +    141913./1920. * GF(f,m,n-7) \
    - 3515771./28800. * GF(f,m,n-6)  +    857909./5760. * GF(f,m,n-5) \
    - 1524691./11520. * GF(f,m,n-4)  +  2872633./34560. * GF(f,m,n-3) \
    - 1206253./34560. * GF(f,m,n-2)  + 1517383./172800. * GF(f,m,n-1)

#define extrapolate_TO8(f,m,n)  GF(f,m,n) = \
    -       35717./17418240. * GF(f,m,n-16) +   25226807./609638400. * GF(f,m,n-15) \
    -    11770133./29030400. * GF(f,m,n-14) +   223675643./87091200. * GF(f,m,n-13) \
    -   205215421./17418240. * GF(f,m,n-12) +    400625593./9676800. * GF(f,m,n-11) \
    -  9982231807./87091200. * GF(f,m,n-10) +  7342823041./29030400. * GF(f,m,n-9) \
    - 30179572973./67737600. * GF(f,m,n-8)  + 10853340863./17418240. * GF(f,m,n-7) \
    - 59628334229./87091200. * GF(f,m,n-6)  + 16888083827./29030400. * GF(f,m,n-5) \
    - 32514179699./87091200. * GF(f,m,n-4)  + 15201876083./87091200. * GF(f,m,n-3) \
    -    108392833./1935360. * GF(f,m,n-2)  + 6714244639./609638400. * GF(f,m,n-1)

/////////////////////////////////////////////////////////////////////////////////////////

#ifdef FIXED_EXTRAPOLATION

    #define extrapolate_R(f,m,n)  extrapolate_TO4(f,m,n)
    #define extrapolate_D(f,m,n)  extrapolate_LIN4(f,m,n)

#elif CFDS_ORDER == 2

    #define extrapolate_R(f,m,n)  extrapolate_TO2(f,m,n)
    #define extrapolate_D(f,m,n)  extrapolate_TO2(f,m,n)

#elif CFDS_ORDER == 4

    #define extrapolate_R(f,m,n)  extrapolate_TO4(f,m,n)
    #define extrapolate_D(f,m,n)  extrapolate_TO4(f,m,n)

#elif CFDS_ORDER == 6

    #define extrapolate_R(f,m,n)  extrapolate_TO6(f,m,n)
    #define extrapolate_D(f,m,n)  extrapolate_TO6(f,m,n)

#elif CFDS_ORDER == 8

    #define extrapolate_R(f,m,n)  extrapolate_TO8(f,m,n)
    #define extrapolate_D(f,m,n)  extrapolate_TO8(f,m,n)

#endif
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g11c Upwind differencing scheme for convection terms.
  * References:
  * [1] arXiv:1810.12346 [gr-qc] (p.5),
  * [2] https://en.wikipedia.org/wiki/Upwind_scheme,
  * [3] https://en.wikipedia.org/wiki/Upwind_differencing_scheme _for _convection,
  * [4] Patankar, S. V. (1980). Numerical Heat Transfer and Fluid Flow. Taylor &
  *     Francis. ISBN 978-0-89116-522-4 (p.83).
  * [5] Versteeg, H.K. and Malalasekera, W., An Introduction to Computational Fluid
  *     Dynamics: The Finite Volume Method, Pearson Education, 2007,
  *     ISBN: 9780131274983 (Chapter 5, p.134)
  */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
/** GF_ud_r(bf,f,m,n) applies either up or down difference scheme based on
 *  the sign of the velocity field grid function `bf`.
 */
#define GF_ud_r(bf,f,m,n)   ( bf(m,n) > 0 ? GF_up_r(f,m,n) : GF_down_r(f,m,n) )

#ifndef _UPWIND
    #define _UPWIND 1     /// @todo Why is _UPWIND=1 by default?
#endif

#if _UPWIND == 1

    /** Convective derivative term, where a velocity field `bf` is multiplying
     *  the spatial derivative of `f`: `bf(m,n) * GF_ud_r(bf,f,m,n)`.
     */
    #define GF_convr(bf,f,m,n)   ( bf(m,n) > 0 ? bf(m,n) * GF_up_r(f,m,n) \
                                               : bf(m,n) * GF_down_r(f,m,n) )
#else

    #define GF_convr(bf,f,m,n)   ( bf(m,n) * GF_r(f,m,n) )

#endif
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g11d Kreiss-Oliger term                                                  */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
/** `KreissOligerTerm` is a macro that gives a Kreiss-Oliger dissipation term.
 *
 *  Coefficients d_k at gridpoints for the centred finite-difference discretisation
 *  of the derivative `partial_x^{2N}` (of the dissipative Kreiss-Oliger operator).
 *  <pre>
 *      ---------------------------------------------
 *       N  -4  -3   -2   -1    0   +1   +2  +3  +4
 *      ---------------------------------------------
 *       1                +1   -2   +1
 *       2           -1   +4   -6   +4   -1
 *       3      +1   -6  +15  -20  +15   -6  +1
 *       4  -1  +8  -28  +56  -70  +56  -28  +8  -1
 *      ---------------------------------------------  </pre>
 *
 *  Here, `N >= 1` is an integer and `D_h^{2N}` is a centered difference operator
 *  approximating a partial spatial derivative of order `2N`, e.g.,
 *  a second- (`N = 1`) or fourth-order (`N = 2`) spatial derivative.
 *
 *  To be subtracted from the evolution equation:
 *
 *      `epsilon (-1)^N  delta_t {delta_x}^{2N-1}  D_h^{2N}`
 *
 *  where:
 *
 *      `D_h^{2N} = (2 delta_x)^{-2N} * ( sum_k d_k f[n+k] )`
 *
 *  For `N = 2`, subtract from the evolution equation `u_t`:
 *
 *      `dissip(f,dt) = epsilon 2^{-2N} * dt/delta_x * ( sum_k d_k f[n+k] )`
 *
 *  @see Rezzolla and Zanotti, Relativistic Hydrodynamics, 2013,
 *       @cite Rezzolla:2013rel
 *
 *  @todo fixme: The order of the Kreiss-Oliger term depends on the order of the TIME
 *        integration? This is not implemented yet. Check [1, p.5].
 *
 */

#if !defined(KO_ORDER) || KO_ORDER <= 0
    /** KO_ORDER defines the order of the Kreiss-Oliger term order.
     *  Defaults to CFDS_ORDER, if it is undefined or set to a non-positive value.
     */
    #define KO_ORDER CFDS_ORDER
#endif

#if KO_ORDER == 2

    #define KreissOligerTerm(f,dt) \
        ( - ( GF(f,m,n-1) - 2* GF(f,m,n) + GF(f,m,n+1) ) * dissip_delta_r * (dt) )

#elif KO_ORDER == 4

    #define KreissOligerTerm(f,dt) \
        ( ( GF(f,m,n-2) - 4* GF(f,m,n-1) + 6* GF(f,m,n) - 4* GF(f,m,n+1) + GF(f,m,n+2) \
           ) * dissip_delta_r * (dt) )

#elif KO_ORDER == 6

    #define KreissOligerTerm(f,dt) \
        ( - (      GF(f,m,n-3) - 6* GF(f,m,n-2) + 15* GF(f,m,n-1) - 20* GF(f,m,n) \
             + 15* GF(f,m,n+1) - 6* GF(f,m,n+2) +     GF(f,m,n+3) \
            ) * dissip_delta_r * (dt) )

#elif KO_ORDER == 8

    #define KreissOligerTerm(f,dt) \
        ( (       GF(f,m,n-4) -  8* GF(f,m,n-3) + 28* GF(f,m,n-2) - 56* GF(f,m,n-1) \
            + 70* GF(f,m,n)   - 56* GF(f,m,n+1) + 28* GF(f,m,n+2) -  8* GF(f,m,n+3) \
            +     GF(f,m,n+4) \
            ) * dissip_delta_r * (dt) )

#endif
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _FINITE_DIFFERENCES_H_INCLUDED
