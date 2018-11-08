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
/** @defgroup g11a Centered stencils                                                   */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
#ifndef CFDS_ORDER
    /** The order of the centered final difference scheme.
     */
    #define CFDS_ORDER 4
#endif

/// @todo A run-time CFDS_ORDER?

/////////////////////////////////////////////////////////////////////////////////////////
// Macros to emit the spatial derivatives (centered in space)
//
// Convention. If f() denotes a function, then f_r() denotes a spatial derivative,
// f_t() denotes a time derivative, and f_rr() denotes the 2nd derivative in r.
//
#if CFDS_ORDER == 2 // Second order central finite difference scheme

    #define GF_r(f)   ( 0.5 * ( f(m,n+1) - f(m,n-1) ) * inv_delta_r )

    #define GF_rr(f)  ( ( f(m,n+1) - 2. * f(m,n) + f(m,n-1) ) * inv_delta_rr )

#elif CFDS_ORDER == 4 // Fourth order central finite difference scheme

    #define GF_r(f,m,n)  \
        ( ( - 1./12. * f(m,n+2) + 2./3.  * f(m,n+1) \
            -  2./3. * f(m,n-1) + 1./12. * f(m,n-2) \
           ) * inv_delta_r )

    #define GF_rr(f,m,n) \
        ( ( - 1./12. * f(m,n+2) +  4./3. * f(m,n+1) - 5./2. * f(m,n) \
            +  4./3. * f(m,n-1) - 1./12. * f(m,n-2) \
           ) * inv_delta_rr )

#elif CFDS_ORDER == 6 // Sixth order central finite difference scheme

    #define GF_r(f,m,n) \
        ( (    1./60. * f(m,n+3) - 9./60. * f(m,n+2) + 45./60. * f(m,n+1) \
            - 45./60. * f(m,n-1) + 9./60. * f(m,n-2) -  1./60. * f(m,n-3) \
           ) * inv_delta_r )

    #define GF_rr(f,m,n) \
        ( (     2./180. * f(m,n+3) -  27./180. * f(m,n+2) + 270./180. * f(m,n+1) \
            - 490./180. * f(m,n)   + 270./180. * f(m,n-1) -  27./180. * f(m,n-2) \
            +   2./180. * f(m,n-3) \
           ) * inv_delta_rr )

#elif CFDS_ORDER == 8 // Eighth order central finite difference scheme

    #define GF_r(f,m,n) \
        ( ( - 1./280. * f(m,n+4) + 4./105. * f(m,n+3) - 1./5. * f(m,n+2) \
            +   4./5. * f(m,n+1) -   4./5. * f(m,n-1) + 1./5. * f(m,n-2) \
            - 4./105. * f(m,n-3) + 1./280. * f(m,n-4) \
           ) * inv_delta_r )

    #define GF_rr(f,m,n) \
        ( ( - 1./560. * f(m,n+4) +  8./135. * f(m,n+3) -   1./5. * f(m,n+2) \
            +   8./5. * f(m,n+1) - 205./72. * f(m,n)   +   8./5. * f(m,n-1) \
            -   1./5. * f(m,n-2) +  8./315. * f(m,n-3) - 1./560. * f(m,n-4) \
           ) * inv_delta_rr )

#else
    #error "CFDS_ORDER must be either 2, 4, 6, or 8"
#endif

/////////////////////////////////////////////////////////////////////////////////////////

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
/** @defgroup g11b Upwind differencing scheme for convection terms.
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
/////////////////////////////////////////////////////////////////////////////////////////
// Macros to emit the convective spatial derivatives (upwind stencils)
//
// Convention. If f() denotes a function, then f_r() denotes a spatial derivative,
// f_t() denotes a time derivative, and f_rr() denotes the 2nd derivative in r.
//

/// @todo Why is _UPWIND=1 by default?

#ifndef _UPWIND
    #define _UPWIND 1
#endif // _UPWIND

#if _UPWIND

    #define UPWIND_ORDER CFDS_ORDER

    #define GF_convr(bf,f,m,n)   ( bf(m,n) > 0 ? bf(m,n) * GF_up_r(f,m,n) \
                                               : bf(m,n) * GF_down_r(f,m,n) )
#else

    #define GF_convr(bf,f,m,n)   ( bf(m,n) * GF_r(f,m,n) )

#endif // _UPWIND

#if UPWIND_ORDER == 2 // Second order upwind/downwind

    #define GF_up_r(f,m,n)   ( ((-3.*f(m,n))/2. + 2.*f(m,1 + n) - f(m,2 + n)/2.) * inv_delta_r )

    #define GF_up_rr(f,m,n)  ( ( f(m,n) - 2.*f(m,1 + n) + f(m,2 + n) ) * inv_delta_rr )

    #define GF_down_r(f,m,n)   ( (f(m,-2 + n)/2. - 2.*f(m,-1 + n) + (3.*f(m,n))/2.) * inv_delta_r )

    #define GF_down_rr(f,m,n)  ( ( f(m,-2 + n) - 2.*f(m,-1 + n) + f(m,n) ) * inv_delta_rr )

#elif UPWIND_ORDER == 4 // Fourth order upwind/downwind

    #define GF_up_r(f,m,n)   (( -f(m,-1 + n)/4. - (5.*f(m,n))/6. + (3.*f(m,1 + n))/2. - f(m,2 + n)/2. + f(m,3 + n)/12. ) * inv_delta_r )

    #define GF_up_rr(f,m,n)  ( ( (11.*f(m,-1 + n))/12. - (5.*f(m,n))/3. + f(m,1 + n)/2. + f(m,2 + n)/3. - f(m,3 + n)/12. ) * inv_delta_rr )

    #define GF_down_r(f,m,n)   (( -f(m,-3 + n)/12. + f(m,-2 + n)/2. - (3.*f(m,-1 + n))/2. + (5.*f(m,n))/6. + f(m,1 + n)/4. ) * inv_delta_r )

    #define GF_down_rr(f,m,n)  ( ( -f(m,-3 + n)/12. + f(m,-2 + n)/3. + f(m,-1 + n)/2. - (5.*f(m,n))/3. + (11.*f(m,1 + n))/12. ) * inv_delta_rr )

#elif UPWIND_ORDER == 6 // Sixth order upwind/downwind

    #define GF_up_r(f,m,n)   (( -f(m,-1 + n)/6. - (77.*f(m,n))/60. + (5.*f(m,1 + n))/2. - (5.*f(m,2 + n))/3. + (5.*f(m,3 + n))/6. - f(m,4 + n)/4. + f(m,5 + n)/30. ) * inv_delta_r )

    #define GF_up_rr(f,m,n)  ( ( (137.*f(m,-1 + n))/180. - (49.*f(m,n))/60. - (17.*f(m,1 + n))/12. + (47.*f(m,2 + n))/18. - (19.*f(m,3 + n))/12. + (31.*f(m,4 + n))/60. - (13.*f(m,5 + n))/180. ) * inv_delta_rr )

    #define GF_down_r(f,m,n)   (( -f(m,-5 + n)/30. + f(m,-4 + n)/4. - (5.*f(m,-3 + n))/6. + (5.*f(m,-2 + n))/3. - (5.*f(m,-1 + n))/2. + (77.*f(m,n))/60. + f(m,1 + n)/6. ) * inv_delta_r )

    #define GF_down_rr(f,m,n)  ( ( (-13.*f(m,-5 + n))/180. + (31.*f(m,-4 + n))/60. - (19.*f(m,-3 + n))/12. + (47.*f(m,-2 + n))/18. - (17.*f(m,-1 + n))/12. - (49.*f(m,n))/60. + (137.*f(m,1 + n))/180. ) * inv_delta_rr )

#elif UPWIND_ORDER == 8 // Eighth order upwind/downwind

    #define GF_up_r(f,m,n)   (( -f(m,-1 + n)/8. - (223.*f(m,n))/140. + (7.*f(m,1 + n))/2. - (7.*f(m,2 + n))/2. + (35.*f(m,3 + n))/12. - (7.*f(m,4 + n))/4. + (7.*f(m,5 + n))/10. - f(m,6 + n)/6. + f(m,7 + n)/56. ) * inv_delta_r )

    #define GF_up_rr(f,m,n)  ( ( (363.*f(m,-1 + n))/560. + (8.*f(m,n))/315. - (83.*f(m,1 + n))/20. + (153.*f(m,2 + n))/20. - (529.*f(m,3 + n))/72. + (47.*f(m,4 + n))/10. - (39.*f(m,5 + n))/20. + (599.*f(m,6 + n))/1260. - (29.*f(m,7 + n))/560. ) * inv_delta_rr )

    #define GF_down_r(f,m,n)   (( -f(m,-7 + n)/56. + f(m,-6 + n)/6. - (7.*f(m,-5 + n))/10. + (7.*f(m,-4 + n))/4. - (35.*f(m,-3 + n))/12. + (7.*f(m,-2 + n))/2. - (7.*f(m,-1 + n))/2. + (223.*f(m,n))/140. + f(m,1 + n)/8. ) * inv_delta_r )

    #define GF_down_rr(f,m,n)  ( ( (-29.*f(m,-7 + n))/560. + (599.*f(m,-6 + n))/1260. - (39.*f(m,-5 + n))/20. + (47.*f(m,-4 + n))/10. - (529.*f(m,-3 + n))/72. + (153.*f(m,-2 + n))/20. - (83.*f(m,-1 + n))/20. + (8.*f(m,n))/315. + (363.*f(m,1 + n))/560. ) * inv_delta_rr )

#endif

/////////////////////////////////////////////////////////////////////////////////////////
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
 *  @fixme The order of the Kreiss-Oliger term depends on the order of the TIME
 *         integration. This is not implemented yet. Check [1, p.5].
 *
 */
#if CFDS_ORDER == 2

    #define KreissOligerTerm(f,dt) \
        ( - ( GF(f,m,n-1) - 2* GF(f,m,n) + GF(f,m,n+1) ) * dissip_delta_r * (dt) )

#elif CFDS_ORDER == 4

    #define KreissOligerTerm(f,dt) \
        ( ( GF(f,m,n-2) - 4* GF(f,m,n-1) + 6* GF(f,m,n) - 4* GF(f,m,n+1) + GF(f,m,n+2) \
           ) * dissip_delta_r * (dt) )

#elif CFDS_ORDER == 6

    #define KreissOligerTerm(f,dt) \
        ( - (      GF(f,m,n-3) - 6* GF(f,m,n-2) + 15* GF(f,m,n-1) - 20* GF(f,m,n) \
             + 15* GF(f,m,n+1) - 6* GF(f,m,n+2) +     GF(f,m,n+3) \
            ) * dissip_delta_r * (dt) )

#else

    #define KreissOligerTerm(f,dt) \
        ( (       GF(f,m,n-4) -  8* GF(f,m,n-3) + 28* GF(f,m,n-2) - 56* GF(f,m,n-1) \
            + 70* GF(f,m,n)   - 56* GF(f,m,n+1) + 28* GF(f,m,n+2) -  8* GF(f,m,n+3) \
            +     GF(f,m,n+4) \
            ) * dissip_delta_r * (dt) )

#endif
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _FINITE_DIFFERENCES_H_INCLUDED
