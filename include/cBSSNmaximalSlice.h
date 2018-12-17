/**
 *  @file      maximalSlice.h
 *  @brief     Maximal slicing algorithm in the covariant BSSN formalism.
 *  @author    Mikica Kocic, Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _MAXIMAL_SLICE_H_INCLUDED
#define _MAXIMAL_SLICE_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////////////////

//#define W(m,n) -eq_gW(m,n) / ( TINY_Real +  eq_fW(m,n))

/** PP is P(x) in the equation: y''(x) = P(x) y'(x) + Q(x) y(x)

#define PP(n) ( gDA(m,n) - 2 * (gDB(m,n) + gconf(m,n) * gDconf(m,n)) - 2 / r(m,n) )

 QQ is Q(x) in the equation: y''(x) = P(x) y'(x) + Q(x) y(x)

#define QQ(n) exp(4 * gconf(m,n)) * (pow2(gK1(m,n)) + 2* pow2(gK2(m,n)))* pow2(gA(m,n)) + k_g * ((exp(4 * gconf(m,n)) * pow2 (gA(m,n)) * (eq_pf_gJ11(m,n) + 2 * eq_pf_gJ22(m,n) + eq_pf_grho(m,n) + 2 * P_1_0(R(m,n)))) / 2. - exp(2 * fconf(m,n) + 2 * gconf(m,n)) * fA(m,n) * gA(m,n) *Lt(m,n) * P_2_1(R(m,n)) + (exp(2 * fconf(m,n) + 2 * gconf(m,n)) * fA(m,n) * gA(m,n) * (2 * P_1_1(R(m,n)) + P_2_1(R(m,n)))) / (2. * Lt(m,n)) + (-(exp(2 * fconf (m,n) + 2 * gconf(m,n)) * fA(m,n) * gA(m,n) * P_1_2(R(m,n))) - (exp (4 * gconf(m,n))* pow2(gA(m,n)) * (2 * P_1_1(R(m,n)) + P_2_1(R(m,n)))) / (2. * Lt(m,n))) * W(m,n))

 QQ for GR

#define QQGR(n) exp(4 * gconf(m,n)) * ( pow2(gK1(m,n)) + 2 * pow2(gK2(m,n)))* pow2(gA(m,n)) + k_g * ((exp(4 * gconf(m,n)) * pow2 (gA(m,n)) * (eq_pf_gJ11(m,n) + 2 * eq_pf_gJ22(m,n) + eq_pf_grho(m,n))) / 2.
*/

/// @todo fixme: Should be  k_g/2 in the above. Removing 1/2, however, speeds up the collapse.
/// @todo fixme: Current implementation uses rho + J of the perfect fluid (should use bimetric).

/** Finds the maximal slice using a second-order accurate finite difference scheme.
 *
 *  The maximal slicing involves solving a differential equation:
 *  <pre>
 *        y''(x) = P(x) y'(x) + Q(x) y(x)   </pre>
 *
 *  with the Neumann BC on the left `y�(0) = 0`,
 *  and the Dirichlet BC on the right `y(r_1) = 1`.
 *
 *  The solution is found using a linear finite-difference method.
 *  The algorithm comprises two major steps:
 *
 *    - Prepare the band-diagonal linear system corresponding to a discretized DE.
 *
 *    - Solve the linear system (using Crout's method).
 *
 *  The coefficients for the discretized DE are generated in the Mathematica notebook:
 *  `Maximal slicing, FD method.nb`. The linear system is solved used Bandec.
 */

/* Real BimetricEvolve::W( Int m, Int n )
{
    if( isGR () )
        {
            return 0;

        } else {

            return  -eq_gW(m,n) / ( TINY_Real +  eq_fW(m,n));

        }
}*/

Real BimetricEvolve::PP( Int m, Int n )
{
    return gDA(m,n) - 2 * (gDB(m,n) + gconf(m,n) * gDconf(m,n)) - 2 / \
r(m,n);
}

Real BimetricEvolve::QQ( Int m, Int n )
    {
        if( isGR () )
        {
            return (k_g * exp(4 * gconf(m,n)) * (eq_pf_gJ11(m,n) + 2 * eq_pf_gJ22(m,n)                   + eq_pf_grho(m,n))
            /*1*/  * pow2(gA(m,n))) / 2. + exp(4 * gconf(m,n)) * (pow2 (                    gK1(m,n))
            /*2*/  + 2 * pow2(gK2(m,n))) * pow2(gA(m,n));
        } else {

            return exp(4 * gconf(m,n)) * (pow2(gK1(m,n)) + 2 * \
                    pow2(gK2(m,n)))
        /*1*/   * pow2(gA(m,n)) + k_g * ((exp(4 * gconf(m,n)) * \
                pow2(gA(m,n))
        /*2*/  * (eq_pf_gJ11(m,n) + 2 * eq_pf_gJ22(m,n) + eq_pf_grho(m,n) + 2 * \
                P_1_0(R(m,n)))) / 2. - exp(2
        /*2*/  * fconf(m,n) + 2 * gconf(m,n)) * fA(m,n) * gA(m,n) * \
                Lt(m,n) * P_2_1(R(m,n))
        /*1*/  + (exp(2 * fconf(m,n) + 2 * gconf(m,n)) * fA(m,n) * \
                gA(m,n) * (2
        /*3*/  * P_1_1(R(m,n)) + P_2_1(R(m,n)))) / (2. * Lt(m,n)) + \
                (-(exp(2 * fconf(m,n)
        /*4*/  + 2 * gconf(m,n)) * fA(m,n) * gA(m,n) * P_1_2(R(m,n))) - \
                (exp(4 * gconf(m,n))
        /*3*/  * pow2(gA(m,n)) * (2 * P_1_1(R(m,n)) + P_2_1(R(m,n)))) / \
                (2. * Lt(m,n)))
        /*1*/  * (-eq_gW(m,n) / ( TINY_Real +  eq_fW(m,n))));

        }
    }


void BimetricEvolve::maximalSlice_2
(
    Int m,          //!< The time slice
    Int N,          //!< The radial coordinate of the right boundary
    Real gAlp_at_N  //!< The right boundary condition
    )
{
    const Int offset = nGhost - 1;  N = N + 1;
    const Real h = delta_r;
    MatReal A( N, 3 );
    VecReal b( N );

    /////////////////////////////////////////////////////////////////////////////////////

    OMP_parallel_for( Int i = 0; i < N; ++i )
    {
        A[i][0] = -2         - h * PP( m, offset + i );
        A[i][1] =  4 + 2 * h * h * QQ( m, offset + i );
        A[i][2] = -2         + h * PP( m, offset + i );
        b[i]    =  0;
    }

    // Override the default row values

    A[0][0] = 0;
    A[0][2] = -4;

    A[N-1][2] = 0;
    b[N-1] = gAlp_at_N * ( 2 - h * PP( m, offset + (N-1) ) ) ;

    /////////////////////////////////////////////////////////////////////////////////////

    VecReal_O x( N );
    BandLUDecomposition( A, 1, 1 ).solve( b, x );

    // Retrieve the results from x

    gAlp( m, offset + N ) = gAlp_at_N;

    for( Int i = 0; i < N; ++i ) {
        gAlp( m, offset + i ) = x[i];
    }

    maximalSlice_PostSteps( m, N );
}

/** Finds the maximal slice using a fourth-order accurate finite difference scheme.
 *  @see maximalSlice_2
 */
void BimetricEvolve::maximalSlice_4
(
    Int m,          //!< The time slice
    Int N,          //!< The radial coordinate of the right boundary
    Real gAlp_at_N  //!< The right boundary condition
    )
{
    const Int offset = nGhost;
    const Real h = delta_r;
    MatReal A( N, 7 );
    VecReal b( N );

    /////////////////////////////////////////////////////////////////////////////////////

    OMP_parallel_for( Int i = 0; i < N; ++i )
    {
        A[i][0] =   0;
        A[i][1] =  -1          - h * PP( m, offset + i );
        A[i][2] =  16      + 8 * h * PP( m, offset + i );
        A[i][3] = -30 - 12 * h * h * QQ( m, offset + i );
        A[i][4] =  16      - 8 * h * PP( m, offset + i );
        A[i][5] =  -1        +   h * PP( m, offset + i );
        A[i][6] =   0;
        b[i]    =   0;
    }

    // Override the default row values

    Int i = 0;
    A[i][1] =   0;
    A[i][2] =   0;
    A[i][3] = -20 + 10 * h * PP( m, offset + i ) - 12 * h * h * QQ( m, offset + i );
    A[i][4] =  17 - 15 * h * PP( m, offset + i );
    A[i][5] =   4 +  6 * h * PP( m, offset + i );
    A[i][6] =  -1 -      h * PP( m, offset + i );

    i = 1;
    A[i][1] = 0;
    A[i][3] = -31 - h * PP( m, offset + i ) - 12 * h * h * QQ( m, offset + i );

    i = N - 2;
    A[i][5] = 0;
    b[i] = gAlp_at_N * ( 1 - 8 * h * PP( m, offset + i ) ) ;

    i = N - 1;
    A[i][0] =  -1 +      h * PP( m, offset + i );
    A[i][1] =   4 -  6 * h * PP( m, offset + i );
    A[i][2] =   6 + 18 * h * PP( m, offset + i );
    A[i][3] = -20 - 10 * h * PP( m, offset + i ) - 12 * h * h * QQ( m, offset + i );
    A[i][4] =   0;
    A[i][5] =   0;
    b[i] = gAlp_at_N * ( -11 + 3 * h * PP( m, offset + i ) ) ;

    /////////////////////////////////////////////////////////////////////////////////////

    VecReal_O x( N );
    BandLUDecomposition( A, 3, 3 ).solve( b, x );

    // Retrieve the results from x

    gAlp( m, offset + N ) = gAlp_at_N;

    for( Int i = 0; i < N; ++i ) {
        gAlp( m, offset + i ) = x[i];
    }

    maximalSlice_PostSteps( m, N );
}

/** Finds the maximal slice using a second-order accurate finite difference scheme.
 *
 *  This is a slightly optimized version where the differential equation is solved
 *  using Algorithm 11.3 from Burden & Faires, Numerical Analysis, 2016
 *  (found in the section Finite-Difference Methods for Linear Problems on p. 700).
 *
 *  @see maximalSlice_2
 */
void BimetricEvolve::maximalSlice_2opt
(
    Int m,          //!< The time slice
    Int N,          //!< The radial coordinate of the right boundary
    Real gAlp_at_N  //!< The right boundary condition
    )
{
    const Int offset = nGhost - 1; N = N + 1;
    const Real h = delta_r;

    // First, allocate memory used to store the intermediate variables.
    //
    static Real* workArea = NULL;
    static Real *a, *b, *c, *d; // Used in steps 1-3
    static Real *l, *u, *z, *w; // Used in steps 4-8
    if ( workArea == NULL ) {
        size_t chunkSize = nLen + 2 * nGhost;
        workArea = new Real[ 8 * chunkSize ];
        a = workArea;      b = a + chunkSize; c = b + chunkSize; d = c + chunkSize;
        l = d + chunkSize; u = l + chunkSize; z = u + chunkSize; w = z + chunkSize;
    }

    // Note: a[] = diagonal, b[] = upper band, c[] = lower band,
    //       and d[] is on the right hand side of the equation A x = d

    /////////////////////////////////////////////////////////////////////////////////////
    // Prepare the tridiagonal linear system corresponding to DE.

    // Steps 1-3

    OMP_parallel_for( Int i = 0; i < N; ++i )
    {
        c[i] = -2         - h * PP( m, offset + i );
        a[i] =  4 + 2 * h * h * QQ( m, offset + i );
        b[i] = -2         + h * PP( m, offset + i );
        d[i] =  0;
    }

    c[0] =  0;
    b[0] = -4;  // Set by Neumann BC: w{-1} = w_{1}

    b[N-1] =   0;
    d[N-1] = ( 2 - h * PP( m, offset + N-1 ) ) * gAlp_at_N;

    /////////////////////////////////////////////////////////////////////////////////////
    // Solve the linear system using Crout's factorization.

    // Step 4
    //
    l[0] = a[0];
    u[0] = b[0] / a[0];
    z[0] = d[0] / l[0];

    // Step 5
    //
    for( Int i = 0; i < N - 1; ++i )
    {
        l[i] = a[i] - c[i] * u[i-1];
        u[i] = b[i] / l[i];
        z[i] = ( d[i] - c[i] * z[i-1] ) / l[i];
    }

    // Step 6-7
    //
    Int i = N - 1;
    l[i] = a[i] - c[i] * u[i-1];
    gAlp( m, offset + i ) = w[i] = z[i] = ( d[i] - c[i] * z[i-1] ) / l[i];

    // Step 8
    //
    for( Int i = N - 2; i >= 1; --i ) {
        gAlp( m, offset + i ) = w[i] = z[i] - u[i] * w[i+1];
    }

    gAlp( m, offset + N ) = w[N] = gAlp_at_N;

    maximalSlice_PostSteps( m, N );
}

/** Additional steps for gAlp and gDAlp (like fixing boundary conditions.)
 */
void BimetricEvolve::maximalSlice_PostSteps( Int m, Int N )
{
    // Fix the ghost cells on the left (using parity BC)
    //
    /*for( Int i = 0; i < nGhost; ++i ) {
        gAlp( m, nGhost - 1 - i  ) = gAlp( m, nGhost + i );
    }*/
    for( Int i = 0; i < nGhost +1; ++i )
    {
        Int n  = nGhost - i - 1;
        Int nR = nGhost + i;

        gAlp (m,n) = gAlp  (m,nR);
        fAlp (m,n) = fAlp  (m,nR);
        //gDAlp(m,n) = -gDAlp(m,nR);
        //fDAlp(m,n) = -fDAlp(m,nR);

    }

    //smoothenGF0 ( m, nGhost, nSmoothUpTo, 32,  fld::gAlp,  fld::tmp,  fld::gAlp,  1 );

    // Determine cells on the right by solving the differential equation.
    //
    Real h = delta_r;
    Int overlap = 10;
    for( Int n = nGhost + N - overlap; n < nGhost + nLen + nGhost - 1; ++n )
    {
        /*Real dbg0  = W(m,n);
        Real dbg1  = eq_gW(m,n);
        Real dbg2  = eq_fW(m,n);

        Real dbg3  = gA1(m,n);
        Real dbg4  = fA1(m,n);
        Real dbg5  = gconf(m,n);
        Real dbg6  = fconf(m,n);
        Real dbg7  = gDB(m,n);
        Real dbg8  = fDB(m,n);
        Real dbg9  = gDconf(m,n);
        Real dbg10 = fDconf(m,n);
        Real dbg11 = eq_pf_gJ11(m,n);
        Real dbg12 = eq_pf_fJ11(m,n);
        Real dbg13 = eq_pf_gJ22(m,n);
        Real dbg14 = eq_pf_fJ22(m,n);
        Real dbg15 = eq_pf_gj(m,n);
        Real dbg16 = eq_pf_fj(m,n);
        Real dbg17 = eq_pf_grho(m,n);
        Real dbg18 = eq_pf_frho(m,n);
        Real dbg19 = gA(m,n);
        Real dbg20 = fA(m,n);
        Real dbg21 = gB(m,n);
        Real dbg22 = fB(m,n);
        Real dbg23 = Lt(m,n);
        Real dbg24 = Lt2(m,n);
        Real dbg25 = P_0_2(R(m,n));
        Real dbg26 = P_0_3(R(m,n));
        Real dbg27 = P_1_0(R(m,n));
        Real dbg28 = P_1_1(R(m,n));
        Real dbg29 = P_1_2(R(m,n));
        Real dbg30 = P_1_3(R(m,n));
        Real dbg31 = P_2_0(R(m,n));
        Real dbg32 = P_2_1(R(m,n));
        Real dbg33 = P_2_2(R(m,n));
        Real dbg34 = p(m,n);
        Real dbg35 = R(m,n);
        Real dbg36 = gtrK(m,n);
        Real dbg37 = ftrK(m,n);

        Real dbg38 = gA2_r(m,n);
        Real dbg39 = fA1_r(m,n);
        Real dbg40 = p_r(m,n);
        Real dbg41 = gtrK_r(m,n);
        Real dbg42 = ftrK_r(m,n);

        Real dbg49 = 3 * gA1(m,n) * (1 + gDB(m,n)
      * r(m,n) + 2 * gconf(m,n) * gDconf(m,n) * r(m,n));

        Real dbg50 = - r(m,n) * (3 * gA2_r(m,n)
      + gtrK_r(m,n));*/

        gAlp( m, n+1 ) =
            ( ( 4  + 2 * h * h * QQ( m,n) ) * gAlp(m,n) - ( 2 + h * PP( m,n) ) * gAlp(m,n-1) )
            / ( 2 - h * PP( m,n) );
    }

    // Fix the left boundary for gAlp using parity BC
    //

    // Calculate gDAlp
    //
    OMP_parallel_for( Int n = 0; n < nGhost + 1; ++n )
    {
        gAlp_r ( m, n ) = GF_right_r( gAlp,  m, n );
        gDAlp(m,n) = gAlp_r(m,n) /*/ (TINY_Real + gAlp(m,n))*/;

    }
    OMP_parallel_for( Int n = nGhost + 1; n < nGhost + nLen + 1; ++n )
    {
        gAlp_r ( m, n ) = GF_up_r( gAlp,  m, n );
        gDAlp(m,n) = gAlp_r(m,n) /*/ (TINY_Real + gAlp(m,n))*/;

    }
    OMP_parallel_for( Int n = nGhost + nLen + 1; n < 2*nGhost + nLen + 1; ++n )
    {
        gAlp_r ( m, n ) = GF_left_r( gAlp,  m, n );
        gDAlp(m,n) = gAlp_r(m,n) /*/ (TINY_Real + gAlp(m,n))*/;

    }

    // The accuracy is very low at low r. Assume that a few first cells are linear,
    // then apply a six-point cubic spline to smooth the region around r = 0.
    //
    cubicSplineSmooth( m, fld::gDAlp, lin2n, cub2n );
    //smoothenGF0 ( m, nGhost, nSmoothUpTo, 32,  fld::gDAlp,  fld::tmp,  fld::gDAlp,  -1 );

    // Fix the left boundary for gDAlp using parity BC
    //
    for( Int i = 0; i < nGhost +1; ++i )
    {
        Int n  = nGhost - i - 1;
        Int nR = nGhost + i;

        gDAlp(m,n) = -gDAlp(m,nR);
        fDAlp(m,n) = -fDAlp(m,nR);

    }
}

#undef PP
#undef QQ

#endif // _MAXIMAL_SLICE_H_INCLUDED
