/**
 *  @file      cubicSpline.h
 *  @brief     The natural cubic spline of degree 3 with continuity C2.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _CUBIC_SPLINE_H_INCLUDED
#define _CUBIC_SPLINE_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g11 Numerical Methods                                                    */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
/** CubicSpline encapsulates the natural cubic spline of degree 3 with continuity C2.
 *  (Natural means that the second derivatives of the spline polynomials are set
 *  equal to zero at the endpoints of the interval of interpolation.)
 *
 *  Input: a set of N = k+1 coordinates.
 *  Output: a spline as a set of polynomial pieces.
 *
 *  @see <a href="https://en.wikipedia.org/wiki/Spline_(mathematics)">
 *       Algorithm for computing natural cubic splines </a>
 *  @see <a href="http://mathfaculty.fullerton.edu/mathews/n2003/CubicSplinesMod.html">
 *       Cubic splines module </a>
 *  @see Mathmatica notebook `Natural cubic spline.nb` used for tests
 */
class CubicSpline
{
    Int  k;              //!< Dimension: k+1 points
    VecReal x;           //!< The intervals: x[0]..x[1], x[1]..x[2], ...
    VecReal a, b, c, d;  //!< The polynomial coefficients, i = 0,...,k-1

    VecReal h, alp, l, z, mu; // Temporary

public:

    /** Constructor. Only allocates memory for the arrays.
     */
    CubicSpline( Int N )
        : k( N - 1 )
        , x(N), a(N), b(N), c(N), d(N)
        , h(N), alp(N), l(N), z(N), mu(N)
    {}

    /** Calculates the spline based on the endpoints of the interval of interpolation.
     */
    void initialize( Real pts[][2] )
    {
        for( Int i = 0; i < k + 1; ++i ) {
            x[i] = pts[i][0];
            a[i] = pts[i][1];
        }

        // Calculate the differences

        for( Int i = 0; i < k; ++i ) {
            h[i] = x[i+1] - x[i];
        }

        for( Int i = 1; i < k; ++i ) {
            alp[i] =  3 * ( a[i+1] - a[i] ) / h[i] - 3 * ( a[i] - a[i-1] ) / h[i-1];
        }

        // Solve the tridiagonal linear equation (Crout's factorization)

        l[0] = 1;  mu[0] = z[0] = 0;

        for( Int i = 1; i < k; ++i ) {
            l[i] = 2 * ( x[i+1] - x[i-1] ) - h[i-1] * mu[i-1];
            mu[i] = h[i] / l[i];
            z[i] = ( alp[i] - h[i-1] * z[i-1] ) / l[i];
        }

        // Compute the polynomial coefficients

        l[k] = 1;  z[k] = c[k] = 0;

        for( Int j = k - 1; j >= 0; --j ) {
            c[j] = z[j] - mu[j] * c[j+1];
            b[j] = ( a[j+1] - a[j] ) / h[j] - h[j] * ( c[j+1] + 2 * c[j] ) / 3;
            d[j] = ( c[j+1] - c[j] ) / ( 3 * h[j] );
        }
    }

    /** Dumps the polynomial coefficients
     */
    void dump ()
    {
        std::cout << "Polynomial coefficients:" << std::endl << std::endl;

        for( Int i = 0; i < k; ++i ) {
            std::cout << a[i] << "\t" << b[i] << "\t" << c[i] << "\t" << d[i]
                 << " x-range " << x[i] << " .. " << x[i+1] << std::endl;
        }
    }
    /** Returns a value approximated by the cubic splice.
     */
    Real operator()( Real v )
    {
        // Find the polynomial. The polynomial index i goes from 0 to k-1.
        // The interval where the polynomial is valid is x[i-1]..x[i]
        //
        Int i = 0;
        while( i < k - 1 && v >= x[i+1] ) {
            ++i;
        }

        // The interval has not been found. Use the last polynomial for approximation.
        if ( i >= k - 1 ) {
            i = k - 1;
        }

        v -= x[i]; // Subtract the x-value from the reference value

        return a[i] + v * ( b[i] + v * ( c[i] + v * d[i] ) );
    }
};

void CubicSpline_sanityCheck ()
{
    Real pts[5][2] =
    {
        {0.005,   4.892186615295462}, {0.045, 44.02967953765916},
        {0.105, 102.7359189212047  }, {0.145, 40.74442799816186},
        {0.195,  10.719645882407347}
    };

    CubicSpline ncs( 5 );
    ncs.initialize( pts );
    ncs.dump ();

    for( Real x = 0; x < 0.2001; x += 0.01 ) {
        std::cout << x << "\t" << ncs( x ) << std::endl;
    }
}
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _CUBIC_SPLINE_H_INCLUDED