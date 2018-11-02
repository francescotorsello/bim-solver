/**
 *  @file      bandSol.h
 *  @brief     Solving a linear equation `A * x = b` for a band-diagonal matrix `A`.
 *  @authors   Mikica Kocic, Wilkinson & Reinsch, Press et al
 *  @copyright TODO
 */

#ifndef _BAND_SOL_H_INCLUDED
#define _BAND_SOL_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g11 Numerical Methods                                                    */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
/** BandLUDecomposition implements a method for solving linear equations
 *  `A * x = b` for a band-diagonal matrix `A` using LU decomposition.
 *
 *  `A` is band-diagonal with `m1` rows below the diagonal and `m2` rows above.
 *  The input vector is `x[0..n-1]` and the output vector is `b[0..n-1]`.
 *
 *  The LU factorization refers to the factorization of `A`, with proper row and/or
 *  column orderings or permutations, into two factors -- a lower triangular matrix
 *  `L` and an upper triangular matrix `U` such that `A = L U`.
 *
 *  The array `A[0..n-1][0..m1+m2]` stores `A` as follows: The diagonal elements are in
 *  `A[0..n-1][m1]`. Subdiagonal elements are in `A[j ..n-1][0..m1-1]` with `j > 0`
 *  appropriate to the number of  elements on each subdiagonal. Superdiagonal elements
 *  are in `A[0..j ][m1+1..m1+m2]` with `j < n-1` appropriate to the number of elements
 *  on each superdiagonal.
 *
 *  The implementation is based on the routines `bandet1` and `bansol1` in
 *  Wilkinson, J.H., and Reinsch, C. 1971, Linear Algebra, vol. II of Handbook for
 *  Automatic Computation (New York: Springer), Chapter I/6. The C++ code is based on
 *  `bandec.h` @cite Press:2007nrec, Section 2.4.2 Band-Diagonal Systems (pp. 58--61).
 */
class BandLUDecomposition
{
    Int nn;        //!< The full matrix size.
    Int m1;        //!< The size of subdiagonal.
    Int m2;        //!< The size of superdiagonal.
    MatReal L;     //!< The lower triangular matrix (stored compactly).
    MatReal U;     //!< The upper triangular matrix (stored compactly).
    VecInt indx;   //!< The row permutation index effected by the partial pivoting
    Real d;        //!< `+/-1` whether the number of row interchanges was even or odd

public:

    /** The constructor which does the LU decomposition. It takes as arguments the
     *  compactly stored matrix `A`, and the integers `m1 >= 0` (subdiagonal size)
     *  and `m2 >= 0` (superdiagonal size). The upper and lower triangular matrices
     *  are stored in `U` and `L`, respectively.
     */
    BandLUDecomposition( MatReal_I& A, const Int mm1, const Int mm2 );

    /** Gets the lower triangular matrix.
     */
    MatReal& getL () {
        return L;
    }

    /** Gets the upper triangular matrix.
     */
    MatReal& getU () {
        return U;
    }

    /** Given a right-hand side vector `b[0..n-1]`, solves the band-diagonal
     *  linear equation `A * x = b`.
     */
    void solve( VecReal_I& b, VecReal_O& x );

    /** Performs a sanity check test with sample data.
     */
    static void sanityCheck ();

    /** Returns the determinant of the matrix A.
     */
    Real det () const
    {
        Real dd = d;
        for( Int i = 0; i < nn; ++i ) {
            dd *= U[i][0];
        }
        return dd;
    }
};

/** Given an `n*n` band-diagonal matrix `A` with `m1` subdiagonal rows and `m2`
 *  superdiagonal rows, compactly stored in the array `A[0..n-1][0..m1+m2]`,
 *  an LU decomposition of a rowwise permutation of A is constructed. The  upper
 *  and lower triangular matrices are stored in `U` and `L`, respectively.
 *  The stored vector `indx[0..n-1]` records the row permutation effected by the partial
 *  pivoting; `d` is `+/-1` depending on whether the number of row interchanges was even
 *  or odd, respectively.
 *
 *  Some pivoting is possible within the storage limitations, and the routine
 *  does take advantage of the opportunity. Also, when `TINY` is returned as
 *  a diagonal element of `U`, then the original matrix (perhaps as modified by roundoff
 *  error) is in fact singular.
 */
BandLUDecomposition::BandLUDecomposition(
    MatReal_I& A,   //!< Compactly stored diagonal-band matrix
    const Int mm1,  //!< Size of subdiagonal
    const Int mm2   //!< Size of superdiagonal
)
    : nn(A.nrows()), m1(mm1), m2(mm2), L(nn,m1), U(A), indx(nn)
{
    const Real TINY = 1e-40;

    Int mm = m1 + m2 + 1;
    Int l = m1;

    for( Int i = 0; i < m1; ++i ) {
        for( Int j = m1 - i; j < mm; ++j ) {
            U[i][j-l] = U[i][j];
        }
        l--;
        for( Int j = mm - l - 1; j < mm; ++j ) {
            U[i][j] = 0;
        }
    }

    Real d = 1;
    l = m1;

    for( Int k = 0; k < nn; ++k )
    {
        Real dum = U[k][0];
        Int i = k;
        if( l < nn ) {
            l++;
        }
        for( Int j = k + 1; j < l; ++j ) {
            if( abs( U[j][0] ) > abs( dum ) ) {
                dum = U[j][0];
                i = j;
            }
        }
        indx[k] = i + 1;
        if( dum == 0 ) {
            U[k][0]=TINY;
        }
        if (i != k) {
            d = -d;
            for( Int j = 0; j < mm; ++j ) {
                swap( U[k][j], U[i][j] );
            }
        }
        for( i = k + 1; i < l; ++i )
        {
            dum = U[i][0] / U[k][0];
            L[k][i-k-1] = dum;
            for( Int j = 1; j < mm; ++j ) {
                U[i][j-1] = U[i][j] - dum * U[k][j];
            }
            U[i][mm-1] = 0;
        }
    }
}

/** Given a right-hand side vector `b[0..n-1]`, solves the band-diagonal
 *  linear equation `A * x = b`.
 */
void BandLUDecomposition::solve( VecReal_I& b, VecReal_O& x )
{
    Int mm = m1 + m2 + 1;
    Int l  = m1;

    for( Int k = 0; k < nn; ++k ) {
        x[k] = b[k];
    }

    for( Int k = 0; k < nn; ++k )
    {
        Int j = indx[k] - 1;
        if( j != k ) {
            swap( x[k], x[j] );
        }
        if( l < nn ) {
            l++;
        }
        for( j = k + 1; j < l; ++j ) {
            x[j] -= L[k][j-k-1] * x[k];
        }
    }

    l = 1;

    for( Int i = nn - 1; i >= 0; --i )
    {
        Real dum = x[i];
        for( Int k = 1; k < l; ++k ) {
            dum -= U[i][k]*x[k+i];
        }

        x[i] = dum / U[i][0];
        if( l < mm ) {
            l++;
        }
    }
}

/** Sanity check for BandLUDecomposition.
 */
void BandLUDecomposition::sanityCheck ()
{
    const Real A_values[] = {
    // -2  -1   0  +1  <--- the matrix diagonal lies at [0]
        0,  0,  3,  1,
        0,  4,  1,  5,
        9,  2,  6,  5,
        3,  5,  8,  9,
        7,  9,  3,  2,
        3,  8,  4,  6,
        2,  4,  4,  0
    };
    const Real b_values[] = {
        7, 6, 5, 4, 3, 2, 1
    };

    MatReal_I A( 7, 4, A_values );
    VecReal_I b( 7, b_values );
    VecReal_O x( 7 );

    BandLUDecomposition( A, 2, 1 ).solve( b, x );

    for( Int i = 0; i < 7; ++i ) {
        std::cout << x[i] << ", ";
    }

    std::cout << std::endl << std::endl << "Compare to:" << std::endl << std::endl
         << "4.66912, -7.00737, -1.13382, -3.24088, 6.29092, 10.616, -13.5114"
         << std::endl << std::endl;
}
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _BAND_SOL_H_INCLUDED