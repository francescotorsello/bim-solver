/**
 *  @file      SavitzkyGolayFilter.h
 *  @brief     Implementation of the algorithm to determine the Savitzky-Golay filter coefficients.
 *  @authors   Francesco Torsello, Mikica Kocic, William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
 *  @copyright GNU General Public License (GPLv3).
 */

/**
 *
 *  References: @cite Press:2007num (William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery, "Numerical Recipes 3rd Edition: The Art of Scientific Computing", Cambridge University Press, 2007, ISBN: 9780521880688)
 *
 *  See also:
 *      <URL: https://github.com/blackstonep/Numerical-Recipes/blob/master/savgol.h>.
 */

/*void savgol( VecReal_O &c, const Int np, const Int nl, const Int nr,
	const Int ld, const Int polOrder )
{
	Int j, k, imj, ipj, kk, mm;
	Real fac, sum;

	if ( np < nl + nr + 1 || nl < 0 || nr < 0 || ld > polOrder || nl + nr < polOrder )

		throw("Bad arguments in savgol");

	MatReal a( polOrder + 1, polOrder + 1 );
	VecReal b( polOrder + 1 );

	for ( ipj = 0; ipj <= (polOrder << 1); ++ipj ) /// Set up the normal equations of the desired least-squared fit
    {
		sum = (ipj ? 0.0 : 1.0);

		for ( k=1; k<=nr; ++k ) sum += pow(Real(k),Real(ipj));

		for ( k=1; k<=nl; ++k ) sum += pow(Real(-k),Real(ipj));

		mm = std::min( ipj, 2 * polOrder - ipj );
		for ( imj = -mm; imj <= mm; imj += 2 ) a[ (ipj+imj) / 2 ][ (ipj-imj) / 2 ] = sum;
	}

	LUdcmp alud(a);                             /// Solve them: LU decomposition

	for ( j=0; j< polOrder + 1; ++j ) b[j]=0.0;

	b[ld]=1.0;

	/// Right-hand side vector is unit vector, depending on which derivative we want.

	alud.solve(b,b);                            /// Get one row of the inverse matrix

	for ( kk = 0; kk < np ; kk++ ) c[kk]=0.0;   /// Zero the output array (it may be bigger
                                                /// than the number of coefficients)

	for (k = -nl;k<=nr;k++) {

		sum = b[0];                             /// Each Savitzky-Golay coefficient is the
                                                /// dot product of powers of an integer with
                                                /// the inverse matrix row.
		fac = 1.0;

		for ( mm = 1; mm<=polOrder; ++m ) sum += b[mm]*( fac *= k );

		kk = ( np - k ) % np;                   /// Stoer in wraparound order.

		c[kk] = sum;
	}

}*/

/** SavGolFilter contains the algorithm to perform the Savitzy-Golay filtering of the grid
 *  functions.
 */
class SavGolFilter : GridUser
{

    VecReal_O coeff;          //!< The array to store the Savitzky-Golay coefficients

public:
    /////////////////////////////////////////////////////////////////////////////////////

    SavGolFilter( UniformGrid& ug, const Int np, const Int nl, const Int nr, const Int ld, const Int polOrder );

    void convlv(VecReal_I &data, VecReal_I &respns, const Int isign,
	VecReal_O &ans);

};

/** Implementation of the algorithm to determine the Savitzky-Golay filter coefficients.
 *
 *  The routine savgol returns in c[0..np-1], in wraparound order (N.B.!) consistent with
 *  the argument respns in routine convlv, a set of Savitzky-Golay filter coefficients.
 *
 *  nl is the number of leftward (past) data points used, while nr is the number of
 *  rightward (future) data points, making the total number of data points used nl +nr +1.
 *  ld is the order of the derivative desired (e.g., ld D 0 for smoothed function. For the
 *  derivative of order k, you must multiply the array c by k!.)
 *
 *  polOrder is the order of the smoothing polynomial, also equal to the
 *  highestnconserved moment; usual values are polOrder = 2 or polOrder = 4.
 */

SavGolFilter::SavGolFilter( UniformGrid& ug, const Int np, const Int nl, const Int nr, const Int ld, const Int polOrder ) : GridUser( ug )
    {
        Int j, k, imj, ipj, kk, mm;
        Real fac, sum;

        if ( np < nl + nr + 1 || nl < 0 || nr < 0 || ld > polOrder || nl + nr < polOrder )

            throw("Bad arguments in savgol");

        MatReal a( polOrder + 1, polOrder + 1 );
        VecReal b( polOrder + 1 );

        for ( ipj = 0; ipj <= (polOrder << 1); ++ipj ) /// Set up the normal equations of
                                                       /// the desired least-squared fit
        {
            sum = (ipj ? 0.0 : 1.0);

            for ( k = 1; k <= nr; ++k ) sum += pow( Real( k ), Real( ipj ) );

            for ( k = 1; k <= nl; ++k ) sum += pow( Real( -k ), Real( ipj ) );

            mm = std::min( ipj, 2 * polOrder - ipj );

            for ( imj = -mm; imj <= mm; imj += 2 ) a[ ( ipj + imj ) / 2 ][ ( ipj - imj ) / 2 ] = sum;
        }

        LUdcmp alud(a);                         /// Solve them: LU decomposition

        for ( j=0; j< polOrder + 1; ++j ) b[j]=0.0;

        b[ld]=1.0;

        /// Right-hand side vector is unit vector, depending on which derivative we want.

        alud.solve(b,b);                        /// Get one row of the inverse matrix

        for ( kk = 0; kk < np ; kk++ )
        {
            coeff[kk] = 0.0;                    /// Zero the output array (it may be bigger
        }                                       /// than the number of coefficients)

        for ( k = -nl; k <= nr; ++k )
        {

            sum = b[0];                         /// Each Savitzky-Golay coefficient is the
                                                /// dot product of powers of an integer with
                                                /// the inverse matrix row.
            fac = 1.0;

            for ( mm = 1; mm <= polOrder; ++mm ) sum += b[ mm ] * ( fac *= k );

            kk = ( np - k ) % np;               /// Store in wraparound order.

            coeff[kk] = sum;

        }

    }

/** Convolves or deconvolves a real data set data[0..n-1] (including any user-supplied zero
 *  padding) with a response function respns[0..m-1], where m is an odd integer <= n. The
 *  response function must be stored in wraparound order: The first half of the array respns
 *  contains the impulse response function at positive times, while the second half of the
 *  array contains the impulse response function at negative times, counting down from the
 *  highest element respns[m-1]. On input isign is +1 for convolution, -1 for deconvolution.
 *  The answer is returned in ans[0..n-1]. n must be an integer power of 2.
 */

 /*
void SavGolFilter::convlv( VecReal_I &data, VecReal_I &respns, const Int isign,
	VecReal_O &ans )
	{
	Int i, no2, n = data.size(), m = respns.size();

	Real mag2, tmp;

	VecReal temp( n );

	temp[0] = respns[0];

	for ( i = 1; i < (m+1)/2; ++i )              /// Put rspns in array of length n.
    {
		temp[i] = respns[i];
		temp[n-i] = respns[m-i];
	}

	for ( i = (m+1)/2; i < n-(m-1)/2; ++i )      /// Pad with zeros
		temp[i] = 0.0;

	for ( i=0; i < n; ++i )
		ans[i] = data[i];

	realft( ans, 1 );                            /// FFT both arrays.
	realft( temp, 1 );

	no2 = n >> 1;

	if ( isign == 1)
    {
		for ( i=2; i < n; i+=2 )
		{
			tmp=ans[i];

			ans[i]=(ans[i]*temp[i]-ans[i+1]*temp[i+1])/no2;

			ans[i+1]=(ans[i+1]*temp[i]+tmp*temp[i+1])/no2;
		}

		ans[0]=ans[0]*temp[0]/no2;
		ans[1]=ans[1]*temp[1]/no2;

	} else if (isign == -1)
	{
		for (i=2;i<n;i+=2)
        {
			if ((mag2=pow2(temp[i])+pow2(temp[i+1])) == 0.0)
				throw("Deconvolving at response zero in convlv");

			tmp=ans[i];

			ans[i]=(ans[i]*temp[i]+ans[i+1]*temp[i+1])/mag2/no2;

			ans[i+1]=(ans[i+1]*temp[i]-tmp*temp[i+1])/mag2/no2;

		}
		if (temp[0] == 0.0 || temp[1] == 0.0)
			throw("Deconvolving at response zero in convlv");

		ans[0]=ans[0]/temp[0]/no2;
		ans[1]=ans[1]/temp[1]/no2;

	} else throw("No meaning for isign in convlv");

	realft(ans,-1);
}
*/
