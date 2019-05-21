/**
 *  @file      ChebyshevCoefficients.h
 *  @brief     Loads the values of the Chebyshev polynomials at the collocation points,
 *             the values of the modified Chebyshev coefficients of the regularized
 *             derivatives, and the evolution matrices needed to compute the time
 *             derivatives of the spectral coefficients.
 *  @authors   Francesco Torsello. Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

/**
 *  'ChebyshevCoefficients' reads the values of the Chebyshev polynomials at the
 *  collocation points and the elements of the evolution matrix needed to evolve the
 *  spectral coefficients.
 *  TODO: merge with bispecInput?
 */

class ChebyshevCoefficients
{
    size_t derivative_order;
    size_t expansion_order;
    size_t number_collocations;
    size_t number_chebycoeffs;
    size_t data_size;
    size_t reg_coeff_size;
    size_t ee_matrix_size;

public:

    Real* coeff;
    Real* reg_der_coeff;
    Real* reg_derr_coeff;
    Real* ee_matrix_even;
    Real* ee_matrix_odd;

    // isOK() returns true if the coefficients were successfully read
    // from a file (so the object is correctly instantiated).
    //
    bool isOK () const {
        return data_size != 0 && coeff != nullptr;
    }

    // Methods to access the dimensions
    //
    size_t orders   () const{ return derivative_order; }
    size_t expOrd   () const{ return expansion_order; }
    size_t chebys   () const{ return number_chebycoeffs; }
    size_t colpoints() const{ return number_collocations; }

    // Method to access the coefficients
    //
    Real operator()( size_t der_order, size_t cheby_index, size_t n )
    {
        return coeff[ der_order * (expansion_order + 1) * number_collocations
                    + cheby_index * number_collocations + n ];
    }

    // Method to access the coefficients of the regularized first derivatives
    //
    Real reg_ders_cheby( size_t cheby_index, size_t n )
    {
        return reg_der_coeff[ cheby_index * number_collocations + n ];
    }

    // Method to access the coefficients of the regularized second derivatives
    //
    Real reg_derrs_cheby( size_t cheby_index, size_t n )
    {
        return reg_derr_coeff[ cheby_index * number_collocations + n ];
    }

    // Method to print the coefficients to a file
    //
    void exportChebyCoeffs()
    {
        FILE* outFile;
        outFile = fopen( "chebyshev_coeffs_cpp.dat", "wb" );
        fwrite( coeff, sizeof(Real), data_size, outFile );
        fclose(outFile);
    }

    // Constructor: Reads the coefficients from a file.
    //
    ChebyshevCoefficients( const std::string fileName )
        : derivative_order(4)
        , expansion_order(0)
        , number_collocations(0)
        , number_chebycoeffs(number_collocations)
        , data_size(0)
        , coeff(nullptr)
    {
        FILE* inf = fopen( fileName.c_str(), "rb" );

        if( ! inf ) {
            std::cerr << "err: CC: Cannot open file" << std::endl;
            return;
        }

        // Read magic number
        //
        Real x = 0;
        if( 1 != fread( &x, sizeof(Real), 1, inf ) || x != 201905021128 )
        {
            std::cerr << "err: CC: Magic number wrong" << std::endl;
            fclose( inf );
            return;
        }

        // Read the expansion order
        //
        if( 1 == fread( &x, sizeof(Real), 1, inf ) ) {
            expansion_order     = size_t(x);
            number_collocations = size_t(( x + 1 ) / 2);
            number_chebycoeffs  = number_collocations;
        }
        else {
            std::cerr << "err: CC: Cannot read expansion_order" << std::endl;
            fclose( inf );
            return;
        }

        // Read the values of the Chebyshev polynomials and their derivatives
        //
        data_size = derivative_order * (expansion_order + 1) * number_collocations;
        coeff = new Real[derivative_order * (expansion_order + 1) * number_collocations];

        if ( data_size != fread( coeff, sizeof(Real), data_size, inf ) )
        {
            std::cerr << "err: CC: Cannot read all the Chebyshev coefficients"
                << std::endl;
            delete [] coeff;
            coeff = nullptr;
            fclose( inf );
            return;
        }

        // Read the coefficients of the regularized first derivatives
        //
        reg_coeff_size = ( expansion_order + 1 ) * number_collocations;
        reg_der_coeff = new Real[ ( expansion_order + 1 ) * number_collocations ];

        if ( reg_coeff_size != fread( reg_der_coeff, sizeof(Real), reg_coeff_size, inf ) )
        {
            std::cerr << "err: CC: Cannot read all the Chebyshev coefficients"
                << std::endl;
            delete [] reg_der_coeff;
            reg_der_coeff = nullptr;
            fclose( inf );
            return;
        }

        // Read the coefficients of the regularized second derivatives
        //
        reg_derr_coeff = new Real[ ( expansion_order + 1 ) * number_collocations ];

        if ( reg_coeff_size != fread( reg_derr_coeff, sizeof(Real), reg_coeff_size, inf ) )
        {
            std::cerr << "err: CC: Cannot read all the Chebyshev coefficients"
                << std::endl;
            delete [] reg_derr_coeff;
            reg_derr_coeff = nullptr;
            fclose( inf );
            return;
        }

        // Read the even evolution matrix
        //
        ee_matrix_size = ( expansion_order + 1 ) / 2 * number_collocations;
        ee_matrix_even = new Real[ ( expansion_order + 1 ) / 2 * number_collocations ];

        if( ee_matrix_size != fread( ee_matrix_even, sizeof(Real), ee_matrix_size, inf ) )
        {
            std::cerr << "err: CC: Cannot read all the elements in the even ev. matrix"
                << std::endl;
            delete [] ee_matrix_even;
            ee_matrix_even = nullptr;
            fclose( inf );
            return;
        }

        // Read the odd evolution matrix
        //
        ee_matrix_odd = new Real[ ( expansion_order + 1 ) / 2 * number_collocations ];

        if( ee_matrix_size != fread( ee_matrix_odd, sizeof(Real), ee_matrix_size, inf ) )
        {
            std::cerr << "err: CC: Cannot read all the elements in the odd ev. matrix"
                << std::endl;
            delete [] ee_matrix_odd;
            ee_matrix_odd = nullptr;
            fclose( inf );
            return;
        }

        fclose( inf );
    }

};
