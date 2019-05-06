/**
 *  @file      bispec.cpp
 *  @brief     The covariant BSSN evolution for spherically symmetric bimetric spacetimes, with the pseudospectral method.
 *  @authors   Mikica Kocic, Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <vector>

#ifndef OBSERVER
    #define OBSERVER 1
#endif // OBSERVER

#ifndef _EVOLVE_DSIG
    #define _EVOLVE_DSIG 0
#endif // _EVOLVE_DSIG

#ifndef _DETECT_NAN
    #define _DETECT_NAN 1
#endif // _DETECT_NAN

/** The following vector contains the ORDERED fields whose spectral coefficients
    have to be set equal to the values read by the class bispecInput.
    This order MUST match the order of exportation in the Mathematica file.
  */

#define FIELDS gconf, fconf, trgK, trfK, gA, fA, gB, fB, gA1, fA1, gL, fL

#define DERS gconf_r, fconf_r, trgK_r, trfK_r, gA_r, fA_r, gB_r, fB_r, gA1_r, fA1_r, gL_r, fL_r

#define DERRS gconf_rr, fconf_rr, trgK_rr, trfK_rr, gA_rr, fA_rr, gB_rr, fB_rr, gA1_rr, fA1_rr, gL_rr, fL_rr

using namespace std;

typedef double Real;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

/** Namespace 'fields' contains the indexing of the fields and their derivatives,
    to be used in the other classes
  */
namespace fields
{
    /// 'fldCheby' contains all the fields that need to be expanded in Chebyshev series
    enum fldCheby { FIELDS };
    static const fldCheby flds[] = { FIELDS }; // this is needed to be able to iterate over the enumeration ( taken from https://stackoverflow.com/questions/261963/how-can-i-iterate-over-an-enum )

    /// 'derCheby' contains all the first radial derivatives that need to be expanded in Chebyshev series
    enum derCheby { DERS };
    static const derCheby ders[] = { DERS };

    /// 'derrCheby' contains all the second radial derivatives that need to be expanded in Chebyshev series
    enum derrCheby { DERRS };
    static const derrCheby derrs[] = { DERRS };

    /// perhaps the vector below is useless
    static const std::vector<int> bispecInput_fields = { FIELDS };

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

class ChebyshevCoefficients
{
    size_t derivative_order;
    size_t chebyshev_index;
    size_t collocation_point;
    size_t data_size;
    size_t ee_matrix_size;
    Real* coeff;
    Real* ee_matrix;

public:

    // isOK() returns true if the coefficients were successfully read
    // from a file (so the object is correctly instantiated).
    //
    bool isOK () const {
        return data_size != 0 && coeff != nullptr;
    }

    // Methods to access the dimensions
    //
    size_t orders () const{ return derivative_order; }
    size_t chebys () const{ return chebyshev_index; }
    size_t points () const{ return collocation_point; }

    // Method to access the coefficients
    //
    Real operator() ( size_t i, size_t j, size_t k )
    {
        return coeff[ i * chebyshev_index * collocation_point
                    + j * collocation_point + k ];
    }

    Real evolutionMatrix ( size_t j, size_t k )
    {
        return ee_matrix[ j * collocation_point + k ];
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
        : derivative_order(0)
        , chebyshev_index(0)
        , collocation_point(0)
        , data_size(0)
        , coeff(nullptr)
    {
        FILE* inf = fopen( fileName.c_str(), "rb" );

        if( ! inf ) {
            std::cerr << "err: Cannot open file" << std::endl;
            return;
        }

        // Read magic number
        //
        Real x = 0;
        if( 1 != fread( &x, sizeof(Real), 1, inf ) || x != 201905021128 )
        {
            std::cerr << "err: Magic number wrong" << std::endl;
            fclose( inf );
            return;
        }

        derivative_order = 4;

        // Read chebyshev_index
        //
        if( 1 == fread( &x, sizeof(Real), 1, inf ) ) {
            chebyshev_index = size_t(x + 1);
            collocation_point = size_t(x + 1);
        }
        else {
            std::cerr << "err: Cannot read expansion_order" << std::endl;
            fclose( inf );
            return;
        }

        // Read the coefficients
        //
        data_size = derivative_order * chebyshev_index * collocation_point;
        coeff = new Real[ derivative_order * chebyshev_index * collocation_point ];

        if ( data_size != fread( coeff, sizeof(Real), data_size, inf ) )
        {
            std::cerr << "err: Cannot read all the Chebyshev coefficients" << std::endl;
            delete [] coeff;
            coeff = nullptr;
            fclose( inf );
            return;
        }

        // Read the evolution matrix
        //
        ee_matrix_size = chebyshev_index * collocation_point;
        ee_matrix = new Real[ chebyshev_index * collocation_point ];

        if ( ee_matrix_size != fread( ee_matrix, sizeof(Real), ee_matrix_size, inf ) )
        {
            std::cerr << "err: Cannot read all the elements in the evolution matrix" << std::endl;
            delete [] ee_matrix;
            ee_matrix = nullptr;
            fclose( inf );
            return;
        }

        fclose( inf );
    }

};

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

class bispecInput
{
    size_t number_fields;
    size_t expansion_order;
    size_t data_size;
    Real* spec_coeff;

public:

    // isOK() returns true if the coefficients were successfully read
    // from a file (so the object is correctly instantiated).
    //
    bool isOK () const {
        return data_size != 0 && spec_coeff != nullptr;
    }

    // Methods to access the dimensions
    //
    size_t n_fields () const{ return number_fields; }
    size_t exp_order () const{ return expansion_order; }
    size_t size_ID () const{ return data_size; }

    // Method to access the initial data
    //
    Real operator() ( size_t field, size_t n )
    {
        return spec_coeff[ ( expansion_order + 1 ) * field + n ];
    }

    // Constructor: Reads the spectral initial data from a file.
    //
    bispecInput( const std::string fileName )
    {
        FILE* specid = fopen( fileName.c_str(), "rb" );

        if( ! specid ) {
            std::cerr << "err: ID: Cannot open file" << std::endl;
            return;
        }

        // Read magic number
        //
        Real x = 0;
        if( 1 != fread( &x, sizeof(Real), 1, specid ) || x != 4478245647 )
        {
            std::cerr << "err: ID: Magic number wrong" << std::endl;
            fclose( specid );
            return;
        }

        // Read the number of evolved fields
        //
        if( 1 == fread( &x, sizeof(Real), 1, specid ) ) {
            number_fields = size_t(x);
        }
        else {
            std::cerr << "err: ID: Cannot read number_fields" << std::endl;
            fclose( specid );
            return;
        }

        // Read the order of Chebyshev expansion
        //
        if( 1 == fread( &x, sizeof(Real), 1, specid ) ) {
            expansion_order = size_t(x);
        }
        else {
            std::cerr << "err: ID: Cannot read expansion_order" << std::endl;
            fclose( specid );
            return;
        }

        // Read the spectral initial data
        //
        data_size = number_fields * ( expansion_order + 1 );
        spec_coeff = new Real[ data_size ];

        if ( data_size != fread( spec_coeff, sizeof(Real), data_size, specid ) )
        {
            std::cerr << "err: ID: Cannot read all the initial data" << std::endl;
            delete [] spec_coeff;
            spec_coeff = nullptr;
            fclose( specid );
            return;
        }

        fclose(specid);
    }

};

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

/** The following class should be inherited from bispecID. However, if I do that, an error occurs...
    The first exp_ord * field + n elements of this array is set equal to the initial data.
  */
class ChebyshevExpansion
{

private:

    /// The size of the cached time steps
    size_t mDim;
    size_t mLen;
    size_t mExtra;
    size_t exp_ord;
    size_t n_flds;
    /// The array containing the spectral coefficients for all the Chebychev series of the fields
    Real *spcoeffs;

public:

    /// Method to access the spectral coefficients. Since we will make use of the integrator from bim-solver, we copy the structure of GF in gridDriver.
    Real specC( size_t m, size_t field, size_t n )
    {
        return spcoeffs[ ( n_flds * ( exp_ord + 1 ) ) * m + ( exp_ord + 1 ) * field + n ];
    }

    /// Constructor: Save the spectral coefficients to their initial values.
    ChebyshevExpansion( bispecInput& bispecID )
    {

        mLen = 5;
        mExtra = 9;
        mDim = mLen + mExtra;
        exp_ord = bispecID.exp_order();
        n_flds = bispecID.n_fields();
        spcoeffs = new Real[ mDim * n_flds * ( exp_ord + 1 ) ];

        //std::cout << std::endl;
        for( size_t field = 0; field < n_flds; ++field )
        {
            //std::cout << "This is the field " << field << std::endl;
            for( size_t n = 0; n < exp_ord + 1; ++n )
            {
                spcoeffs[ ( exp_ord + 1 ) * field + n ] = bispecID( field, n );
                 //std::cout << "(" << field << "," << n << ") : "
                 //     << spcoeffs[ ( exp_ord + 1 ) * field + n ] << " = " << bispecID( field, n ) << std::endl;
            }
        }
        //std::cout << std::endl;

    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

/** TODO: learn how to access the namespace variables and put them as arguments of the functions
    defined here.
  */
class bispecEvolve
{
    size_t m; /// The time step
    size_t n; /// The collocation point index

    Real gconf  ( size_t m, size_t n );
    Real gtrK   ( size_t m, size_t n );
    Real gA     ( size_t m, size_t n );
    Real gB     ( size_t m, size_t n );
    Real gA1    ( size_t m, size_t n );
    Real gL     ( size_t m, size_t n );

    Real fconf  ( size_t m, size_t n );
    Real ftrK   ( size_t m, size_t n );
    Real fA     ( size_t m, size_t n );
    Real fB     ( size_t m, size_t n );
    Real fA1    ( size_t m, size_t n );
    Real fL     ( size_t m, size_t n );

    /** The constructor computs the values of the fields and the evolution equations on the initial hypersurface
     */
    bispecEvolve( bispecInput& bispecID, ChebyshevCoefficients& chebyC, ChebyshevExpansion& chebyExp )
    {
        /// Here we can probably define a macro and run over all the fields by using the enumeration.
        for( n = 0; n <= bispecID.exp_order(); ++n )
        {
            Real sum = 0;
            for( size_t i = 0; i <= bispecID.exp_order(); ++i )
            {
                /// The idea is to put gA from the enum fields into specC as the index referring to the field, similar to GF in bim-solver.
                sum += chebyExp.specC( fields::gA, m, i ) * chebyC( 0, i, n );
            }
            //gA( m, n ) = sum;
        }
    }

    /// Running over the enumeration (i)
    /*for ( const auto e : fields::flds )
    {
        cout << e << endl;
    }*/
};

/** Implementation of the constructor
     */
/*bispecEvolve::bispecEvolve( bispecInput& bispecID, ChebyshevCoefficients& chebyC, ChebyshevExpansion& chebyExp )
{

    gA( m, n )
    {
        Real sum = 0;
        for( size_t i = 0; i <= bispecID.exp_order(); ++i )
        {
            /// The idea is to put gA from the enum fields into specC as the index referring to the field, similar to GF in bim-solver.
            sum += chebyExp.specC( fields::gA, m, i ) * chebyC( 0, i, n );

        }
        return sum;
    }
}*/

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{

    /** Import the values of the Chebyshev polynomials at the collocation points.
        The are constant in time, hence we can compute them in Mathematica with high precision and load them once here.
     */

    ChebyshevCoefficients cc( "include/chebyshev-values/testBin.dat" );
    if( ! cc.isOK () ) {
        return -1;
    }

    /*std::cout << cc.orders() << " x " << cc.chebys()
        << " x " << cc.points() << std::endl;

    for( size_t i = 0; i < cc.orders(); ++i )
    {
        for( size_t j = 0; j < cc.chebys(); ++j )
        {
            for( size_t k = 0; k < cc.points(); ++k )
            {
                std::cout << "(" << i << "," << j << "," << k << ") = "
                          << cc(i,j,k) << std::endl;
            }
        }
    }

    cc.exportChebyCoeffs();

    std::cout << "The following is the evolution matrix," << std::endl;

    for( size_t j = 0; j < cc.chebys(); ++j )
    {
        for( size_t k = 0; k < cc.points(); ++k )
        {
            std::cout << "(" << j << "," << k << ") = "
                        << cc.evolutionMatrix(j,k) << std::endl;
        }
    }*/

    /** Import the initial data, i.e., the values of the Chebyshev spectral coefficients in the Chebyshev series of the fields on the initial hypersurface.
     */

    bispecInput ID("C:/Users/Francesco/Dropbox/Dottorato/Research/3+1_Numerical_bimetric_relativity/BSSN_formalism/C++/bimetric-ss-20181026/run/specInput.dat");
    if( ! ID.isOK () ) {
        return -1;
    }

    /*std::cout << "The following is the spectral initial data," << std::endl;


    std::cout << ID.n_fields() << " x " << ID.exp_order() + 1 << std::endl;

    for( size_t field = 0; field < ID.n_fields(); ++field )
    {
        std::cout << "This is the field " << field << std::endl;
        for( size_t n = 0; n < ID.exp_order() + 1; ++n )
        {
            std::cout << "(" << field << "," << n << ") = "
                      << ID(field,n) << std::endl;
        }

    }

    std::cout << std::endl;*/

    /** Define the array containing the spectral coefficients as a function of the fields defined in the namespace fields, the time step m and the collocation point index n.
     */

    ChebyshevExpansion chebySeries( ID );

    //chebySeries.instantiateSpecID();

    /*std::cout << "The following is the spectral initial data assigned to the coefficients," << std::endl;
    for( size_t field = 0; field < ID.n_fields(); ++field )
    {
        std::cout << "This is the field " << field << std::endl;
        for( size_t n = 0; n < ID.exp_order() + 1; ++n )
        {
            std::cout << "(" << field << "," << n << ") = "
                      << chebySeries.specC( 0, field, n ) << std::endl;
        }
    }*/

    /// Running over the enumeration (i)
    for ( const auto e : fields::flds )
    {
        cout << e << endl;
    }

    return 0;
}
