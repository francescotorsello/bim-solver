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

using namespace std;

typedef double Real;

/// The following namespace contains the indexing of the fields, to be used in the other classes
namespace fields
{
    /// This enumeration contains all the fields that need to be expanded in Chebyshev series
    enum fldCheby
    {
        gconf, fconf,
        trgK, trfK,
        gA, fA,
        gB, fB,
        gA1, fA1,
        gL, fL,

        gconf_r, fconf_r,
        trgK_r, trfK_r,
        gA_r, fA_r,
        gB_r, fB_r,
        gA1_r, fA1_r,
        gL_r, fL_r,

        gconf_rr, fconf_rr,
        trgK_rr, trfK_rr,
        gA_rr, fA_rr,
        gB_rr, fB_rr,
        gA1_rr, fA1_rr,
        gL_rr, fL_rr,

    };

    /// The following vector contains the ORDERED fields whose spectral coefficients
    /// have to be set equal to the values read by the class bispecInput.
    /// This order MUST match the order of exportation in the Mathematica file.

    static const std::vector<int> bispecInput_fields =
    {
        gconf,
        fconf,
        trgK,
        trfK,
        gA,
        gB,
        fA,
        fB,
        gA1,
        fA1,
        gL,
        fL

    };

}

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

    // Method to access the initial data
    //
    Real operator() ( size_t i, size_t j )
    {
        return spec_coeff[ i * ( expansion_order + 1 ) + j ];
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

class bispecEvolve
{
    int m; /// The time step
    int n; /// The collocation point index

    Real gconf( int m, int n );
    Real gtrK( int m, int n );
    Real gA( int m, int n );
    Real gB( int m, int n );
    Real gA1( int m, int n );
    Real gL( int m, int n );

    Real fconf( int m, int n );
    Real ftrK( int m, int n );
    Real fA( int m, int n );
    Real fB( int m, int n );
    Real fA1( int m, int n );
    Real fL( int m, int n );

    /** Declaration of the constructor
     */
    bispecEvolve( bispecInput& bispecID, ChebyshevCoefficients& chebyC );

};

/** Implementation of the constructor
     */
/*bispecEvolve::bispecEvolve( bispecInput& bispecID, ChebyshevCoefficients& chebyC )
{
    gA ( m , n ) = 0;
    for( i = 0; i <= ID.exp_order(), ++i )
    {
        /// The idea is to put gA from the enum fields into specC as the index referring to the field, similar to GF in bim-solver.
       gA ( m, n ) += specC( gA, m, i ) * chebyC( 0, i, n );
    }
}*/

int main()
{
    ChebyshevCoefficients cc( "include/chebyshev-values/testBin.dat" );
    if( ! cc.isOK () ) {
        return -1;
    }

    std::cout << cc.orders() << " x " << cc.chebys()
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
    }

    bispecInput ID("C:/Users/Francesco/Dropbox/Dottorato/Research/3+1_Numerical_bimetric_relativity/BSSN_formalism/C++/bimetric-ss-20181026/run/specInput.dat");
    if( ! ID.isOK () ) {
        return -1;
    }

    std::cout << "The following is the spectral initial data," << std::endl;


    std::cout << ID.n_fields() << " x " << ID.exp_order() << std::endl;

    for( size_t i = 0; i < ID.n_fields(); ++i )
    {
        for( size_t j = 0; j < ID.exp_order() + 1; ++j )
        {
            std::cout << "(" << i << "," << j << ") = "
                        << ID(i,j) << std::endl;
        }
    }

    return 0;
}
