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

    return 0;
}
