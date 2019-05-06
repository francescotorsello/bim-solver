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

#define FIELDS gconf, fconf, gtrK, ftrK, gA, fA, gB, fB, gA1, fA1, gL, fL

#define DERS   gconf_r, fconf_r, gtrK_r, ftrK_r, gA_r, fA_r, gB_r, fB_r, gA1_r, fA1_r, gL_r, fL_r

#define DERRS  gconf_rr, fconf_rr, gtrK_rr, ftrK_rr, gA_rr, fA_rr, gB_rr, fB_rr, gA1_rr, fA1_rr, gL_rr, fL_rr


/** Note that the macro 'setfield' below, at present, only works inside bispecEvolve.
    It defines the functions that are included in the evolution equations exported by Mathematica
  */
#define setfield( field ) \
        inline Real field( int m, int n ) \
        { \
            return values_fields[ ( fields::field * ( exp_ord + 1 ) ) * m + ( exp_ord + 1 ) * fields::field + n ]; \
        }

#define setder( field ) \
        inline Real field( int m, int n ) \
        { \
            return values_ders[ ( fields::field * ( exp_ord + 1 ) ) * m + ( exp_ord + 1 ) * fields::field + n ]; \
        }

#define setderr( field ) \
        inline Real field( int m, int n ) \
        { \
            return values_derrs[ ( fields::field * ( exp_ord + 1 ) ) * m + ( exp_ord + 1 ) * fields::field + n ]; \
        }

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

/** 'ChebyshevCoefficients' reads the values of the Chebyshev polynomials at the collocation points and the elements of the evolution matrix needed to evolve the spectral coefficients.
    TODO: merge with bispecInput.
  */
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

/** 'bispecInput' reads the initial data, i.e., the values of the spectral coefficients on the initial hypersurface.
    TODO: merge with ChebyshevCoefficients.
  */
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

protected:

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
            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; ++cheby_index )
            {
                spcoeffs[ ( exp_ord + 1 ) * field + cheby_index ] = bispecID( field, cheby_index );
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

    /// The size of the cached time steps
    size_t mDim;
    size_t mLen;
    size_t mExtra;
    size_t exp_ord;
    size_t n_flds;

    Real *values_fields;
    Real *values_ders;
    Real *values_derrs;

    Real *fields_t;

public:

    /// The loop above would define automatically all the fields in terms of the array initialized by the constructor. However, apparently a loop cannot be 'alone' inside the body of a class. In addition, this is a loop of declarations, and declarations cannot be inside the body of a function. Hnece, encapsulating the loop in a void method in the class does not work. For now, we just type all the declarations by hand, as we do in bim-solver with emitfield().
    /*for ( const auto fld : fields::flds )
    {
        setfield( fld )
    }*/

    setfield( gconf )   setfield( gtrK )    setfield( gA )
    setfield( gB )      setfield( gA1 )     setfield( gL )
    setfield( fconf )   setfield( ftrK )    setfield( fA )
    setfield( fB )      setfield( fA1 )     setfield( fL )

    setder( gconf_r )   setder( gtrK_r )    setder( gA_r )
    setder( gB_r )      setder( gA1_r )     setder( gL_r )
    setder( fconf_r )   setder( ftrK_r )    setder( fA_r )
    setder( fB_r )      setder( fA1_r )     setder( fL_r )

    setderr( gconf_rr ) setderr( gtrK_rr )  setderr( gA_rr )
    setderr( gB_rr )    setderr( gA1_rr )   setderr( gL_rr )
    setderr( fconf_rr ) setderr( ftrK_rr )  setderr( fA_rr )
    setderr( fB_rr )    setderr( fA1_rr )   setderr( fL_rr )

    /** The constructor computes the values of the fields, the derivatives and the evolution equations on the initial hypersurface
     */
    bispecEvolve( bispecInput& bispecID, ChebyshevCoefficients& chebyC, ChebyshevExpansion& chebyExp )
    {

        /// These definitions below are temporary. This class should inherit from ChebyshevExpansion, but I get errors from the inheritance.
        mLen    = 5;
        mExtra  = 9;
        mDim    = mLen + mExtra;
        exp_ord = bispecID.exp_order();
        n_flds  = bispecID.n_fields();

        /// These arrays contain the values of the fields and their spatial first and second derivatives at the collocation points at each time step
        values_fields   = new Real[ mDim * n_flds * ( exp_ord + 1 ) ];
        values_ders     = new Real[ mDim * n_flds * ( exp_ord + 1 ) ];
        values_derrs    = new Real[ mDim * n_flds * ( exp_ord + 1 ) ];

        /// Computation of the values of the fields at the collocation points (CPs) on the initial hypersurface
        for( const auto field : fields::flds )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; ++cheby_index )
                {
                    sum += chebyExp.specC( 0, field, cheby_index ) * chebyC( 0, cheby_index, n );
                }
                values_fields[ ( exp_ord + 1 ) * field + n ] = sum;
            }
        }

        /// Computation of the values of the first radial derivatives at the collocation points (CPs) on the initial hypersurface
        for( const auto der : fields::ders )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; ++cheby_index )
                {
                    sum += chebyExp.specC( 0, der, cheby_index ) * chebyC( 1, cheby_index, n );
                }
                values_ders[ ( exp_ord + 1 ) * der + n ] = sum;
            }
        }

        /// Computation of the values of the second radial derivatives at the collocation points (CPs) on the initial hypersurface
        for( const auto derr : fields::derrs )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; ++cheby_index )
                {
                    sum += chebyExp.specC( 0, derr, cheby_index ) * chebyC( 2, cheby_index, n );
                }
                values_derrs[ ( exp_ord + 1 ) * derr + n ] = sum;
            }
        }

        /// The printouts below print the values of the fields at the collocation points on the initial hypersurface. They are compared against the values in Mathematica and they coincide (up to MachinePrecision).
        /*std::cout << "The following are the values of the fields on the collocation points on the initial hypersurface (g-sector)," << std::endl;

        std::cout << std::endl;

        std::cout << "Conformal factor," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gconf( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gtrK," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gtrK( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gA( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gB," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gB( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA1," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gA1( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gL," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gL( 0, n ) << std::endl;
        }*/

        //#include "eom-BSSN/eomBSSNEvolutionCompEul.h"

        /// This array contains the right-hand sides of the evolution equations for the fields.
        fields_t   = new Real[ n_flds ];

        for( size_t field = 0; field < n_flds; ++field )
        {
            fields_t[ field ] = 0;
        }

    }
};



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

    /** Evolve the coefficients
     */

    bispecEvolve bispecEvolution( ID, cc, chebySeries );

    return 0;
}
