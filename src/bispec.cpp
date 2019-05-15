/**
 *  @file      bispec.cpp
 *  @brief     The covariant BSSN evolution for spherically symmetric bimetric spacetimes,
               with the pseudospectral method.
 *  @authors   Francesco Torsello, Mikica Kocic
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
#include<functional>

#include "../include/numMethods/dataTypes.h"

#include "../include/sys/slog.h"
#include "../include/sys/paramsHolder.h"
#include "../include/bimetricModel.h"

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

#define EVOLVED_FIELDS  gconf, fconf, gtrK, ftrK, gA, fA, gB, fB, gA1, fA1, gL, fL, \
                        gsig, fsig, gAsig, fAsig, pfD, pfS, pftau

#define EVEN_FIELDS  gconf, fconf, gtrK, ftrK, gA, fA, gB, fB, gA1, fA1, gL, fL, \
                     gsig, fsig, gAsig, fAsig, pfD, pfS, pftau

#define ODD_FIELDS  gL, fL, pfS

#define FIELDS_T    gconf_t, fconf_t, gtrK_t, ftrK_t, gA_t, fA_t, gB_t, fB_t, gA1_t,  \
                    fA1_t, gsig_t, fsig_t, gAsig_t, fAsig_t,  gD_t, gS_t, gtau_t

#define DERS    gconf_r, fconf_r, gtrK_r, ftrK_r, gA_r, fA_r, gB_r, fB_r, gA1_r, fA1_r, \
                gL_r, fL_r, gsig_r, fsig_r, gAsig_r, fAsig_r, pfD_r, pfS_r, pftau_r

#define REG_DERS    gBetr_r, fBetr_r, gLr_r, fLr_r, gderAlpr_r, fderAlpr_r, \
                    gderconfr_r, fderconfr_r

#define EVEN_DERS   gconf_r, fconf_r, gtrK_r, ftrK_r, gA_r, fA_r, gB_r, fB_r, gA1_r, \
                    fA1_r, gL_r, fL_r, gsig_r, fsig_r, gAsig_r, fAsig_r, pfD_r, pfS_r, \
                    pftau_r

#define ODD_DERS  gL_r, fL_r, pfS_r

#define DERRS   gconf_rr, fconf_rr, gtrK_rr, ftrK_rr, gA_rr, fA_rr, gB_rr, fB_rr, \
                gA1_rr, fA1_rr, gL_rr, fL_rr, gsig_rr, fsig_rr, gAsig_rr, fAsig_rr, \
                pfD_rr, pfS_rr, pftau_rr

#define EVEN_DERRS  gconf_rr, fconf_rr, gtrK_rr, ftrK_rr, gA_rr, fA_rr, gB_rr, fB_rr, \
                    gA1_rr, fA1_rr, gL_rr, fL_rr, gsig_rr, fsig_rr, gAsig_rr, fAsig_rr, \
                    pfD_rr, pfS_rr, pftau_rr

#define ODD_DERRS   gL_rr, fL_rr, pfS_rr

#define GAUGE       gAlp, fAlp, hAlp, gBet, fBet, hBet, p

#define GAUGE_R     gAlp_r, fAlp_r, hAlp_r, gBet_r, fBet_r, hBet_r, p_r

#define GAUGE_RR    gAlp_rr, fAlp_rr, hAlp_rr, gBet_rr, fBet_rr, hBet_rr, p_rr

#define VALENCIA    pfD, pfS, pftau

#define VALENCIA_R  pfD_r, pfS_r, pftau_r

#define VALENCIA_RR pfD_rr, pfS_rr, pftau_rr

#define ALL_FIELDS  gconf, fconf, gtrK, ftrK, gA, fA, gB, fB, gA1, fA1, gL, fL, \
                    gsig, fsig, gAsig, fAsig, gAlp, fAlp, hAlp, gBet, fBet, hBet, p, \
                    pfD, pfS, pftau

#define ALL_DERS    gconf_r, fconf_r, gtrK_r, ftrK_r, gA_r, fA_r, gB_r, fB_r, gA1_r, \
                    fA1_r, gL_r, fL_r, gsig_r, fsig_r, gAsig_r, fAsig_r, gAlp_r, fAlp_r, \
                    hAlp_r, gBet_r, fBet_r, hBet_r, p_r, pfD_r, pfS_r, pftau_r

#define ALL_DERRS   gconf_rr, fconf_rr, gtrK_rr, ftrK_rr, gA_rr, fA_rr, gB_rr, fB_rr, \
                    gA1_rr, fA1_rr, gL_rr, fL_rr, gsig_rr, fsig_rr, gAsig_rr, fAsig_rr, \
                    gAlp_rr, fAlp_rr, hAlp_rr, gBet_rr, fBet_rr, hBet_rr, p_rr, \
                    pfD_rr, pfS_rr, pftau_rr


/** Note that the macro 'setfield' below, at present, only works inside bispecEvolve.
    It defines the functions that are included in the evolution equations exported by
    Mathematica
  */
#define setfield( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_fields[ ( n_all_flds * ( exp_ord + 1 ) ) * m \
            + ( exp_ord + 1 ) * fields::field + n ]; \
        }

#define setder( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_ders[ ( n_all_flds * ( exp_ord + 1 ) ) * m \
            + ( exp_ord + 1 ) * fields::field + n ]; \
        }

#define setderr( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_derrs[ ( n_all_flds * ( exp_ord + 1 ) ) * m \
            + ( exp_ord + 1 ) * fields::field + n ]; \
        }

#define setregder( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_reg_ders[ ( 8 * ( exp_ord + 1 ) ) * m \
            + ( exp_ord + 1 ) * fields::field + n ]; \
        }

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/** Namespace 'fields' contains the indexing of the fields and their derivatives,
    to be used in the other classes
  */
namespace fields
{
    /*/// 'fldCheby' contains all the fields that need to be expanded in Chebyshev series
    enum fldCheby { FIELDS };
    static const fldCheby flds[] = { FIELDS }; // this is needed to be able to iterate over the enumeration ( taken from https://stackoverflow.com/questions/261963/how-can-i-iterate-over-an-enum )

    /// 'derCheby' contains all the first radial derivatives that need to be expanded in Chebyshev series
    enum derCheby { DERS };
    static const derCheby ders[] = { DERS };

    /// 'derrCheby' contains all the second radial derivatives that need to be expanded in Chebyshev series
    enum derrCheby { DERRS };
    static const derrCheby derrs[] = { DERRS };

    /// 'gaugeVars' contains all the gauge variables
    enum gaugeVars { GAUGE };
    static const gaugeVars gauges[] = { GAUGE };

    /// 'gaugeDers' contains the first radial derivatives of the gauge variables
    enum gaugeDers { GAUGE_R };
    static const gaugeDers gauges_r[] = { GAUGE_R };

    /// 'gaugeDerrs' contains the second radial derivatives of the gauge variables
    enum gaugeDerrs { GAUGE_RR };
    static const gaugeDerrs gauges_rr[] = { GAUGE_RR };

    /// 'valenciaVars' contains all the gauge variables
    enum valenciaVars { VALENCIA };
    static const valenciaVars valencia[] = { VALENCIA };

    /// 'valenciaDers' contains the first radial derivatives of the gauge variables
    enum valenciaDers { VALENCIA_R };
    static const valenciaDers valencia_r[] = { VALENCIA_R };

    /// 'valenciaDerrs' contains the second radial derivatives of the gauge variables
    enum valenciaDerrs { VALENCIA_RR };
    static const valenciaDerrs valencia_rr[] = { VALENCIA_RR };*/

    /// TODO: instead of defining different enumerations, define different vectors over which you can run over.

    /// 'allFields' contains all the gauge variables
    enum allFields { ALL_FIELDS };
    static const allFields all_flds[]       = { ALL_FIELDS };
    static const allFields evolved_flds[]   = { EVOLVED_FIELDS };
    static const allFields even_flds[]      = { EVEN_FIELDS };
    static const allFields odd_flds[]       = { ODD_FIELDS };

    /// 'allDers' contains the first radial derivatives of the gauge variables
    enum allDers { ALL_DERS };
    static const allDers all_ders[]  = { ALL_DERS };
    static const allDers even_ders[] = { EVEN_DERS };
    static const allDers odd_ders[]  = { ODD_DERS };

    /// 'regDers' contains the regularizing radial derivatives
    enum regDers { REG_DERS };
    static const regDers reg_ders[] = { REG_DERS };

    /// 'allDerrs' contains the second radial derivatives of the gauge variables
    enum allDerrs { ALL_DERRS };
    static const allDerrs all_derrs[]   = { ALL_DERRS };
    static const allDerrs even_derrs[]  = { EVEN_DERRS };
    static const allDerrs odd_derrs[]   = { ODD_DERRS };

    /// perhaps the vector below is useless
    //static const std::vector<Int> bispecInput_fields = { FIELDS };

}



//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/** 'ChebyshevCoefficients' reads the values of the Chebyshev polynomials at the
    collocation points and the elements of the evolution matrix needed to evolve the
    spectral coefficients.
    TODO: merge with bispecInput.
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
    size_t orders () const{ return derivative_order; }
    size_t chebys () const{ return expansion_order; }
    size_t points () const{ return number_collocations; }

    // Method to access the coefficients
    //
    Real operator()( size_t der_order, size_t cheby_index, size_t n )
    {
        return coeff[ der_order * ( expansion_order + 1 ) * number_collocations
                    + cheby_index * number_collocations + n ];
    }

    // Method to access the coefficients of the regularized first derivatives
    //
    Real reg_ders_cheby( size_t cheby_index, size_t n )
    {
        return reg_der_coeff[ ( expansion_order + 1 ) * n + n ];
    }

    // Method to access the coefficients of the regularized second derivatives
    //
    Real reg_derrs_cheby( size_t cheby_index, size_t n )
    {
        return reg_derr_coeff[ ( expansion_order + 1 ) * n + n ];
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

        if ( reg_coeff_size != fread( reg_der_coeff, sizeof(Real), data_size, inf ) )
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

        if ( reg_coeff_size != fread( reg_derr_coeff, sizeof(Real), data_size, inf ) )
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

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/** 'bispecInput' reads the initial data, i.e., the values of the spectral coefficients
    on the initial hypersurface.
    TODO: merge with ChebyshevCoefficients.
  */
class BispecInput
{
    size_t input_fields;
    size_t number_fields;
    size_t number_gauge_vars = 6;
    size_t number_valencia   = 3;
    size_t expansion_order;
    size_t data_size;
    size_t number_collocations;
    size_t number_chebycoeffs;
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
    size_t n_allfields      () const{ return input_fields; }
    size_t n_fields         () const{ return number_fields; }
    size_t n_gauges         () const{ return number_gauge_vars; }
    size_t n_valencia       () const{ return number_valencia; }
    size_t exp_order        () const{ return expansion_order; }
    size_t n_collocations   () const{ return number_collocations; }
    size_t n_chebycoeffs    () const{ return number_chebycoeffs; }
    size_t size_ID          () const{ return data_size; }

    // Method to access the initial data
    //
    Real operator() ( size_t field, size_t n )
    {
        return spec_coeff[ ( expansion_order + 1 ) * field + n ];
    }

    // Constructor: Reads the spectral initial data from a file.
    //
    BispecInput( const std::string fileName )
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
            input_fields    = size_t(x);
            /// The number of evolved fields is the number of input fields, \
                minus the number of gauge variables ( at present, 6 ) minus p minus r
            number_fields   = input_fields - number_gauge_vars - 1 - 1;
        }
        else {
            std::cerr << "err: ID: Cannot read number_fields" << std::endl;
            fclose( specid );
            return;
        }

        // Read the order of Chebyshev expansion
        //
        if( 1 == fread( &x, sizeof(Real), 1, specid ) ) {
            expansion_order     = size_t(x);
            number_collocations = size_t(( x + 1 ) / 2);
            number_chebycoeffs  = number_collocations;
        }
        else {
            std::cerr << "err: ID: Cannot read expansion_order" << std::endl;
            fclose( specid );
            return;
        }

        // Read the spectral initial data
        //
        data_size   = input_fields * number_chebycoeffs;
        spec_coeff  = new Real[ data_size ];

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



//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/** 'ChebyshevExpansion' defines the array containing the spectral coefficients and a
    method to access them. It also stores the initial data into this array.
  */
class ChebyshevExpansion
{

protected:

    /// The size of the cached time steps
    size_t mDim;
    size_t mLen;
    size_t mExtra;
    size_t exp_ord;
    size_t n_collocs;
    size_t n_chebycffs;
    size_t n_allflds;

public:

    /// The array containing the spectral coefficients for all the Chebyshev series \
        of the fields
    Real *spcoeffs;

    /** The constructor defines the array containing the spectral coefficients for all
        the fields and assigns their initial values
      */
    ChebyshevExpansion( BispecInput& bispecID )
    {

        /// The parameter mLen should be given by the user in config.ini
        mLen        = 5;
        mExtra      = 9;
        mDim        = mLen + mExtra;

        exp_ord     = bispecID.exp_order();
        n_collocs   = bispecID.n_collocations();
        n_chebycffs = bispecID.n_chebycoeffs();
        n_allflds   = bispecID.n_allfields();

        spcoeffs    = new Real[ mDim * n_allflds * n_chebycffs ];

        //std::cout << std::endl;
        for( size_t field = 0; field < n_allflds; ++field )
        {
            //std::cout << "This is the field " << field << std::endl;
            for( size_t cheby_index = 0; cheby_index < n_chebycffs; ++cheby_index )
            {
                spcoeffs[ n_chebycffs * field + cheby_index ] =
                    bispecID( field, cheby_index );
                 //std::cout << "(" << field << "," << n << ") : "
                 //     << spcoeffs[ ( exp_ord + 1 ) * field + n ] << " = "
                 //     << bispecID( field, n ) << std::endl;
            }
        }
        //std::cout << std::endl;

    }
};



//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/** 'primaryFields' defines the methods to access the values of the fields, (the gauges
    variables), the Valencia variables and their first and second radial derivatives. A
    'primary field' is a field which needs to access the values of the Chebyshev
    coefficients to be computed.

    The format is: field( m, n ), with m being the time step and n being the index of the
    collocation point. This is the same format exported by Mathematica.

    The constructor assigns the initial data to the fields, i.e., it assignes values to
    field( 0, n ) for all n.

    'primaryFields' was part of 'bispecEvolve', but it is now a separate class, since the
    classes gaugeVariables and and dependentFields need the primary fields to be defined.
    In turn, bispecEvolve needs gaugeVariables and dependentFields in order to use the
    evolution equations. Hence, bispecEvolve is now relegated strictly to the computation
    of the time derivative of the spectral coefficients and their evolution through
    'integrator.h'.

    Note that bispecEvolve is also responsible to assign values to the arrays
    values_fields, values_fields, values_fields at the next time step. The primary
    fields are called through methods accessing these arrays. Hence, if these arrays
    have a meaningful value, the primary fields will access them. This is necessary for
    the dependent fields in the derived class dependentFields to have a meaningful value.
    Hence, in the end what we need to do in bispecEvolve is to evolve the coefficients
    and to assign values to values_fields, values_fields, values_fields. The rest will
    follow.
  */
class PrimaryFields
    : //ChebyshevCoefficients,
      ChebyshevExpansion
{

protected:

    size_t m; /// The time step
    size_t n; /// The collocation point index

    /// The size of the cached time steps
    size_t mDim;
    size_t mLen;
    size_t mExtra;
    size_t exp_ord;
    size_t n_flds;
    size_t n_gauge;
    size_t n_vlc;
    size_t n_all_flds;
    size_t n_collocs;
    size_t n_chebycffs;

    Real *values_fields;
    Real *values_ders;
    Real *values_reg_ders;
    Real *values_derrs;
    Real *values_colpoints;

    Real *fields_t;
    Real *chebyCoeff;
    Real *regDerCheby;
    Real *regDerrCheby;

public:

    /// Method to access the spectral coefficients. Since we will make use of the \
        integrator from bim-solver, we copy the structure of GF in gridDriver.
    Real specC( size_t m, size_t field, size_t n )
    {
        return spcoeffs[ ( n_flds * n_chebycffs ) * m + n_chebycffs * field + n ];
    }

    /** Method to access the values of the collocation points.
      */

    inline Real r( Int m, Int n )
    {
        return values_colpoints[ n ];
    }

    inline Real r_minus( Int m, Int n )
    {
        return pow2( values_colpoints[ n ] ) - 1;
    }

    inline Real r_plus( Int m, Int n )
    {
        return pow2( values_colpoints[ n ] ) + 1;
    }

    /** Methods to access the values of the fields at the collocation points. Everyone of
        these functions depend on the time step m and the collocation point index n.
      */

    setfield( gconf )   setfield( gtrK )    setfield( gA )
    setfield( gB )      setfield( gA1 )     setfield( gL )
    setfield( fconf )   setfield( ftrK )    setfield( fA )
    setfield( fB )      setfield( fA1 )     setfield( fL )

    setfield( gsig )    setfield( fsig )
    setfield( gAsig )   setfield( fAsig )

    setfield( gAlp )    setfield( fAlp )    setfield( hAlp )
    setfield( gBet )    setfield( fBet )    setfield( hBet )
    setfield( p )

    setfield( pfD )     setfield( pfS )     setfield( pftau )

    /** Define methods to access the values of the first radial derivatives of the fields
        at the collocation points. Everyone of these functions depend on the time step m
        and the collocation point index n.
      */

    setder( gconf_r )   setder( gtrK_r )    setder( gA_r )
    setder( gB_r )      setder( gA1_r )     setder( gL_r )
    setder( fconf_r )   setder( ftrK_r )    setder( fA_r )
    setder( fB_r )      setder( fA1_r )     setder( fL_r )

    setder( gsig_r )    setder( fsig_r )
    setder( gAsig_r )   setder( fAsig_r )

    setder( gAlp_r )    setder( fAlp_r )    setder( hAlp_r )
    setder( gBet_r )    setder( fBet_r )    setder( hBet_r )
    setder( p_r )

    setder( pfD_r )     setder( pfS_r )     setder( pftau_r )

    /** Define methods to access the values of the second radial derivatives of the
        fields at the collocation points. Everyone of these functions depend on the
        time step m and the collocation point index n.
      */

    setderr( gconf_rr ) setderr( gtrK_rr )  setderr( gA_rr )
    setderr( gB_rr )    setderr( gA1_rr )   setderr( gL_rr )
    setderr( fconf_rr ) setderr( ftrK_rr )  setderr( fA_rr )
    setderr( fB_rr )    setderr( fA1_rr )   setderr( fL_rr )

    setderr( gsig_rr )  setderr( fsig_rr )
    setderr( gAsig_rr ) setderr( fAsig_rr )

    setderr( gAlp_rr )  setderr( fAlp_rr )  setderr( hAlp_rr )
    setderr( gBet_rr )  setderr( fBet_rr )  setderr( hBet_rr )
    setderr( p_rr )

    setderr( pfD_rr )   setderr( pfS_rr )   setderr( pftau_rr )

    /** Define methods to access the values of the first and second regularized
        radial derivatives of some of the fields at the collocation points.
        Everyone of these functions depend on the time step m and the collocation
        point index n.
      */

    setregder( gBetr_r )        setregder( fBetr_r )
    setregder( gLr_r )          setregder( fLr_r )
    setregder( gderAlpr_r )     setregder( fderAlpr_r )
    setregder( gderconfr_r )    setregder( fderconfr_r )

    inline void computeFields( size_t m )
    {
        /** TODO: In the expansion below, one should include only the appropriate parity
            in the Chebyshev series of the fields
          */

        /// Computation of the values of the fields at the collocation points (CPs) on \
            the initial hypersurface
        for( const auto field : fields::even_flds )
        {
            for( size_t n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += specC( m, field, cheby_index )
                            * chebyCoeff[ cheby_index * n_collocs + n ];
                }
                values_fields[ m * n_collocs * n_flds + n_collocs * field
                    + n ] = sum;
            }
        }
        for( const auto field : fields::odd_flds )
        {
            for( size_t n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += specC( m, field, cheby_index )
                            * chebyCoeff[ cheby_index * n_collocs + n ];
                }
                values_fields[ m * n_collocs * n_flds + n_collocs * field
                    + n ] = sum;
            }
        }

        /// Computation of the values of the first radial derivatives at the \
            collocation points (CPs) on the initial hypersurface
        for( const auto der : fields::even_ders )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += specC( m, der, cheby_index )
                            * chebyCoeff[ ( exp_ord + 1 ) * ( exp_ord + 1 )
                                + cheby_index * ( exp_ord + 1 ) + n ];
                }
                values_ders[ m * n_collocs * n_flds + n_collocs * der
                    + n ] = sum;
            }
        }
        for( const auto der : fields::odd_ders )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += specC( m, der, cheby_index )
                            * chebyCoeff[ ( exp_ord + 1 ) * ( exp_ord + 1 )
                                + cheby_index * ( exp_ord + 1 ) + n ];
                }
                values_ders[ m * n_collocs * n_flds + n_collocs * der
                    + n ] = sum;
            }
        }


        /// Computation of the values of the second radial derivatives at the \
            collocation points (CPs) on the initial hypersurface
        for( const auto derr : fields::even_derrs )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += specC( m, derr, cheby_index )
                            * chebyCoeff[ 2*( exp_ord + 1 ) * ( exp_ord + 1 )
                                + cheby_index * ( exp_ord + 1 ) + n ];
                }
                values_derrs[ m * n_collocs * n_flds + n_collocs * derr
                    + n ] = sum;
            }
        }
        for( const auto derr : fields::odd_derrs )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += specC( m, derr, cheby_index )
                            * chebyCoeff[ 2*( exp_ord + 1 ) * ( exp_ord + 1 )
                                + cheby_index * ( exp_ord + 1 ) + n ];
                }
                values_derrs[ m * n_collocs * n_flds + n_collocs * derr
                    + n ] = sum;
            }
        }
    }

    /** The following fields needs to access the spectral coefficients even if they
        are not primary. Hence, they need to be computed in here. TODO: check the parity conditions for these fields.
      */

    inline void computeRegDers( size_t m )
    {
        /*for( const auto regder : fields::reg_ders )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; ++cheby_index )
                {
                    sum += specC( m, derr, cheby_index )
                            * chebyCoeff[ 2*( exp_ord + 1 ) * ( exp_ord + 1 )
                                + cheby_index * ( exp_ord + 1 ) + n ];
                }
                values_derrs[ ( exp_ord + 1 ) * derr + n ] = sum;
            }
        }*/

        for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += specC( m, fields::gBet, cheby_index )
                        * regDerCheby[ cheby_index * ( exp_ord + 1 ) + n ];
            }
            values_reg_ders[ m * ( exp_ord + 1 ) * n_flds +
                ( exp_ord + 1 ) * 0 + n ] = sum;
        }

        for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += specC( m, fields::fBet, cheby_index )
                        * regDerCheby[ cheby_index * ( exp_ord + 1 ) + n ];
            }
            values_reg_ders[ m * ( exp_ord + 1 ) * n_flds +
                ( exp_ord + 1 ) * 1 + n ] = sum;
        }

        for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += specC( m, fields::gL, cheby_index )
                        * regDerCheby[ cheby_index * ( exp_ord + 1 ) + n ];
            }
            values_reg_ders[ m * ( exp_ord + 1 ) * n_flds +
                ( exp_ord + 1 ) * 2 + n ] = sum;
        }

        for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += specC( m, fields::fL, cheby_index )
                        * regDerCheby[ cheby_index * ( exp_ord + 1 ) + n ];
            }
            values_reg_ders[ m * ( exp_ord + 1 ) * n_flds +
                ( exp_ord + 1 ) * 3 + n ] = sum;
        }

        for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += specC( m, fields::gAlp, cheby_index )
                        * regDerrCheby[ cheby_index * ( exp_ord + 1 ) + n ];
            }
            values_reg_ders[ m * ( exp_ord + 1 ) * n_flds +
                ( exp_ord + 1 ) * 4 + n ] = sum;
        }

        for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += specC( m, fields::fAlp, cheby_index )
                        * regDerrCheby[ cheby_index * ( exp_ord + 1 ) + n ];
            }
            values_reg_ders[ m * ( exp_ord + 1 ) * n_flds +
                ( exp_ord + 1 ) * 5 + n ] = sum;
        }

        for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += specC( m, fields::gconf, cheby_index )
                        * regDerrCheby[ cheby_index * ( exp_ord + 1 ) + n ];
            }
            values_reg_ders[ m * ( exp_ord + 1 ) * n_flds +
                ( exp_ord + 1 ) * 6 + n ] = sum;
        }

        for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += specC( m, fields::fconf, cheby_index )
                        * regDerrCheby[ cheby_index * ( exp_ord + 1 ) + n ];
            }
            values_reg_ders[ m * ( exp_ord + 1 ) * n_flds +
                ( exp_ord + 1 ) * 7 + n ] = sum;
        }

    }


    /** The constructor computes the values of the fields, the derivatives and the
        evolution equations on the initial hypersurface
     */
    PrimaryFields(
        BispecInput&            bispecID,
        //const std::string   fileName
        ChebyshevCoefficients&  chebyC
        //ChebyshevExpansion&     chebyExp
    ) :
        //ChebyshevCoefficients( fileName ),
        ChebyshevExpansion( bispecID )
    {
        mLen    = 5;
        mExtra  = 9;
        mDim    = mLen + mExtra;
        exp_ord = bispecID.exp_order();
        n_collocs   = bispecID.n_collocations();
        n_chebycffs = bispecID.n_chebycoeffs();
        n_all_flds = bispecID.n_allfields();
        n_flds  = bispecID.n_fields();
        n_gauge = bispecID.n_gauges();
        n_vlc   = bispecID.n_valencia();

        /// These arrays contain the values of the fields and their spatial first and  \
        second derivatives at the collocation points at each time step

        values_colpoints = new Real[ n_collocs ];
        values_fields    = new Real[ mDim * n_all_flds * n_collocs ];
        values_ders      = new Real[ mDim * n_all_flds * n_collocs ];
        values_reg_ders  = new Real[ mDim * 8          * n_collocs ];
        values_derrs     = new Real[ mDim * n_all_flds * n_collocs ];

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            values_colpoints[ n ] = bispecID( n_all_flds - 1, n );
        }

        chebyCoeff  = new Real[ 4 * ( exp_ord + 1 ) * n_collocs ];
        regDerCheby = new Real[ ( exp_ord + 1 ) * n_collocs ];

        for( size_t der_order = 0; der_order < 4; ++der_order )
        {
            for( size_t cheby_index = 0; cheby_index < ( exp_ord + 1 ); ++cheby_index )
            {
                for( size_t n = 0; n < n_collocs; ++n )
                {
                    chebyCoeff[ der_order * ( exp_ord + 1 ) * n_collocs
                        + cheby_index * n_collocs + n ] =
                            chebyC( der_order, cheby_index, n );
                }
            }
        }

        for( size_t cheby_index = 0; cheby_index < ( exp_ord + 1 ); ++cheby_index )
        {
            for( size_t n = 0; n < ( exp_ord + 1 ); ++n )
            {
                regDerCheby[ cheby_index * n_collocs + n ] =
                        chebyC.reg_derrs_cheby( cheby_index, n );
            }
        }

        computeFields( 0 );
        computeRegDers( 0 );

        /** TODO: In the expansion below, one should include only the appropriate parity
            in the Chebyshev series of the fields
          */

        /// Computation of the values of the fields at the collocation points (CPs) on \
            the initial hypersurface
        /*for( const auto field : fields::all_flds )
        {
            for( size_t n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; ++cheby_index )
                {
                    sum += chebyExp.specC( 0, field, cheby_index )
                            * chebyC( 0, cheby_index, n );
                }
                values_fields[ ( exp_ord + 1 ) * field + n ] = sum;
            }
        }

        /// Computation of the values of the first radial derivatives at the \
            collocation points (CPs) on the initial hypersurface
        for( const auto der : fields::all_ders )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; ++cheby_index )
                {
                    sum += chebyExp.specC( 0, der, cheby_index )
                            * chebyC( 1, cheby_index, n );
                }
                values_ders[ ( exp_ord + 1 ) * der + n ] = sum;
            }
        }

        /// Computation of the values of the second radial derivatives at the \
            collocation points (CPs) on the initial hypersurface
        for( const auto derr : fields::all_derrs )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; ++cheby_index )
                {
                    sum += chebyExp.specC( 0, derr, cheby_index )
                            * chebyC( 2, cheby_index, n );
                }
                values_derrs[ ( exp_ord + 1 ) * derr + n ] = sum;
            }
        }*/

        /// The printouts below print the values of the fields at the collocation points \
            on the initial hypersurface. They are compared against the values in \
            Mathematica and they coincide (up to MachinePrecision).

        /*std::cout << "The following are the values of the fields on the collocation
            points on the initial hypersurface (g-sector)," << std::endl;

        std::cout << std::endl;

        std::cout << "Collocation points," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << r( 0, n ) << std::endl;
        }
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
        }
        std::cout << std::endl;
        std::cout << "gAlp," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gAlp( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gBet," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gBet( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "p," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << p( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfD," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << pfD( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfS," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << pfS( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pftau," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << pftau( 0, n ) << std::endl;
        }*/

    }
};


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/** 'gaugeVariables' defines and compute the gauge variables (maybe not needed, or maybe
    it can be used for the maximal slicing BVP)
  */
/*class gaugeVariables
    : defineFields
{

protected:

    Real Kdiff = 0.01;
    Real Kelas = 0.03;

    //#include "../include/eom-BSSN/eomBSSNKDGaugeComp.h"

public:

    gaugeVariables(
        bispecInput& bispecID,
        ChebyshevCoefficients& chebyC,
        ChebyshevExpansion& chebyExp
    ) :
        defineFields( bispecID, chebyC, chebyExp )
    {

    }

};*/


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/** 'dependentFields' defines the methods to access the dependent (nonprimary) fields
    appearing in the evolution equations. If the primary fields have an assigned value
    at the time step m, then the dependent fields will also have a defined value.
  */
class DependentFields
    : virtual PrimaryFields,
      BimetricModel
{

protected:

    /** Method to access the time step dimension
      */

    inline get_mDim()
    {
        return mDim;
    }

    inline get_m()
    {
        return m;
    }

    inline get_n()
    {
        return n;
    }

    /// Dependent fields depending only on the primary fields

    #include "../include/eom-BSSN/eomBSSNRicciComp.h"

    /// Dependent fields which other dependent fields depend on

    inline Real R  ( Int m, Int n )
    {
        Real x;
        isGR() ? x = 1 : x = ( fB(m,n)*exp(2*fconf(m,n))) / (gB(m,n)*exp(2*gconf(m,n)));
        return x;
    }

    inline Real Lt ( Int m, Int n )
    {
        Real x;
        isGR() ? x = 1 : x = sqrt( 1 + p(m,n) * p(m,n) );
        return x;
    }

    inline Real Lt2 ( Int m, Int n )
    {
        Real x;
        isGR() ? x = 1 : x = 1 + p(m,n) * p(m,n);
        return x;
    }

    inline Real gA2 ( Int m, Int n )    { return - gA1(m,n) / 2; }

    inline Real fA2 ( Int m, Int n )    { return - fA1(m,n) / 2; }

    inline Real gK1 ( Int m, Int n )    { return gA1(m,n) + 1/3 *gtrK(m,n); }

    inline Real gK2 ( Int m, Int n )    { return gA2(m,n) + 1/3 *gtrK(m,n); }

    inline Real fK1 ( Int m, Int n )    { return fA1(m,n) + 1/3 *ftrK(m,n); }

    inline Real fK2 ( Int m, Int n )    { return fA2(m,n) + 1/3 *ftrK(m,n); }

    /// Regularizing fields

    inline Real gBetr( Int m, Int n )   { return -2 * r_minus(m,n) * gBet(m,n)
        / ( r(m,n) + TINY_Real ); }

    inline Real fBetr( Int m, Int n )   { return -2 * r_minus(m,n) * fBet(m,n)
        / ( r(m,n) + TINY_Real ); }

    inline Real gLr( Int m, Int n )     { return -2 * r_minus(m,n) * gL(m,n)
        / ( r(m,n) + TINY_Real ); }

    inline Real fLr( Int m, Int n )     { return -2 * r_minus(m,n) * fL(m,n)
        / ( r(m,n) + TINY_Real ); }

    /// Last dependent fields

    #include "../include/eom-BSSN/eomBSSNMatterComp.h"
    #include "../include/eom-BSSN/eomBSSNSourcesCompEul.h"
    #include "../include/eom-BSSN/eomBSSNEvolutionCompEul.h"

public:

    DependentFields(
        BispecInput&            bispecID,
        //const std::string       fileName,
        ChebyshevCoefficients&  chebyC,
        //ChebyshevExpansion&     chebyExp,
        Parameters&             params
    ) :
        PrimaryFields( bispecID, chebyC ),
        BimetricModel( params )
    {

        /// The printouts below print the values of the dependent fields at the collocation points on the initial hypersurface.

        /*std::cout << "The following are the values of the dependent fields on the collocation points on the initial hypersurface," << std::endl << std::endl;

        std::cout << "The collocation points, "<< std::endl << std::endl;
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << r( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "The Ricci terms," << std::endl << std::endl;

        std::cout << "gRicci:   ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gRicci( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fRicci:   ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fRicci( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Lorentz factor, extrinsic curvatures and auxiliary field," << std::endl << std::endl;

        std::cout << "R:   ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << R( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "Lt:  ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << Lt( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "Lt2: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << Lt2( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gA2: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gA2( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gK1: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gK1( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gK2: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gK2( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fA2: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fA2( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fK1: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fK1( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fK2: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fK2( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "The matter fields," << std::endl << std::endl;

        std::cout << "pfv:   ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << pfv( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "pfv_r: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << pfv_r( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gD_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gD_t( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gS_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gS_t( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gtau_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gtau_t( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gW: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gW( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "grhobar: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << grhobar( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "grho: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << grho( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gj: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gj( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gJ11: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gJ11( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gJ22: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gJ22( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "frho: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << grho( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fj: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fj( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fJ11: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fJ11( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fJ22: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fJ22( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;;

        std::cout << "The sources in the evolution equations at the collocation points,"
            << std::endl << std::endl;

        std::cout << "gJK: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gJK( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gJA1: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gJA1( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gJL: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gJL( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fJK: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fJK( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fJA1: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fJA1( 0, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fJL: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fJL( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "The evolution equations," << std::endl << std::endl;

        std::cout << "gconf_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gconf_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gtrK_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gtrK_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gA_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gA_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gB_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gB_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gA1_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gA1_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gL_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gL_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fconf_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fconf_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "ftrK_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << ftrK_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fA_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fA_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fB_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fB_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fA1_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fA1_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fL_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << fL_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gD_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gD_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gS_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gS_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gtau_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gtau_t( 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Other fields," << std::endl << std::endl;

        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << (gBet(0,n) * (2 - 2 * r(0,n))) / (1 + TINY_Real
                + r(0,n)) << ", ";
        }
        std::cout << std::endl << std::endl;

        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << gBet(0,n) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << std::endl;*/

    }
};


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/** The evolution of the gauge variables is not included yet in 'bispecEvolve'. It must be handled separately because it depends on the choice of the user.
  */
class BispecEvolve
    : virtual PrimaryFields,
      DependentFields,
      ChebyshevExpansion
{

    Int m;
    Int n;
    Real *spec_t;

    size_t exp_ord;

    typedef Real (DependentFields::*FP)( Int m, Int n );
    std::vector<FP> fields_t = { FIELDS_T };

    Real *b;
    Real *A_even;
    Real *A_odd;
    //ChebyshevCoefficients* cheby_pointer;

    /*Real evolutionMatrix( size_t cheby_index, size_t n )
    {
        return ee_matrix[ cheby_index * collocation_point + n ];
    }*/

public:

    /** Method to access the vector b
      */

    inline Real evolution_eqs( Int m, Int field, Int n )
    {
        return b[ fields_t.size() * ( exp_ord + 1 ) * m
            + ( exp_ord + 1 ) * field + n ];
    }

    /** Method to access the time derivatives of the spectral coefficients
      */

    inline Real get_spec_t( size_t m, size_t field, size_t n )
    {
        return spec_t[ ( exp_ord + 1 * fields_t.size() ) * m
            + ( exp_ord + 1 ) * field + n ];
    }

    /** Method to fill the vector b
      */

    inline void arrange_fields_t( size_t m )
    {
        for( size_t field = 0; field < fields_t.size(); ++field )
        {
            for( size_t n = 0; n < exp_ord + 1; ++n )
            {
                b[ ( ( exp_ord + 1 ) * fields_t.size() ) * m
                    + ( exp_ord + 1 ) * field + n ] =
                        ( this ->* fields_t[ field ] )( m, n );
            }
        }
    }

    /** Method solve for the time derivatives of the spectral coefficients
      */

    inline void solveDerivatives( size_t m )
    {
        /// For each field
        for( size_t field = 0; field < fields_t.size(); ++field )
        { /// Solve the linear algebraic system, i.e.,
            /// for each Chebyshev index (row)
            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; ++cheby_index )
            {
                spec_t[ ( ( exp_ord + 1 ) * fields_t.size() ) * m
                    + ( exp_ord + 1 ) * field + cheby_index ] = 0;
                /// For each collocation point (column)
                for( size_t n = 0; n < exp_ord +1; ++n )
                {
                    spec_t[ ( ( exp_ord + 1 ) * fields_t.size() ) * m
                        + ( exp_ord + 1 ) * field + cheby_index ] +=
                            A_even[ cheby_index * ( exp_ord + 1 ) + n ]
                                * evolution_eqs( m, field, n );
                            /*(cheby_pointer -> evolutionMatrix( cheby_index, n ))
                            * evolution_eqs( m, field, n );*/
                            /*ChebyshevCoefficients::evolutionMatrix( cheby_index, n )
                                * evolution_eqs( m, field, n );*/
                            /*ChebyshevCoefficients::evolutionMatrix( cheby_index, n )
                                * evolution_eqs( m, field, n );*/

                }
            }
        }
    }

    /** The constructor sets up the linear algebraic system on the initial hypersurface
      */

    BispecEvolve(
        BispecInput&            bispecID,
        //const std::string       fileName,
        ChebyshevCoefficients&  chebyC,
        //ChebyshevExpansion&     chebyExp,
        Parameters&             params
    ) :
        PrimaryFields  ( bispecID, chebyC ),
        DependentFields( bispecID, chebyC, params ),//chebyC, chebyExp, params ),
        ChebyshevExpansion( bispecID )
        //cheby_pointer( &chebyC )
    {

        exp_ord  = bispecID.exp_order();

        A_even  = new Real[ ( bispecID.exp_order() + 1 ) / 2
                    * ( bispecID.exp_order() + 1 ) / 2 ];
        A_odd   = new Real[ ( bispecID.exp_order() + 1 ) / 2
                    * ( bispecID.exp_order() + 1 ) / 2 ];
        b       = new Real[ get_mDim() * fields_t.size()
                    * ( bispecID.exp_order() + 1 ) / 2 ];
        spec_t  = new Real[ get_mDim() * fields_t.size()
                    * ( bispecID.exp_order() + 1 ) / 2 ];

        if( bispecID.n_fields() != fields_t.size() )
        {
            std::cerr << "err: bispecEvolve: mismatch in the number of fields"
                << std::endl << std::endl;
            std::cout << bispecID.n_fields() << ", " << fields_t.size() << std::endl
                 << std::endl;

            /// To stop the run from here, an exception must be thrown
            return;
        }

        for( size_t cheby_index = 0; cheby_index < ( bispecID.exp_order() + 1 );
            ++cheby_index )
        {
            for( size_t n = 0; n < ( bispecID.exp_order() + 1 ); ++n )
            {
                A_even[ cheby_index * ( bispecID.exp_order() + 1 ) + n ] =
                    chebyC.ee_matrix_even[ cheby_index * ( bispecID.exp_order() + 1 )
                        + n ];
            }
        }

        for( size_t cheby_index = 0; cheby_index < ( bispecID.exp_order() + 1 );
            ++cheby_index )
        {
            for( size_t n = 0; n < ( bispecID.exp_order() + 1 ); ++n )
            {
                A_odd[ cheby_index * ( bispecID.exp_order() + 1 ) + n ] =
                    chebyC.ee_matrix_odd[ cheby_index * ( bispecID.exp_order() + 1 )
                        + n ];
            }
        }


        arrange_fields_t( 0 );
        solveDerivatives( 0 );
        specC( 0, 1, 3 );
        computeFields( 1 );

        /*std::cout << "The evolution equations stored in the vector," << std::endl
            << std::endl;

        std::cout << "gconf_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << evolution_eqs( 0, 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gtrK_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << evolution_eqs( 0, 2, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gA_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << evolution_eqs( 0, 4, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gB_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << evolution_eqs( 0, 5, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gA1_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << evolution_eqs( 0, 8, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gL_t: ";
        for( size_t n = 0; n < bispecID.exp_order() + 1; ++n )
        {
            std::cout << evolution_eqs( 0, 10, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << std::endl;*/

    }

};



//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int main()
{

    /** Import the values of the Chebyshev polynomials at the collocation points.
        The are constant in time, hence we can compute them in Mathematica with high
        precision and load them once here.
      */

    ChebyshevCoefficients cc( "../bim-solver/include/chebyshev-values/chebyshev-values.dat" );
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

    cc.exportChebyCoeffs();*/

    /*std::cout << "The following is the evolution matrix," << std::endl;

    for( size_t j = 0; j < cc.chebys(); ++j )
    {
        for( size_t k = 0; k < cc.points(); ++k )
        {
            std::cout << "(" << j << "," << k << ") = "
                        << cc.evolutionMatrix(j,k) << std::endl;
        }
    }*/

    /** Import the initial data, i.e., the values of the Chebyshev spectral coefficients
        in the Chebyshev series of the fields on the initial hypersurface.
      */

    BispecInput ID("../run/specInput.dat");
    if( ! ID.isOK () ) {
        return -1;
    }

    /*std::cout << "The following is the spectral initial data," << std::endl;


    std::cout << ID.n_fields() << " x " << ID.exp_order() + 1 << std::endl;

    for( size_t field = 0; field < ID.n_allfields(); ++field )
    {
        std::cout << "This is the field " << field << std::endl;
        for( size_t n = 0; n < ID.exp_order() + 1; ++n )
        {
            std::cout << "(" << field << "," << n << ") = "
                      << ID(field,n) << std::endl;
        }

    }

    std::cout << std::endl;*/

    /// - Read the run-time configuration parameters (taken from bim-solver)
    ///
    //Parameters params( argc >= 2 ? argv[1] : "config.ini" );
    Parameters params( "../run/config.ini" );

    /** Define the array containing the spectral coefficients as a function of the fields
        defined in the namespace fields, the time step m and the
        collocation point index n.
      */

    ChebyshevExpansion chebySeries( ID );

    //chebySeries.instantiateSpecID();

    /*std::cout << "The following is the spectral initial data assigned to the
        coefficients," << std::endl;
    for( size_t field = 0; field < ID.n_fields(); ++field )
    {
        std::cout << "This is the field " << field << std::endl;
        for( size_t n = 0; n < ID.exp_order() + 1; ++n )
        {
            std::cout << "(" << field << "," << n << ") = "
                      << chebySeries.specC( 0, field, n ) << std::endl;
        }
    }*/

    /** Define and compute the values of the fields at the collocation points
      */

     //primaryFields defineFlds( ID, cc );

    /** Define and compute everything on the initial hypersurface. The constructor of
        dependentFields calls primaryFields' one.
      */

     //dependentFields setupIH( ID, cc, params );

    /** Evolve the spectral coefficients
      */

    BispecEvolve evolution( ID, cc, params );

    //evolution.solveDerivatives( 0 );

    /*std::cout <<
        "The time derivatives of the fields arranged in the vector b,"
            << std::endl << std::endl;
    for( size_t field = 0; field < ID.n_fields(); ++field )
    {
        std::cout << "These are the derivatives of the fields "
            << "(field,collocation_point), " << std::endl << std::endl;
        for( size_t n = 0; n < ID.exp_order() + 1; ++n )
        {
            std::cout << "(" << field << "," << n << ") = "
                      << evolution.evolution_eqs( 0, field, n ) << std::endl;
        }
        std::cout << std::endl;
    }*/

    std::cout <<
        "The time derivatives of the spectral coefficients on the initial hypersurface,"
            << std::endl << std::endl;
    for( size_t field = 0; field < ID.n_fields(); ++field )
    {
        std::cout << "These are the derivatives of the spectral coefficients of "
            << "(field,cheby_index), " << std::endl << std::endl;
        for( size_t cheby_index = 0; cheby_index < ID.exp_order() + 1; ++cheby_index )
        {
            std::cout << "(" << field << "," << cheby_index << ") = "
                      << evolution.get_spec_t( 0, field, cheby_index ) << std::endl;
        }
        std::cout << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    /// The following is the last line in bim-solver. It should be the same here.
    /// - Evolve the equations of motion
    ///
    //return integrator.evolveEquations () ? 0 : -1;

    //////////////////////////////////////////////////////////////////////////////////////

    return 0;
}
