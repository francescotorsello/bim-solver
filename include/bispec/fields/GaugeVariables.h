/**
 *  @file      GaugeVariables.h
 *  @brief     Computes the lapses and shifts at each time slice.
 *  @authors   Francesco Torsello, Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

/**
 *  'GaugeVariables' defines and compute the gauge variables (maybe not needed, or maybe
 *  it can be used for the maximal slicing BVP)
 */

class GaugeVariables
{

protected:

    BispecEvolve& be;

public:

    virtual void applyGaugeCondition( Int m ) = 0;

    GaugeVariables(
        BispecEvolve& BE
    ) :
        be( BE )
    {}

};

class GeodesicSlicing
    : public GaugeVariables
{

public:

    inline void applyGaugeCondition( Int m )
    {
        for( Int n = 0; n < be.n_colpoints(); ++n )
        {
            be.values_fields[ ( be.get_n_all_flds() * be.n_colpoints() ) * m
                + be.n_colpoints() * fields::gAlp + n ] = 1;

            be.values_fields[ ( be.get_n_all_flds() * be.n_colpoints() ) * m
                + be.n_colpoints() * fields::gBet + n ] = 0;

            be.values_fields[ ( be.get_n_all_flds() * be.n_colpoints() ) * m
                + be.n_colpoints() * fields::fAlp + n ] = 1;

            be.values_fields[ ( be.get_n_all_flds() * be.n_colpoints() ) * m
                + be.n_colpoints() * fields::fBet + n ] = 0;
        }
    }

    //Real Kdiff = 0.01;
    //Real Kelas = 0.03;

    //#include "../include/eom-BSSN/eomBSSNKDGaugeComp.h"

    GeodesicSlicing(
        BispecEvolve& BE
    ) :
        // initialize the pointer to the object of type BispecEvolve, to the object passed
        // as argument to the constructor.
        GaugeVariables( BE )
    {}

};

void BispecEvolve::develop( GaugeVariables& gaugeSetter )
{
    // Open the file for text output (one text line per one time step)
    FILE* outf = fopen( "output.dat", "w" );

    Real t = 0;
    Real delta_t = 0.5;
    size_t next_m;

    //cout << "mDim = " << mDim << endl;

    for( size_t m = 0; m <= 1 + 0*(mDim - 1) && t < 10; ++m, t += delta_t )
    {
        //cout << endl << "Time step m = " << m << endl << endl;
        next_m = m + 1;

        if( m == mDim - 1 )
        {
            next_m = 0;
        }

        // Perform the integration step

            arrange_fields_t( m );
            solveSpectralDerivatives( m ); // integStep_Prepare ends

            evolveEvenFields( &EulerMethod, m, next_m, delta_t );
            evolveOddFields ( &EulerMethod, m, next_m, delta_t );

            computeFields   ( next_m );
            gaugeSetter.applyGaugeCondition( next_m );
            computeRegDers  ( next_m ); // Up to here, it could be integStep_Finalize

            dumpPrimaryFields( m );
            //dumpDependentFields( next_m );

        // Over spectral coefficients of the even fields
        for( size_t field = 0; field <= 0; ++field )//const auto field : fields::even_flds )
        {
            for( size_t cheby_index = 0; cheby_index < n_chebycffs; ++cheby_index )
            {
                if( cheby_index + field > 0 )
                { // Separate values using tabs
                    fprintf( outf, "\t");
                    //cout << "\t";
                }
                // here comes output of the coefficients
                fprintf( outf, "%g", spcoeffs[ ( n_all_flds * n_chebycffs ) * m
                                                + n_chebycffs * field
                                                    + cheby_index ] );
                //cout << spcoeffs[ ( n_all_flds * n_chebycffs ) * m
                 //                               + n_chebycffs * field
                 //                                   + cheby_index ];
            }
            //cout << "\n";
        }
        // Over spectral coefficients of the odd fields
        /*for( const auto field : fields::odd_flds )
        {
            for( size_t cheby_index = 0; cheby_index < n_chebycffs; ++cheby_index )
            {
                if( cheby_index + field > 0 )
                { // Separate values using tabs
                    fprintf( outf, "\t");
                    cout << "\t";
                }
                // here comes output of the coefficients
                fprintf( outf, "%g", spcoeffs[ ( n_all_flds * n_chebycffs ) * m
                                                + n_chebycffs * field
                                                    + cheby_index ] );
                cout << spcoeffs[ ( n_all_flds * n_chebycffs ) * m
                                                + n_chebycffs * field
                                                    + cheby_index ];
            }
            cout << "\n";
        }*/

        fprintf( outf, "\n\n" ); // Marks end of line
        //cout << "\n\n";

        if( m == mDim - 1 )
        {
            m = -1;
        }
    }

    // Close the output
    fclose( outf );
}
