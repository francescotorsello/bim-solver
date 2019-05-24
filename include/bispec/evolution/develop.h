/**
 *  @file      develop.h
 *  @brief     Method in BispecEvolve to perform the time integration.
 *  @authors   Francesco Torsello, Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

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
        cout << endl << "---------------------Time step m = " << m << endl << endl;
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

        /*std::cout << std:: endl << "Inside develop()" << std::endl << std::endl;

        for( Int n = 0; n < n_collocs; ++n )
        {
            std::cout << "values_fields[ " << next_m << ", " << "fields::gAlp, "
                      << n << " ] = ";
            std::cout << values_fields[ ( n_all_flds * n_collocs ) * next_m
                    + n_collocs * fields::gAlp + n ] << ", " << "\t";
            std::cout << "gAlp(" << next_m << "," << n << ") = ";
            std::cout << gAlp( next_m, n ) << std::endl;
        }*/

        //dumpPrimaryFields( m );
        dumpDependentFields( m );

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
