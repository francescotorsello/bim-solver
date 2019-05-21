/**
 *  @file      BispecInput.h
 *  @brief     Loads the spectral coefficients on the initial hypersurface from a file.
 *  @authors   Francesco Torsello, Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

/**
 *  'bispecInput' reads the initial data, i.e., the values of the spectral coefficients
 *  on the initial hypersurface.
 *  TODO: merge with ChebyshevCoefficients?
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
        return spec_coeff[ number_chebycoeffs * field + n ];
    }

    // Constructor: Reads the spectral initial data from a file.
    //
    BispecInput( const std::string fileName )
    {
        FILE* specID = fopen( fileName.c_str(), "rb" );

        if( ! specID ) {
            std::cerr << "err: ID: Cannot open file" << std::endl;
            return;
        }

        // Read magic number
        //
        Real x = 0;
        if( 1 != fread( &x, sizeof(Real), 1, specID ) || x != 4478245647 )
        {
            std::cerr << "err: ID: Magic number wrong" << std::endl;
            fclose( specID );
            return;
        }

        // Read the order of Chebyshev expansion
        //
        if( 1 == fread( &x, sizeof(Real), 1, specID ) ) {
            expansion_order     = size_t(x);
            number_collocations = size_t(( x + 1 ) / 2);
            number_chebycoeffs  = number_collocations;
        }
        else {
            std::cerr << "err: ID: Cannot read expansion_order" << std::endl;
            fclose( specID );
            return;
        }

        // Read the number of evolved fields
        //
        if( 1 == fread( &x, sizeof(Real), 1, specID ) ) {
            input_fields    = size_t(x);
            /// The number of evolved fields is the number of input fields,
            /// minus the number of gauge variables ( at present, 6 ) minus p minus r
            number_fields   = input_fields - number_gauge_vars - 1 - 1;
        }
        else {
            std::cerr << "err: ID: Cannot read number_fields" << std::endl;
            fclose( specID );
            return;
        }

        // Read the spectral initial data
        //
        data_size   = input_fields * number_chebycoeffs;
        spec_coeff  = new Real[ data_size ];

        if ( data_size != fread( spec_coeff, sizeof(Real), data_size, specID ) )
        {
            std::cerr << "err: ID: Cannot read all the initial data" << std::endl;
            delete [] spec_coeff;
            spec_coeff = nullptr;
            fclose( specID );
            return;
        }

        fclose(specID);
    }

};
