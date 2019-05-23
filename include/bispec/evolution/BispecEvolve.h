/**
 *  @file      BispecEvolve.h
 *  @brief     Sets up and performs the time integration of the spectral coefficients.
 *  @authors   Francesco Torsello, Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

/**
 *  The evolution of the gauge variables is not included yet in 'BispecEvolve'.
 *  It must be handled separately because it depends on the choice of the user.
 */

class BispecEvolve
    : virtual public PrimaryFields,
      DependentFields
{

    Int m;
    Int n;
    Real *even_spec_t;
    Real *odd_spec_t;

    size_t exp_ord;
    size_t n_collocs;
    size_t n_chebycffs;

    typedef Real (DependentFields::*FP)( Int m, Int n );
    std::vector<FP> even_fields_t   = { EVEN_FIELDS_T };
    std::vector<FP> odd_fields_t    = { ODD_FIELDS_T };

    Real *b_even;
    Real *b_odd;
    const Real *A_even;
    const Real *A_odd;

public:

    inline size_t get_n_all_flds()  { return n_all_flds; }
    inline size_t n_evenflds()  { return even_fields_t.size(); }
    inline size_t n_oddflds()   { return odd_fields_t.size(); }
    inline size_t n_colpoints() { return n_collocs; }
    inline size_t n_chebys()    { return n_chebycffs; }

    //////////////////////////////////////////////////////////////////////////////////////
    /////   Methods to access values
    //////////////////////////////////////////////////////////////////////////////////////

    /** Method to access the vector b_even
      */

    inline Real even_evolution_eqs( size_t m, size_t field, size_t n )
    {
        return b_even[ even_fields_t.size() * n_collocs * m
            + n_collocs * field + n ];
    }

    /** Method to access the vector b_odd
      */

    inline Real odd_evolution_eqs( size_t m, size_t field, size_t n )
    {
        return b_odd[ odd_fields_t.size() * n_collocs * m
            + n_collocs * field + n ];
    }

    /** Method to access the time derivatives of the even spectral coefficients
      */

    inline Real get_even_spec_t( size_t m, size_t field, size_t cheby_index )
    {
        return even_spec_t[ ( n_collocs * even_fields_t.size() ) * m
            + n_chebycffs * field + cheby_index ];
    }

    /** Method to access the time derivatives of the odd spectral coefficients
      */

    inline Real get_odd_spec_t( size_t m, size_t field, size_t cheby_index )
    {
        return odd_spec_t[ ( n_collocs * odd_fields_t.size() ) * m
            + n_chebycffs * field + cheby_index ];
    }

    //////////////////////////////////////////////////////////////////////////////////////
    /////   Methods to compute the derivatives of the spectral coefficients
    //////////////////////////////////////////////////////////////////////////////////////

    /** Method to fill the vectors b_even and b_odd
      */

    inline void arrange_fields_t( size_t m )
    {
        for( size_t field = 0; field < even_fields_t.size(); ++field )
        {
            for( size_t n = 0; n < n_collocs; ++n )
            {
                b_even[ ( n_collocs * even_fields_t.size() ) * m
                    + n_collocs * field + n ] =
                        ( this ->* even_fields_t[ field ] )( m, n );
            }
        }

        for( size_t field = 0; field < odd_fields_t.size(); ++field )
        {
            for( size_t n = 0; n < n_collocs; ++n )
            {
                b_odd[ ( n_collocs * odd_fields_t.size() ) * m
                    + n_collocs * field + n ] =
                        ( this ->* odd_fields_t[ field ] )( m, n );
            }
        }
    }

    /** Method to solve for the time derivatives of the spectral coefficients
      */

    inline void solveSpectralDerivatives( size_t m )
    {
        // For each even field
        for( size_t field = 0; field < even_fields_t.size(); ++field )
        { // Solve the linear algebraic system, i.e.,
            // for each even Chebyshev index (row)
            for( size_t cheby_index = 0; cheby_index < n_chebycffs; ++cheby_index )
            {
                even_spec_t[ ( n_chebycffs * even_fields_t.size() ) * m
                    + n_chebycffs * field + cheby_index ] = 0;
                // For each collocation point (column)
                for( size_t n = 0; n < n_collocs; ++n )
                {
                    even_spec_t[ ( n_chebycffs * even_fields_t.size() ) * m
                    + n_chebycffs * field + cheby_index ] +=
                            A_even[ cheby_index * n_collocs + n ]
                                * even_evolution_eqs( m, field, n );
                }
            }
        }

        // For each odd field
        for( size_t field = 0; field < odd_fields_t.size(); ++field )
        { // Solve the linear algebraic system, i.e.,
            // for each odd Chebyshev index (row)
            for( size_t cheby_index = 0; cheby_index < n_chebycffs; ++cheby_index )
            {
                odd_spec_t[ ( n_chebycffs * odd_fields_t.size() ) * m
                    + n_chebycffs * field + cheby_index ] = 0;
                // For each collocation point (column)
                for( size_t n = 0; n < n_collocs; ++n )
                {
                    odd_spec_t[ ( n_chebycffs * odd_fields_t.size() ) * m
                    + n_chebycffs * field + cheby_index ] +=
                            A_odd[ cheby_index * n_collocs + n ]
                                * odd_evolution_eqs( m, field, n );
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    /////   Methods to evolve the spectral coefficients
    //////////////////////////////////////////////////////////////////////////////////////

    inline Real EulerMethod( Real value, Real der, Real delta_t )
    {
        return value + der * delta_t;
    }

    inline void evolveEvenFields( Real (BispecEvolve::*method)( Real, Real, Real ),
                                    size_t m, size_t next_m, Real delta_t )
    {
        //cout << "Inside method evolveEvenFields." << endl << endl;
        // counter select the field in the vector containing the derivatives
        Real counter = 0;
        for( const auto field : fields::even_flds ) // for each even field
        {
            for( size_t n = 0; n < n_collocs; ++n ) // for each collocation point
            {
                for( size_t cheby_index; cheby_index < n_chebycffs; ++cheby_index )
                        // for each spectral coefficients
                {
                    spcoeffs[ ( n_all_flds * n_chebycffs ) * next_m
                                    + n_chebycffs * field
                                        + cheby_index ]

                    = ( this ->* method )(

                        spcoeffs[ ( n_all_flds * n_chebycffs ) * m
                                    + n_chebycffs * field
                                        + cheby_index ], //value

                        get_even_spec_t( m, counter, cheby_index ), //derivative

                        delta_t ); // time step

                    //cout << spcoeffs[ ( n_all_flds * n_chebycffs ) * next_m
                      //              + n_chebycffs * field
                      //                  + cheby_index ];
                }
            }
            ++counter;
        }
    }

    inline void evolveOddFields( Real (BispecEvolve::*method)( Real, Real, Real ),
                                    size_t m, size_t next_m, Real delta_t )
    {
        // counter select the field in the vector containing the derivatives
        Real counter = 0;
        for( const auto field : fields::odd_flds ) // for each odd field
        {
            for( size_t n = 0; n < n_collocs; ++n ) // for each collocation point
            {
                for( size_t cheby_index; cheby_index < n_chebycffs; ++cheby_index )
                    // for each spectral coefficients
                {
                    spcoeffs[ ( n_all_flds * n_chebycffs ) * next_m
                                    + n_chebycffs * field
                                        + cheby_index ]

                    = ( this ->* method )(

                        spcoeffs[ ( n_all_flds * n_chebycffs ) * m
                                    + n_chebycffs * field
                                        + cheby_index ], //value

                        get_odd_spec_t( m, counter, cheby_index ), //derivative

                        delta_t ); // time step
                }
            }
            ++counter;
        }
    }

    void develop( class GaugeVariables& );

    //////////////////////////////////////////////////////////////////////////////////////
    /////   Constructor
    /////////////////////////////////////////////////////////////////////////////////////

    /** The constructor perform the first step
      */

    BispecEvolve(
        BispecInput&            bispecID,
        ChebyshevCoefficients&  chebyC,
        Parameters&             params
    ) :
        PrimaryFields       ( bispecID, chebyC ),
        DependentFields     ( bispecID, chebyC, params )
    {

        exp_ord         = bispecID.exp_order();
        n_collocs       = bispecID.n_collocations();
        n_chebycffs     = bispecID.n_chebycoeffs();

        b_even      = new Real[ mDim * even_fields_t.size() * n_collocs ];
        b_odd       = new Real[ mDim * odd_fields_t.size()  * n_collocs ];
        even_spec_t = new Real[ mDim * even_fields_t.size() * n_collocs ];
        odd_spec_t  = new Real[ mDim * odd_fields_t.size()  * n_collocs ];

        // Access to the evolution matrices via pointers

        A_even      = chebyC.ee_matrix_even;
        A_odd       = chebyC.ee_matrix_odd;

        // Perform the first integration step (testing)

        /*Real t = 0;
        Real delta_t = 0.5;

        arrange_fields_t( 0 );
        solveSpectralDerivatives( 0 ); // integStep_Prepare ends

        evolveEvenFields( &EulerMethod, 0, 1, delta_t );
        evolveOddFields ( &EulerMethod, 0, 1, delta_t );

        computeFields   ( 1 );
        computeRegDers  ( 1 );*/

        /*
        std::cout << std::endl << "-------This is the constructor of BispecEvolve."
                  << std::endl << std::endl << std::endl;

        std::cout << "The evolution equations stored in the vectors," << std::endl
            << std::endl;

        std::cout << "gconf_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << even_evolution_eqs( 0, 0, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gtrK_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << even_evolution_eqs( 0, 2, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gA_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << even_evolution_eqs( 0, 4, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gB_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << even_evolution_eqs( 0, 5, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gA1_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << even_evolution_eqs( 0, 8, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gL_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << odd_evolution_eqs( 0, 1, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << std::endl;

        Int m1 = 1;

        std::cout << "The following are the values of the fields at the collocation "
                  << "points on the initial hypersurface (g-sector)," << std::endl;

        std::cout << std::endl;

        std::cout << "Collocation points," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << r( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Conformal factor," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gconf( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gtrK," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gtrK( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gB," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gB( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA1," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA1( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gL," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gL( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gsig," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gsig( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAsig," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAsig( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAlp," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAlp( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gBet," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gBet( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "p," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << p( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfD," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfD( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfS," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfS( m1, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pftau," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pftau( m1, n ) << std::endl;
        }

        std::cout << std::endl << "-------End of the constructor of BispecEvolve."
                  << std::endl << std::endl << std::endl;*/

    }

};
