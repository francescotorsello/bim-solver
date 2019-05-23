/**
 *  @file      PrimaryFields.h
 *  @brief     Defines and computes the fields which access the values of Chebyshev pols.
 *  @authors   Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

/**
 *  The macros define methods to access the values of the fields.
 */

#define setfield( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_fields[ ( n_flds * n_collocs ) * m \
            + n_collocs * fields::field + n ]; \
        }

#define setder( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_ders[ ( n_flds * n_collocs ) * m \
            + n_collocs * fields::field + n ]; \
        }

#define setderr( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_derrs[ ( n_flds * n_collocs ) * m \
            + n_collocs * fields::field + n ]; \
        }

#define setregder( field ) \
        inline Real field( Int m, Int n ) \
        { \
            return values_reg_ders[ ( n_reg_ders * n_collocs ) * m \
            + n_collocs * fields::field + n ]; \
        }

/**
 *  'PrimaryFields' defines the methods to access the values of the fields accessing the
 *  values of the Chebyshev polynomials, and their first and second radial derivatives.
 *
 *  The notation is: field( m, n ), with m being the time step and n being the index of
 *  the collocation point.
 *
 *  The constructor assigns the initial data to the fields, i.e., it assigns values to
 *  field( 0, n ) for all n.
 */

class PrimaryFields
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
    size_t n_lapses = 3;
    size_t n_shifts = 3;
    size_t n_reg_ders = 8;

    Int ini_cheby;

    Real *values_reg_ders;

    Real *values_colpoints;

    Real *fields_t;

    Real *chebyValues;
    Real *regDerCheby;
    Real *regDerrCheby;

    /// The array containing the spectral coefficients for all the Chebyshev series
    /// of the fields
    Real *spcoeffs;

public:

    Real *values_fields;
    //Real *values_lapses;
    //Real *values_shifts;

    Real *values_ders;
    //Real *values_lapses_r;
    //Real *values_shifts_r;

    Real *values_derrs;
    //Real *values_lapses_rr;
    //Real *values_shifts_rr;

    /// Method to access the spectral coefficients. Since we will make use of the
    /// integrator from bim-solver, we copy the structure of GF in gridDriver.
    Real spectral_coeffs( size_t m, size_t field, size_t cheby_index )
    {
        return spcoeffs[ ( n_all_flds * n_chebycffs ) * m + n_chebycffs * field
            + cheby_index ];
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

    /** Methods to compute the values of the fields via parity-restricted Chebyshev series
      */

    inline void computeFields( size_t m )
    {

        /// Computation of the values of the fields at the collocation points (CPs)

        for( const auto field : fields::even_flds )
        {
            for( size_t n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_fields[ m * n_collocs * n_all_flds + n_collocs * field + n ] = 0;

                for( size_t cheby_index = 0; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_fields[ m * n_collocs * n_all_flds + n_collocs * field + n ]
                        += spectral_coeffs( m, field, cheby_index / 2 )
                            * chebyValues[ n_collocs * cheby_index + n ];
                }
            }
        } // TODO: spectral_coeffs is a function. Why not calling the array directly?

        for( const auto field : fields::odd_flds )
        {
            for( size_t n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_fields[ m * n_collocs * n_all_flds + n_collocs * field + n ] = 0;

                for( size_t cheby_index = 1; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_fields[ m * n_collocs * n_all_flds + n_collocs * field + n ]
                        += spectral_coeffs( m, field, ( cheby_index - 1 ) / 2 )
                            * chebyValues[ cheby_index * n_collocs + n ];
                }
            }
        }

        for( const auto field : fields::lapses )
        {
            for( size_t n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_fields[ m * n_collocs * n_all_flds + n_collocs * field + n ] = 0;

                for( size_t cheby_index = 0; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_fields[ m * n_collocs * n_all_flds + n_collocs * field + n ]
                        += spectral_coeffs( m, field, cheby_index / 2 )
                            * chebyValues[ n_collocs * cheby_index + n ];
                }
            }
        }

        for( const auto field : fields::shifts )
        {
            for( size_t n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_fields[ m * n_collocs * n_all_flds + n_collocs * field + n ] = 0;

                for( size_t cheby_index = 1; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_fields[ m * n_collocs * n_all_flds + n_collocs * field + n ]
                        += spectral_coeffs( m, field, ( cheby_index - 1 ) / 2 )
                            * chebyValues[ cheby_index * n_collocs + n ];
                }
            }
        }

        /// Computation of the values of the first radial derivatives at the
        /// collocation points (CPs) on the initial hypersurface

        for( const auto der : fields::even_ders )
        {
            for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_ders[ m * n_collocs * n_all_flds + n_collocs * der + n ] = 0;

                for( size_t cheby_index = 0; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_ders[ m * n_collocs * n_all_flds + n_collocs * der + n ]
                        += spectral_coeffs( m, der, cheby_index / 2 )
                            * chebyValues[ ( exp_ord + 1 ) * n_collocs
                                + cheby_index * n_collocs + n ];
                }
            }
        }
        for( const auto der : fields::odd_ders )
        {
            for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_ders[ m * n_collocs * n_all_flds + n_collocs * der + n ] = 0;

                for( size_t cheby_index = 1; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_ders[ m * n_collocs * n_all_flds + n_collocs * der + n ]
                        += spectral_coeffs( m, der, ( cheby_index - 1 ) / 2 )
                            * chebyValues[ ( exp_ord + 1 ) * n_collocs
                                + cheby_index * n_collocs + n ];
                }
            }
        }
        for( const auto field : fields::lapses_r )
        {
            for( size_t n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_ders[ m * n_collocs * n_all_flds + n_collocs * field + n ] = 0;

                for( size_t cheby_index = 0; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_ders[ m * n_collocs * n_all_flds + n_collocs * field + n ]
                        += spectral_coeffs( m, field, cheby_index / 2 )
                            * chebyValues[ ( exp_ord + 1 ) * n_collocs
                                + n_collocs * cheby_index + n ];
                }
            }
        }
        for( const auto field : fields::shifts_r )
        {
            for( size_t n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_ders[ m * n_collocs * n_all_flds + n_collocs * field + n ] = 0;

                for( size_t cheby_index = 1; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_ders[ m * n_collocs * n_all_flds + n_collocs * field + n ]
                        += spectral_coeffs( m, field, ( cheby_index - 1 ) / 2 )
                            * chebyValues[ ( exp_ord + 1 ) * n_collocs
                                + cheby_index * n_collocs + n ];
                }
            }
        }


        /// Computation of the values of the second radial derivatives at the
        /// collocation points (CPs) on the initial hypersurface

        for( const auto derr : fields::even_derrs )
        {
            for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_derrs[ m * n_collocs * n_all_flds + n_collocs * derr + n ] = 0;

                for( size_t cheby_index = 0; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_derrs[ m * n_collocs * n_all_flds + n_collocs * derr + n ]
                        += spectral_coeffs( m, derr, cheby_index / 2 )
                            * chebyValues[ 2*( exp_ord + 1 ) * n_collocs
                                + cheby_index * n_collocs + n ];
                }

            }
        }
        for( const auto derr : fields::odd_derrs )
        {
            for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_derrs[ m * n_collocs * n_all_flds + n_collocs * derr + n ] = 0;

                for( size_t cheby_index = 1; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_derrs[ m * n_collocs * n_all_flds + n_collocs * derr + n ]
                        += spectral_coeffs( m, derr, ( cheby_index - 1 ) / 2 )
                            * chebyValues[ 2*( exp_ord + 1 ) * n_collocs
                                + cheby_index * n_collocs + n ];
                }
            }
        }
         for( const auto field : fields::lapses_rr )
        {
            for( size_t n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_derrs[ m * n_collocs * n_all_flds + n_collocs * field + n ] = 0;

                for( size_t cheby_index = 0; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_derrs[ m * n_collocs * n_all_flds + n_collocs * field + n ]
                        += spectral_coeffs( m, field, cheby_index / 2 )
                            * chebyValues[ 2*( exp_ord + 1 ) * n_collocs
                                + n_collocs * cheby_index + n ];
                }
            }
        }
        for( const auto field : fields::shifts_rr )
        {
            for( size_t n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_ders[ m * n_collocs * n_all_flds + n_collocs * field + n ] = 0;

                for( size_t cheby_index = 1; cheby_index < exp_ord; cheby_index += 2 )
                {
                    values_derrs[ m * n_collocs * n_all_flds + n_collocs * field + n ]
                        += spectral_coeffs( m, field, ( cheby_index - 1 ) / 2 )
                            * chebyValues[ 2*( exp_ord + 1 ) * n_collocs
                                + cheby_index * n_collocs + n ];
                }
            }
        }

    }

    /** Methods to access the regularized derivatives
      */

    inline void computeRegDers( size_t m )
    {
        for( const auto regder : fields::reg_ders ) // for each regularized derivative
            // regder goes from 0 to 3
        {
            if( /*// If the field fields::flds_reg_ders[ regder ] is even,
                std::find(
                          std::begin(fields::even_flds),
                          std::end  (fields::even_flds),
                          fields::flds_reg_ders[ regder ]
                          // fields::flds_reg_ders[ regder ] contains the fields
                          // whose derivatives are regularized and in fields::reg_ders
                          ) != std::end  (fields::even_flds)*/
                fields::parity_fields_reg_ders[ regder ]
               )
            {
                // do this,
                ini_cheby = 0;

            } else { // otherwise (i.e., if it's odd) do this,

                ini_cheby = 1;
            }

            for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_reg_ders[ m * n_collocs * n_reg_ders
                                 + n_collocs * regder + n ] = 0;

                for( size_t cheby_index = ini_cheby; cheby_index < exp_ord;
                    cheby_index += 2 )
                {
                    values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * regder + n ]

                        += spectral_coeffs( m,
                                            fields::flds_reg_ders[ regder ],
                                            ( cheby_index - ini_cheby ) / 2
                                          )

                            * regDerCheby[ n_collocs * cheby_index + n ];

                }
            }
        }

        for( const auto regder : fields::reg_derrs ) // for each regularized derivative
            // regder goes from 4 to 7
        {
            if(
                fields::parity_fields_reg_ders[ regder ]
               )
            {
                // do this,
                ini_cheby = 0;

            } else { // otherwise (i.e., if it's odd) do this,

                ini_cheby = 1;
            }

            for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                values_reg_ders[ m * n_collocs * n_reg_ders
                                 + n_collocs * regder + n ] = 0;

                for( size_t cheby_index = ini_cheby; cheby_index < exp_ord;
                    cheby_index += 2 )
                {
                    values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * regder + n ]

                        += spectral_coeffs( m,
                                            fields::flds_reg_ders[ regder ],
                                            ( cheby_index - ini_cheby ) / 2
                                          )

                            * regDerrCheby[ n_collocs * cheby_index + n ];

                }
            }
        }

        /// TODO:make a unique loop for all the computations in computeRegDers

        /*for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 0 + n ] = 0;

            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 0 + n ]
                    += spectral_coeffs( m, fields::gBet, ( cheby_index + 1 ) / 2 )
                        * regDerCheby[ n +
                            cheby_index * n_collocs ];
            }
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 1 + n ] = 0;

            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 1 + n ]
                    += spectral_coeffs( m, fields::fBet, ( cheby_index + 1 ) / 2 )
                        * regDerCheby[ n +
                            cheby_index * n_collocs ];
            }
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 2 + n ] = 0;

            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 2 + n ]
                    += spectral_coeffs( m, fields::gL, ( cheby_index + 1 ) / 2 )
                        * regDerCheby[ n +
                            cheby_index * n_collocs ];
            }
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 3 + n ] = 0;

            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 3 + n ]
                    += spectral_coeffs( m, fields::fL, ( cheby_index + 1 ) / 2 )
                        * regDerCheby[ n +
                            cheby_index * n_collocs ];
            }
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 4 + n ] = 0;

            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 4 + n ]
                    += spectral_coeffs( m, fields::gAlp, cheby_index / 2 )
                        * regDerrCheby[ n +
                            cheby_index * n_collocs ];
            }
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 5 + n ] = 0;

            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 5 + n ]
                    += spectral_coeffs( m, fields::fAlp, cheby_index / 2 )
                        * regDerrCheby[ n +
                            cheby_index * n_collocs ];
            }
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 6 + n ] = 0;

            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 6 + n ]
                    += spectral_coeffs( m, fields::gconf, cheby_index / 2 )
                        * regDerrCheby[ n +
                            cheby_index * n_collocs ];
            }
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 7 + n ] = 0;

            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                values_reg_ders[ m * n_collocs * n_reg_ders + n_collocs * 7 + n ]
                    += spectral_coeffs( m, fields::fconf, cheby_index / 2 )
                        * regDerrCheby[ n +
                            cheby_index * n_collocs ];
            }
        }*/

    }

    inline void dumpPrimaryFields( Int m )
    {

        std::cout << std::endl<< std::endl
                  << "-------This is the method dumpPrimaryFields."
                  << std::endl << std::endl << std::endl;

        std::cout << std::endl << std::endl
                  << "The following are the values of the primary fields at the "
                  << "collocation points at time step " << m << "." << std::endl;

        std::cout << std::endl;

        std::cout << "r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gconf," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gconf( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gtrK," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gtrK( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gB," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gB( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA1," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA1( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gL," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gL( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gsig," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gsig( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAsig," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAsig( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAlp," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAlp( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gBet," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gBet( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "p," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << p( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfD," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfD( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfS," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfS( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pftau," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pftau( m, n ) << std::endl;
        }

        std::cout << std::endl;
        std::cout << "The following are the values of the first derivatives of the fields"
                  << " at the collocation points at time step " << m << "."
                  << std::endl;

        std::cout << std::endl;

        std::cout << "gconf_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gconf_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gtrK_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gtrK_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gB_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gB_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA1_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA1_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gL_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gL_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gsig_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gsig_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAsig_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAsig_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAlp_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAlp_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gBet_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gBet_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "p_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << p_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfD_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfD_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfS_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfS_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pftau_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pftau_r( m, n ) << std::endl;
        }

        std::cout << std::endl;
        std::cout << "The following are the values of the second ders of the fields"
                  << " at the collocation points at time step " << m << "."
                  << std::endl;

        std::cout << std::endl;

        std::cout << "gconf_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gconf_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gtrK_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gtrK_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gB_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gB_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA1_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA1_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gL_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gL_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gsig_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gsig_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAsig_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAsig_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAlp_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAlp_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gBet_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gBet_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "p_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << p_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfD_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfD_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfS_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfS_rr( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pftau_rr," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pftau_rr( m, n ) << std::endl;
        }

        std::cout << std::endl;
        std::cout << "The following are the values of the regularized ders of the fields"
                  << " at the collocation points at time step " << m << "."
                  << std::endl;

        std::cout << std::endl;

        std::cout << "gBetr_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gBetr_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gLr_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gLr_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gderAlpr_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gderAlpr_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gderconfr_r," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gderconfr_r( m, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAlp," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAlp( m, n ) << std::endl;
        }
        std::cout << std::endl;

        std::cout << std::endl << "-------End of the method dumpPrimaryFields."
                  << std::endl << std::endl << std::endl;
    }

    /** The constructor computes the values of the fields, the derivatives and the
        evolution equations on the initial hypersurface
     */
    PrimaryFields(
        BispecInput&            bispecID,
        ChebyshevCoefficients&  chebyC
    )
    {
        mLen        = 5;
        mExtra      = 0;
        mDim        = mLen + mExtra;
        exp_ord     = bispecID.exp_order();
        n_collocs   = bispecID.n_collocations();
        n_chebycffs = bispecID.n_chebycoeffs();
        n_all_flds  = bispecID.n_allfields();
        n_flds      = bispecID.n_fields();
        n_gauge     = bispecID.n_gauges();
        n_vlc       = bispecID.n_valencia();

        /// These arrays contain the values of the fields and their spatial first and
        /// second derivatives at the collocation points at each time step

        /// TODO: split even and odd fields at this level.

        values_colpoints = new Real[ n_collocs ];

        spcoeffs         = new Real[ mDim * n_all_flds * n_chebycffs ];

        values_fields    = new Real[ mDim * n_all_flds * n_collocs ];
        //values_lapses    = new Real[ mDim * 3      * n_collocs ];
        //values_shifts    = new Real[ mDim * 3      * n_collocs ];

        values_ders      = new Real[ mDim * n_all_flds * n_collocs ];
        //values_lapses_r  = new Real[ mDim * 3      * n_collocs ];
        //values_shifts_r  = new Real[ mDim * 3      * n_collocs ];

        values_derrs     = new Real[ mDim * n_all_flds * n_collocs ];
        //values_lapses_rr = new Real[ mDim * 3      * n_collocs ];
        //values_shifts_rr = new Real[ mDim * 3      * n_collocs ];

        values_reg_ders  = new Real[ mDim * n_reg_ders * n_collocs ];

        /// Assign initial values to the spectral coefficients

        for( const auto field : fields::all_flds )
        {
            for( size_t cheby_index = 0; cheby_index < n_chebycffs; ++cheby_index )
            {
                spcoeffs[ n_chebycffs * field + cheby_index ] =
                    bispecID( field, cheby_index );
            }
        }

        /*std::cout << "The following is the spectral initial data in PrimaryFields"
                    << std::endl;

        std::cout << bispecID.n_allfields() << " x " << bispecID.n_chebycoeffs()
                  << std::endl;

        for( size_t field = 0; field < bispecID.n_allfields(); ++field )
        {
            std::cout << std::endl << "These are the spectral coefficients of the field "
                      << field << std::endl << std::endl;
            for( size_t cheby_index = 0; cheby_index < bispecID.n_chebycoeffs();
                ++cheby_index )
            {
                std::cout << "(" << field << "," << cheby_index << ") = "
                          << bispecID( field, cheby_index ) << ", "
                          << spcoeffs[ n_chebycffs * field + cheby_index ] << ", "
                          << spectral_coeffs( 0, field, cheby_index )
                          << std::endl;
            }

        }
        std::cout << std::endl;*/

        /// Assign values to the collocation points

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            values_colpoints[ n ] = bispecID( n_all_flds - 1, n );
        }

        /** Define pointers to access the values of the Chebyshev polynomials at the
            collocation points
          */

        chebyValues  = chebyC.coeff;
        regDerCheby  = chebyC.reg_der_coeff;
        regDerrCheby = chebyC.reg_derr_coeff;

        /** Compute values of all the fields on the initial hypersurface
          */

        computeFields( 0 );
        computeRegDers( 0 );

    }
};
