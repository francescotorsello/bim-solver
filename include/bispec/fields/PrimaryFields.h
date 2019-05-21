/**
 *  @file      PrimaryFields.h
 *  @brief     Defines and computes the fields which access the values of Chebyshev pols.
 *  @authors   Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

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

    Real *values_fields;
    Real *values_ders;
    Real *values_reg_ders;
    Real *values_derrs;
    Real *values_colpoints;

    Real *fields_t;
    Real *chebyValues;
    Real *regDerCheby;
    Real *regDerrCheby;

    /// The array containing the spectral coefficients for all the Chebyshev series
    /// of the fields
    Real *spcoeffs;

public:

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

    /** Method to compute the values of the fields via parity-restricted Chebyshev series
      */

    inline void computeFields( size_t m )
    {

        /// Computation of the values of the fields at the collocation points (CPs)

        for( const auto field : fields::even_flds )
        {
            for( size_t n = 0; n < n_collocs; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += spectral_coeffs( m, field, cheby_index / 2 )
                            * chebyValues[ n_collocs * cheby_index + n ];
                }
                values_fields[ m * n_collocs * n_flds + n_collocs * field + n ] = sum;
            }
        } // spectral_coeffs is a function. Why not calling the array directly?
        for( const auto field : fields::odd_flds )
        {
            for( size_t n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += spectral_coeffs( m, field, ( cheby_index - 1 ) / 2 )
                            * chebyValues[ cheby_index * n_collocs + n ];
                }
                values_fields[ m * n_collocs * n_flds + n_collocs * field + n ] = sum;
            }
        }

        /// Computation of the values of the first radial derivatives at the
        /// collocation points (CPs) on the initial hypersurface
        for( const auto der : fields::even_ders )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += spectral_coeffs( m, der, cheby_index / 2 )
                            * chebyValues[ ( exp_ord + 1 ) * n_collocs
                                + cheby_index * n_collocs + n ];
                }
                values_ders[ m * n_collocs * n_flds + n_collocs * der + n ] = sum;
            }
        }
        for( const auto der : fields::odd_ders )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += spectral_coeffs( m, der, ( cheby_index - 1 ) / 2 )
                            * chebyValues[ ( exp_ord + 1 ) * n_collocs
                                + cheby_index * n_collocs + n ];
                }
                values_ders[ m * n_collocs * n_flds + n_collocs * der + n ] = sum;
            }
        }


        /// Computation of the values of the second radial derivatives at the
        /// collocation points (CPs) on the initial hypersurface
        for( const auto derr : fields::even_derrs )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += spectral_coeffs( m, derr, cheby_index / 2 )
                            * chebyValues[ 2*( exp_ord + 1 ) * n_collocs
                                + cheby_index * n_collocs + n ];
                }
                values_derrs[ m * n_collocs * n_flds + n_collocs * derr + n ] = sum;
            }
        }
        for( const auto derr : fields::odd_derrs )
        {
            for( n = 0; n < exp_ord + 1; ++n ) // loop over the collocation points
            {
                Real sum = 0;
                for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
                {
                    sum += spectral_coeffs( m, derr, ( cheby_index - 1 ) / 2 )
                            * chebyValues[ 2*( exp_ord + 1 ) * n_collocs
                                + cheby_index * n_collocs + n ];
                }
                values_derrs[ m * n_collocs * n_flds + n_collocs * derr + n ] = sum;
            }
        }
    }

    /** The following fields needs to access the spectral coefficients even if they
        are not primary. Hence, they need to be computed in here.
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
                    sum += spectral_coeffs( m, derr, cheby_index )
                            * chebyValues[ 2*( exp_ord + 1 ) * ( exp_ord + 1 )
                                + cheby_index * ( exp_ord + 1 ) + n ];
                }
                values_derrs[ ( exp_ord + 1 ) * derr + n ] = sum;
            }
        }*/

        /// TODO:make a unique loop for all the computations in computeRegDers

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += spectral_coeffs( m, fields::gBet, ( cheby_index + 1 ) / 2 )
                        * regDerCheby[ n +
                            cheby_index * n_collocs ];
            }
            values_reg_ders[ m * n_collocs * n_flds +
                n_collocs * 0 + n ] = sum;
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += spectral_coeffs( m, fields::fBet, ( cheby_index + 1 ) / 2 )
                        * regDerCheby[ n +
                            cheby_index * n_collocs ];
            }
            values_reg_ders[ m * n_collocs * n_flds +
                n_collocs * 1 + n ] = sum;
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += spectral_coeffs( m, fields::gL, ( cheby_index + 1 ) / 2 )
                        * regDerCheby[ n +
                            cheby_index * n_collocs ];
            }
            values_reg_ders[ m * n_collocs * n_flds +
                n_collocs * 2 + n ] = sum;
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 1; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += spectral_coeffs( m, fields::fL, ( cheby_index + 1 ) / 2 )
                        * regDerCheby[ n +
                            cheby_index * n_collocs ];
            }
            values_reg_ders[ m * n_collocs * n_flds +
                n_collocs * 3 + n ] = sum;
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += spectral_coeffs( m, fields::gAlp, cheby_index / 2 )
                        * regDerrCheby[ n +
                            cheby_index * n_collocs ];
            }
            values_reg_ders[ m * n_collocs * n_flds +
                n_collocs * 4 + n ] = sum;
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += spectral_coeffs( m, fields::fAlp, cheby_index / 2 )
                        * regDerrCheby[ n +
                            cheby_index * n_collocs ];
            }
            values_reg_ders[ m * n_collocs * n_flds +
                n_collocs * 5 + n ] = sum;
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += spectral_coeffs( m, fields::gconf, cheby_index / 2 )
                        * regDerrCheby[ n +
                            cheby_index * n_collocs ];
            }
            values_reg_ders[ m * n_collocs * n_flds +
                n_collocs * 6 + n ] = sum;
        }

        for( n = 0; n < n_collocs; ++n ) // loop over the collocation points
        {
            Real sum = 0;
            for( size_t cheby_index = 0; cheby_index < exp_ord + 1; cheby_index += 2 )
            {
                sum += spectral_coeffs( m, fields::fconf, cheby_index / 2 )
                        * regDerrCheby[ n +
                            cheby_index * n_collocs ];
            }
            values_reg_ders[ m * n_collocs * n_flds +
                n_collocs * 7 + n ] = sum;
        }

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

        values_colpoints = new Real[ n_collocs ];
        values_fields    = new Real[ mDim * n_flds * n_collocs ];
        values_ders      = new Real[ mDim * n_flds * n_collocs ];
        values_reg_ders  = new Real[ mDim * 8      * n_collocs ];
        values_derrs     = new Real[ mDim * n_flds * n_collocs ];
        spcoeffs         = new Real[ mDim * n_all_flds * n_chebycffs ];

        /// Assign initial values to the spectral coefficients

        for( size_t field = 0; field < n_all_flds; ++field )
        {
            for( size_t cheby_index = 0; cheby_index < n_chebycffs; ++cheby_index )
            {
                spcoeffs[ n_chebycffs * field + cheby_index ] =
                    bispecID( field, cheby_index );
            }
        }

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

        /*chebyValues  = new Real[ 4 * ( exp_ord + 1 ) * n_collocs ];
        regDerCheby  = new Real[ ( exp_ord + 1 ) * n_collocs ];
        regDerrCheby = new Real[ ( exp_ord + 1 ) * n_collocs ];

        for( size_t der_order = 0; der_order < 4; ++der_order )
        {
            for( size_t cheby_index = 0; cheby_index < ( exp_ord + 1 ); ++cheby_index )
            {
                for( size_t n = 0; n < n_collocs; ++n )
                {
                    chebyValues[ der_order * ( exp_ord + 1 ) * n_collocs
                        + cheby_index * n_collocs + n ] =
                            chebyC( der_order, cheby_index, n );
                }
            }
        }

        for( size_t cheby_index = 0; cheby_index < ( exp_ord + 1 ); ++cheby_index )
        {
            for( size_t n = 0; n < n_collocs; ++n )
            {
                regDerCheby[ cheby_index * n_collocs + n ] =
                        chebyC.reg_ders_cheby( cheby_index, n );
            }
        }

        for( size_t cheby_index = 0; cheby_index < ( exp_ord + 1 ); ++cheby_index )
        {
            for( size_t n = 0; n < n_collocs; ++n )
            {
                regDerrCheby[ cheby_index * n_collocs + n ] =
                        chebyC.reg_derrs_cheby( cheby_index, n );
            }
        }*/

        computeFields( 0 );
        computeRegDers( 0 );

        /// The printouts below print the values of the fields at the collocation points
        /// on the initial hypersurface. They are compared against the values in
        /// Mathematica and they coincide (up to MachinePrecision).

        std::cout << "The following are the values of the fields at the collocation "
                  << "points on the initial hypersurface (g-sector)," << std::endl;

        std::cout << std::endl;

        std::cout << "Collocation points," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << r( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Conformal factor," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gconf( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gtrK," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gtrK( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gB," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gB( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gA1," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA1( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gL," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gL( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gsig," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gsig( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAsig," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAsig( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gAlp," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAlp( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "gBet," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gBet( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "p," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << p( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfD," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfD( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pfS," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfS( 0, n ) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "pftau," << std::endl;
        std::cout << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pftau( 0, n ) << std::endl;
        }

    }
};
