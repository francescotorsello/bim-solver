/**
 *  @file      DependentFields.h
 *  @brief     Defines the fields which do not access the values of the Chebyshev pols.
 *  @authors   Francesco Torsello, Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

/**
 *  'dependentFields' defines the methods to access the dependent (nonprimary) fields
 *  appearing in the evolution equations. If the primary fields have an assigned value
 *  at the time step m, then the dependent fields will also have a defined value.
 */

class DependentFields
    : virtual PrimaryFields,
      BimetricModel
{

public:

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

    #include "../eom-BSSN/eomBSSNRicciComp.h"

    /// Dependent fields which other dependent fields depend on

    inline Real R  ( Int m, Int n )
    {
        Real x;
        isGR() ? x = 1 : x = (fB(m,n)*exp(2*fconf(m,n))) / (gB(m,n)*exp(2*gconf(m,n)));
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

    inline Real gK1 ( Int m, Int n )    { return gA1(m,n) + 1/3 * gtrK(m,n); }

    inline Real gK2 ( Int m, Int n )    { return gA2(m,n) + 1/3 * gtrK(m,n); }

    inline Real fK1 ( Int m, Int n )    { return fA1(m,n) + 1/3 * ftrK(m,n); }

    inline Real fK2 ( Int m, Int n )    { return fA2(m,n) + 1/3 * ftrK(m,n); }

    /// Regularizing fields

    inline Real gBetr( Int m, Int n )   { return - r_minus(m,n) * gBet(m,n)
        / ( r(m,n) + TINY_Real ); }

    inline Real fBetr( Int m, Int n )   { return - r_minus(m,n) * fBet(m,n)
        / ( r(m,n) + TINY_Real ); }

    /// Last dependent fields

    #include "../eom-BSSN/eomBSSNMatterComp.h"
    #include "../eom-BSSN/eomBSSNSourcesCompEul.h"
    #include "../eom-BSSN/eomBSSNEvolutionCompEul.h"

    inline void dumpDependentFields ( Int m )
    {
        std::cout << std::endl<< std::endl
                  << "-------This is the method dumpDependentFields."
                  << std::endl << std::endl << std::endl;

        /// The printouts below print the values of the dependent fields at the collocation points on the initial hypersurface.

        std::cout << "The following are the values of the dependent fields at the collocation points at time step " << m << "." << std::endl << std::endl;

        std::cout << "The collocation points, "<< std::endl << std::endl;
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << r( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "The Ricci terms," << std::endl << std::endl;

        std::cout << "gRicci:   ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gRicci( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fRicci:   ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fRicci( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Lorentz factor, extrinsic curvatures and auxiliary field," << std::endl << std::endl;

        std::cout << "R:   ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << R( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "Lt:  ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << Lt( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "Lt2: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << Lt2( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gA2: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA2( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gK1: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gK1( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gK2: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gK2( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fA2: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fA2( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fK1: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fK1( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fK2: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fK2( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "The matter fields," << std::endl << std::endl;

        std::cout << "pfv:   ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfv( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "pfv_r: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << pfv_r( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gD_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gD_t( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gS_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gS_t( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gtau_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gtau_t( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gW: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gW( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "grhobar: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << grhobar( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "grho: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << grho( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gj: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gj( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gJ11: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gJ11( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gJ22: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gJ22( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "frho: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << frho( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fj: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fj( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fJ11: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fJ11( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fJ22: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fJ22( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;;

        std::cout << "The sources in the evolution equations at the collocation points,"
            << std::endl << std::endl;

        std::cout << "gJK: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gJK( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gJA1: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gJA1( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "gJL: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gJL( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fJK: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fJK( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fJA1: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fJA1( m, n ) << ", ";
        }
        std::cout << std::endl;

        std::cout << "fJL: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fJL( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "The evolution equations," << std::endl << std::endl;

        std::cout << "gconf_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gconf_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gtrK_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gtrK_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gA_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gB_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gB_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gA1_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gA1_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gL_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gL_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gsig_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gsig_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gAsig_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gAsig_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fconf_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fconf_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "ftrK_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << ftrK_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fA_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fA_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fB_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fB_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fA1_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fA1_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fL_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fL_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fsig_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fsig_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "fAsig_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << fAsig_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gD_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gD_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gS_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gS_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "gtau_t: ";
        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gtau_t( m, n ) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Other fields," << std::endl << std::endl;

        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << (gBet(m,n) * (2 - 2 * r(m,n))) / (1 + TINY_Real
                + r(m,n)) << ", ";
        }
        std::cout << std::endl << std::endl;

        for( size_t n = 0; n < n_collocs; ++n )
        {
            std::cout << gBet(m,n) << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << std::endl;

        std::cout << std::endl << "-------End of the method dumpDependentFields."
                  << std::endl << std::endl << std::endl;
    }

    DependentFields(
        BispecInput&            bispecID,
        ChebyshevCoefficients&  chebyC,
        Parameters&             params
    ) :
        PrimaryFields( bispecID, chebyC ),
        BimetricModel( params )
    {

        dumpDependentFields( 0 );

    }
};
