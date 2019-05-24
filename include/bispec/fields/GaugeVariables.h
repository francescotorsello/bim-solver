/**
 *  @file      GaugeVariables.h
 *  @brief     Computes the lapses and shifts at each time slice.
 *  @authors   Francesco Torsello, Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

/**
 *  'GaugeVariables' is an abstract class serving as base for the different gauges
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

/**
 *  'GeodesicSlicing' inherits from GaugeVariables and defines the method
 *  applyGaugeCondition to compute the geodesic slicing.
 */

class GeodesicSlicing
    : public GaugeVariables
{

public:

    inline void applyGaugeCondition( Int m )
    {
        std::cout << "Inside GeodesicSlicing" << std:: endl << std::endl;
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

            std::cout << "gAlp(" << m << "," << n << ") = "
                      << be.values_fields[ ( be.get_n_all_flds() * be.n_colpoints() ) * m
                + be.n_colpoints() * fields::gAlp + n ] << std::endl;
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
