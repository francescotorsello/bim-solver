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

/*class GaugeVariables
{

    BispecEvolve *be; //

public:

    inline void geodesicSlicing_GR( Int m )
    {
        for( Int n = 0; n < (*be).n_colpoints(); ++n )
        {
            (*be).values_fields[ ( n_flds * n_colpoints() ) * m \
            + n_collocs * fields::gAlp + n ];
            //be.gBet( m, n ) = 0;

            //be.fAlp( m, n ) = 1;
            //be.fBet( m, n ) = 0;
        }
    }

    //Real Kdiff = 0.01;
    //Real Kelas = 0.03;

    //#include "../include/eom-BSSN/eomBSSNKDGaugeComp.h"

    GaugeVariables(
        BispecEvolve BE
    ) :
        // initialize the pointer to the object of type BispecEvolve, to the object passed
        // as argument to the constructor.
        be( &BE )
    {}

};*/
