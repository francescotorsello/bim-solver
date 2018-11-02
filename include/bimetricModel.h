/**
 *  @file      bimetricModel.h
 *  @brief     Implements the parameters which define bimetric models.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _BIMETRIC_MODEL_H_INCLUDED
#define _BIMETRIC_MODEL_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g2 The bimetric model                                                    */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
#ifndef V_SIGN
    #define V_SIGN -   //!< V_SIGN should be defined either as empty or as `-`.
#endif

/** A structure that holds parameters which define a specific bimetric model.
 */
struct BimetricModel
{
    std::string title; //!< Description

    Real k_g;  //!< kappa_g parameter
    Real k_f;  //!< kappa_f parameter
    Real b_0;  //!< beta_0 parameter
    Real b_1;  //!< beta_1 parameter
    Real b_2;  //!< beta_2 parameter
    Real b_3;  //!< beta_3 parameter
    Real b_4;  //!< beta_4 parameter

    /** Constructs the bimetric model based on values from the parameter file.
     */
    BimetricModel( Parameters& params )
    {
        params.get( "model.title", title, "" );
        params.get( "model.k_g",   k_g,   1.0 );
        params.get( "model.k_f",   k_f,   1.0 );
        params.get( "model.b_0",   b_0,   0.0 );
        params.get( "model.b_1",   b_1,   0.0 );
        params.get( "model.b_2",   b_2,   0.0 );
        params.get( "model.b_3",   b_3,   0.0 );
        params.get( "model.b_4",   b_4,   0.0 );

        slog << "Bimetric Model: " << std::endl << std::endl
            << "    title = " << title << std::endl
            << "    k_g = " << k_g << ",  k_f = " << k_f
            << ",  b_0 = " << b_0 << ",  b_1 = " << b_1
            << ",  b_2 = " << b_2 << ",  b_3 = " << b_3
            << ",  b_4 = " << b_4 << std::endl << std::endl;
    }

    /**  Returns true if the bimetric model comprises two decoupled GR sectors.
     */
    bool isGR ()
    {
        return b_1 == 0 && b_2 == 0 && b_3 == 0;
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /** @defgroup g3 Shifted elementary symmetric polynomials                          */
    /////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
    inline Real P_0_0( Real R ) { return V_SIGN b_0; }
    inline Real P_0_1( Real R ) { return V_SIGN b_1; }
    inline Real P_0_2( Real R ) { return V_SIGN b_2; }
    inline Real P_0_3( Real R ) { return V_SIGN b_3; }
    inline Real P_0_4( Real R ) { return V_SIGN b_4; }

    inline Real P_1_0( Real R ) { return V_SIGN ( b_0 + b_1 * R ); }
    inline Real P_1_1( Real R ) { return V_SIGN ( b_1 + b_2 * R ); }
    inline Real P_1_2( Real R ) { return V_SIGN ( b_2 + b_3 * R ); }
    inline Real P_1_3( Real R ) { return V_SIGN ( b_3 + b_4 * R ); }

    inline Real P_2_0( Real R ) { return V_SIGN ( b_0 + ( 2 * b_1 + b_2 * R ) * R ); }
    inline Real P_2_1( Real R ) { return V_SIGN ( b_1 + ( 2 * b_2 + b_3 * R ) * R ); }
    inline Real P_2_2( Real R ) { return V_SIGN ( b_2 + ( 2 * b_3 + b_4 * R ) * R ); }
};
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _BIMETRIC_MODEL_H_INCLUDED
