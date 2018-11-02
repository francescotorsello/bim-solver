/**
 *  @file      gridFunctions.h
 *  @brief     Contains definitions of the variables that are tracked on a grid point.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 *
 *  Here we configure the grid functions we are dealing with on a grid point.
 *  The customization should be done in the @ref fld namespace below.
 *
 *  @copydoc newgf
 *
 *  @page newgf Adding new grid functions
 *  @par
 *  A **grid function** is a discretization of a variable (a "field") which is defined
 *  on every point of a given grid. Examples are the components of the metric tensor
 *  in GR, or the rest-mass density of a fluid which is defined in all cells
 *  within a given domain.
 *
 *  ## Instructions to add a new grid function for the field *gf* ##
 *
 *  - add the identifier **gf** to an enumeration in the @ref fld namespace
 *
 *  - (optionally) add **gf** to a list passed to GridInitialData::gridFunction()
 *
 *  - (optionally) add **gf** to a list passed to GridOutputWriter::gridFunction()
 *
 *  - (optionally) add **gf** to a list passed to MoLIntegrator::keepConstant()
 *    or to MoLIntegrator::keepEvolved()
 *
 *  - to be able to access **gf**'s values using **gf(m,n)**, add @ref emitField(gf)
 *    inside the declaration.
 *
 *  You will also need to:
 *
 *  - (optionally) add @ref emitDerivative_r(gf) and/or @ref emitDerivative_rr(gf)
 *    (in case you need the spatial derivatives of the field).
 *
 *  - (optionally) define the RHS of the evolution equation `Real eq_gf_t( Int m, Int n )`
 *    and then assign the RHS to **gf_t** in BimetricEvolve::integStep_Prepare.
 *
 *  - setup the behavior of 'gf' at boundaries in the @ref IntegFace interface methods
 *    IntegFace::applyLeftBoundaryCondition and IntegFace::applyRightBoundaryCondition.
 */

#ifndef _GRID_FUNCTIONS_H_INCLUDED
#define _GRID_FUNCTIONS_H_INCLUDED

#include <vector>

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g4 Grid driver                                                           */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
/** Verbose descriptor of a grid function.
 */
struct GF_Descriptor
{
    Int gf;            //!< Grid function index
    std::string name;  //!< Name (in Mathematica, char '_' is translated to '$')
    std::string tex;   //!< TeX equation
};

/** @namespace fld Contains the localized variable names for all known grid functions.
 */
namespace fld
{
    /** Identifiers of all known grid functions (the variables or fields on a grid).
     */
    enum sysIndex
    {
        t = 0, r,     //!< Coordinates
        sysStat,      //!< Grid function used to store various statistics
        sysLast       //!< Used as the size marker
    };

    #define GFCNT fld::sysLast

    /** A pair of two grid functions where the first is evolved in time by the other.
     */
    struct EvolvedBy { Int f; Int f_t; };
}
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _GRID_FUNCTIONS_H_INCLUDED
