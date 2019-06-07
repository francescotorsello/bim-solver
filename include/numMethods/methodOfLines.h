/**
 *  @file      methodOfLines.h
 *  @brief     Interface and descriptors for various methods of lines (MoL).
 *  @authors   Mikica Kocic, Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _METHOD_OF_LINES_H_INCLUDED
#define _METHOD_OF_LINES_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g11 Numerical Methods                                                    */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
/** Interface used for MoL integration that prepares and finalizes integration steps.
  */
class IntegFace
{
    friend class MoL;

    /** Calculate the dependent variables from the prime state variables
     *  which are needed for the integration.
     */
    virtual void integStep_Prepare( Int m ) = 0;

    /** Perform the additional steps after each integration step.
     */
    virtual void integStep_Finalize( Int mNext, Int mPrev ) = 0;

    /** Calculate diagnostics variables (also checking for eventual NaNs).
     */
    virtual bool integStep_Diagnostics( Int m, Int chkNaNs_nFrom, Int chkNaNs_nTo ) = 0;

    /** Prepare the right-hand side of the evolution equations.
     */
    virtual void integStep_CalcEvolutionRHS( Int m ) = 0;

    /** Fix the state variables at the left boundary.
     */
    virtual void applyLeftBoundaryCondition( Int m ) = 0;

    /** Fix the state variables at the right boundary.
     */
    virtual void applyRightBoundaryCondition( Int m ) = 0;
};

/** @struct MoLDescriptor
 *  A method of lines (MoL) descriptor.
 *
 *  <pre>
 *       partial_t Y = L[Y],  </pre>
 *
 *  using the following algorithm:
 *  <pre>
 *      Y^{(0)}   = Y^{m},
 *      Y^{(i+1)} = Sum_{j=0}^{i} ( alpha_{ij} Y^{(j)} ) +            <- integStep_MoL
 *                  Delta t beta_i L[Y^{(i)}],   for i = 0, ..., N-1, <- integStep_MoLInit
 *      Y^{m+1}   = Y^{(N)}.    </pre>
 *
 *  The variables `Y^{(i)}`, `i = 0, ..., N`, are intermediate.
 *
 *  The method is completely specified by `N`, and arrays alpha and beta
 *  in the structure MoLDescriptor.
 *
 *  The alpha and beta coefficients are not Butcher's tableaus.
 *  See the accompanying Mathematica notebook `MoL descriptor.nb` for the conversion.
 *
 *  The beta2 coefficients are used for the adaptive step method. When beta2 is used,
 *  the row alpha[N] is used instead of the row alpha[N-1], that is,
 *  the matrix consists of alpha[0], ..., alpha[N-2], alpha[N].
 *
 *  @see <a href="http://cactuscode.org/documentation/thorns/CactusBase-MoL.pdf">
 *       Method of Lines - Cactus Code</a>
 *  @see Butcher, Numerical Methods for Ordinary Differential Equations, 2008,
 *       @cite Butcher:2008num
 */
struct MoLDescriptor
{
    const char* name;    //!< Name of the MoL method
    int N;               //!< Number of intermediate steps (N <= 16)
    int pOrder;          //!< The highest adaptive order (will be used by adaptive method)
    Real alpha[16][16];  //!< Alpha coefficients
    Real beta[16];       //!< Beta coefficients
    Real beta2[16];      //!< Alternative beta coefficients for the adaptive step
};

// Time evolution methods provided by MoL

MoLDescriptor ICN2 = //!< Iterative Crank-Nicolson, 2-iterations
{
    "Iterative Crank-Nicolson", 2, 0,
    { {  1., 0. },
      {  0., 1. } },
    { 1./2., 1. }
};

MoLDescriptor ICN3 = //!< Iterative Crank-Nicolson, 3-iterations
{
    "Iterative Crank-Nicolson", 3, 0,
    { {  1., 0., 0. },
      {  0., 1., 0. },
      {  0., 0., 1. } },
    { 1./2., 1./2., 1. }
};

MoLDescriptor RK1 = //!< Classical Runge-Kutta, 1st order (the Euler method)
{
    "Runge-Kutta", 1, 0,
    { { 1. } },
    { 1. }
};

MoLDescriptor RK2 = //!< Classical Runge-Kutta, 2nd order
{
    "Runge-Kutta", 2, 0,
    { {     1.,    0. },
      {  1./2., 1./2. } },
    { 1., 1./2. }
};

MoLDescriptor RK3 = //!< Classical Runge-Kutta, 3rd order
{
    "Runge-Kutta", 3, 0,
    { {     1.,    0.,    0. },
      {  3./4., 1./4.,    0. },
      {  1./3.,    0., 2./3. } },
    { 1., 1./4., 2./3. }
};

MoLDescriptor RK4 = //!< Classical Runge-Kutta, 4th order
{
    "Runge-Kutta", 4, 0,
    { {     1.,    0.,    0.,    0. },
      {     1.,    0.,    0.,    0. },
      {     1.,    0.,    0.,    0. },
      { -1./3., 1./3., 2./3., 1./3. } },
    { 1./2., 1./2., 1., 1./6. }
};

// The adaptive step methods are in a separate file:
//
#include "embeddedMoL.h"

// The diagonally implicit (DIRK) step methods are in a separate file:
//
#include "DIRK_MoL.h"
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _METHOD_OF_LINES_H_INCLUDED
