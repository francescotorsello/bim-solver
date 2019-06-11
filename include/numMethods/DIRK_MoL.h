/**
 *  @file      DIRK_MoL.h
 *  @brief     Descriptors for diagonally implicit Runge_Kutta methods.
 *  @authors   Francesco Torsello, Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

 /**    References
  *     [1] Kværnø, Singly Diagonally Implicit Runge\[Dash]Kutta Methods with an Explicit
  *         First Stage, A.BIT Numerical Mathematics (2004) 44 : 489.
  *         https : // doi.org/10.1023/B : BITN .0000046811 .70614 .38
  */

#ifndef _DIRK_MOL_H_INCLUDED
#define _DIRK_MOL_H_INCLUDED

struct ButcherTable
{
    const char* name;    //!< Name of the DIRK method
    int s;               //!< Number of stages (s <= 16 ?)
    int pOrder;          //!< The highest adaptive order (will be used by adaptive method)
    Real A[16][16];      //!< Matrix A in the Butcher table
    Real c[16];          //!< Vector c in the Butcher table
    Real b2[16];         //!< Alternative beta coefficients for the adaptive step
};

// Singly diagonally implicit Runge-Kutta methods with an explicit first stage (ESDIRK)

ButcherTable ESDIRK32 =  //!< ESDIRK 3/2 with 4 stages [1, p.497]
{
    "ESDIRK32", 4, 0,

    {
        {0, 0, 0, 0},
        {0.435866521508458999416019451194, 0.435866521508458999416019451194, 0, 0},
        {0.49056338842178057062846795847, 0.073570090069760429955512590340,
            0.435866521508458999416019451194, 0},
        {0.308809969976746523348162469887, 1.4905633884217805706284679585,
            -1.2352398799069860933926498795, 0.435866521508458999416019451194}},

    { 0, 0.871733043016917998832038902387, 1, 1 }
};

#endif // _DIRK_MOL_H_INCLUDED
