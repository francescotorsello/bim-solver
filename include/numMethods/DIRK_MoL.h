/**
 *  @file      DIRK_MoL.h
 *  @brief     Descriptors for diagonally implicit Runge_Kutta methods.
 *  @authors   Francesco Torsello, Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _DIRK_MOL_H_INCLUDED
#define _DIRK_MOL_H_INCLUDED

/** DIRK, third order, L-stable method
 */

MoLDescriptor DIRK3L =
{
    "DIRK3L", 3, 0,
    { { 1.00000000000000000000000000000, 0., 0.},
      { 0.352859819838216309645204976060, \
0.64714018016178369035479502394, 0.},
      { -1.25097989510569161029602026045, \
3.72932966245123233779360198613, -1.47834976734554072749758172569 }
    },
    { 0.435866521500000031076953809640, \
0.435866521500000031076953809640, 0.435866521500000031076953809640}
};

#endif // _DIRK_MOL_H_INCLUDED
