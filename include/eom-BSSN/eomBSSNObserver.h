/** @file  eomBSSNObserver.h
 *  @author Francesco Torsello
 *  @brief The BSSN observer and the trace of the conformal extrinsic curvature.
 *  @version 2018-09-03T11:11:01
 *  @image html BSSNobserver.png
 */

 /////////////////////////////////////////////////////////////////////////////////////
    /** @defgroup Determinants
     *  @note These definitions are not used in the Eulerian case
    /////////////////////////////////////////////////////////////////////////////////////*/

                                                                                    /*@{*/
    Real BimetricEvolve::gdet(Int m, Int n){
        return pow2(gA(m, n))*pow(gB(m, n),4);
    }
    Real BimetricEvolve::fdet(Int m, Int n){
        return pow2(fA(m, n))*pow(fB(m, n),4);
    }
                                                                                   /*@}*/

/////////////////////////////////////////////////////////////////////////////////////
    /** @defgroup The choice of the observer
    /////////////////////////////////////////////////////////////////////////////////////*/
    // The choice of the traces of the modified extrinsic curvatures (MECs) and of the determinant of the metrics (Eulerian or Lagrangian observers)
    // Here: The traces are 0, so the MECs are the traceless parts of the extrinsic curvatures, and the observer is Eulerian, i.e the Pfaffian of the determinant is zero


#if OBSERVER==1

    Real BimetricEvolve::gdet_pff(Int m, Int n ) {
        return  0;
    }

    Real BimetricEvolve::gdet_pff_r( Int m, Int n ) {
        return  0;
    }

    Real BimetricEvolve::gdet_r( Int m, Int n ) {
        return  0;
    }// For the Eulerian observer, this term is not important since it is multiplied by detg_pf (m,n)

    Real BimetricEvolve::fdet_pff( Int m, Int n ) {
        return  0;
    }

    Real BimetricEvolve::fdet_pff_r( Int m, Int n ) {
        return  0;
    }

    Real BimetricEvolve::fdet_r( Int m, Int n ) {
        return  0;
    } // For the Eulerian observer, this term is not important since it is multiplied by detg_pf (m,n)

#elif OBSERVER==2

    Real BimetricEvolve::gdet_pff( Int m, Int n ) {
        return  0;
    }

    Real BimetricEvolve::fdet_pff( Int m, Int n ) {
        return  0;
    }

#endif // OBSERVER

/////////////////////////////////////////////////////////////////////////////////////
    // Definition of the traces of the shifted extrinsic curvatures

#define gAvalue 0

#if gAvalue==0

    Real BimetricEvolve::gtrA( Int m, Int n ) {
        return  1e-100;
    }

    Real BimetricEvolve::gtrA_pff( Int m, Int n ) {
        return  0;
    }

    Real BimetricEvolve::gtrA_r( Int m, Int n ) {
        return  0;
    }

    Real BimetricEvolve::gtrAv( Int m, Int n ) {
        return  gA1 (m,n) + 2*gA2 (m,n);
    }

#else

    Real BimetricEvolve::gtrAv( Int m, Int n ) {
        return  0;
    }

#endif

#define fAvalue 0

#if fAvalue==0

    Real BimetricEvolve::ftrA( Int m, Int n ) {
        return  1e-100;
    }

    Real BimetricEvolve::ftrA_pff( Int m, Int n ) {
        return  0;
    }

    Real BimetricEvolve::ftrA_r( Int m, Int n ) {
        return  0;
    }

    Real BimetricEvolve::ftrAv( Int m, Int n ) {
        return  fA1 (m,n) + 2*fA2 (m,n);
    }

#else

    Real BimetricEvolve::ftrAv( Int m, Int n ) {
        return  0;
    }

#endif
