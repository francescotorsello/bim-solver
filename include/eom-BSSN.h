/**
 *  @file      eom-BSSN.h
 *  @brief     Equations of Motion.
 *  @author    Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

// The following files comprise EoM which are generated in Mathematica

#include "eom-BSSN/eomGauge.h"
#include "eom-BSSN/eomBSSNStdGaugeReg.h"
#include "eom-BSSN/eomBSSNKDGaugeReg.h"
#include "eom-BSSN/eomBSSNMatterReg.h"
#include "eom-BSSN/eomBSSNBimetricSourcesReg.h"
#include "eom-BSSN/eomBSSNRicciReg.h"

#if _EVOLVE_DSIG

    #include "eom-BSSN/eomBSSNEvolutionDR.h"
    #include "eom-BSSN/eomBSSNSourcesReg.h"
    #include "eom-BSSN/eomBSSNConstraintsReg.h"

#else

    #if OBSERVER==1

        #include "eom-BSSN/eomBSSNEvolutionEul.h"
        #include "eom-BSSN/eomBSSNConstraintsEul.h"
        #include "eom-BSSN/eomBSSNSourcesEul.h"

    #else

        #include "eom-BSSN/eomBSSNEvolutionConv.h"
        #include "eom-BSSN/eomBSSNSourcesReg.h"
        #include "eom-BSSN/eomBSSNConstraintsReg.h"

    #endif // OBSERVER

#endif // _EVOLVE_DSIG

#include "eom-BSSN/eomBSSNLapseRatiosReg.h"
#include "eom-BSSN/eomBSSNMiscEquationsReg.h"
#include "eom-BSSN/eomBSSNObserver.h"
