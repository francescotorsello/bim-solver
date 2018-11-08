/**
 *  @file      eom-BSSN.h
 *  @brief     Equations of Motion.
 *  @author    Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

// The following files comprise EoM which are generated in Mathematica

#include "eom-BSSN/eomGauge.h"
#include "eom-BSSN/eomBSSNStdGaugeReg.h"
#include "eom-BSSN/eomBSSNMatterReg.h"
#include "eom-BSSN/eomBSSNSourcesReg.h"
#include "eom-BSSN/eomBSSNBimetricSourcesReg.h"
#include "eom-BSSN/eomBSSNRicciReg.h"

#if _EVOLVE_DSIG

    #include "eom-BSSN/eomBSSNEvolutionDR.h"

#else

    #include "eom-BSSN/eomBSSNEvolutionReg.h"

#endif

#include "eom-BSSN/eomBSSNConstraintsReg.h"
#include "eom-BSSN/eomBSSNLapseRatiosReg.h"
#include "eom-BSSN/eomBSSNMiscEquationsReg.h"
#include "eom-BSSN/eomBSSNObserver.h"
