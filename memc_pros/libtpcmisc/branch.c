/// @file branch.c
/// @author Vesa Oikonen
/// @brief Isotope branching ratio correction.
///
/******************************************************************************/
#include "libtpcmisc.h"
/******************************************************************************/

/******************************************************************************/
/** Branching fraction for specified isotope.
 *
\return Returns the branching factor, or 0 in case branching fraction is unknown.
 */
float branchingFraction(
  /** Isotope code; see hlIsotopeFromHalflife() */
  int isotope
) {
  float bf=0.0;
  switch(isotope) {
    case TPCISOT_CU_64: bf=BRANCHING_Cu64; break;
    case TPCISOT_GA_68: bf=BRANCHING_Ga; break;
    case TPCISOT_GE_68: bf=BRANCHING_Ge; break;
    case TPCISOT_RB_82: bf=BRANCHING_Rb; break;
    case TPCISOT_F_18:  bf=BRANCHING_F; break;
    case TPCISOT_C_11:  bf=BRANCHING_C; break;
    case TPCISOT_N_13:  bf=BRANCHING_N; break;
    case TPCISOT_O_15:  bf=BRANCHING_O; break;
    case TPCISOT_BR_75:
    case TPCISOT_BR_76:
    case TPCISOT_CU_62:
    case TPCISOT_FE_52:
    case TPCISOT_NA_22:
    case TPCISOT_O_14:
    case TPCISOT_I_124:
    case TPCISOT_ZN_62:
    case TPCISOT_UNKNOWN:
    default:            bf=0.0;
  }
  return(bf);
}
/******************************************************************************/

/******************************************************************************/
