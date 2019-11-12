/// @file point.c
/// @brief A three dimensional point.
/// @author Riku Klén
///
/***********************************************************************/

/***********************************************************************/
#include "libtpcimgp.h"
/***********************************************************************/

/***********************************************************************/
/** Round float to int.
    @return Rounded value.
 */
int pRound(
  /** Any number */
  float number
) {
  if(number - floor(number) < 0.5){
    return floor(number);
  }
  return ceil(number);
}
/***********************************************************************/

/***********************************************************************/
/** Calculate distance between points
    @return Distance.
 */
float getDistance(
  /** First point */
  point begin, 
  /** Second point */
  point end
) {
  float dx, dy, dz;

  dx = begin.x - end.x;	// Difference in x cooridnates.
  dy = begin.y - end.y;	// Difference in y cooridnates.
  dz = begin.z - end.z;	// Difference in z cooridnates.
  return sqrt(dx*dx + dy*dy + dz*dz);
}
/***********************************************************************/

/***********************************************************************/
/** Calculates xy-projection of angle FCX in degrees, where F=first point, 
    C=centre point and X=point with higher x coordinate (y and z coordinate remain the same).

    This is used to calculate polar angle with two dimensional points (z=constant).
    @return Angle first - centre - x in degrees (0-360) and -360.0 if first point is 
     equal to centre point.
 */
float getAngle(
  /** First point. */
  point begin, 
  /** Centre point. */
  point center
) {
  float dx, dy;
  dx = begin.x - center.x;	// Difference in x coordinates.
  dy = begin.y - center.y;	// Difference in y coordinates.

  // Check trivial cases.
  if(dx == 0.0){
    if(dy == 0.0) return -360.0;
    if(dy > 0.0) return 90.0;
    else return 270.0;
  }
  if(dy==0.0){
    if(dx > 0.0) return 0.0;
    else return 180.0;
  }

  // If it was not a trivial case, then calculate angle in...
  // ...the second and the third quarter.
  if(dx<0.0) return 180.0 + atan(dy/dx)*180.0/M_PI;

  // ...the first quarter.
  if(dy>0.0) return atan(dy/dx)*180.0/M_PI;

  // ...the fourth quarter.
  return 360.0 + atan(dy/dx)*180.0/M_PI;
}
/***********************************************************************/

/***********************************************************************/
