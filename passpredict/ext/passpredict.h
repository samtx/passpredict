#ifndef PASSPREDICT_H
#define PASSPREDICT_H

extern "C"
{
#include "sofa.h"
}
#include "SGP4.h"

#define PASSPREDICT_DEG2RAD 1.745329251994329576923691e-2 // sofam.m DD2R
#define PASSPREDICT_RAD2DEG 57.29577951308232087679815 // sofam.h DR2D
#define PASSPREDICT_PI 3.141592653589793238462643 // sofam.h DPI
#define PASSPREDICT_2PI 6.283185307179586476925287        // sofam.m D2PI
#define PASSPREDICT_R_EARTH 6378.137                      //  Earth mean equatorial radius [km]
#define PASSPREDICT_e2_EARTH 0.006694385000               // Earth eccentricity squared

#include "omm.h"
#include "orbit.h"
#include "satellite.h"
#include "location.h"
#include "observer.h"

namespace passpredict {

double ComputeElevationAngle(double, Location, Satellite);
double FindAOS(double, double, Location, Satellite);

} // passpredict


#endif