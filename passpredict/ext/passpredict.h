#ifndef PASSPREDICT_H
#define PASSPREDICT_H

#include "sofa.h"

#define PASSPREDICT_R_EARTH 6378.137
#define PASSPREDICT_DEG2RAD 1.745329251994329576923691e-2  // sofam.m DD2R
#define PASSPREDICT_2PI 6.283185307179586476925287      // sofam.m D2PI


namespace passpredict {

}

#define D_R_EARTH 6378.137  //  Earth mean equatorial radius [km]

/*  crotations.c  */
void c_mod2ecef(double *jd, double *r, int n);
void c_teme2ecef(double *jd, double *r, int n);
void c_ecef2sez(double *r, double phi, double lmda, int n);

/*  csolar.c  */
void c_sun_pos_ecef(double *jd, double *r, int n);
void c_sun_pos(double *jd, double *r, int n);
void c_sun_sat_illumination_distance(double *rsat, double *rsun, double *illum_dist, int n);

/*  coverpass.c  */
void c_sez2rngel(double *r_view, double *rng_view, double *el_view, int n);

#endif