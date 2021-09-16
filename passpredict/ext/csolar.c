#include <stdlib.h>
#include <stdio.h>
#include "passpredict.h"

/* Compute the Sun position vector

    Inputs:
       jd, double[n]: julian date array
       r, double[n, 3]: position vector to return in-place
       n, int: number of time steps, size of jd

    References:
        Vallado, p. 279, Alg. 29
        Vallado software, AST2BODY.FOR, subroutine SUN
    */
void c_sun_pos(double *jd, double *r, int n){
    
    int i;
    double t_ut1, t_tdb, lmda_Msun, M_sun, lmda_eclp;
    double r_sun_mag, eps, coslmda, sinlmda, coseps, sineps;

    for(i=0; i<n; i++){
        t_ut1 = (jd[i] - DJ00) / DJC; // DJ00=2451545.0, DJC=36525.0
        t_tdb = t_ut1;
        lmda_Msun = fmod(280.4606184 + 36000.77005361*t_tdb, 360);
        // M_sun = fmod(357.5291092 + 35999.05034*t_tdb, 360);
        M_sun = fmod(357.5277233 + 35999.05034*t_tdb, 360);
        lmda_eclp = lmda_Msun + 1.914666471 * sin(M_sun * DD2R);
        lmda_eclp += 0.019994643 * sin(2 * M_sun * DD2R);
        r_sun_mag = 1.000140612 - 0.016708617 * cos(M_sun * DD2R);
        r_sun_mag -= 0.000139589 * cos(2 * M_sun * DD2R) * (DAU/1000);
        eps = 23.439291 - 0.0130042 * t_tdb;
        coslmda = cos(lmda_eclp * DD2R);
        sinlmda = sin(lmda_eclp * DD2R);
        coseps = cos(eps * DD2R);
        sineps = sin(eps * DD2R);
        r[i*3] = r_sun_mag * coslmda;
        r[i*3 + 1] = r_sun_mag * coseps * sinlmda;
        r[i*3 + 2] = r_sun_mag * sineps * sinlmda;
    }
}


// /* Compute the Sun position vector in ECEF coordinates

//     Inputs:
//        jd, double[n]: julian date array
//        r, double[n, 3]: position vector to return in-place
//        n, int: number of time steps, size of jd

//     References:
//         Vallado, p. 279, Alg. 29
//         Vallado software, AST2BODY.FOR, subroutine SUN
//     */
// void c_sun_pos_ecef(double *jd, double *r, int n){

//     int i, j;
//     double t_ut1, t_tdb, lmda_Msun, M_sun, lmda_eclp;
//     double r_sun_mag, eps, coslmda, sinlmda, coseps, sineps;

//     double dp80, de80, dpsi, deps, epsa, rn[3][3], ee, gst;

//     // UTC1 correction set to zero
//     // const double dut1 = 0.0;

//     // EOP corrections set to zero
//     // const double xp = 0.0, yp = 0.0;
//     const double ddp80 = 0.0;
//     const double dde80 = 0.0;

//     // Set partial jd to zero
//     const double tt = 0.0;
//     const double tut = 0.0;

//     // temporary position vector
//     double p[3] = {0, 0, 0};
//     double rnp[3] = {0, 0, 0};

//     for(i=0; i<n; i++){

//         /*  from c_sun_pos  */
//         t_ut1 = (jd[i] - DJ00) / DJC; // DJ00=2451545.0, DJC=36525.0
//         t_tdb = t_ut1;
//         lmda_Msun = fmod(280.4606184 + 36000.77005361*t_tdb, 360);
//         // M_sun = fmod(357.5291092 + 35999.05034*t_tdb, 360);
//         M_sun = fmod(357.5277233 + 35999.05034*t_tdb, 360);
//         lmda_eclp = lmda_Msun + 1.914666471 * sin(M_sun * DD2R);
//         lmda_eclp += 0.019994643 * sin(2 * M_sun * DD2R);
//         r_sun_mag = 1.000140612 - 0.016708617 * cos(M_sun * DD2R);
//         r_sun_mag -= 0.000139589 * cos(2 * M_sun * DD2R) * (DAU/1000);
//         eps = 23.439291 - 0.0130042 * t_tdb;
//         coslmda = cos(lmda_eclp * DD2R);
//         sinlmda = sin(lmda_eclp * DD2R);
//         coseps = cos(eps * DD2R);
//         sineps = sin(eps * DD2R);
//         r[i*3] = r_sun_mag * coslmda;
//         r[i*3 + 1] = r_sun_mag * coseps * sinlmda;
//         r[i*3 + 2] = r_sun_mag * sineps * sinlmda;

//         /*  from c_mod2ecef  */
//         /* IAU 1980 Nutation */
//         iauNut80(jd[i], tt, &dp80, &de80); 
//         dpsi = dp80 + ddp80;
//         deps = de80 + dde80;
//         /* Mean obliquity. */
//         epsa = iauObl80(jd[i], tt);
//         /* Build the rotation matrix. */
//         iauNumat(epsa, dpsi, deps, rn);
//         /* Equation of the equinoxes, including nutation correction. */
//         ee = iauEqeq94(jd[i], tt) + ddp80 * cos(epsa);
//         /* Greenwich apparent sidereal time (IAU 1982/1994). */
//         gst = iauAnp(iauGmst82(jd[i], tut) + ee);
//         /* Form celestial-terrestrial matrix (no polar motion yet). */
//         iauRz(gst, rn);

//         /* Rotate the i position vector in-place */
//         for(j=0; j<3; j++) 
//             p[j] = r[i*3 + j];
//         iauRxp(rn, p, rnp);
//         for(j=0; j<3; j++) 
//             r[i*3 + j] = rnp[j];   
//     }

// }


void c_sun_pos_ecef(double *jd, double *r, int n){

    c_sun_pos(jd, r, n);
    c_mod2ecef(jd, r, n);

}


void c_sun_sat_illumination_distance(double *rsat, double *rsun, double *illum_dist, int n){
    int i, j;
    double zeta, d;
    double sat_i[3] = {0, 0, 0};
    double sun_i[3] = {0, 0, 0};

    for(i=0; i<n; i++)
    {
        for(j=0; j<3; j++)
        {
            sat_i[j] = rsat[i*3 + j];
            sun_i[j] = rsun[i*3 + j];
        }
        /*  compute sun-sat angle  */
        zeta = iauSepp(sat_i, sun_i);

        /*  compute normal distance from satellite vector  */
        d = iauPm(sat_i) * cos(zeta - DPI*0.5);

        /*  store difference from radius of Earth  */
        illum_dist[i] = d - D_R_EARTH;
    }


}

int main(){

    int i;
    int n = 1;
    double jd[] = {2459003.0};
    double rmod[] = {-2792.501172, -8082.480612, -3503.739214};
    double *r;
    r = (double *)malloc(3 * sizeof(double));
    c_sun_pos_ecef(jd, r, n);
    
    printf("Before mod2ecef\n");
    for(i=0; i<3; i++)
      printf("%f ", rmod[i]);
    printf("\n");
    

    c_mod2ecef(jd, rmod, 1);

    printf("After mod2ecef\n");
    for(i=0; i<3; i++)
      printf("%f ", rmod[i]);
    printf("\n");
    free(r);
    
}