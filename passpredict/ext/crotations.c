#include "passpredict.h"

/*
Rotate MOD position vector r[n, 3] using 1980 Nutation theory and
1976/FK5 sidereal time 

Rotates in-place

Rotates MOD -> ITRF (ECEF)
*/
void c_mod2ecef(double *jd, double *r, int n){
    int i, j;
    int err;

    double dp80, de80, dpsi, deps, epsa, rn[3][3], ee, gst;
    int year, month, day;
    double day_fraction, delta_at;


    // UTC1 correction set to zero
    // const double dut1 = 0.0;

    // EOP corrections set to zero
    // const double xp = 0.0, yp = 0.0;
    const double ddp80 = 0.0;
    const double dde80 = 0.0;

    double tai, tt, tt2000, tt_date, tt_time, tt1, tt2;
    double mjd2000, mjd;
    double jd1, jd2;
    
    // temporary position vector
    double p[3] = {0, 0, 0};
    double rnp[3] = {0, 0, 0};

    for(i=0; i<n; i++){
        /*  Break julian date into date and time components */
        jd1 = fmod(jd[i], 1.0);
        jd2 = jd[i] - jd1;
        
        /* Find terrestial time, ignore delta_UT1*/
        err = iauJd2cal(jd1, jd2, &year, &month, &day, &day_fraction);
        // Have some error checking here, if err != 0
        err = iauDat(year, month, day, day_fraction, &delta_at);
        tai = (jd1 + jd2) + delta_at/DAYSEC;
        tt = tai + 32.184/DAYSEC;
        tt -= DJ00;   // Leave terrestial time in J2000 format

        /* IAU 1980 Nutation */
        iauNut80(DJ00, tt, &dp80, &de80); 
        dpsi = dp80 + ddp80;
        deps = de80 + dde80;
        /* Mean obliquity. */
        epsa = iauObl80(DJ00, tt);
        /* Build the rotation matrix. */
        iauNumat(epsa, dpsi, deps, rn);
        /* Equation of the equinoxes, including nutation correction. */
        ee = iauEqeq94(DJ00, tt) + ddp80 * cos(epsa);
        /* Greenwich apparent sidereal time (IAU 1982/1994). */
        /* Use date and time method for highest accuracy */
        tt1 = fmod(tt, 1);
        tt2 = tt - tt1;
        gst = iauAnp(iauGmst82(DJ00 + tt1, tt2) + ee);
        /* Form celestial-terrestrial matrix (no polar motion yet). */
        iauRz(gst, rn);

        /* Rotate the i position vector in-place */
        for(j=0; j<3; j++) 
            p[j] = r[i*3 + j];
        iauRxp(rn, p, rnp);
        for(j=0; j<3; j++) 
            r[i*3 + j] = rnp[j];   
    }
}


/*
Rotate TEME position vector r[n, 3] using FK5/1976/1980 theory and 

Rotates in-place

Rotates TEME -> ECEF (ITRF)
*/
void c_teme2ecef(double *jd, double *r, int n){
    int i, j;

    double dp80, de80, dpsi, deps, epsa, rn[3][3], ee, gst;

    int year, month, day, err;
    double day_fraction, delta_at;

    // UTC1 correction set to zero
    // const double dut1 = 0.0;

    // EOP corrections set to zero
    // const double xp = 0.0, yp = 0.0;
    const double ddp80 = 0.0, dde80 = 0.0;

    // Set partial jd to zero
    double tai, tt, tt2000, tt_date, tt_time, tt1, tt2;
    double mjd2000, mjd;
    double jd1, jd2;
    

    // temporary position vector
    double p[3] = {0, 0, 0};
    double rnp[3] = {0, 0, 0};

    for(i=0; i<n; i++){

        /*  Break julian date into date and time components */
        jd1 = fmod(jd[i], 1.0);
        jd2 = jd[i] - jd1;
        
        /* Find terrestial time, ignore delta_UT1*/
        err = iauJd2cal(jd1, jd2, &year, &month, &day, &day_fraction);
        // Have some error checking here, if err != 0
        err = iauDat(year, month, day, day_fraction, &delta_at);
        tai = (jd1 + jd2) + delta_at/DAYSEC;
        tt = tai + 32.184/DAYSEC;
        // tt -= DJ00;   // Leave terrestial time in J2000 format

        /* IAU 1980 Nutation */
        iauNut80(DJ00, tt, &dp80, &de80); 
        dpsi = dp80 + ddp80;
        /* Mean obliquity. */
        epsa = iauObl80(DJ00, tt);
        /* Equation of the equinoxes, including nutation correction. */
        ee = iauEqeq94(DJ00, tt) + ddp80 * cos(epsa);
        /* Greenwich apparent sidereal time (IAU 1982/1994). */
        /* Use date and time method for highest accuracy */
        tt1 = fmod(tt, 1);
        tt2 = tt - tt1;
        gst = iauAnp(iauGmst82(DJ00 + tt1, tt2) + ee);
        // gst = iauAnp(iauGmst82(tt1, tt2));
        /* Form celestial-terrestrial matrix (no polar motion yet). */
        iauIr(rn); // Initialize rn to identity matrix
        iauRz(gst, rn);

        /* Rotate the i position vector in-place */
        for(j=0; j<3; j++) 
            p[j] = r[i*3 + j];
        iauRxp(rn, p, rnp);
        // iauTrxp(rn, p, rnp);
        for(j=0; j<3; j++) 
            r[i*3 + j] = rnp[j];   
    }
}


/*
Rotate a position vector r[n, 3] from ECEF to SEZ 

Rotates in-place
*/
void c_ecef2sez(double *r, double phi, double lmda, int n){
    int i, j;
    double phi_rad, lmda_rad;
    double ang1, cosang1, sinang1, cosang2, sinang2;

    // temporary position vector
    double p[3] = {0, 0, 0};
    double mp[3] = {0, 0, 0};
    double m[3][3] = 
        {
            {0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}
        };

    phi_rad = phi * DD2R;
    lmda_rad = lmda * DD2R;
    ang1 = (90 - phi) * DD2R;
    cosang1 = cos(ang1);
    sinang1 = sin(ang1);
    cosang2 = cos(lmda_rad);
    sinang2 = sin(lmda_rad);
    m[0][0] = cosang1 * cosang2;
    m[0][1] = cosang1 * sinang2;
    m[0][2] = -sinang1;
    m[1][0] = -sinang2;
    m[1][1] = cosang2;
    m[2][0] = sinang1 * cosang2;
    m[2][1] = sinang1 * sinang2;
    m[2][2] = cosang1;
    for(i=0; i<n; i++){
        /* Rotate the i position vector in-place */
        for(j=0; j<3; j++) 
            p[j] = r[i*3 + j];
        iauRxp(m, p, mp);
        for(j=0; j<3; j++) 
            r[i*3 + j] = mp[j];   
    }
}