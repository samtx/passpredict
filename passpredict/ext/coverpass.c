#include "passpredict.h"

/*  compute range and elevation from topocentric SEZ position vector  */
void c_sez2rngel(double *rSEZ, double *rng_ary, double *el_ary, int n){
    int i, j;
    double rng, el;
    double p[3] = {0, 0, 0};

    for(i=0; i<n; i++)
    {
        for(j=0; j<3; j++)
            p[j] = rSEZ[i*3 + j];
        
        /*  compute range  */
        rng = iauPm(p);

        /*  compute elevation  */
        el = asin(p[2]/rng) * DR2D;

        rng_ary[i] = rng;
        el_ary[i] = el;
    }

}