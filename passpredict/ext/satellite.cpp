#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "SGP4.h"
#include "orbit.h"
#include "satellite.h"
#include "passpredict.h"
#include "sofa.h"

namespace passpredict {

Satellite::Satellite(Orbit orbit) : orbit_(orbit)
{
    epoch_ = orbit_.satrec_.jdsatepoch + orbit_.satrec_.jdsatepochF;
};

void Satellite::Sgp4() {
    SGP4Funcs::sgp4(orbit_.satrec_, tsince_, rteme_, vteme_);
};

void Satellite::PropagateTSince(double tsince)
{
    /*
        inputs:
            t_tsince    time since epoch in minutes
        */
    // find time since julian date in minutes
    tsince_ = tsince;
    jd_ = epoch_ + (tsince_ / 1440.0);
    Satellite::Sgp4();
    Satellite::Teme2Ecef();
};

void Satellite::PropagateJd(double jd)
{
    /*
        inputs:
            t_jd    julian date
        */
    // find time of julian date since epoch in minutes
    jd_ = jd;
    tsince_ = (epoch_ - jd_) / 1440.0;
    Satellite::Sgp4();
    Satellite::Teme2Ecef();
};


int Satellite::Teme2Ecef()
{
    int err;
    err = ComputeTeme2Ecef(jd_, rteme_, recef_);
    return err;
}


void Satellite::ComputeAltitude() {
    // Compute latitude, longitude, and altitude of satellite
    // Vallado, Alg 12, p 172
    double rdelta_sat, C;
    double phi0, phi, sin_phi0;
    double tol = 1E-6;
    u_int i = 0;
    using namespace std;
    rdelta_sat = sqrt(recef_[0]*recef_[0] + recef_[1]*recef_[1]);
    phi = atan2(recef_[2], rdelta_sat);
    do {
        phi0 = phi;
        sin_phi0 = sin(phi0);
        C = PASSPREDICT_R_EARTH / sqrt(1 - PASSPREDICT_e2_EARTH*sin_phi0*sin_phi0);
        phi = atan2(recef_[2] + C*PASSPREDICT_e2_EARTH*sin_phi0, rdelta_sat);
        i++;
    } while (fabs(phi - phi0) >= tol && i < 100);
    alt_ = rdelta_sat / cos(phi0) - C;
}

Subpoint Satellite::ComputeSubpoint() {
    // Compute latitude, longitude, and altitude of satellite
    // Vallado, Alg 12, p 172
    double rdelta_sat, delta, tan_delta, C, lmbda, alpha;
    double phi0, phi, sin_phi0;
    double tol = 1E-6;
    int i = 0;
    Subpoint subpoint;
    using namespace std;
    rdelta_sat = sqrt(recef_[0]*recef_[0] + recef_[1]*recef_[1]);
    alpha = asin(recef_[1] / rdelta_sat);
    subpoint.lon = alpha * PASSPREDICT_RAD2DEG;
    delta = atan2(recef_[2], rdelta_sat);
    phi = delta; // longitude
    do {
        phi0 = phi;
        sin_phi0 = sin(phi0);
        C = PASSPREDICT_R_EARTH / sqrt(1 - PASSPREDICT_e2_EARTH*sin_phi0*sin_phi0);
        phi = atan2(recef_[2] + C*PASSPREDICT_e2_EARTH*sin_phi0, rdelta_sat);
        // cout << "i=" << i << "  phi0=" << phi0 * PASSPREDICT_RAD2DEG << "  C=" << C << endl;
        i++;
    } while (fabs(phi - phi0) >= tol && i < 100);
    subpoint.lat = phi0 * PASSPREDICT_RAD2DEG; // latitude
    alt_ = rdelta_sat / cos(phi0) - C;
    subpoint.alt = alt_;
    return subpoint;
};


void Satellite::Print()
{
    int i;
    using namespace std;
    cout.precision(10);
    cout << "jd    =  " << fixed << jd_ << endl;
    cout.precision(4);
    cout << "rteme = [" << fixed;
    for (i = 0; i < 3; i++)
    {
        cout << rteme_[i] << ", ";
    }
    cout << "]" << endl;
    cout << "vteme = [" << fixed;
    for (i = 0; i < 3; i++)
        cout << vteme_[i] << ", ";
    cout << "]" << endl;
};

void Satellite::PrintOneline()
{
    int i;
    using namespace std;
    cout.precision(8);
    cout << "   " << setw(14) << tsince_;
    for (i = 0; i < 3; i++)
        cout << "   " << setw(14) << rteme_[i];
    for (i = 0; i < 3; i++)
        cout << "   " << setw(12) << vteme_[i];
}


std::vector<double> PropagateSatelliteJd(double jd, Satellite satellite){
    // Function to return ECEF vector of Satellite at time jd
    double tsince;
    std::vector<double> recef(3, 0.0);
    double rteme[3], vteme[3], recef_ary[3];
    int i, err;

    // Propagate satellite using SGP4
    tsince = (satellite.epoch_ - jd) / 1440.0;
    SGP4Funcs::sgp4(satellite.orbit_.satrec_, tsince, rteme, vteme);

    // Rotate TEME vector to ECEF
    err = ComputeTeme2Ecef(jd, rteme, recef_ary);

    // populate vector object and return
    for (i = 0; i < 3; i++)
        recef[i] = recef_ary[i];

    return recef;
};

int Utc2tt(double jd1, double jd2, double &tt1, double &tt2){
    // Julian date UTC to IAU terrestial time
    int j = 0;
    double tai1, tai2;
    j = iauUtctai(jd1, jd2, &tai1, &tai2);
    if ( j ) return 1;
    j = iauTaitt(tai1, tai2, &tt1, &tt2);
    if ( j ) return 1;
    return 0;
};


int ComputeTeme2Ecef(double jd, double rteme[3], double recef[3]){
    /*
    Rotate rteme and vteme to respective ecef vectors
    */
    int i, j, err;
    double gmst;
    double rotz[3][3];
    double costheta, sintheta;

    // find GMST
    gmst = iauGmst82(jd, 0.0);

    // Rotate around z axis
    // Create z-rotation matrix and transpose it
    costheta = std::cos(gmst);
    sintheta = std::sin(gmst);
    rotz[0][0] = costheta;
    rotz[0][1] = sintheta;
    rotz[0][2] = 0;
    rotz[1][0] = -sintheta;
    rotz[1][1] = costheta;
    rotz[1][2] = 0;
    rotz[2][0] = 0;
    rotz[2][1] = 0;
    rotz[2][2] = 1;

    // rotate position vector
    // recef = rotz * rteme
    for (i=0; i<3; i++) {
        recef[i] = 0.0;
        for (j=0; j<3; j++) {
            recef[i] += rotz[i][j] * rteme[j];
        }
    }
    return 0;
}


} // namespace passpredict
