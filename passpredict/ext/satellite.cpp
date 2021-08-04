#include <iostream>
#include <iomanip>
#include <math.h>

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
    /*
    Rotate rteme and vteme to respective ecef vectors
    */
    int i, j, err;
    double gmst, epsa, dp80, de80;
    // double dpsi, , deps;
    double ddp80 = 0.0, dde80 = 0.0;
    double ee, gst, w;
    double tai_1, tai_2, tt_1, tt_2;
    double rotz[3][3], cosw, sinw, tmp[3] = {0, 0, 0};

    // find GMST
    gmst = iauGmst82(jd_, 0.0);

    // Find terrestial time
    // Convert UTC to TAI
    err = iauUtctai(jd_, 0.0, &tai_1, &tai_2);
    if (err) return 1;
    // Convert TAI to TT
    err = iauTaitt(tai_1, tai_2, &tt_1, &tt_2);
    if (err) return 1;

    // Find omega from nutation 1980 theory
    iauNut80(tt_1, tt_2, &dp80, &de80);

    // Add adjustments, in this case assume zero
    // dpsi = dp80 + ddp80;
    // deps = de80 + dde80;

    // Mean obliquity
    epsa = iauObl80(tt_1, tt_2);

    // Equation of equinoxes
    ee = iauEqeq94(tt_1, tt_2) + ddp80 * std::cos(epsa);

    // find GAST with kinematic terms
    gst = gmst + ee;

    // normalize to between 0 and 2pi
    w = fmod(gst, PASSPREDICT_2PI);
    if (w < 0) w += PASSPREDICT_2PI;

    // Rotate around z axis
    // Create z-rotation matrix and transpose it
    cosw = std::cos(w);
    sinw = std::sin(w);
    rotz[0][0] = cosw;
    rotz[0][1] = sinw;
    rotz[0][2] = 0;
    rotz[1][0] = -sinw;
    rotz[1][1] = cosw;
    rotz[1][2] = 0;
    rotz[2][0] = 0;
    rotz[2][1] = 0;
    rotz[2][2] = 1;

    // rotate position vector
    // recef = rotz * rteme
    for (i=0; i<3; i++) {
        recef_[i] = 0.0;
        for (j=0; j<3; j++) {
            recef_[i] += rotz[i][j] * rteme_[j];
        }
    }
    return 0;
}

void Satellite::ComputeAltitude() {
    // Compute latitude, longitude, and altitude of satellite
    // Vallado, Alg 12, p 172
    double rdelta_sat, C;
    double phi0, phi, sin_phi0;
    double tol = 1E-6;
    u_int i;
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

} // namespace passpredict
