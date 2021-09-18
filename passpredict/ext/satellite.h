#ifndef SATELLITE_H
#define SATELLITE_H

#include <vector>

#include "orbit.h"

namespace passpredict {

typedef struct Subpoint {
    double lat;
    double lon;
    double alt;
} Subpoint;

class Satellite {
    public:
        double rteme_[3] = {0, 0, 0}; // TEME position vector at time jd
        double vteme_[3] = {0, 0, 0}; // TEME velocity vector at time jd
        double jd_;                   // julian date
        double tsince_ = 0;               // minutes since epoch
        double epoch_ = 0;                // epoch, julian date
        double alt_ = 0;                  // altitude [km]

        Orbit orbit_;      // orbit data, OMM/TLE
        std::string name_; // satellite name
        int satid_;        // NORAD satellite ID number
        double recef_[3];  // ECEF position vector at time jd
        double vecef_[3];  // ECEF velocity vector at time jd

        Satellite(){};
        Satellite(Orbit);
        int Teme2Ecef();
        void Sgp4();
        void PropagateTSince(double t_tsince);
        void PropagateJd(double t_jd);
        Subpoint ComputeSubpoint();
        void ComputeAltitude();

        void Print();
        void PrintOneline();
};

std::vector<double> PropagateSatelliteJd(double, Satellite);
int ComputeTeme2Ecef(double, double[3], double[3]);
int Utc2tt(double jd1, double jd2, double &tt1, double &tt2);

} // namespace passpredict

#endif