#ifndef SATELLITE_H
#define SATELLITE_H

#include "orbit.h"

namespace passpredict
{
    class Satellite
    {
    public:
        double rteme[3] = {0, 0, 0}; // TEME position vector at time jd
        double vteme[3] = {0, 0, 0}; // TEME velocity vector at time jd
        double jd;                   // julian date
        double tsince;               // minutes since epoch
        double epoch;                // epoch, julian date

        Orbit orbit;      // orbit data, OMM/TLE
        std::string name; // satellite name
        int satid;        // NORAD satellite ID number
        double recef[3];  // ECEF position vector at time jd
        double vecef[3];  // ECEF velocity vector at time jd

        Satellite(){};
        Satellite(Orbit);
        int teme2ecef();
        void sgp4();
        void propagate_tsince(double t_tsince);
        void propagate_jd(double t_jd);
        void print();
        void print_oneline();
    };

}

#endif