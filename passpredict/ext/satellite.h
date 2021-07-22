#ifndef SATELLITE_H
#define SATELLITE_H

#include "orbit.h"

namespace passpredict
{
    class Satellite
    {
    private:
        double m_rteme[3] = {0, 0, 0}; // TEME position vector at time jd
        double m_vteme[3] = {0, 0, 0}; // TEME velocity vector at time jd
        double m_jd;                   // julian date
        double m_tsince;               // minutes since epoch
        double m_epoch;                // epoch, julian date
    public:
        Orbit orbit;      // orbit data, OMM/TLE
        std::string name; // satellite name
        int satid;        // NORAD satellite ID number
        double recef[3];  // ECEF position vector at time jd
        double vecef[3];  // ECEF velocity vector at time jd

        Satellite(){};
        Satellite(Orbit);
        void propagate_tsince(double t_tsince);
        void propagate_jd(double t_jd);
        void teme2ecef();
        void print();
        void print_oneline();
    };
}

#endif