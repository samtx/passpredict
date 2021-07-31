#ifndef SATELLITE_H
#define SATELLITE_H

#include "orbit.h"

namespace passpredict {
    
class Satellite {
    public:
        double rteme_[3] = {0, 0, 0}; // TEME position vector at time jd
        double vteme_[3] = {0, 0, 0}; // TEME velocity vector at time jd
        double jd_;                   // julian date
        double tsince_;               // minutes since epoch
        double epoch_;                // epoch, julian date

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
        void Print();
        void PrintOneline();
};

} // namespace passpredict

#endif