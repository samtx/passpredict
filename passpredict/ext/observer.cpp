#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <list>
#include <memory>
#include <array>

#include "observer.h"
#include "passpredict.h"

namespace passpredict {

Observer::Observer(Location location, Satellite satellite){
    loc_ptr_ = std::make_shared<Location>(location);
    sat_ptr_ = std::make_shared<Satellite>(satellite);
};


void Observer::ComputeRazel(){
    // Compute observation angle, azimuth, and range
    int i;
    std::array<double, 3> rho;
    std::array<double, 3> rsez;

    // Propagate satellite to time jd_
    sat_ptr_->PropagateJd(jd_);

    // get vector rho from location to satellite
    for (i=0; i<3; i++)
        rho[i] = sat_ptr_->recef_[i] - loc_ptr_->recef_[i];

    // Rotate rho ECEF vector to topographic SEZ vector
    Observer::Ecef2Sez(rho, rsez);

    // Compute elevation angle, range, and azimuth
    Observer::Sez2Razel(rsez);
}

// double Observer::ComputeElevationAngle(double jd, Satellite sat){
//     // Compute observation elevation angle at time jd
//     int i;
//     std::array<double, 3> rho;
//     std::array<double, 3> rsez;

//     // Propagate satellite to time jd_
//     sat_ptr_->PropagateJd(jd);

//     // get vector rho from location to satellite
//     for (i=0; i<3; i++)
//         rho[i] = sat_ptr_->recef_[i] - loc_ptr_->recef_[i];

//     // Rotate rho ECEF vector to topographic SEZ vector
//     Observer::Ecef2Sez(rho, rsez);

//     // Compute elevation angle, range, and azimuth
//     Observer::Sez2Razel(rsez);

//     return el_;

// };

double Observer::FindAos(double t0, double tmax) {
    // find time of acquiring of signal
    double t = t0;
    double aos_time = 0.0;
    double dt, told, t_step;

    std::cout << "t0 = " << std::fixed << std::setprecision(8) << t0 << std::endl;

    Observer::UpdateToJd(t0);

    if (el_ > 0) {
        t = Observer::FindLos(t, tmax);
    }

    if (tmax > 0.0) {

        int i = 0;
        // coarse time steps
        while ((el_ < -1.0) && (t <= tmax)) {
            told = t;
            // t_step = 0.00035 * (el_ * ((sat_ptr_->alt_ / 8400.0) + 0.46) - 2.0);
            t_step = 20. / 86400.0;  // add 20 seconds
            t = t + t_step;
            Observer::UpdateToJd(t);
            // print iterations
            dt = (t - told) * 86400;
            std::cout << "i=" << i << "  dt=" << dt << "  t=" << t << "  el=" << el_ << "  alt=" << sat_ptr_->alt_<< std::endl;
            i++;
        }

        // fine time steps
        int j = 0;
        while ((aos_time == 0.0) && (t <= tmax)) {
            told = t;
            if (std::fabs(el_) < 0.005) {
                aos_time = t;
            }
            else {
                // t -= el_ * std::sqrt(sat_ptr_->alt_) / (5300.0 * 5);
                t += 1.0 / 86400.0;  // add 1 second
                Observer::UpdateToJd(t);
            }
            // print iterations
            dt = (t - told) * 86400;
            std::cout << "j=" << j << "  dt=" << dt <<  "  t=" << t << "  el=" << el_ << "  alt=" << sat_ptr_->alt_<< std::endl;
            j++;
            if (j > 100) break;
        }
    }

    return aos_time;
};

double Observer::FindLos(double t0, double tmax) {
    // Find time of loss of signal
    double los_time = 0.0;
    double el_tmp;
    double t;

    Observer::UpdateToJd(t0);

    if (el_ < 0.0) {
        t = Observer::FindAos(t0, tmax);
    }

    if (t < 0.01) {
        // invalid time, potentially returned by FindAos
        return 0.0;
    }

    if (tmax > 0.0) {

        // coarse time steps
        while ((el_ >= 5.0) && (t <= tmax)) {
            // t += std::cos((el_ - 1.0) * PASSPREDICT_DEG2RAD) * std::sqrt(sat_ptr_->alt_) / 25000.0;
            t += 5.0 / 86400.0;  // add 5 seconds
            Observer::UpdateToJd(t);
        }

        // fine time steps
        while ((los_time == 0.0) && (t < tmax)) {
            // t += el_ * std::sqrt(sat_ptr_->alt_) / 502500.0;
            t += 0.25 / 86400.0;  // add 0.25 seconds
            Observer::UpdateToJd(t);

            if (std::fabs(el_) < 0.1) {
                el_tmp = el_;
                // check elevation 1 second earlier
                Observer::UpdateToJd(t - 1.0/86400.0);

                if (el_ > el_tmp) {
                    los_time = t;
                }
            }
        }
    }

    return los_time;
};


void Observer::UpdateToJd(double jd) {
    sat_ptr_->PropagateJd(jd);
    sat_ptr_->ComputeAltitude();
    Observer::ComputeRazel();
}


void Observer::Ecef2Sez(std::array<double, 3> recef, std::array<double, 3>& rsez) {
    // Rotate vector from ECEF to SEZ frame based on
    // location geodetic coordinates
    double lmda_rad, ang1, cosang1, sinang1;
    double cosang2, sinang2;
    lmda_rad = loc_ptr_->lon_ * PASSPREDICT_DEG2RAD;
    ang1 = (90 - loc_ptr_->lat_) * PASSPREDICT_DEG2RAD;
    cosang1 = std::cos(ang1);
    sinang1 = std::sin(ang1);
    cosang2 = std::cos(lmda_rad);
    sinang2 = std::sin(lmda_rad);
    rsez[0] = cosang1*cosang2*recef[0] + cosang1*sinang2*recef[1] - sinang1*recef[2];
    rsez[1] = -sinang2*recef[0] + cosang2*recef[1] + 0;
    rsez[2] = sinang1*cosang2*recef[0] + sinang1*sinang2*recef[1] + cosang1*recef[2];
};

void Observer::Sez2Razel(std::array<double, 3> rsez) {
    // Get Range, Elevation, and Azimuth from SEZ vector
    range_ = std::sqrt(
        rsez[0]*rsez[0] + rsez[1]*rsez[1] + rsez[2]*rsez[2]
    );
    el_ = std::asin(rsez[2] / range_) * PASSPREDICT_RAD2DEG;
    az_ = std::atan2(rsez[0], rsez[1]) * PASSPREDICT_RAD2DEG;
    if (rsez[0] < 0 && rsez[1] < 0)
        az_ = std::fmod(az_, 360.0);
};

std::shared_ptr<Overpass> Observer::GetNextOverpass(double t0, double tmax) {
    // Find next overpass from time t0

    double t_aos, t_los, t_max;
    auto overpass = std::make_shared<Overpass>();
    Point point_aos, point_los;
    t_aos = Observer::FindAos(t0, tmax);
    t_los = Observer::FindLos(t_aos + 1.0/86400.0, tmax);
    point_aos.jd = t_aos;
    point_los.jd = t_los;
    overpass->aos = point_aos;
    overpass->los = point_los;
    return overpass;
};

std::list<std::shared_ptr<Overpass>> Observer::GetOverpasses(double t0, double tmax) {
    // Return a singly-linked list of pointers to overpass structs

    double t;
    int k;
    // std::shared_ptr<Overpass> overpass_ptr;
    std::list<std::shared_ptr<Overpass>> overpasses;

    // update observer and satellite to time t0
    Observer::UpdateToJd(t0);
    t = t0;
    while ((t < tmax) && (k < 25)) {
        auto overpass_ptr = std::make_shared<Overpass>();
        overpass_ptr = Observer::GetNextOverpass(t, tmax);
        overpasses.push_back(overpass_ptr);
        t = overpass_ptr->los.jd  + 10.0/1440.0; // last LOS time plus 10 minutes
        ++k;
    }

    // reverse the list to have earliest overpass start first
    // overpasses.reverse();

    return overpasses;
};

Observer MakeObserver(Location location, Satellite satellite) {
    Observer observer = Observer(location, satellite);
    return observer;
};

std::list<std::shared_ptr<Overpass>> Predict(Location location, Satellite satellite, double t0, double tmax) {
    Observer observer = Observer(location, satellite);
    std::list<std::shared_ptr<Overpass>> overpasses;
    overpasses = observer.GetOverpasses(t0, tmax);
    return overpasses;
};


void Overpass::PrintLn() {
    using namespace std;
    cout << fixed;
    cout << setprecision(6) << aos.jd << setw(6) << aos.el << endl;
};



}; // namespace passpredict