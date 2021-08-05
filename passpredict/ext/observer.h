#ifndef OBSERVER_H
#define OBSERVER_H

#include <forward_list>
#include <memory>
#include <array>

#include "location.h"
#include "satellite.h"
#include "passpredict.h"

namespace passpredict {

enum class PassType {
    visible = 0,
    unlit,
    daylight,
};

struct Point {
    double az;  // azimuth [degrees]
    double el;  // elevation [degrees]
    double range; // range [km]
    double jd;  // time in jd
};


struct Overpass {
    Point aos;  // time of acquisition of signal in jd
    Point los;  // time of loss of signal in jd
    Point max;  // time of max elevation in jd
    PassType pass_type;  // visible, unlit, daylight
    double brightness;  // brightness magnitude
    double altitude;  // altitude at max elevation in km
    double duration;  // duration of overpass in seconds
    double duration_vis;  // visible duration of overpass in seconds
    Point aos_vis;
    Point los_vis;
    // Overpass(Point a_aos, Point a_los, Point a_max,  ) : aos(Point) {}
};

class Observer
{
private:
    std::shared_ptr<Location> loc_ptr_;
    std::shared_ptr<Satellite> sat_ptr_;
public:
    std::array<double, 3> r_; // position in ECI coordinates (km)
    std::array<double, 3> v_; // velocity in ECI coordinates (km/s)
    double jd_;               // julian date for time
    double el_;               // elevation (degrees)
    double az_;               // azimuth   (degrees)
    double range_;            // range (km)

    Observer(Location location, Satellite satellite);
    void UpdateToJd(double);
    double ComputeElevationAngle();
    void Ecef2Sez(std::array<double, 3>, std::array<double, 3>&);
    void Sez2Razel(std::array<double, 3>);
    double FindLos(double, double);
    double FindAos(double, double);
    std::shared_ptr<Location> GetLocationPtr();
    std::shared_ptr<Satellite> GetSatellitePtr();
    std::shared_ptr<Overpass> GetNextOverpass(double, double);
    std::forward_list<std::shared_ptr<Overpass>> GetOverpasses(double, double);


    // void FindAos(double t0, double tmax);
    // void FindLos(double t0, double tmax);

};

} // passpredict
#endif