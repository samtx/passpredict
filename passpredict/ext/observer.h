#ifndef OBSERVER_H
#define OBSERVER_H

#include <iostream>
#include <list>
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
    double rng; // range [km]
    double jd;  // time in jd
};


struct Overpass {
    PassType pass_type;  // visible, unlit, daylight
    double brightness;  // brightness magnitude
    double altitude;  // altitude at max elevation in km
    double duration;  // duration of overpass in seconds
    double duration_vis;  // visible duration of overpass in seconds
    Point aos;  // time of acquisition of signal in jd
    Point los;  // time of loss of signal in jd
    Point max;  // time of max elevation in jd
    Point aos_vis;
    Point los_vis;
    // Overpass(Point a_aos, Point a_los, Point a_max,  ) : aos(Point) {}

    void PrintLn();
};

class Observer
{
private:
    std::shared_ptr<Location> loc_ptr_;
    std::shared_ptr<Satellite> sat_ptr_;
public:
    double jd_;               // julian date for time
    double el_;               // elevation (degrees)
    double az_;               // azimuth   (degrees)
    double range_;            // range (km)
    std::array<double, 3> r_; // position in ECI coordinates (km)
    std::array<double, 3> v_; // velocity in ECI coordinates (km/s)

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
    std::list<std::shared_ptr<Overpass>> GetOverpasses(double, double);
};

Observer MakeObserver(Location, Satellite);
std::list<std::shared_ptr<Overpass>> Predict(Location, Satellite, double, double);

} // passpredict
#endif