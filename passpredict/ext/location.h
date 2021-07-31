#ifndef LOCATION_H
#define LOCATION_H

namespace passpredict {

class Location
{
public:
    double lat_;      // latitude (degrees)
    double lon_;      // longitude (degrees)
    double h_;        // height above MSL (m)
    double recef_[3]; // ECEF position vector (km)

    Location(){};
    Location(double lat, double lon, double h);
    void SiteEcef();
};

} // namespace passpredict
#endif