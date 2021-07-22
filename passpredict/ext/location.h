#ifndef LOCATION_H
#define LOCATION_H

namespace passpredict
{
    class Location
    {
    public:
        double lat;      // latitude (degrees)
        double lon;      // longitude (degrees)
        double h;        // height above MSL (m)
        double recef[3]; // ECEF position vector (km)

        Location(){};
        Location(double alat, double alon, double ah);
        void site_ecef();
    };
}
#endif