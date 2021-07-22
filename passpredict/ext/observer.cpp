#include "observer.h"
#include "passpredict.h"

namespace passpredict
{
    Observer::Observer(Location location, Satellite satellite)
    {
        this->location = location;
        this->satellite = satellite;
    }
};