#ifndef ORBIT_H
#define ORBIT_H

#include "SGP4.h"
#include "omm.h"

namespace passpredict
{
    class Orbit
    {
    private:
        gravconsttype m_whichconst = wgs84;

    public:
        const char *tle1;
        const char *tle2;
        elsetrec satrec;

        Orbit(){};
        Orbit(const Omm &omm);
        Orbit(char *atle1, char *atle2);
    };
}

#endif