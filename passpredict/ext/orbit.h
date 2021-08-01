#ifndef ORBIT_H
#define ORBIT_H

#include "SGP4.h"
#include "omm.h"

namespace passpredict {

class Orbit {
    gravconsttype whichconst_;

    public:
        const char *tle1_;
        const char *tle2_;
        elsetrec satrec_;

        Orbit(){};
        Orbit(const Omm &omm) : Orbit(omm, wgs84) {};
        Orbit(const Omm &omm, gravconsttype whichconst);
        Orbit(char *atle1, char *atle2) : Orbit(atle1, atle2, wgs84) {};
        Orbit(char *atle1, char *atle2, gravconsttype whichconst);
};

}

#endif