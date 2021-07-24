
#include "SGP4.h"
#include "orbit.h"
#include "passpredict.h"

namespace passpredict
{

    Orbit::Orbit(const Omm &omm, gravconsttype whichconst)
        : whichconst(whichconst)
    {
        bool err;
        double ndot, nddot; // These aren't used in sgp4 so they are dummy variables
        const double deg2rad = PASSPREDICT_DEG2RAD;
        const double xpdotp = 1440.0 / PASSPREDICT_2PI;
        strcpy(satrec.satnum, omm.satnum);
        satrec.jdsatepoch = omm.jdsatepoch;
        satrec.jdsatepochF = omm.jdsatepochF;
        satrec.bstar = omm.bstar;
        satrec.ecco = omm.ecco;
        satrec.argpo = omm.argpo * deg2rad;
        satrec.inclo = omm.inclo * deg2rad,
        satrec.mo = omm.mo * deg2rad;
        satrec.no_kozai = omm.no_kozai / xpdotp;
        satrec.nodeo = omm.nodeo * deg2rad;
        satrec.revnum = omm.revnum;
        err = SGP4Funcs::sgp4init(
            whichconst, 'i', satrec.satnum, satrec.jdsatepoch + satrec.jdsatepochF - 2433281.5,
            satrec.bstar, ndot, nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo,
            satrec.no_kozai, satrec.nodeo, satrec);
    };

    Orbit::Orbit(char *atle1, char *atle2, gravconsttype whichconst)
        : whichconst(whichconst)
    {
        double dummy;
        SGP4Funcs::twoline2rv(
            atle1, atle2, ' ', ' ', 'i', whichconst, dummy, dummy, dummy, satrec);
    };
}