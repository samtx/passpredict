
#include "SGP4.h"
#include "orbit.h"
#include "passpredict.h"

namespace passpredict
{

    Orbit::Orbit(const Omm &omm, gravconsttype whichconst)
        : whichconst_(whichconst)
    {
        bool err;
        double ndot, nddot; // These aren't used in sgp4 so they are dummy variables
        const double deg2rad = PASSPREDICT_DEG2RAD;
        const double xpdotp = 1440.0 / PASSPREDICT_2PI;
        strcpy(satrec_.satnum, omm.satnum);
        satrec_.jdsatepoch = omm.jdsatepoch;
        satrec_.jdsatepochF = omm.jdsatepochF;
        satrec_.bstar = omm.bstar;
        satrec_.ecco = omm.ecco;
        satrec_.argpo = omm.argpo * deg2rad;
        satrec_.inclo = omm.inclo * deg2rad,
        satrec_.mo = omm.mo * deg2rad;
        satrec_.no_kozai = omm.no_kozai / xpdotp;
        satrec_.nodeo = omm.nodeo * deg2rad;
        satrec_.revnum = omm.revnum;
        err = SGP4Funcs::sgp4init(
            whichconst_, 'i', satrec_.satnum, satrec_.jdsatepoch + satrec_.jdsatepochF - 2433281.5,
            satrec_.bstar, ndot, nddot, satrec_.ecco, satrec_.argpo, satrec_.inclo, satrec_.mo,
            satrec_.no_kozai, satrec_.nodeo, satrec_);
    };

    Orbit::Orbit(char *tle1, char *tle2, gravconsttype whichconst)
        : whichconst_(whichconst)
    {
        double dummy;
        SGP4Funcs::twoline2rv(
            tle1, tle2, ' ', ' ', 'i', whichconst_, dummy, dummy, dummy, satrec_);
    };
}