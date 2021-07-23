#include <iostream>
#include <iomanip>
#include <math.h>

#include "SGP4.h"
#include "orbit.h"
#include "satellite.h"
#include "passpredict.h"
#include "sofa.h"

namespace passpredict
{
    Satellite::Satellite(Orbit orbit)
    {
        this->orbit = orbit;
        m_epoch = orbit.satrec.jdsatepoch + orbit.satrec.jdsatepochF;
    };

    void Satellite::propagate_tsince(double t_tsince)
    {
        /*
            inputs:
                t_tsince    time since epoch in minutes
            */
        // find time since julian date in minutes
        m_jd = m_epoch + (t_tsince / 1440.0);
        m_tsince = t_tsince;
        SGP4Funcs::sgp4(orbit.satrec, m_tsince, m_rteme, m_vteme);

    };

    void Satellite::propagate_jd(double t_jd)
    {
        /*
            inputs:
                t_jd    julian date
            */
        // find time of julian date since epoch in minutes
        m_jd = t_jd;
        m_tsince = (m_epoch - m_jd) / 1440.0;
        SGP4Funcs::sgp4(orbit.satrec, m_tsince, m_rteme, m_vteme);
    };

    int Satellite::teme2ecef()
    {
        /*
        Rotate rteme and vteme to respective ecef vectors
        */
        int i, j, err;
        double gmst, dpsi, deps, dp80, de80, epsa;
        double ddp80 = 0.0, dde80 = 0.0;
        double ee, gst, w;
        double jd_tt, tai_1, tai_2, tt_1, tt_2;
        double rotz[3][3], cosw, sinw, tmp[3] = {0, 0, 0};

        // find GMST
        gmst = iauGmst82(m_jd, 0.0);

        // Find terrestial time
        // Convert UTC to TAI
        err = iauUtctai(m_jd, 0.0, &tai_1, &tai_2);
        if (err) return 1;
        // Convert TAI to TT
        err = iauTaitt(tai_1, tai_2, &tt_1, &tt_2);
        if (err) return 1;

        // Find omega from nutation 1980 theory
        iauNut80(tt_1, tt_2, &dp80, &de80);

        // Add adjustments, in this case assume zero
        dpsi = dp80 + ddp80;
        deps = de80 + dde80;

        // Mean obliquity
        epsa = iauObl80(tt_1, tt_2);

        // Equation of equinoxes
        ee = iauEqeq94(tt_1, tt_2) + ddp80 * std::cos(epsa);

        // find GAST with kinematic terms
        gst = gmst + ee;

        // normalize to between 0 and 2pi
        w = fmod(gst, PASSPREDICT_2PI);
        if (w < 0) w += PASSPREDICT_2PI;

        // Rotate around z axis
        // Create z-rotation matrix and transpose it
        cosw = std::cos(w);
        sinw = std::sin(w);
        rotz[0][0] = cosw;
        rotz[0][1] = -sinw;
        rotz[0][2] = 0;
        rotz[1][0] = sinw;
        rotz[1][1] = cosw;
        rotz[1][2] = 0;
        rotz[2][0] = 0;
        rotz[2][1] = 0;
        rotz[2][2] = 1;

        // rotate position vector
        // recef = rotz * rteme
        for (i=0; i<3; i++) {
            recef[i] = 0.0;
            for (j=0; j<3; j++) {
                recef[i] += rotz[i][j] * m_rteme[j];
            }
        }
    }

    void Satellite::print()
    {
        int i;
        using namespace std;
        cout.precision(10);
        cout << "jd    =  " << fixed << m_jd << endl;
        cout.precision(4);
        cout << "rteme = [" << fixed;
        for (i = 0; i < 3; i++)
        {
            cout << m_rteme[i] << ", ";
        }
        cout << "]" << endl;
        cout << "vteme = [" << fixed;
        for (i = 0; i < 3; i++)
            cout << m_vteme[i] << ", ";
        cout << "]" << endl;
    };

    void Satellite::print_oneline()
    {
        int i;
        using namespace std;
        cout.precision(8);
        cout << "   " << setw(14) << m_tsince;
        for (i = 0; i < 3; i++)
            cout << "   " << setw(14) << m_rteme[i];
        for (i = 0; i < 3; i++)
            cout << "   " << setw(12) << m_vteme[i];
    }
}
