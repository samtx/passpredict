#include <iostream>
#include <iomanip>

#include "SGP4.h"
#include "orbit.h"
#include "satellite.h"
#include "passpredict.h"

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

    void Satellite::teme2ecef()
    {
        /*
                Rotate rteme and vteme to respective ecef vectors
            */
        int i;
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
