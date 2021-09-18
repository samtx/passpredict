/*     -------------------------------------------------------------------------
*
*                                ast2Body.cpp
*
*   this file contains fundamental astrodynamic procedures and functions
*   using 2-body dynamics. the routines span a wide range of material, and
*   they come from chapters 2, 3, 5, and 11.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2013
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              11 jan 18  david vallado
*                           misc cleanup
*    changes :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*               4 may 09  david vallado
*                           misc updates
*              23 feb 07  david vallado
*                           3rd edition baseline
*              21 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */

#include "ast2Body.h"


namespace ast2Body
{

	/* -----------------------------------------------------------------------------
	*
	*                           function rv2coe
	*
	*  this function finds the classical orbital elements given the geocentric
	*    equatorial position and velocity vectors.
	*
	*  author        : david vallado                  719-573-2600   21 jun 2002
	*
	*  revisions
	*    vallado     - fix special cases                              5 sep 2002
	*    vallado     - delete extra check in inclination code        16 oct 2002
	*    vallado     - add constant file use                         29 jun 2003
	*
	*  inputs          description                    range / units
	*    r           - ijk position vector            km
	*    v           - ijk velocity vector            km / s
	*
	*  outputs       :
	*    p           - semilatus rectum               km
	*    a           - semimajor axis                 km
	*    ecc         - eccentricity
	*    incl        - inclination                    0.0  to pi rad
	*    raan       - longitude of ascending node    0.0  to 2pi rad
	*    argp        - argument of perigee            0.0  to 2pi rad
	*    nu          - true anomaly                   0.0  to 2pi rad
	*    m           - mean anomaly                   0.0  to 2pi rad
	*    eccanom     - eccentric, parabolic,
	*                  hyperbolic anomaly             rad
	*    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
	*    truelon     - true longitude            (ce) 0.0  to 2pi rad
	*    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
	*
	*  locals        :
	*    hbar        - angular momentum h vector      km2 / s
	*    ebar        - eccentricity     e vector
	*    nbar        - line of nodes    n vector
	*    c1          - v**2 - u/r
	*    rdotv       - r dot v
	*    hk          - hk unit vector
	*    sme         - specfic mechanical energy      km2 / s2
	*    i           - index
	*    temp        - temporary variable
	*    typeorbit   - type of orbit                  ee, ei, ce, ci
	*
	*  coupling      :
	*    mag         - magnitude of a vector
	*    cross       - cross product of two vectors
	*    angle       - find the angle between two vectors
	*    newtonnu    - find the mean anomaly
	*
	*  references    :
	*    vallado       2013, 113, alg 9, ex 2-5
	* --------------------------------------------------------------------------- */

	void rv2coe
		(
		double r[3], double v[3], const double mu,
		double& p, double& a, double& ecc, double& incl, double& raan, double& argp,
		double& nu, double& m, double& eccanom, double& arglat, double& truelon, double& lonper
		)
	{
		double small, hbar[3], nbar[3], magr, magv, magn, ebar[3], sme, rdotv, temp, c1, hk, magh, halfpi;

		int i;
		// switch this to an integer msvs seems to have problems with this and strncpy_s
		//char typeorbit[2];
		int typeorbit;
		// here
		// typeorbit = 1 = 'ei'
		// typeorbit = 2 = 'ce'
		// typeorbit = 3 = 'ci'
		// typeorbit = 4 = 'ee'

		halfpi = 0.5 * pi;
		small = 0.00000001;
		eccanom = 0.0;

		// -------------------------  implementation   -----------------
		magr = astMath::mag(r);
		magv = astMath::mag(v);

		// ------------------  find h n and e vectors   ----------------
		astMath::cross(r, v, hbar);
		magh = astMath::mag(hbar);
		if (magh > small)
		{
			nbar[0] = -hbar[1];
			nbar[1] = hbar[0];
			nbar[2] = 0.0;
			magn = astMath::mag(nbar);
			c1 = magv * magv - mu / magr;
			rdotv = astMath::dot(r, v);
			temp = 1.0 / mu;
			for (i = 0; i <= 2; i++)
				ebar[i] = (c1 * r[i] - rdotv * v[i]) * temp;
			ecc = astMath::mag(ebar);

			// ------------  find a e and semi-latus rectum   ----------
			sme = (magv * magv * 0.5) - (mu / magr);
			if (fabs(sme) > small)
				a = -mu / (2.0 * sme);
			else
				a = infinite;
			p = magh * magh * temp;

			// -----------------  find inclination   -------------------
			hk = hbar[2] / magh;
			incl = acos(hk);

			// --------  determine type of orbit for later use  --------
			// ------ elliptical, parabolic, hyperbolic inclined -------
			//#ifdef _MSC_VER  // chk if compiling under MSVS
			//		   strcpy_s(typeorbit, 2 * sizeof(char), "ei");
			//#else
			//		   strcpy(typeorbit, "ei");
			//#endif
			typeorbit = 1;

			if (ecc < small)
			{
				// ----------------  circular equatorial ---------------
				if ((incl < small) | (fabs(incl - pi) < small))
				{
					//#ifdef _MSC_VER
					//				   strcpy_s(typeorbit, sizeof(typeorbit), "ce");
					//#else
					//				   strcpy(typeorbit, "ce");
					//#endif
					typeorbit = 2;
				}
				else
				{
					// --------------  circular inclined ---------------
					//#ifdef _MSC_VER
					//				   strcpy_s(typeorbit, sizeof(typeorbit), "ci");
					//#else
					//				   strcpy(typeorbit, "ci");
					//#endif
					typeorbit = 3;
				}
			}
			else
			{
				// - elliptical, parabolic, hyperbolic equatorial --
				if ((incl < small) | (fabs(incl - pi) < small))
				{
					//#ifdef _MSC_VER
					//				   strcpy_s(typeorbit, sizeof(typeorbit), "ee");
					//#else
					//				   strcpy(typeorbit, "ee");
					//#endif
					typeorbit = 4;
				}
			}

			// ----------  find right ascension of the ascending node ------------
			if (magn > small)
			{
				temp = nbar[0] / magn;
				if (fabs(temp) > 1.0)
					temp = astMath::sgn(temp);
				raan = acos(temp);
				if (nbar[1] < 0.0)
					raan = twopi - raan;
			}
			else
				raan = undefined;

			// ---------------- find argument of perigee ---------------
			//if (strcmp(typeorbit, "ei") == 0)
			if (typeorbit == 1)
			{
				argp = astMath::angle(nbar, ebar);
				if (ebar[2] < 0.0)
					argp = twopi - argp;
			}
			else
				argp = undefined;

			// ------------  find true anomaly at epoch    -------------
			//if (typeorbit[0] == 'e')
			if ((typeorbit == 1) || (typeorbit == 4))
			{
				nu = astMath::angle(ebar, r);
				if (rdotv < 0.0)
					nu = twopi - nu;
			}
			else
				nu = undefined;

			// ----  find argument of latitude - circular inclined -----
			//if (strcmp(typeorbit, "ci") == 0)
			if (typeorbit == 3)
			{
				arglat = astMath::angle(nbar, r);
				if (r[2] < 0.0)
					arglat = twopi - arglat;
				m = arglat;
			}
			else
				arglat = undefined;

			// -- find longitude of perigee - elliptical equatorial ----
			//if ((ecc>small) && (strcmp(typeorbit, "ee") == 0))
			if ((ecc>small) && (typeorbit == 4))
			{
				temp = ebar[0] / ecc;
				if (fabs(temp) > 1.0)
					temp = astMath::sgn(temp);
				lonper = acos(temp);
				if (ebar[1] < 0.0)
					lonper = twopi - lonper;
				if (incl > halfpi)
					lonper = twopi - lonper;
			}
			else
				lonper = undefined;

			// -------- find true longitude - circular equatorial ------
			//if ((magr>small) && (strcmp(typeorbit, "ce") == 0))
			if ((magr > small) && (typeorbit == 2))
			{
				temp = r[0] / magr;
				if (fabs(temp) > 1.0)
					temp = astMath::sgn(temp);
				truelon = acos(temp);
				if (r[1] < 0.0)
					truelon = twopi - truelon;
				if (incl > halfpi)
					truelon = twopi - truelon;
				m = truelon;
			}
			else
				truelon = undefined;

			// ------------ find mean anomaly for all orbits -----------
			//if (typeorbit[0] == 'e')
			if ((typeorbit == 1) || (typeorbit == 4))
				newtonnu(ecc, nu, eccanom, m);
		}
		else
		{
			p = undefined;
			a = undefined;
			ecc = undefined;
			incl = undefined;
			raan = undefined;
			argp = undefined;
			nu = undefined;
			m = undefined;
			arglat = undefined;
			truelon = undefined;
			lonper = undefined;
		}
	}  // rv2coe


	/* ------------------------------------------------------------------------------
	*
	*                           function coe2rv
	*
	*  this function finds the position and velocity vectors in geocentric
	*    equatorial (ijk) system given the classical orbit elements.
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    p           - semilatus rectum               km
	*    ecc         - eccentricity
	*    incl        - inclination                    0.0 to pi rad
	*    raan       - longitude of ascending node    0.0 to 2pi rad
	*    argp        - argument of perigee            0.0 to 2pi rad
	*    nu          - true anomaly                   0.0 to 2pi rad
	*    arglat      - argument of latitude      (ci) 0.0 to 2pi rad
	*    lamtrue     - true longitude            (ce) 0.0 to 2pi rad
	*    lonper      - longitude of periapsis    (ee) 0.0 to 2pi rad
	*
	*  outputs       :
	*    r           - ijk position vector            km
	*    v           - ijk velocity vector            km / s
	*
	*  locals        :
	*    temp        - temporary real*8 value
	*    rpqw        - pqw position vector            km
	*    vpqw        - pqw velocity vector            km / s
	*    sinnu       - sine of nu
	*    cosnu       - cosine of nu
	*    tempvec     - pqw velocity vector
	*
	*  coupling      :
	*    rot3        - rotation about the 3rd axis
	*    rot1        - rotation about the 1st axis
	*
	*  references    :
	*    vallado       2013, 118, alg 10, ex 2-5
	* --------------------------------------------------------------------------- */

	void coe2rv
		(
		double p, double ecc, double incl, double raan, double argp, double nu,
		double arglat, double truelon, double lonper,
		double r[3], double v[3]
		)
	{
		double rpqw[3], vpqw[3], tempvec[3], temp, sinnu, cosnu, mu, small;

		small = 0.0000001;
		mu = 398600.4418;

		// --------------------  implementation   ----------------------
		//       determine what type of orbit is involved and set up the
		//       set up angles for the special cases.
		// -------------------------------------------------------------
		if (ecc < small)
		{
			// ----------------  circular equatorial  ------------------
			if ((incl < small) | (fabs(incl - pi) < small))
			{
				argp = 0.0;
				raan = 0.0;
				nu = truelon;
			}
			else
			{
				// --------------  circular inclined  ------------------
				argp = 0.0;
				nu = arglat;
			}
		}
		else
		{
			// ---------------  elliptical equatorial  -----------------
			if ((incl < small) | (fabs(incl - pi) < small))
			{
				argp = lonper;
				raan = 0.0;
			}
		}

		// ----------  form pqw position and velocity vectors ----------
		cosnu = cos(nu);
		sinnu = sin(nu);
		temp = p / (1.0 + ecc*cosnu);
		rpqw[0] = temp*cosnu;
		rpqw[1] = temp*sinnu;
		rpqw[2] = 0.0;
		if (fabs(p) < 0.00000001)
			p = 0.00000001;

		vpqw[0] = -sinnu    * sqrt(mu / p);
		vpqw[1] = (ecc + cosnu) * sqrt(mu / p);
		vpqw[2] = 0.0;

		// ----------------  perform transformation to ijk  ------------
		astMath::rot3(rpqw, -argp, tempvec);
		astMath::rot1(tempvec, -incl, tempvec);
		astMath::rot3(tempvec, -raan, r);

		astMath::rot3(vpqw, -argp, tempvec);
		astMath::rot1(tempvec, -incl, tempvec);
		astMath::rot3(tempvec, -raan, v);
	}  // coe2rv


	/* ----------------------------------------------------------------------------
	*
	*                           function rv2eq.m
	*
	*  this function transforms a position and velocity vector into the flight
	*    elements - latgc, lon, fpa, az, position and velocity magnitude.
	*
	*  author        : david vallado                  719 - 573 - 2600    7 jun 2002
	*
	*  inputs          description                    range / units
	*    r           - eci position vector            km
	*    v           - eci velocity vector            km / s
	*
	*  outputs       :
	*    n           - mean motion                    rad
	*    a           - semi major axis                km
	*    af          - component of ecc vector
	*    ag          - component of ecc vector
	*    chi         - component of node vector in eqw
	*    psi         - component of node vector in eqw
	*    meanlon     - mean longitude                 rad
	*    truelon     - true longitude                 rad
	*
	*  locals        :
	*    none -
	*
	*  coupling :
	*    none -
	*
	*  references :
	*    vallado       2013, 108
	*    chobotov            30
	---------------------------------------------------------------------------- */

	void rv2eq
		(
		double r[3], double v[3],
		double& a, double& n, double& af, double& ag, double& chi, double& psi, double& meanlonM, double& meanlonNu, int& fr
		)
	{
		double p, ecc, incl, raan, argp, nu, m, eccanom, arglat, truelon, lonper, mu, small;
		mu = 3.986004418e5;
		small = 0.00000001;

		// --------convert to classical elements----------------------
		ast2Body::rv2coe(r, v, mu, p, a, ecc, incl, raan, argp, nu, m, eccanom, arglat, truelon, lonper);

		// --------setup retrograde factor----------------------------
		fr = 1;
		// ----------set this so it for orbits over 90 deg !!-------- -
		if (incl > pi * 0.5)
			fr = -1;

		if (ecc < small)
		{
			// ----------------circular equatorial------------------
			if ((incl < small) || (abs(incl - pi) < small))
			{
				argp = 0.0;
				raan = 0.0;
				// nu = truelon;
			}
			else
			{
				// --------------circular inclined------------------
				argp = 0.0;
			}
			//  nu = arglat;
		}
		else
		{
			// -------------- - elliptical equatorial---------------- -
			if ((incl < small) || (abs(incl - pi) < small))
			{
				argp = lonper;
				raan = 0.0;
			}
		}

		af = ecc * cos(fr*raan + argp);
		ag = ecc * sin(fr*raan + argp);

		if (fr > 0)
		{
			chi = tan(incl*0.5) * sin(raan);
			psi = tan(incl*0.5) * cos(raan);
		}
		else
		{
			chi = astMath::cot(incl*0.5) * sin(raan);
			psi = astMath::cot(incl*0.5) * cos(raan);
		}

		n = sqrt(mu / (a * a * a));

		meanlonM = fr * raan + argp + m;
		meanlonM = fmod(meanlonM, twopi);

		meanlonNu = fr * raan + argp + nu;
		meanlonNu = fmod(meanlonNu, twopi);
	} // rv2eq


	/* ------------------------------------------------------------------------------
	*
	*                           function eq2rv
	*
	*  this function finds the classical orbital elements given the equinoctial
	*   elements.
	*
	*  author        : david vallado                  719 - 573 - 2600    9 jun 2002
	*
	*  revisions
	*    vallado - fix elliptical equatorial orbits case         19 oct 2002
	* vallado - add constant file use                         29 jun 2003
	*
	*  inputs          description                    range / units
	*    n           - mean motion                    rad
	*    a           - semi major axis                km
	*    af          - component of ecc vector
	*    ag          - component of ecc vector
	*    chi         - component of node vector in eqw
	*    psi         - component of node vector in eqw
	*    meanlon     - mean longitude                 rad
	*    truelon     - true longitude                 rad
	*
	*  outputs       :
	*    r           - position vector                km
	*    v           - velocity vector                km / s
	*
	*  locals :
	*    temp - temporary variable
	*    p - semilatus rectum               km
	*    ecc - eccentricity
	*    incl - inclination                    0.0  to pi rad
	*    raan - longitude of ascending node    0.0  to 2pi rad
	*    argp - argument of perigee            0.0  to 2pi rad
	*    nu - true anomaly                   0.0  to 2pi rad
	*    m - mean anomaly                   0.0  to 2pi rad
	*    arglat - argument of latitude(ci) 0.0  to 2pi rad
	*    truelon - true longitude(ce) 0.0  to 2pi rad
	*    lonper - longitude of periapsis(ee) 0.0  to 2pi rad
	*
	*  coupling      :
	*
	*  references :
	*    vallado 2013 : 108
	------------------------------------------------------------------------------ */

	void eq2rv
		(
		double a, double af, double ag, double chi, double psi, double meanlon, int fr,
		double r[3], double v[3]
		)
	{
		double p, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper, mu, small, e0;
		mu = 3.986004418e5;
		small = 0.00000001;

		// ------------------------ - implementation----------------
		arglat = 999999.1;
		lonper = 999999.1;
		truelon = 999999.1;

		// ---- if n is input----
		//a = (mu / n ^ 2) ^ (1.0 / 3.0);

		ecc = sqrt(af * af + ag * ag);
		p = a * (1.0 - ecc * ecc);
		incl = pi * ((1.0 - fr) * 0.5) + 2.0 * fr * atan(sqrt(chi * chi + psi * psi));
		raan = atan2(chi, psi);
		argp = atan2(ag, af) - fr * atan2(chi, psi);

		if (ecc < small)
		{
			// ----------------circular equatorial------------------
			if ((incl < small) || (abs(incl - pi) < small))
			{
				argp = 0.0;
				raan = 0.0;
				// truelon = nu;
			}
			else
			{
				// --------------circular inclined------------------
				argp = 0.0;
				// arglat = nu;
			}
		}
		else
		{
			// -------------- - elliptical equatorial----------------
			if ((incl < small) || (abs(incl - pi) < small))
			{
				// argp = lonper;
				raan = 0.0;
			}
		}

		m = meanlon - fr * raan - argp;
		m = fmod(m + twopi, twopi);

		newtonm(ecc, m, e0, nu);

		// ----------fix for elliptical equatorial orbits------------
		if (ecc < small)
		{
			// ----------------circular equatorial------------------
			if ((incl < small) || (abs(incl - pi) < small))
			{
				argp = undefined;
				raan = undefined;
				truelon = nu;
			}
			else
			{
				// --------------circular inclined------------------
				argp = undefined;
				arglat = nu;
			}
			nu = undefined;
		}
		else
		{
			// -------------- - elliptical equatorial---------------- -
			if ((incl < small) || (abs(incl - pi) < small))
			{
				lonper = argp;
				argp = undefined;
				raan = undefined;
			}
		}

		// --------now convert back to position and velocity vectors
		coe2rv(p, ecc, incl, raan, argp, nu, arglat, truelon, lonper, r, v);
	}  // eq2rv


	/* -----------------------------------------------------------------------------
	*
	*                           function findc2c3
	*
	*  this function calculates the c2 and c3 functions for use in the universal
	*    variable calculation of z.
	*
	*  author        : david vallado                  719-573-2600   27 may 2002
	*
	*  revisions
	*                -
	*
	*  inputs          description                    range / units
	*    znew        - z variable                     rad2
	*
	*  outputs       :
	*    c2new       - c2 function value
	*    c3new       - c3 function value
	*
	*  locals        :
	*    sqrtz       - square root of znew
	*
	*  coupling      :
	*    sinh        - hyperbolic sine
	*    cosh        - hyperbolic cosine
	*
	*  references    :
	*    vallado       2013, 63, alg 1
	* --------------------------------------------------------------------------- */

	void findc2c3
		(
		double znew,
		double& c2new, double& c3new
		)
	{
		double small, sqrtz;
		small = 0.00000001;

		// -------------------------  implementation   -----------------
		if (znew > small)
		{
			sqrtz = sqrt(znew);
			c2new = (1.0 - cos(sqrtz)) / znew;
			c3new = (sqrtz - sin(sqrtz)) / (sqrtz*sqrtz*sqrtz);
		}
		else
		{
			if (znew < -small)
			{
				sqrtz = sqrt(-znew);
				c2new = (1.0 - cosh(sqrtz)) / znew;
				c3new = (sinh(sqrtz) - sqrtz) / (sqrtz*sqrtz*sqrtz);
			}
			else
			{
				c2new = 0.5;
				c3new = 1.0 / 6.0;
			}
		}
	}  // findc2c3

	/* -----------------------------------------------------------------------------
	*
	*                           function kepler
	*
	*  this function solves keplers problem for orbit determination and returns a
	*    future geocentric equatorial (ijk) position and velocity vector.  the
	*    solution uses universal variables.
	*
	*  author        : david vallado                  719-573-2600   22 jun 2002
	*
	*  revisions
	*    vallado     - fix some mistakes                             13 apr 2004
	*
	*  inputs          description                    range / units
	*    ro          - ijk position vector - initial  km
	*    vo          - ijk velocity vector - initial  km / s
	*    dtsec       - length of time to propagate    s
	*
	*  outputs       :
	*    r           - ijk position vector            km
	*    v           - ijk velocity vector            km / s
	*    error       - error flag                     'ok', ...
	*
	*  locals        :
	*    f           - f expression
	*    g           - g expression
	*    fdot        - f dot expression
	*    gdot        - g dot expression
	*    xold        - old universal variable x
	*    xoldsqrd    - xold squared
	*    xnew        - new universal variable x
	*    xnewsqrd    - xnew squared
	*    znew        - new value of z
	*    c2new       - c2(psi) function
	*    c3new       - c3(psi) function
	*    dtsec       - change in time                 s
	*    timenew     - new time                       s
	*    rdotv       - result of ro dot vo
	*    a           - semi or axis                   km
	*    alpha       - reciprocol  1/a
	*    sme         - specific mech energy           km2 / s2
	*    period      - time period for satellite      s
	*    s           - variable for parabolic case
	*    w           - variable for parabolic case
	*    h           - angular momentum vector
	*    temp        - temporary real*8 value
	*    i           - index
	*
	*  coupling      :
	*    mag         - magnitude of a vector
	*    findc2c3    - find c2 and c3 functions
	*
	*  references    :
	*    vallado       2013, 93, alg 8, ex 2-4
	---------------------------------------------------------------------------- */

	void kepler
		(
		double ro[3], double vo[3], double dtseco, double r[3], double v[3]
		)
	{
		int ktr, i, numiter, mulrev;
		double h[3], f, g, fdot, gdot, rval, xold, xoldsqrd, xnew,
			xnewsqrd, znew, p, c2new, c3new, dtnew, rdotv, a, dtsec,
			alpha, sme, period, s, w, temp,
			magro, magvo, magh, magr, magv;
		char show, errork[10];
		show = 'n';
		double re, velkmps, small, mu, halfpi;

		re = 6378.137;
		velkmps = 7.905365719014;
		mu = 398600.4418;
		small = 0.00000001;
		halfpi = pi * 0.5;

		// -------------------------  implementation   -----------------
		// set constants and intermediate printouts
		numiter = 50;

		if (show == 'y')
		{
			printf(" ro %16.8f %16.8f %16.8f ER \n", ro[0] / re, ro[1] / re, ro[2] / re);
			printf(" vo %16.8f %16.8f %16.8f ER/TU \n", vo[0] / velkmps, vo[1] / velkmps, vo[2] / velkmps);
		}

		// --------------------  initialize values   -------------------
		ktr = 0;
		xold = 0.0;
		znew = 0.0;
		strcpy(errork, "      ok");
		dtsec = dtseco;
		mulrev = 0;

		if (fabs(dtseco) > small)
		{
			magro = astMath::mag(ro);
			magvo = astMath::mag(vo);
			rdotv = astMath::dot(ro, vo);

			// -------------  find sme, alpha, and a  ------------------
			sme = ((magvo * magvo) * 0.5) - (mu / magro);
			alpha = -sme * 2.0 / mu;

			if (fabs(sme) > small)
				a = -mu / (2.0 * sme);
			else
				a = infinite;
			if (fabs(alpha) < small)   // parabola
				alpha = 0.0;

			if (show == 'y')
			{
				printf(" sme %16.8f  a %16.8f alp  %16.8f ER \n", sme / (mu / re), a / re, alpha * re);
				printf(" sme %16.8f  a %16.8f alp  %16.8f km \n", sme, a, alpha);
				printf(" ktr      xn        psi           r          xn+1        dtn \n");
			}

			// ------------   setup initial guess for x  ---------------
			// -----------------  circle and ellipse -------------------
			if (alpha >= small)
			{
				period = twopi * sqrt(fabs(a * a * a) / mu);
				// ------- next if needed for 2body multi-rev ----------
				if (fabs(dtseco) > fabs(period))
					// including the truncation will produce vertical lines that are parallel
					// (plotting chi vs time)
					//                    dtsec = rem( dtseco,period );
					mulrev = int(dtseco / period);
				if (fabs(alpha - 1.0) > small)
					xold = sqrt(mu) * dtsec * alpha;
				else
					// - first guess can't be too close. ie a circle, r=a
					xold = sqrt(mu) * dtsec * alpha * 0.97;
			}
			else
			{
				// --------------------  parabola  ---------------------
				if (fabs(alpha) < small)
				{
					astMath::cross(ro, vo, h);
					magh = astMath::mag(h);
					p = magh * magh / mu;
					s = 0.5  * (halfpi - atan(3.0 * sqrt(mu / (p * p * p)) * dtsec));
					w = atan(pow(tan(s), (1.0 / 3.0)));
					xold = sqrt(p) * (2.0 * astMath::cot(2.0 * w));
					alpha = 0.0;
				}
				else
				{
					// ------------------  hyperbola  ------------------
					temp = -2.0 * mu * dtsec /
						(a * (rdotv + astMath::sgn(dtsec) * sqrt(-mu * a) *	(1.0 - magro * alpha)));
					xold = astMath::sgn(dtsec) * sqrt(-a) * log(temp);
				}
			} // if alpha

			ktr = 1;
			dtnew = -10.0;
			// conv for dtsec to x units
			double tmp = 1.0 / sqrt(mu);

			while ((fabs(dtnew * tmp - dtsec) >= small) && (ktr < numiter))
			{
				xoldsqrd = xold * xold;
				znew = xoldsqrd * alpha;

				// ------------- find c2 and c3 functions --------------
				findc2c3(znew, c2new, c3new);

				// ------- use a newton iteration for new values -------
				rval = xoldsqrd * c2new + rdotv * tmp * xold * (1.0 - znew * c3new) +
					magro * (1.0 - znew * c2new);
				dtnew = xoldsqrd * xold * c3new + rdotv * tmp * xoldsqrd * c2new +
					magro * xold * (1.0 - znew * c3new);

				// ------------- calculate new value for x -------------
				xnew = xold + (dtsec * sqrt(mu) - dtnew) / rval;

				// ----- check if the univ param goes negative. if so, use bissection
				if (xnew < 0.0)
					xnew = xold*0.5;

				if (show == 'y')
				{
					printf("%3i %11.7f %11.7f %11.7f %11.7f %11.7f \n", ktr, xold, znew, rval, xnew, dtnew);
					printf("%3i %11.7f %11.7f %11.7f %11.7f %11.7f \n", ktr, xold / sqrt(re), znew, rval / re, xnew / sqrt(re), dtnew / sqrt(mu));
				}

				ktr = ktr + 1;
				xold = xnew;
			}  // while

			if (ktr >= numiter)
			{
				strcpy(errork, "knotconv");
				printf("not converged in %2i iterations \n", numiter);
				for (i = 0; i < 3; i++)
				{
					v[i] = 0.0;
					r[i] = v[i];
				}
			}
			else
			{
				// --- find position and velocity vectors at new time --
				xnewsqrd = xnew * xnew;
				f = 1.0 - (xnewsqrd * c2new / magro);
				g = dtsec - xnewsqrd * xnew * c3new / sqrt(mu);

				for (i = 0; i < 3; i++)
					r[i] = f * ro[i] + g * vo[i];
				magr = astMath::mag(r);
				gdot = 1.0 - (xnewsqrd * c2new / magr);
				fdot = (sqrt(mu) * xnew / (magro * magr)) * (znew * c3new - 1.0);
				for (i = 0; i < 3; i++)
					v[i] = fdot * ro[i] + gdot * vo[i];
				magv = astMath::mag(v);
				temp = f * gdot - fdot * g;
				if (fabs(temp - 1.0) > 0.00001)
					strcpy(errork, "fandg");

				if (show == 'y')
				{
					printf("f %16.8f g %16.8f fdot %16.8f gdot %16.8f \n", f, g, fdot, gdot);
					printf("f %16.8f g %16.8f fdot %16.8f gdot %16.8f \n", f, g, fdot, gdot);
					printf("r1 %16.8f %16.8f %16.8f ER \n", r[0] / re, r[1] / re, r[2] / re);
					printf("v1 %16.8f %16.8f %16.8f ER/TU \n", v[0] / velkmps, v[1] / velkmps, v[2] / velkmps);
				}
			}
		} // if fabs
		else
			// ----------- set vectors to incoming since 0 time --------
		for (i = 0; i < 3; i++)
		{
			r[i] = ro[i];
			v[i] = vo[i];
		}

		//       fprintf( fid,"%11.5f  %11.5f %11.5f  %5i %3i ",znew, dtseco/60.0, xold/(rad), ktr, mulrev );
	}   // kepler


	// actually for astPert, but leave here for now
	/* ----------------------------------------------------------------------------
	*
	*                           function pkepler
	*
	*  this function propagates a satellite's position and velocity vector over
	*    a given time period accounting for perturbations caused by j2.
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    ro          - original position vector       km
	*    vo          - original velocity vector       km/sec
	*    ndot        - time rate of change of n       rad/sec
	*    nddot       - time accel of change of n      rad/sec2
	*    dtsec       - change in time                 sec
	*
	*  outputs       :
	*    r           - updated position vector        km
	*    v           - updated velocity vector        km/sec
	*
	*  locals        :
	*    p           - semi-paramter                  km
	*    a           - semior axis                    km
	*    ecc         - eccentricity
	*    incl        - inclination                    rad
	*    argp        - argument of periapsis          rad
	*    argpdot     - change in argument of perigee  rad/sec
	*    raan       - longitude of the asc node      rad
	*    raandot    - change in raan                rad
	*    e0          - eccentric anomaly              rad
	*    e1          - eccentric anomaly              rad
	*    m           - mean anomaly                   rad/sec
	*    mdot        - change in mean anomaly         rad/sec
	*    arglat      - argument of latitude           rad
	*    arglatdot   - change in argument of latitude rad/sec
	*    truelon     - true longitude of vehicle      rad
	*    truelondot  - change in the true longitude   rad/sec
	*    lonper     - longitude of periapsis         rad
	*    lonperodot  - longitude of periapsis change  rad/sec
	*    n           - mean angular motion            rad/sec
	*    nuo         - true anomaly                   rad
	*    j2op2       - j2 over p sqyared
	*    sinv,cosv   - sine and cosine of nu
	*
	*  coupling:
	*    rv2coe      - orbit elements from position and velocity vectors
	*    coe2rv      - position and velocity vectors from orbit elements
	*    newtonm     - newton rhapson to find nu and eccentric anomaly
	*
	*  references    :
	*    vallado       2013, 691, alg 65
	* -------------------------------------------------------------------------- - */

	void pkepler
		(
		double r1[3], double v1[3], double& dtsec, double& ndot, double& nddot, double r2[3], double v2[3]
		)
	{
		double  truelondot, arglatdot, lonperdot, e0;
		double p, a, ecc, incl, raan, argp, nu, m, eccanom, arglat, truelon, lonper;

		double re = 6378.137;
		double velkmps = 7.905365719014;
		double mu = 398600.4418;
		double small = 0.00000001;
		double halfpi = pi * 0.5;
		double j2 = 0.00108263;

		ast2Body::rv2coe(r1, v1, mu, p, a, ecc, incl, raan, argp, nu, m, eccanom, arglat, truelon, lonper);

		//%     fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
		//%     fprintf(1,'coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n',...
		//%             p,a,ecc,incl*rad,raan*rad,argp*rad,nu*rad,m*rad, ...
		//%             arglat*rad,truelon*rad,lonper*rad );

		double n1 = (mu / (a * a * a));
		double n = sqrt(n1);

		//     % ------------- find the value of j2 perturbations -------------
		double j2op2 = (n * 1.5 * re * re * j2) / (p * p);
		//%     nbar    = n*( 1.0 + j2op2*sqrt(1.0-ecc*ecc)* (1.0 - 1.5*sin(incl)*sin(incl)) );
		double raandot = -j2op2 * cos(incl);
		double argpdot = j2op2 * (2.0 - 2.5 * sin(incl) * sin(incl));
		double mdot = n;

		a = a - 2.0 * ndot * dtsec * a / (3.0 * n);
		ecc = ecc - 2.0 * (1.0 - ecc) * ndot * dtsec / (3.0 * n);
		p = a * (1.0 - ecc * ecc);

		// ----- update the orbital elements for each orbit type --------
		if (ecc < small)
		{
			//  -------------  circular equatorial  ----------------
			if ((incl < small) | (fabs(incl - pi) < small))
			{
				truelondot = raandot + argpdot + mdot;
				truelon = truelon + truelondot * dtsec;
				truelon = fmod(truelon, twopi);
			}
			else
			{
				//  -------------  circular inclined    --------------
				raan = raan + raandot * dtsec;
				raan = fmod(raan, twopi);
				arglatdot = argpdot + mdot;
				arglat = arglat + arglatdot * dtsec;
				arglat = fmod(arglat, twopi);
			}
		}
		else
		{
			//  ---- elliptical, parabolic, hyperbolic equatorial ---
			if ((incl < small) | (fabs(incl - pi) < small))
			{
				lonperdot = raandot + argpdot;
				lonper = lonper + lonperdot * dtsec;
				lonper = fmod(lonper, twopi);
				m = m + mdot*dtsec + ndot * dtsec * dtsec + nddot * pow(dtsec, 3);
				m = fmod(m, twopi);
				ast2Body::newtonm(ecc, m, e0, nu);
			}

			else
			{
				//  ---- elliptical, parabolic, hyperbolic inclined --
				raan = raan + raandot * dtsec;
				raan = fmod(raan, twopi);
				argp = argp + argpdot  * dtsec;
				argp = fmod(argp, twopi);
				m = m + mdot * dtsec + ndot * dtsec * dtsec + nddot * dtsec * dtsec * dtsec;
				m = fmod(m, twopi);
				ast2Body::newtonm(ecc, m, e0, nu);
			}
		}

		// ------------- use coe2rv to find new vectors ---------------

		ast2Body::coe2rv(p, ecc, incl, raan, argp, nu, arglat, truelon, lonper, r2, v2);

		//%        fprintf(1,'r    %15.9f%15.9f%15.9f',r );
		//%        fprintf(1,' v %15.10f%15.10f%15.10f\n',v );
	}   // pkepler



	/* -----------------------------------------------------------------------------
	*
	*                           function rv2rsw
	*
	*  this function converts position and velocity vectors into radial, along-
	*    track, and cross-track coordinates. note that sometimes the middle vector
	*    is called in-track.
	*
	*  author        : david vallado                  719-573-2600    9 jun 2002
	*
	*  revisions
	*                -
	*
	*  inputs          description                    range / units
	*    r           - position vector                km
	*    v           - velocity vector                km/s
	*
	*  outputs       :
	*    rrsw        - position vector                km
	*    vrsw        - velocity vector                km/s
	*    transmat    - transformation matrix
	*
	*  locals        :
	*    tempvec     - temporary vector
	*    rvec,svec,wvec - direction cosines
	*
	*  coupling      :
	*
	*  references    :
	*    vallado       2013, 164
	* --------------------------------------------------------------------------- */

	void rv2rsw
		(
		double r[3], double v[3],
		double rrsw[3], double vrsw[3], std::vector< std::vector<double> > &transmat
		)
	{
		double rvec[3], svec[3], wvec[3], tempvec[3];

		// --------------------  Implementation   ----------------------
		// in order to work correctly each of the components must be
		// unit vectors
		// radial component
		astMath::norm(r, rvec);

		// ncross-track component
		astMath::cross(r, v, tempvec);
		astMath::norm(tempvec, wvec);

		// along-track component
		astMath::cross(wvec, rvec, tempvec);
		astMath::norm(tempvec, svec);

		// assemble transformation matrix from to rsw frame (individual
		//  components arranged in row vectors)
		transmat[0][0] = rvec[0];
		transmat[0][1] = rvec[1];
		transmat[0][2] = rvec[2];
		transmat[1][0] = svec[0];
		transmat[1][1] = svec[1];
		transmat[1][2] = svec[2];
		transmat[2][0] = wvec[0];
		transmat[2][1] = wvec[1];
		transmat[2][2] = wvec[2];

		astMath::matvecmult(transmat, r, rrsw);
		astMath::matvecmult(transmat, v, vrsw);
		/*
		*   alt approach
		*       rrsw[0] = mag(r)
		*       rrsw[2] = 0.0
		*       rrsw[3] = 0.0
		*       vrsw[0] = dot(r, v)/rrsw(0)
		*       vrsw[1] = sqrt(v(0)**2 + v(1)**2 + v(2)**2 - vrsw(0)**2)
		*       vrsw[2] = 0.0
		*/
	}  //  rv2rsw


	/* -----------------------------------------------------------------------------
	*
	*                           function rv2ntw
	*
	*  this function converts position and velocity vectors into normal,
	*    velocity, and cross-track coordinates. this is the ntw system in vallado.
	*
	*  author        : david vallado                  719-573-2600    5 jul 2002
	*
	*  revisions
	*                -
	*
	*  inputs          description                    range / units
	*    r           - position vector                km
	*    v           - velocity vector                km/s
	*
	*  outputs       :
	*    rntw        - position vector                km
	*    vntw        - velocity vector                km/s
	*    transmat    - transformation matrix
	*
	*  locals        :
	*    tempvec     - temporary vector
	*    nvec,tvec,wvec - direction cosines
	*
	*  coupling      :
	*
	*  references    :
	*    vallado       2013, 164
	* --------------------------------------------------------------------------- */

	void  rv2ntw
		(
		double r[3], double v[3],
		double rntw[3], double vntw[3], std::vector< std::vector<double> > &transmat
		)
	{
		double tvec[3], nvec[3], wvec[3], tempvec[3];

		// --------------------  Implementation   ----------------------
		// in order to work correctly each of the components must be
		// unit vectors
		// in-velocity component
		astMath::norm(v, tvec);

		// cross-track component
		astMath::cross(r, v, tempvec);
		astMath::norm(tempvec, wvec);

		// along-radial component
		astMath::cross(tvec, wvec, tempvec);
		astMath::norm(tempvec, nvec);

		// assemble transformation matrix from to ntw frame (individual
		//  components arranged in row vectors)
		transmat[0][0] = nvec[0];
		transmat[0][1] = nvec[1];
		transmat[0][2] = nvec[2];
		transmat[1][0] = tvec[0];
		transmat[1][1] = tvec[1];
		transmat[1][2] = tvec[2];
		transmat[2][0] = wvec[0];
		transmat[2][1] = wvec[1];
		transmat[2][2] = wvec[2];

		astMath::matvecmult(transmat, r, rntw);
		astMath::matvecmult(transmat, v, vntw);
	}  // rv2ntw


	/* ----------------------------------------------------------------------------
	*
	*                           procedure newtonm
	*
	*  this procedure performs the newton rhapson iteration to find the
	*    eccentric anomaly given the mean anomaly.  the true anomaly is also
	*    calculated.
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    ecc         - eccentricity                   0.0 to
	*    m           - mean anomaly                   -2pi to 2pi rad
	*
	*  outputs       :
	*    e0          - eccentric anomaly              0.0 to 2pi rad
	*    nu          - true anomaly                   0.0 to 2pi rad
	*
	*  locals        :
	*    e1          - eccentric anomaly, next value  rad
	*    sinv        - sine of nu
	*    cosv        - cosine of nu
	*    ktr         - index
	*    r1r         - cubic roots - 1 to 3
	*    r1i         - imaginary component
	*    r2r         -
	*    r2i         -
	*    r3r         -
	*    r3i         -
	*    s           - variables for parabolic solution
	*    w           - variables for parabolic solution
	*
	*  coupling      :
	*    atan2       - arc tangent function which also resloves quadrants
	*    cubic       - solves a cubic polynomial
	*    power       - raises a base number to an arbitrary power
	*    sinh        - hyperbolic sine
	*    cosh        - hyperbolic cosine
	*    sgn         - returns the sign of an argument
	*
	*  references    :
	*    vallado       2013, 65, alg 2, ex 2-1
	* --------------------------------------------------------------------------- */

	void newtonm
		(
		double ecc, double m, double& e0, double& nu
		)
	{
		const int numiter = 50;
		const double small = 0.00000001;       // small value for tolerances

		double e1, sinv, cosv, cose1, coshe1, temp, r1r = 0.0;
		int ktr;

		/* -------------------------- hyperbolic  ----------------------- */
		if ((ecc - 1.0) > small)
		{
			/* ------------  initial guess ------------- */
			if (ecc < 1.6)
			if (((m < 0.0) && (m > -pi)) || (m > pi))
				e0 = m - ecc;
			else
				e0 = m + ecc;
			else
			if ((ecc < 3.6) && (fabs(m) > pi)) // just edges)
				e0 = m - astMath::sgn(m) * ecc;
			else
				e0 = m / (ecc - 1.0); // best over 1.8 in middle
			ktr = 1;
			e1 = e0 + ((m - ecc * sinh(e0) + e0) / (ecc * cosh(e0) - 1.0));
			while ((fabs(e1 - e0) > small) && (ktr <= numiter))
			{
				e0 = e1;
				e1 = e0 + ((m - ecc * sinh(e0) + e0) / (ecc * cosh(e0) - 1.0));
				ktr++;
			}
			// ---------  find true anomaly  -----------
			coshe1 = cosh(e1);
	    	sinv = -(sqrt(ecc * ecc - 1.0) * sinh(e1)) / (1.0 - ecc * coshe1);
			cosv = (coshe1 - ecc) / (1.0 - ecc * coshe1);
			nu = atan2(sinv, cosv);
		}
		else
		{
			/* ---------------------- parabolic ------------------------- */
			if (fabs(ecc - 1.0) < small)
			{
				//kbn      cubic(1.0 / 3.0, 0.0, 1.0, -m, r1r, r1i, r2r, r2i, r3r, r3i);
				e0 = r1r;
				//kbn      if (fileout != null)
				//        fprintf(fileout, "roots %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f\n",
				//                          r1r, r1i, r2r, r2i, r3r, r3i);
				/*
					 s  = 0.5 * (halfpi - atan(1.5 * m));
					 w  = atan(power(tan(s), 1.0 / 3.0));
					 e0 = 2.0 * cot(2.0* w );
					 */
				ktr = 1;
				nu = 2.0 * atan(e0);
			}
			else
			{
				/* --------------------- elliptical --------------------- */
				if (ecc > small)
				{
					/* ------------  initial guess ------------- */
					if (((m < 0.0) && (m > -pi)) || (m > pi))
						e0 = m - ecc;
					else
						e0 = m + ecc;
					ktr = 1;
					e1 = e0 + (m - e0 + ecc * sin(e0)) / (1.0 - ecc * cos(e0));
					while ((fabs(e1 - e0) > small) && (ktr <= numiter))
					{
						ktr++;
						e0 = e1;
						e1 = e0 + (m - e0 + ecc * sin(e0)) / (1.0 - ecc * cos(e0));
					}
					/* ---------  find true anomaly  ----------- */
					cose1 = cos(e1);
					temp = 1.0 / (1.0 - ecc * cose1);
					sinv = (sqrt(1.0 - ecc * ecc) * sin(e1)) * temp;
					cosv = (cose1 - ecc) * temp;
					nu = atan2(sinv, cosv);
				}
				else
				{
					/* --------------------- circular --------------------- */
					ktr = 0;
					nu = m;
					e0 = m;
				}
			}
		}
		if (ktr > numiter)
			printf("newtonrhapson not converged in %3d iterations\n", numiter);
	}    // procedure newtonm



	/* -----------------------------------------------------------------------------
	*
	*                           function newtonnu
	*
	*  this function solves keplers equation when the true anomaly is known.
	*    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
	*    the parabolic limit at 168� is arbitrary. the hyperbolic anomaly is also
	*    limited. the hyperbolic sine is used because it's not double valued.
	*
	*  author        : david vallado                  719-573-2600   27 may 2002
	*
	*  revisions
	*    vallado     - fix small                                     24 sep 2002
	*
	*  inputs          description                    range / units
	*    ecc         - eccentricity                   0.0  to
	*    nu          - true anomaly                   -2pi to 2pi rad
	*
	*  outputs       :
	*    e0          - eccentric anomaly              0.0  to 2pi rad       153.02 deg
	*    m           - mean anomaly                   0.0  to 2pi rad       151.7425 deg
	*
	*  locals        :
	*    e1          - eccentric anomaly, next value  rad
	*    sine        - sine of e
	*    cose        - cosine of e
	*    ktr         - index
	*
	*  coupling      :
	*    arcsinh     - arc hyperbolic sine
	*    sinh        - hyperbolic sine
	*
	*  references    :
	*    vallado       2013, 77, alg 5
	* --------------------------------------------------------------------------- */

	void newtonnu
		(
		double ecc, double nu,
		double& e0, double& m
		)
	{
		double small, sine, cose, cosnu, temp;

		// ---------------------  implementation   ---------------------
		e0 = 999999.9;
		m = 999999.9;
		small = 0.00000001;

		// --------------------------- circular ------------------------
		if (fabs(ecc) < small)
		{
			m = nu;
			e0 = nu;
		}
		else
			// ---------------------- elliptical -----------------------
		if (ecc < 1.0 - small)
		{
			cosnu = cos(nu);
			temp = 1.0 / (1.0 + ecc * cosnu);
			sine = (sqrt(1.0 - ecc * ecc) * sin(nu)) * temp;
			cose = (ecc + cosnu) * temp;
			e0 = atan2(sine, cose);
			m = e0 - ecc*sin(e0);
		}
		else
			// -------------------- hyperbolic  --------------------
		if (ecc > 1.0 + small)
		{
			if ((ecc > 1.0) && (fabs(nu) + 0.00001 < pi - acos(1.0 / ecc)))
			{
				sine = (sqrt(ecc * ecc - 1.0) * sin(nu)) / (1.0 + ecc * cos(nu));
				e0 = astMath::asinh(sine);
				m = ecc*sinh(e0) - e0;
			}
		}
		else
			// ----------------- parabolic ---------------------
		if (fabs(nu) < 168.0 * pi / 180.0)
		{
			e0 = tan(nu * 0.5);
			m = e0 + (e0 * e0 * e0) / 3.0;
		}

		if (ecc < 1.0)
		{
			m = fmod(m, 2.0 * pi);
			if (m < 0.0)
				m = m + 2.0 * pi;
			e0 = fmod(e0, 2.0 *pi);
		}
	}  // newtonnu



	/* -----------------------------------------------------------------------------
	*
	*                           function gc_gd
	*
	*  this function converts from geodetic to geocentric latitude for positions
	*    on the surface of the earth.  notice that (1-f) squared = 1-esqrd.
	*
	*  author        : david vallado                  719-573-2600    6 dec 2005
	*
	*  revisions
	*
	*  inputs          description                    range / units
	*    latgd       - geodetic latitude              -pi to pi rad
	*
	*  outputs       :
	*    latgc       - geocentric latitude            -pi to pi rad
	*
	*  locals        :
	*    none.
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 140, eq 3-11
	* --------------------------------------------------------------------------- */

	void gc_gd
		(
		double&    latgc,
		edirection direct,
		double&    latgd
		)
	{
		const double eesqrd = 0.006694385000;     // eccentricity of earth sqrd

		if (direct == eTo)
			latgd = atan(tan(latgc) / (1.0 - eesqrd));
		else
			latgc = atan((1.0 - eesqrd) * tan(latgd));
	}   // gc_gd


	/* -----------------------------------------------------------------------------
	*
	*                           function ijk2ll
	*
	*  these subroutines convert a geocentric equatorial position vector into
	*    latitude and longitude.  geodetic and geocentric latitude are found. the
	*    inputs must be ecef.
	*
	*  author        : david vallado                  719-573-2600    6 dec 2005
	*
	*  revisions
	*
	*  inputs          description                    range / units
	*    recef       - ecef position vector           km
	*
	*  outputs       :
	*    latgc       - geocentric latitude            -pi to pi rad
	*    latgd       - geodetic latitude              -pi to pi rad
	*    lon         - longitude (west -)             -2pi to 2pi rad
	*    hellp       - height above the ellipsoid     km
	*
	*  locals        :
	*    temp        - diff between geocentric/
	*                  geodetic lat                   rad
	*    sintemp     - sine of temp                   rad
	*    olddelta    - previous value of deltalat     rad
	*    rtasc       - right ascension                rad
	*    decl        - declination                    rad
	*    i           - index
	*
	*  coupling      :
	*    mag         - magnitude of a vector
	*    gcgd        - converts between geocentric and geodetic latitude
	*
	*  references    :
	*    vallado       2013, 173, alg 12 and alg 13, ex 3-3
	* --------------------------------------------------------------------------- */

	void ijk2ll
		(
		double recef[3],
		double& latgc, double& latgd, double& lon, double& hellp
		)
	{
		const double small = 0.00000001;         // small value for tolerances
		const double re = 6378.137;
		const double eesqrd = 0.006694385000;     // eccentricity of earth sqrd
		double magr, decl, rtasc, olddelta, temp, sintemp, s, c = 0.0;
		int i;

		// ---------------------------  implementation   -----------------------
		magr = astMath::mag(recef);

		// ---------------------- find longitude value  ------------------------
		temp = sqrt(recef[0] * recef[0] + recef[1] * recef[1]);
		if (fabs(temp) < small)
			rtasc = astMath::sgn(recef[2]) * pi * 0.5;
		else
			rtasc = atan2(recef[1], recef[0]);

		lon = rtasc;
		if (fabs(lon) >= pi)   // mod it ?
		{
			if (lon < 0.0)
				lon = twopi + lon;
			else
				lon = lon - twopi;

		}
		decl = asin(recef[2] / magr);
		latgd = decl;

		// ----------------- iterate to find geodetic latitude -----------------
		i = 1;
		olddelta = latgd + 10.0;

		while ((fabs(olddelta - latgd) >= small) && (i<10))
		{
			olddelta = latgd;
			sintemp = sin(latgd);
			c = re / (sqrt(1.0 - eesqrd * sintemp * sintemp));
			latgd = atan((recef[2] + c*eesqrd*sintemp) / temp);
			i = i + 1;
		}

		if ((pi * 0.5 - fabs(latgd)) > pi / 180.0)  // 1 deg
			hellp = (temp / cos(latgd)) - c;
		else
		{
			s = c * (1.0 - eesqrd);
			hellp = recef[2] / sin(latgd) - s;
		}

		gc_gd(latgc, eFrom, latgd);
	}   // ijk2ll


	/*------------------------------------------------------------------------------
	*
	*                           procedure rvsez_razel
	*
	*  this procedure converts range, azimuth, and elevation values with slant
	*    range and velocity vectors for a satellite from a radar site in the
	*    topocentric horizon (sez) system.
	*
	*  author        : david vallado                  719-573-2600   22 jun 2002
	*
	*  inputs          description                    range / units
	*    rhovec      - sez satellite range vector     km
	*    drhovec     - sez satellite velocity vector  km/s
	*    direct      -  direction to convert          eFrom  eTo
	*
	*  outputs       :
	*    rho         - satellite range from site      mk
	*    az          - azimuth                        0.0 to 2pi rad
	*    el          - elevation                      -Math.PI/2 to Math.PI/2 rad
	*    drho        - range rate                     km/s
	*    daz         - azimuth rate                   rad/s
	*    del         - elevation rate                 rad/s
	*
	*  locals        :
	*    Math.Sinel       - variable for Math.Sin( el )
	*    Math.Cosel       - variable for Math.Cos( el )
	*    Math.Sinaz       - variable for Math.Sin( az )
	*    Math.Cosaz       - variable for Math.Cos( az )
	*    temp        -
	*    temp1       -
	*
	*  coupling      :
	*    astMath::mag         - astMath::magnitude of a vector
	*    Math.Sign         - returns the sign of a variable
	*    dot         - dot product
	*    arcsin      - arc Math.Sine function
	*    Math.Atan2       - arc tangent function that resolves quadrant ambiguites
	*
	*  references    :
	*    vallado       2013, 261, eq 4-4, eq 4-5
	-----------------------------------------------------------------------------*/

	void rvsez_razel
		(
		double rhosez[3], double drhosez[3],
		edirection direct,
		double& rho, double& az, double& el, double& drho, double& daz, double& del
		)
	{
		const double small = 0.00000001;
		const double halfpi = pi / 2.0;

		double temp1, temp, sinel, cosel, sinaz, cosaz;

		if (direct == eFrom)
		{
			sinel = sin(el);
			cosel = cos(el);
			sinaz = sin(az);
			cosaz = cos(az);

			/* ----------------- form sez range vector ------------------ */
			rhosez[0] = (-rho * cosel * cosaz);
			rhosez[1] = (rho * cosel * sinaz);
			rhosez[2] = (rho * sinel);

			/* --------------- form sez velocity vector ----------------- */
			drhosez[0] = (-drho * cosel * cosaz +
				rhosez[2] * del * cosaz + rhosez[1] * daz);
			drhosez[1] = (drho * cosel * sinaz -
				rhosez[2] * del * sinaz - rhosez[0] * daz);
			drhosez[2] = (drho * sinel + rho * del * cosel);
		}
		else
		{
			/* ------------ calculate azimuth and elevation ------------- */
			temp = sqrt(rhosez[0] * rhosez[0] + rhosez[1] * rhosez[1]);
			if (fabs(rhosez[1]) < small)
			if (temp < small)
			{
				temp1 = sqrt(drhosez[0] * drhosez[0] +
					drhosez[1] * drhosez[1]);
				az = atan2(drhosez[1] / temp1, drhosez[0] / temp1);
			}
			else
			if (drhosez[0] > 0.0)
				az = pi;
			else
				az = 0.0;
			else
				az = atan2(rhosez[1] / temp, rhosez[0] / temp);

			if (temp < small)   // directly over the north pole
				el = astMath::sgn(rhosez[2]) * halfpi;  // +- 90
			else
				el = asin(rhosez[2] / astMath::mag(rhosez));

			/* ------  calculate range, azimuth and elevation rates ----- */
			drho = astMath::dot(rhosez, drhosez) / rho;
			if (fabs(temp * temp) > small)
				daz = (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) /
				(temp * temp);
			else
				daz = 0.0;

			if (fabs(temp) > small)
				del = (drhosez[2] - drho * sin(el)) / temp;
			else
				del = 0.0;
		}
	}   // rvsez_razel



	/* ------------------------------------------------------------------------------
	*
	*                            function rv2radec
	*
	*   this function converts the right ascension and declination values with
	*     position and velocity vectors of a satellite. uses velocity vector to
	*     find the solution of singular cases.
	*
	*   author        : david vallado                  719-573-2600   25 jun 2002
	*
	*   revisions
	*     vallado     - fix rtasc tests                               29 sep 2002
	*
	*   inputs          description                    range / units
	*     r           -  position vector eci           km
	*     v           -  velocity vector eci           km/s
	*
	*   outputs       :
	*     rr          - radius of the satellite        km
	*     rtasc       - right ascension                rad
	*     decl        - declination                    rad
	*     drr         - radius of the satellite rate   km/s
	*     drtasc      - right ascension rate           rad/s
	*     ddecl       - declination rate               rad/s
	*
	*   locals        :
	*     temp        - temporary position vector
	*     temp1       - temporary variable
	*
	*   coupling      :
	*     none
	*
	*   references    :
	*     vallado       2001, 246-248, alg 25
	*
	*   rv2radec( r, v, out rr, out rtasc, out decl, out drr, out drtasc, out ddecl);
	------------------------------------------------------------------------------  */

	void rv2radec
		(
		double r[3], double v[3],
		double& rr, double& rtasc, double& decl, double& drr, double& drtasc, double& ddecl
		)
	{
		double small = 0.00000001;
		double temp, temp1;

		rr = 0.0;
		rtasc = 0.0;
		decl = 0.0;
		drr = 0.0;
		drtasc = 0.0;
		ddecl = 0.0;

		// -------------------------  implementation   -------------------------
		// ------------- calculate angles and rates ----------------
		rr = astMath::mag(r);
		temp = sqrt(r[0] * r[0] + r[1] * r[1]);
		if (temp < small)
			rtasc = atan2(v[1], v[0]);
		else
			rtasc = atan2(r[1], r[0]);
		decl = asin(r[2] / rr);

		temp1 = -r[1] * r[1] - r[0] * r[0];  // different now
		drr = astMath::dot(r, v) / rr;
		if (fabs(temp1) > small)
			drtasc = (v[0] * r[1] - v[1] * r[0]) / temp1;
		else
			drtasc = 0.0;
		if (fabs(temp) > small)
			ddecl = (v[2] - drr * sin(decl)) / temp;
		else
			ddecl = 0.0;

	}  // rv2radec

	/*------------------------------------------------------------------------------
	*
	*                           procedure rv_razel
	*
	*  this procedure converts range, azimuth, and elevation and their rates with
	*    the geocentric equatorial (ecef) position and velocity vectors.  notice the
	*    value of small as it can affect rate term calculations. uses velocity
	*    vector to find the solution of Math.Singular cases.
	*
	*  author        : david vallado                  719-573-2600   22 jun 2002
	*
	*  inputs          description                    range / units
	*    recef       - ecef position vector           km
	*    vecef       - ecef velocity vector           km/s
	*    rsecef      - ecef site position vector      km
	*    latgd       - geodetic latitude              -Math.PI/2 to Math.PI/2 rad
	*    lon         - geodetic longitude             -2pi to Math.PI rad
	*    direct      -  direction to convert          eFrom  eTo
	*
	*  outputs       :
	*    rho         - satellite range from site      km
	*    az          - azimuth                        0.0 to 2pi rad
	*    el          - elevation                      -Math.PI/2 to Math.PI/2 rad
	*    drho        - range rate                     km/s
	*    daz         - azimuth rate                   rad/s
	*    del         - elevation rate                 rad/s
	*
	*  locals        :
	*    rhovecef    - ecef range vector from site    km
	*    drhovecef   - ecef velocity vector from site km/s
	*    rhosez      - sez range vector from site     km
	*    drhosez     - sez velocity vector from site  km
	*    tempvec     - temporary vector
	*    temp        - temporary extended value
	*    temp1       - temporary extended value
	*    i           - index
	*
	*  coupling      :
	*    astMath::mag         - astMath::magnitude of a vector
	*    addvec      - add two vectors
	*    rot3        - rotation about the 3rd axis
	*    rot2        - rotation about the 2nd axis
	*    Math.Atan2       - arc tangent function which also resloves quadrants
	*    dot         - dot product of two vectors
	*    rvsez_razel - find r and v from site in topocentric horizon (sez) system
	*    lncom2      - combine two vectors and constants
	*    arcsin      - arc Math.Sine function
	*    Math.Sign         - returns the sign of a variable
	*
	*  references    :
	*    vallado       2013, 265, alg 27
	-----------------------------------------------------------------------------*/

	void rv_razel
		(
		double recef[3], double vecef[3], double rsecef[3], double latgd, double lon,
		edirection direct,
		double& rho, double& az, double& el, double& drho, double& daz, double& del
		)
	{
		const double halfpi = pi / 2.0;
		const double small = 0.0000001;

		double temp, temp1;
		double rhoecef[3];
		double drhoecef[3];
		double rhosez[3];
		double drhosez[3];
		double tempvec[3];

		if (direct == eFrom)
		{
			/* ---------  find sez range and velocity vectors ----------- */
			rvsez_razel(rhosez, drhosez, direct, rho, az, el, drho, daz, del);

			/* ----------  perform sez to ecef transformation ------------ */
			astMath::rot2(rhosez, latgd - halfpi, tempvec);
			astMath::rot3(tempvec, -lon, rhoecef);
			astMath::rot2(drhosez, latgd - halfpi, tempvec);
			astMath::rot3(tempvec, -lon, drhoecef);

			/* ---------  find ecef range and velocity vectors -----------*/
			astMath::addvec(1.0, rhoecef, 1.0, rsecef, recef);
			vecef[0] = drhoecef[0];
			vecef[1] = drhoecef[1];
			vecef[2] = drhoecef[2];
		}
		else
		{
			/* ------- find ecef range vector from site to satellite ----- */
			astMath::addvec(1.0, recef, -1.0, rsecef, rhoecef);
			drhoecef[0] = vecef[0];
			drhoecef[1] = vecef[1];
			drhoecef[2] = vecef[2];
			rho = astMath::mag(rhoecef);

			/* ------------ convert to sez for calculations ------------- */
			astMath::rot3(rhoecef, lon, tempvec);
			astMath::rot2(tempvec, halfpi - latgd, rhosez);
			astMath::rot3(drhoecef, lon, tempvec);
			astMath::rot2(tempvec, halfpi - latgd, drhosez);

			/* ------------ calculate azimuth and elevation ------------- */
			temp = sqrt(rhosez[0] * rhosez[0] + rhosez[1] * rhosez[1]);
			if (fabs(rhosez[1]) < small)
			if (temp < small)
			{
				temp1 = sqrt(drhosez[0] * drhosez[0] +
					drhosez[1] * drhosez[1]);
				az = atan2(drhosez[1] / temp1, -drhosez[0] / temp1);
			}
			else
			if (rhosez[0] > 0.0)
				az = pi;
			else
				az = 0.0;
			else
				az = atan2(rhosez[1] / temp, -rhosez[0] / temp);

			if (temp < small)  // directly over the north pole
				el = astMath::sgn(rhosez[2]) * halfpi; // +- 90
			else
				el = asin(rhosez[2] / astMath::mag(rhosez));

			/* ----- calculate range, azimuth and elevation rates ------- */
			drho = astMath::dot(rhosez, drhosez) / rho;
			if (fabs(temp * temp) > small)
				daz = (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) /
				(temp * temp);
			else
				daz = 0.0;

			if (fabs(temp) > 0.00000001)
				del = (drhosez[2] - drho * sin(el)) / temp;
			else
				del = 0.0;
		}
	}  // rv_razel


	/* -----------------------------------------------------------------------------
	*
	*                           function initjplde
	*
	*  this function initializes the jpl planetary ephemeris data. the input
	*  data files are from processing the ascii files into a text file of sun
	*  and moon positions.
	*
	*  author        : david vallado                  719-573-2600   22 jan 2018
	*
	*  revisions
	*
	*  inputs          description                    range / units
	*
	*
	*
	*
	*  outputs       :
	*    jpldearr    - array of jplde data records
	*    jdjpldestart- julian date of the start of the jpldearr data
	*
	*  locals        :
	*                -
	*
	*  coupling      :
	*
	*  references    :
	*
	*  -------------------------------------------------------------------------- */

	void initjplde
		(
		std::vector<jpldedata> &jpldearr,
		char infilename[140],
		double& jdjpldestart, double& jdjpldestartFrac
		)
	{
		jpldearr.resize(jpldesize);

		FILE *infile;
		double rs2, jdtdb, jdtdbf;
		char longstr[170];
		long i;

		// ---- open files select compatible files!!
#ifdef _MSC_VER
		errno_t jpldeInFileErrors = fopen_s(&infile, infilename, "r");
#else
		infile = fopen(infilename, "r");
#endif

		// ---- start finding of data
		// ---- find epoch date
		fgets(longstr, 170, infile);

		i = 0;
		// ---- first record is in the string above
#ifdef _MSC_VER
		sscanf_s(longstr, "%i %i %i %lf %lf %lf %lf %lf %lf %lf %lf ",
			&jpldearr[i].year, &jpldearr[i].mon, &jpldearr[i].day, &jpldearr[i].rsun[0],
			&jpldearr[i].rsun[1], &jpldearr[i].rsun[2], &jpldearr[i].rsmag, &rs2,
			&jpldearr[i].rmoon[0], &jpldearr[i].rmoon[1], &jpldearr[i].rmoon[2]);
#else
		sscanf(longstr, "%i %i %i %lf %lf %lf %lf %lf %lf %lf %lf ",
			&jpldearr[i].year, &jpldearr[i].mon, &jpldearr[i].day, &jpldearr[i].rsun[0],
			&jpldearr[i].rsun[1], &jpldearr[i].rsun[2], &jpldearr[i].rsmag, &rs2,
			&jpldearr[i].rmoon[0], &jpldearr[i].rmoon[1], &jpldearr[i].rmoon[2]);
#endif
		astTime::jday(jpldearr[i].year, jpldearr[i].mon, jpldearr[i].day, 0, 0, 0.0, jdjpldestart, jdjpldestartFrac);
		jpldearr[i].mjd = jdjpldestart + jdjpldestartFrac;

		// ---- process observed records
		for (i = 1; i <= 51830 - 1; i++)  // number of days 1 jan 1957 - 2100
		{
			// use d format for integers with leading 0's
#ifdef _MSC_VER
			fscanf_s(infile, "%i %i %i %lf %lf %lf %lf %lf %lf %lf %lf ",
				&jpldearr[i].year, &jpldearr[i].mon, &jpldearr[i].day, &jpldearr[i].rsun[0],
				&jpldearr[i].rsun[1], &jpldearr[i].rsun[2], &jpldearr[i].rsmag, &rs2,
				&jpldearr[i].rmoon[0], &jpldearr[i].rmoon[1], &jpldearr[i].rmoon[2]);
#else
			fscanf(infile, "%i %i %i %lf %lf %lf %lf %lf %lf %lf %lf ",
				&jpldearr[i].year, &jpldearr[i].mon, &jpldearr[i].day, &jpldearr[i].rsun[0],
				&jpldearr[i].rsun[1], &jpldearr[i].rsun[2], &jpldearr[i].rsmag, &rs2,
				&jpldearr[i].rmoon[0], &jpldearr[i].rmoon[1], &jpldearr[i].rmoon[2]);
#endif
			astTime::jday(jpldearr[i].year, jpldearr[i].mon, jpldearr[i].day, 0, 0, 0.0, jdtdb, jdtdbf);
			jpldearr[i].mjd = jdtdb + jdtdbf - 2400000.5;
		}

		fclose(infile);
	}   //  initjplde


	/* -----------------------------------------------------------------------------
	*
	*                           function findjpldeparam
	*
	*  this routine finds the jplde parameters for a given time. several types of
	*  interpolation are available. the cio and iau76 nutation parameters are also
	*  read for optimizing the speeed of calculations.
	*
	*  author        : david vallado                      719-573-2600   12 dec 2005
	*
	*  inputs          description                          range / units
	*    jdtdb         - epoch julian date              days from 4713 BC
	*    jdtdbF        - epoch julian date fraction     day fraction from jdutc
	*    interp        - interpolation                        n-none, l-linear, s-spline
	*    jpldearr      - array of jplde data records
	*    jdjpldestart  - julian date of the start of the jpldearr data (set in initjplde)
	*
	*  outputs       :
	*    dut1        - delta ut1 (ut1-utc)                  sec
	*    dat         - number of leap seconds               sec
	*    lod         - excess length of day                 sec
	*    xp          - x component of polar motion          rad
	*    yp          - y component of polar motion          rad
	*    ddpsi       - correction to delta psi (iau-76 theory) rad
	*    ddeps       - correction to delta eps (iau-76 theory) rad
	*    dx          - correction to x (cio theory)         rad
	*    dy          - correction to y (cio theory)         rad
	*    x           - x component of cio                   rad
	*    y           - y component of cio                   rad
	*    s           -                                      rad
	*    deltapsi    - nutation longitude angle             rad
	*    deltaeps    - obliquity of the ecliptic correction rad
	*
	*  locals        :
	*                -
	*
	*  coupling      :
	*    none        -
	*
	*  references    :
	*    vallado       2004,
	* --------------------------------------------------------------------------- */

	void findjpldeparam
		(
		double  jdtdb, double jdtdbF, char interp,
		const std::vector<jpldedata> &jpldearr,  // pass by reference, not modify
		double jdjpldestart,
		double rsun[3], double& rsmag,
		double rmoon[3]
		)
	{
		long recnum;
		int  off1, off2;
		jpldedata jplderec, nextjplderec;
		double  fixf, jdjpldestarto, mjd, jdb, mfme;

		// the ephemerides are centered on jdtdb, but it turns out to be 0.5, or 0000 hrs.
		// check if any whole days in jdF
		jdb = floor(jdtdb + jdtdbF) + 0.5;  // want jd at 0 hr
		mfme = (jdtdb + jdtdbF - jdb) * 1440.0;
		if (mfme < 0.0)
			mfme = 1440.0 + mfme;

		//printf("jdtdb %lf  %lf  %lf  %lf \n ", jdtdb, jdtdbF, jdb, mfme);

		// ---- read data for day of interest
		jdjpldestarto = floor(jdtdb + jdtdbF - jdjpldestart);
		recnum = int(jdjpldestarto);

		// check for out of bound values
		if ((recnum >= 1) && (recnum <= jpldesize))
		{
			jplderec = jpldearr[recnum];

			// ---- set non-interpolated values
			rsun[0] = jplderec.rsun[0];
			rsun[1] = jplderec.rsun[1];
			rsun[2] = jplderec.rsun[2];
			rsmag = jplderec.rsmag;
			mjd = jplderec.mjd;
			rmoon[0] = jplderec.rmoon[0];
			rmoon[1] = jplderec.rmoon[1];
			rmoon[2] = jplderec.rmoon[2];

			// ---- find nutation parameters for use in optimizing speed

			// ---- do linear interpolation
			if (interp == 'l')
			{
				nextjplderec = jpldearr[recnum + 1];
				fixf = mfme / 1440.0;

				rsun[0] = jplderec.rsun[0] + (nextjplderec.rsun[0] - jplderec.rsun[0]) * fixf;
				rsun[1] = jplderec.rsun[1] + (nextjplderec.rsun[1] - jplderec.rsun[1]) * fixf;
				rsun[2] = jplderec.rsun[2] + (nextjplderec.rsun[2] - jplderec.rsun[2]) * fixf;
				rsmag = jplderec.rsmag + (nextjplderec.rsmag - jplderec.rsmag) * fixf;
				rmoon[0] = jplderec.rmoon[0] + (nextjplderec.rmoon[0] - jplderec.rmoon[0]) * fixf;
				rmoon[1] = jplderec.rmoon[1] + (nextjplderec.rmoon[1] - jplderec.rmoon[1]) * fixf;
				rmoon[2] = jplderec.rmoon[2] + (nextjplderec.rmoon[2] - jplderec.rmoon[2]) * fixf;
				//printf("sunm %i rsmag %lf fixf %lf n %lf nxt %lf \n", recnum, rsmag, fixf, jplderec.rsun[0], nextjplderec.rsun[0]);
				//printf("recnum l %i fixf %lf %lf rsun %lf %lf %lf \n", recnum, fixf, jplderec.rsun[0], rsun[0], rsun[1], rsun[2]);
			}

			// ---- do spline interpolations
			if (interp == 's')
			{
				off2 = 2;   // every 5 days data...
				off1 = 1;
				mfme = mfme / 1440.0; // get back to days for this since each step is in days
				// setup so the interval is in between points 2 and 3
				rsun[0] = astMath::cubicinterp(jpldearr[recnum - off1].rsun[0], jpldearr[recnum].rsun[0], jpldearr[recnum + off1].rsun[0], jpldearr[recnum + off2].rsun[0],
					jpldearr[recnum - off1].mjd, jpldearr[recnum].mjd, jpldearr[recnum + off1].mjd, jpldearr[recnum + off2].mjd,
					jpldearr[recnum].mjd + mfme);
				rsun[1] = astMath::cubicinterp(jpldearr[recnum - off1].rsun[1], jpldearr[recnum].rsun[1], jpldearr[recnum + off1].rsun[1], jpldearr[recnum + off2].rsun[1],
					jpldearr[recnum - off1].mjd, jpldearr[recnum].mjd, jpldearr[recnum + off1].mjd, jpldearr[recnum + off2].mjd,
					jpldearr[recnum].mjd + mfme);
				rsun[2] = astMath::cubicinterp(jpldearr[recnum - off1].rsun[2], jpldearr[recnum].rsun[2], jpldearr[recnum + off1].rsun[2], jpldearr[recnum + off2].rsun[2],
					jpldearr[recnum - off1].mjd, jpldearr[recnum].mjd, jpldearr[recnum + off1].mjd, jpldearr[recnum + off2].mjd,
					jpldearr[recnum].mjd + mfme);
				rsmag = astMath::cubicinterp(jpldearr[recnum - off1].rsmag, jpldearr[recnum].rsmag, jpldearr[recnum + off1].rsmag, jpldearr[recnum + off2].rsmag,
					jpldearr[recnum - off1].mjd, jpldearr[recnum].mjd, jpldearr[recnum + off1].mjd, jpldearr[recnum + off2].mjd,
					jpldearr[recnum].mjd + mfme);
				rmoon[0] = int(astMath::cubicinterp(jpldearr[recnum - off1].rmoon[0], jpldearr[recnum].rmoon[0], jpldearr[recnum + off1].rmoon[0], jpldearr[recnum + off2].rmoon[0],
					jpldearr[recnum - off1].mjd, jpldearr[recnum].mjd, jpldearr[recnum + off1].mjd, jpldearr[recnum + off2].mjd,
					jpldearr[recnum].mjd + mfme));
				rmoon[1] = int(astMath::cubicinterp(jpldearr[recnum - off1].rmoon[1], jpldearr[recnum].rmoon[1], jpldearr[recnum + off1].rmoon[1], jpldearr[recnum + off2].rmoon[1],
					jpldearr[recnum - off1].mjd, jpldearr[recnum].mjd, jpldearr[recnum + off1].mjd, jpldearr[recnum + off2].mjd,
					jpldearr[recnum].mjd + mfme));
				rmoon[2] = int(astMath::cubicinterp(jpldearr[recnum - off1].rmoon[2], jpldearr[recnum].rmoon[2], jpldearr[recnum + off1].rmoon[2], jpldearr[recnum + off2].rmoon[2],
					jpldearr[recnum - off1].mjd, jpldearr[recnum].mjd, jpldearr[recnum + off1].mjd, jpldearr[recnum + off2].mjd,
					jpldearr[recnum].mjd + mfme));
				//printf("recnum s %i mfme %lf days rsun %lf %lf %lf \n", recnum, mfme, rsun[0], rsun[1], rsun[2]);
				//printf(" %lf %lf %lf %lf \n", jpldearr[recnum - off2].mjd, jpldearr[recnum - off1.mjd, jpldearr[recnum].mjd, jpldearr[recnum + off1].mjd);
			}
		}
		// set default values
		else
		{
			rsun[0] = 0.0;
			rsun[1] = 0.0;
			rsun[2] = 0.0;
			rsmag = 33.0;
			rmoon[0] = 0.0;
			rmoon[1] = 0.0;
			rmoon[2] = 0.0;
		}
	}  //  findjpldeparam



	/* ------------------------------------------------------------------------------
	*
	*                           function sun
	*
	*  this function calculates the geocentric equatorial position vector
	*    the sun given the julian date.  this is the low precision formula and
	*    is valid for years from 1950 to 2050.  accuaracy of apparent coordinates
	*    is 0.01  degrees.  notice many of the calculations are performed in
	*    degrees, and are not changed until later.  this is due to the fact that
	*    the almanac uses degrees exclusively in their formulations.
	*
	*  author        : david vallado                  719-573-2600   27 may 2002
	*
	*  revisions
	*    vallado     - fix mean lon of sun                            7 mat 2004
	*
	*  inputs          description                    range / units
	*    jdtdb         - epoch julian date              days from 4713 BC
	*    jdtdbF        - epoch julian date fraction     day fraction from jdutc
	*
	*  outputs       :
	*    rsun        - ijk position vector of the sun au
	*    rtasc       - right ascension                rad
	*    decl        - declination                    rad
	*
	*  locals        :
	*    meanlong    - mean longitude
	*    meananomaly - mean anomaly
	*    eclplong    - ecliptic longitude
	*    obliquity   - mean obliquity of the ecliptic
	*    tut1        - julian centuries of ut1 from
	*                  jan 1, 2000 12h
	*    ttdb        - julian centuries of tdb from
	*                  jan 1, 2000 12h
	*    hr          - hours                          0 .. 24              10
	*    min         - minutes                        0 .. 59              15
	*    sec         - seconds                        0.0  .. 59.99          30.00
	*    temp        - temporary variable
	*    deg         - degrees
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 279, alg 29, ex 5-1
	* --------------------------------------------------------------------------- */

	void sun
		(
		double jdtdb, double jdtdbF,
		double rsun[3], double& rtasc, double& decl
		)
	{
		double deg2rad;
		double tut1, meanlong, ttdb, meananomaly, eclplong, obliquity, magr;

		deg2rad = pi / 180.0;

		// -------------------------  implementation   -----------------
		// -------------------  initialize values   --------------------
		tut1 = (jdtdb + jdtdbF - 2451545.0) / 36525.0;

		meanlong = 280.460 + 36000.77 * tut1;
		meanlong = fmod(meanlong, 360.0);  //deg

		ttdb = tut1;
		meananomaly = 357.5277233 + 35999.05034 * ttdb;
		meananomaly = fmod(meananomaly * deg2rad, twopi);  //rad
		if (meananomaly < 0.0)
		{
			meananomaly = twopi + meananomaly;
		}
		eclplong = meanlong + 1.914666471 * sin(meananomaly)
			+ 0.019994643 * sin(2.0 * meananomaly); //deg
		obliquity = 23.439291 - 0.0130042 * ttdb;  //deg
		meanlong = meanlong * deg2rad;
		if (meanlong < 0.0)
		{
			meanlong = twopi + meanlong;
		}
		eclplong = eclplong * deg2rad;
		obliquity = obliquity * deg2rad;

		// --------- find magnitude of sun vector, ) components ------
		magr = 1.000140612 - 0.016708617 * cos(meananomaly)
			- 0.000139589 * cos(2.0 *meananomaly);    // in au's

		rsun[0] = magr * cos(eclplong);
		rsun[1] = magr * cos(obliquity) * sin(eclplong);
		rsun[2] = magr * sin(obliquity) * sin(eclplong);

		rtasc = atan(cos(obliquity) * tan(eclplong));

		// --- check that rtasc is in the same quadrant as eclplong ----
		if (eclplong < 0.0)
		{
			eclplong = eclplong + twopi;    // make sure it's in 0 to 2pi range
		}
		if (fabs(eclplong - rtasc) > pi * 0.5)
		{
			rtasc = rtasc + 0.5 * pi * astMath::round((eclplong - rtasc) / (0.5 * pi));
		}
		decl = asin(sin(obliquity) * sin(eclplong));

	}  // sun


	/* ------------------------------------------------------------------------------
	*
	*                           function sunmoonjpl
	*
	*  this function calculates the geocentric equatorial position vector
	*    the sun given the julian date. these are the jpl de ephemerides.
	*
	*  author        : david vallado                  719-573-2600   27 may 2002
	*
	*  revisions
	*
	*  inputs          description                    range / units
	*    jdtdb         - epoch julian date              days from 4713 BC
	*    jdtdbF        - epoch julian date fraction     day fraction from jdutc
	*    interp        - interpolation                        n-none, l-linear, s-spline
	*    jpldearr      - array of jplde data records
	*    jdjpldestart  - julian date of the start of the jpldearr data (set in initjplde)
	*
	*  outputs       :
	*    rsun        - ijk position vector of the sun au
	*    rtasc       - right ascension                rad
	*    decl        - declination                    rad
	*
	*  locals        :
	*    meanlong    - mean longitude
	*    meananomaly - mean anomaly
	*    eclplong    - ecliptic longitude
	*    obliquity   - mean obliquity of the ecliptic
	*    tut1        - julian centuries of ut1 from
	*                  jan 1, 2000 12h
	*    ttdb        - julian centuries of tdb from
	*                  jan 1, 2000 12h
	*    hr          - hours                          0 .. 24              10
	*    min         - minutes                        0 .. 59              15
	*    sec         - seconds                        0.0  .. 59.99          30.00
	*    temp        - temporary variable
	*    deg         - degrees
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 279, alg 29, ex 5-1
	* --------------------------------------------------------------------------- */

	void sunmoonjpl
		(
		double jdtdb, double jdtdbF,
		char interp,
		const std::vector<jpldedata> &jpldearr,
		double jdjpldestart,
		double rsun[3], double& rtascs, double& decls,
		double rmoon[3], double& rtascm, double& declm
		)
	{
		double deg2rad;
		double rsmag, rr, temp;
		double small = 0.000001;

		deg2rad = pi / 180.0;

		// -------------------------  implementation   -----------------
		// -------------------  initialize values   --------------------
		ast2Body::findjpldeparam(jdtdb, jdtdbF, interp, jpldearr, jdjpldestart, rsun, rsmag, rmoon);

		rr = astMath::mag(rsun);
		temp = sqrt(rsun[0] * rsun[0] + rsun[1] * rsun[1]);
		if (temp < small)
			// rtascs = atan2(v[1], v[0]);
			rtascs = 0.0;
		else
			rtascs = atan2(rsun[1], rsun[0]);
		decls = asin(rsun[2] / rr);

		rr = astMath::mag(rmoon);
		temp = sqrt(rmoon[0] * rmoon[0] + rmoon[1] * rmoon[1]);
		if (temp < small)
			// rtascm = atan2(v[1], v[0]);
			rtascm = 0.0;
		else
			rtascm = atan2(rmoon[1], rmoon[0]);
		declm = asin(rmoon[2] / rr);


		//// evaluate DE series
		//const char *ephfile_name = "D:/Dataorig/Planet/DE200/jpleph";  // 1940-2100 or so
		//      #define JPL_MAX_N_CONSTANTS 1018
		//char nams[JPL_MAX_N_CONSTANTS][6], buff[102];
		//double vals[JPL_MAX_N_CONSTANTS];
		//void *ephem = jpl_init_ephemeris(ephfile_name, nams, vals);
		//double jdarr;
		//jdarr = 2453238.5;
		//// jdarr = 0.349572634765;
		//int calcopt[13];
		//int ntarg = 3; // earth
		//int ncent = 11;  // sun 10 = moon
		////   jpl_state(jdarr, calcopt, ntarg, ncent, rrd[], calc_velocity);
		//int calc_velocity = 1;  // true, 0 if not to calc vel
		//double rrd[6];
		//double au = 149597870.0;  // au 2 km
		//int err_code = jpl_pleph(ephem, jdarr, ntarg, ncent, rrd, calc_velocity);
		//printf("%i sun - earth %d %d %d %d %d %d \n", err_code, rrd[0] * au, rrd[1] * au, rrd[2] * au, rrd[3], rrd[4], rrd[5]);

		//err_code = jpl_pleph(ephem, jdarr, 5, 12, rrd, calc_velocity);
		//printf("%i sol sys bary - jup %d %d %d %d %d %d \n", err_code, rrd[0] * au, rrd[1] * au, rrd[2] * au, rrd[3], rrd[4], rrd[5]);

		//err_code = jpl_pleph(ephem, jdarr, 13, 11, rrd, calc_velocity);
		//printf("%i sun - earthm bary %d %d %d %d %d %d \n", err_code, rrd[0] * au, rrd[1] * au, rrd[2] * au, rrd[3], rrd[4], rrd[5]);

		//err_code = jpl_pleph(ephem, jdarr, 10, 3, rrd, calc_velocity);
		//printf("%i earth - moon %d %d %d %d %d %d \n", err_code, rrd[0] * au, rrd[1] * au, rrd[2] * au, rrd[3], rrd[4], rrd[5]);
	}  // sunmoonjpl


	/* -----------------------------------------------------------------------------
	*
	*                           function moon
	*
	*  this function calculates the geocentric equatorial (ijk) position vector
	*    for the moon given the julian date.
	*
	*  author        : david vallado                  719-573-2600   27 may 2002
	*
	*  revisions
	*                -
	*
	*  inputs          description                    range / units
	*    jdtdb         - epoch julian date              days from 4713 BC
	*    jdtdbF        - epoch julian date fraction     day fraction from jdutc
	*
	*  outputs       :
	*    rmoon       - ijk position vector of moon    km
	*    rtasc       - right ascension                rad
	*    decl        - declination                    rad
	*
	*  locals        :
	*    eclplong    - ecliptic longitude
	*    eclplat     - eclpitic latitude
	*    hzparal     - horizontal parallax
	*    l           - geocentric direction cosines
	*    m           -             "     "
	*    n           -             "     "
	*    ttdb        - julian centuries of tdb from
	*                  jan 1, 2000 12h
	*    hr          - hours                          0 .. 24
	*    min         - minutes                        0 .. 59
	*    sec         - seconds                        0.0  .. 59.99
	*    deg         - degrees
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 288, alg 31, ex 5-3
	* --------------------------------------------------------------------------- */

	void moon
		(
		double jdtdb, double jdtdbF,
		double rmoon[3], double& rtasc, double& decl
		)
	{
		double deg2rad, magr;
		double ttdb, l, m, n, eclplong, eclplat, hzparal, obliquity;
		double re = 6378.137;

		deg2rad = pi / 180.0;

		// -------------------------  implementation   -----------------
		ttdb = (jdtdb + jdtdbF - 2451545.0) / 36525.0;

		eclplong = 218.32 + 481267.8813 * ttdb
			+ 6.29 * sin((134.9 + 477198.85 * ttdb) * deg2rad)
			- 1.27 * sin((259.2 - 413335.38 * ttdb) * deg2rad)
			+ 0.66 * sin((235.7 + 890534.23 * ttdb) * deg2rad)
			+ 0.21 * sin((269.9 + 954397.70 * ttdb) * deg2rad)
			- 0.19 * sin((357.5 + 35999.05 * ttdb) * deg2rad)
			- 0.11 * sin((186.6 + 966404.05 * ttdb) * deg2rad);      // deg

		eclplat = 5.13 * sin((93.3 + 483202.03 * ttdb) * deg2rad)
			+ 0.28 * sin((228.2 + 960400.87 * ttdb) * deg2rad)
			- 0.28 * sin((318.3 + 6003.18 * ttdb) * deg2rad)
			- 0.17 * sin((217.6 - 407332.20 * ttdb) * deg2rad);      // deg

		hzparal = 0.9508 + 0.0518 * cos((134.9 + 477198.85 * ttdb)
			* deg2rad)
			+ 0.0095 * cos((259.2 - 413335.38 * ttdb) * deg2rad)
			+ 0.0078 * cos((235.7 + 890534.23 * ttdb) * deg2rad)
			+ 0.0028 * cos((269.9 + 954397.70 * ttdb) * deg2rad);    // deg

		eclplong = fmod(eclplong * deg2rad, twopi);
		eclplat = fmod(eclplat * deg2rad, twopi);
		hzparal = fmod(hzparal * deg2rad, twopi);

		obliquity = 23.439291 - 0.0130042 * ttdb;  //deg
		obliquity = obliquity * deg2rad;

		// ------------ find the geocentric direction cosines ----------
		l = cos(eclplat) * cos(eclplong);
		m = cos(obliquity) * cos(eclplat) * sin(eclplong) - sin(obliquity) * sin(eclplat);
		n = sin(obliquity) * cos(eclplat) * sin(eclplong) + cos(obliquity) * sin(eclplat);

		// ------------- calculate moon position vector ----------------
		magr = re / sin(hzparal);  // km
		rmoon[0] = magr * l;
		rmoon[1] = magr * m;
		rmoon[2] = magr * n;

		// -------------- find rt ascension and declination ------------
		rtasc = atan2(m, l);
		decl = asin(n);
	}  // moon


}   // namespace
