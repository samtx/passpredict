#ifndef _ast2Body_h_
#define _ast2Body_h_
/* --------------------------------------------------------------------
*
*                                ast2body.h
*
*    this file contains miscallaneous two-body motion functions.
*
*    current :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*               4 may 09  david vallado
*                           misc updates
*              23 feb 07  david vallado
*                           3rd edition baseline
*               6 dec 05  david vallado
*                           add ijk2ll
*              20 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
  ----------------------------------------------------------------------      */

#include <vector>
#include <stdint.h>

// be sure to update to your specific paths!!
// " " tells the compiler to look in this directory first, usually the parent directory
// you can leave generic as astMath.h, but you have to set the include directory in the property pages
#include "astMath.h"  // pi, infinite, undefined
#include "astTime.h"  // pi, twopi, edirection
// #include "D:/Dataorig/Planet/jpl_eph-master/jpleph.h"  // several
// #include "D:/Dataorig/Planet/jpl_eph-master/jpl_int.h"  // several

#pragma once



const int jpldesize = 60000; // 60000 if from 1957-2100

typedef struct jpldedata
{
	double rsun[3], rmoon[3];
	int    year, mon, day;
	double rsmag, mjd;
} jpldedata;


namespace ast2Body

{
	void rv2coe
		(
		double r[3], double v[3], const double mu,
		double& p, double& a, double& ecc, double& incl, double& raan, double& argp,
		double& nu, double& m, double& eccanom, double& arglat, double& truelon, double& lonper
		);

	void coe2rv
		(
		double p, double ecc, double incl, double raan, double argp, double nu,
		double arglat, double truelon, double lonper,
		double r[3], double v[3]
		);

	void rv2eq
		(
		double r[3], double v[3],
		double& a, double& n, double& af, double& ag, double& chi, double& psi, double& meanlonM, double& meanlonNu, int& fr
		);

	void eq2rv
		(
		double a, double af, double ag, double chi, double psi, double meanlon, int fr,
		double r[3], double v[3]
		);

	void findc2c3
		(
		double znew,
		double& c2new, double& c3new
		);

	void kepler
		(
		double ro[3], double vo[3], double dtseco, double r[3], double v[3]
		);

	void pkepler
		(
		double r1[3], double v1[3], double& dtsec, double& ndot, double& nddot, double r2[3], double v2[3]
		);

    void rv2rsw
		(
		double r[3], double v[3],
		double rrsw[3], double vrsw[3], std::vector< std::vector<double> > &transmat
		);

	void rv2ntw
		(
		double r[3], double v[3],
		double rntw[3], double vntw[3], std::vector< std::vector<double> > &transmat
		);

	void newtonm
		(
		double ecc, double m, double& e0, double& nu
		);

	void newtonnu
		(
		double ecc, double nu,
		double& e0, double& m
		);

	void gc_gd
		(
		double&    latgc,
		edirection direct,
		double&    latgd
		);

	void ijk2ll
		(
		double recef[3], double jdut1,
		double& latgc, double& latgd, double& lon, double& hellp
		);

	void rvsez_razel
		(
		double rhosez[3], double drhosez[3],
		edirection direct,
		double& rho, double& az, double& el, double& drho, double& daz, double& del
		);

	void rv2radec
		(
		double r[3], double v[3],
		double& rr, double& rtasc, double& decl, double& drr, double& drtasc, double& ddecl
		);

	void rv_razel
		(
		double recef[3], double vecef[3], double rsecef[3], double latgd, double lon,
		edirection direct,
		double& rho, double& az, double& el, double& drho, double& daz, double& del
		);

	void initjplde
		(
		std::vector<jpldedata> &jpldearr,
		char infilename[140],
		double& jdjpldestart, double& jdjpldestartFrac
		);

	void findjpldeparam
		(
		double  jdtdb, double jdtdbF, char interp,
		const std::vector<jpldedata> &jpldearr,
		double jdjpldestart,
		double rsun[3], double& rsmag,
		double rmoon[3]
		);

	void sun
		(
		double jdtdb, double jdtdbF,
		double rsun[3], double& rtasc, double& decl
		);

	void sunmoonjpl
		(
		double jdtdb, double jdtdbF,
		char interp,
		const std::vector<jpldedata> &jpldearr,
		double jdjpldestart,
		double rsun[3], double& rtascs, double& decls,
		double rmoon[3], double& rtascm, double& declm
		);

	void moon
		(
		double jdtdb, double jdtdbF,
		double rmoon[3], double& rtasc, double& decl
		);


	//	};  // Class

};  // namespace

#endif