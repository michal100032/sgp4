#pragma once

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "sgp4/sgp4.hpp"

namespace SGP4Funcs_mod
{
	typedef struct elsetrec
	{
	  char      satnum[6];
	  int       epochyr = 0, epochtynumrev = 0;
	  int       error = 0;
	  char      operationmode = 0;
	  char      init = 0, method = 0;

	  /* Near Earth */
	  int    isimp;
	  double aycof  = 0.0, con41  = 0.0, cc1     = 0.0, cc4    = 0.0, cc5    = 0.0, d2      = 0.0, d3    = 0.0, d4    = 0.0,
			 delmo  = 0.0, eta    = 0.0, argpdot = 0.0, omgcof = 0.0, sinmao = 0.0, t       = 0.0, t2cof = 0.0, t3cof = 0.0,
			 t4cof  = 0.0, t5cof  = 0.0, x1mth2  = 0.0, x7thm1 = 0.0, mdot   = 0.0, nodedot = 0.0, xlcof = 0.0, xmcof = 0.0,
			 nodecf = 0;

	  /* Deep Space */
	  int    irez = 0;
	  double d2201 = 0.0, d2211 = 0.0, d3210 = 0.0, d3222 = 0.0, d4410 = 0.0, d4422 = 0.0, d5220 = 0.0, d5232 = 0.0,
			 d5421 = 0.0, d5433 = 0.0, dedt  = 0.0, del1  = 0.0, del2  = 0.0, del3  = 0.0, didt  = 0.0, dmdt  = 0.0,
			 dnodt = 0.0, domdt = 0.0, e3    = 0.0, ee2   = 0.0, peo   = 0.0, pgho  = 0.0, pho   = 0.0, pinco = 0.0,
			 plo   = 0.0, se2   = 0.0, se3   = 0.0, sgh2  = 0.0, sgh3  = 0.0, sgh4  = 0.0, sh2   = 0.0, sh3   = 0.0,
			 si2   = 0.0, si3   = 0.0, sl2   = 0.0, sl3   = 0.0, sl4   = 0.0, gsto  = 0.0, xfact = 0.0, xgh2  = 0.0,
			 xgh3  = 0.0, xgh4  = 0.0, xh2   = 0.0, xh3   = 0.0, xi2   = 0.0, xi3   = 0.0, xl2   = 0.0, xl3   = 0.0,
			 xl4   = 0.0, xlamo = 0.0, zmol  = 0.0, zmos  = 0.0, atime = 0.0, xli   = 0.0, xni = 0.0;
	  //
	  // from TLE: 
	  // epochdays - day fraction (exactly as in tle),
	  // epochyr - last two digits of the year (tle)
	  // 
	  // 
	  // 
	  //  mo - mean anomaly [rad]
	  //  no_kozai - mean motion  [rad/min]
	  //  ndot - first derivative of mean motion  [rad/min^2]
	  //  nddot - second derivative of mean motion [rad/min^3]
	  //  
	  // jdsatepoch - julian date of epoch (tle type), midnight 0.0 
	  // jdsatepochF - day fraction
	  // inclo - inclination [rad]
	  // nodeo - right ascension of the ascending node [rad]
	  // argpo - argument of perigee [rad]
	  // ecco - eccentricity 


	  double a = 0.0, altp = 0.0, alta = 0.0, epochdays = 0.0, 
		  jdsatepoch = 0.0, jdsatepochF = 0.0, nddot = 0.0, ndot = 0.0,
		  bstar = 0.0, rcse = 0.0, inclo = 0.0, nodeo = 0.0, ecco = 0.0,
		  argpo = 0.0, mo = 0.0, no_kozai = 0.0;
	  // sgp4fix add new variables from tle
	  char  classification = 0, intldesg[11];
	  int   ephtype = 0;
	  long  elnum = 0, revnum = 0;
	  // sgp4fix add unkozai'd variable
	  double no_unkozai = 0.0;
	  // sgp4fix add singly averaged variables
	  double am = 0.0, em = 0.0, im = 0.0, Om = 0.0, om = 0.0, mm = 0.0, nm = 0.0;
	  // sgp4fix add constant parameters to eliminate mutliple calls during execution
	  
	  //       Additional elements to capture relevant TLE and object information:       
	  long dia_mm = 0; // RSO dia in mm
	  double period_sec = 0; // Period in seconds
	  unsigned char active = 0; // "Active S/C" flag (0=n, 1=y) 
	  unsigned char not_orbital = 0; // "Orbiting S/C" flag (0=n, 1=y)  
	  double rcs_m2 = 0; // "RCS (m^2)" storage  
	} elsetrec;
	//	public class SGP4Class
	//	{

	bool sgp4init
		(
		char opsmode, const char satn[9], const double epoch,
		const double xbstar, const double xndot, const double xnddot, const double xecco, const double xargpo,
		const double xinclo, const double xmo, const double xno,
		const double xnodeo, elsetrec& satrec
		);

	bool sgp4(elsetrec& satrec, double tsince, double r[3], double v[3]);

	elsetrec twoline2rv(const sgp4::tle_set& set, char opsmode);

	// older sgp4ext methods
	double  gstime_SGP4
		(
		double jdut1
		);

	void    jday_SGP4
		(
		int year, int mon, int day, int hr, int minute, double sec,
		double& jd, double& jdFrac
		);

	void    days2mdhms_SGP4
		(
		int year, double days,
		int& mon, int& day, int& hr, int& minute, double& sec
		);

	void    invjday_SGP4
		(
		double jd, double jdFrac,
		int& year, int& mon, int& day,
		int& hr, int& minute, double& sec
		);


}  // namespace