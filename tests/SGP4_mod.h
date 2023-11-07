#pragma once

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "sgp4/sgp4.hpp"

namespace SGP4Funcs_mod
{

	// -------------------------- structure declarations ----------------------------
	typedef enum
	{
	  wgs72old,
	  wgs72,
	  wgs84
	} gravconsttype;

	typedef struct elsetrec
	{
	  char      satnum[6];
	  int       epochyr, epochtynumrev;
	  int       error;
	  char      operationmode;
	  char      init, method;

	  /* Near Earth */
	  int    isimp;
	  double aycof  , con41  , cc1    , cc4      , cc5    , d2      , d3   , d4    ,
			 delmo  , eta    , argpdot, omgcof   , sinmao , t       , t2cof, t3cof ,
			 t4cof  , t5cof  , x1mth2 , x7thm1   , mdot   , nodedot, xlcof , xmcof ,
			 nodecf;

	  /* Deep Space */
	  int    irez;
	  double d2201  , d2211  , d3210  , d3222    , d4410  , d4422   , d5220 , d5232 ,
			 d5421  , d5433  , dedt   , del1     , del2   , del3    , didt  , dmdt  ,
			 dnodt  , domdt  , e3     , ee2      , peo    , pgho    , pho   , pinco ,
			 plo    , se2    , se3    , sgh2     , sgh3   , sgh4    , sh2   , sh3   ,
			 si2    , si3    , sl2    , sl3      , sl4    , gsto    , xfact , xgh2  ,
			 xgh3   , xgh4   , xh2    , xh3      , xi2    , xi3     , xl2   , xl3   ,
			 xl4    , xlamo  , zmol   , zmos     , atime  , xli     , xni;
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


	  double a, altp, alta, epochdays, jdsatepoch, jdsatepochF, nddot, ndot,
			 bstar, rcse, inclo, nodeo, ecco, argpo, mo, no_kozai;
	  // sgp4fix add new variables from tle
	  char  classification, intldesg[11];
	  int   ephtype;
	  long  elnum    , revnum;
	  // sgp4fix add unkozai'd variable
	  double no_unkozai;
	  // sgp4fix add singly averaged variables
	  double am     , em     , im     , Om       , om     , mm      , nm;
	  // sgp4fix add constant parameters to eliminate mutliple calls during execution
	  double tumin, mus, radiusearthkm, xke, j2, j3, j4, j3oj2;

	  //       Additional elements to capture relevant TLE and object information:       
	  long dia_mm; // RSO dia in mm
	  double period_sec; // Period in seconds
	  unsigned char active; // "Active S/C" flag (0=n, 1=y) 
	  unsigned char not_orbital; // "Orbiting S/C" flag (0=n, 1=y)  
	  double rcs_m2; // "RCS (m^2)" storage  

	} elsetrec;
	//	public class SGP4Class
	//	{

	bool sgp4init
		(
		gravconsttype whichconst, char opsmode, const char satn[9], const double epoch,
		const double xbstar, const double xndot, const double xnddot, const double xecco, const double xargpo,
		const double xinclo, const double xmo, const double xno,
		const double xnodeo, elsetrec& satrec
		);

	bool sgp4
		(
		// no longer need gravconsttype whichconst, all data contained in satrec
		elsetrec& satrec, double tsince,
		double r[3], double v[3]
		);

	void getgravconst
		(
		gravconsttype whichconst,
		double& tumin,
		double& mus,
		double& radiusearthkm,
		double& xke,
		double& j2,
		double& j3,
		double& j4,
		double& j3oj2
		);

	// older sgp4io methods
	// removed typerun, typeinput, startmfe, stopmfe, deltamin
	void twoline2rv
	(
		const sgp4::tle_set& set,
		char opsmode,
		gravconsttype whichconst,
		elsetrec& satrec
	);

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