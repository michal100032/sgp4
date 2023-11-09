#include "SGP4_mod_near.h"

#include <iostream>
#include <iomanip>

static const double PI = 3.14159265358979323846;
static const double TWO_PI = 2 * PI;
static const double DEG_TO_RAD = PI / 180.0;

// gravitational parameter of the Earth [km^3/s^2]
static const double EARTH_GRV = 398600.5;

// radius of the Earth [km]
static const double EARTH_RAD = 6378.137;

// ??????
// ???? angular speed in an orbit of radius = EARTH_RAD? [rad/min]
static const double X_KE = 60.0 / sqrt(EARTH_RAD * EARTH_RAD * EARTH_RAD / EARTH_GRV);

static const double TUMIN = 1.0 / X_KE;

// earth zonal harmonic model
static const double EARTH_J2 = 0.00108262998905;
static const double EARTH_J3 = -0.00000253215306;
static const double EARTH_J4 = -0.00000161098761;
static const double EARTH_J3_TO_J2 = EARTH_J3 / EARTH_J2;

namespace SGP4Funcs_mod_near
{
	/*-----------------------------------------------------------------------------
	*
	*                           procedure initl
	*
	*  this procedure initializes the spg4 propagator. all the initialization is
	*    consolidated here instead of having multiple loops inside other routines.
	*
	*  author        : david vallado                  719-573-2600   28 jun 2005
	*
	*  inputs        :
	*    satn        - satellite number - not needed, placed in sat
	*    xke         - reciprocal of tumin
	*    j2          - j2 zonal harmonic
	*    ecco        - eccentricity                           0.0 - 1.0
	*    epoch       - epoch time in days from jan 0, 1950. 0 hr
	*    inclo       - inclination of satellite
	*    no          - mean motion of satellite
	*
	*  outputs       :
	*    ainv        - 1.0 / ao
	*    ao          - semi major axis
	*    con41       -
	*    con42       - 1.0 - 5.0 cos(i)
	*    cosio       - cosine of inclination
	*    cosio2      - cosio squared
	*    eccsq       - eccentricity squared
	*    method      - flag for deep space                    'd', 'n'
	*    omeosq      - 1.0 - ecco * ecco
	*    posq        - semi-parameter squared
	*    rp          - radius of perigee
	*    rteosq      - square root of (1.0 - ecco*ecco)
	*    sinio       - sine of inclination
	*    gsto        - gst at time of observation               rad
	*    no          - mean motion of satellite
	*
	*  locals        :
	*    ak          -
	*    d1          -
	*    del         -
	*    adel        -
	*    po          -
	*
	*  coupling      :
	*    getgravconst- no longer used
	*    gstime      - find greenwich sidereal time from the julian date
	*
	*  references    :
	*    hoots, roehrich, norad spacetrack report #3 1980
	*    hoots, norad spacetrack report #6 1986
	*    hoots, schumacher and glover 2004
	*    vallado, crawford, hujsak, kelso  2006
	----------------------------------------------------------------------------*/

	static void initl
	(
		double ecco, double epoch, double inclo, double no_kozai,
		char& method, double& con41, double& con42, double& cosio,
		double& cosio2, double& eccsq, double& omeosq, double& posq,
		double& rp, double& rteosq, double& sinio, double& no_unkozai
	)
	{
		/* --------------------- local variables ------------------------ */
		double ak, d1, del, adel, po;
	
		/* ------------- calculate auxillary epoch quantities ---------- */
		eccsq = ecco * ecco;
		omeosq = 1.0 - eccsq;
		rteosq = sqrt(omeosq);
		cosio = cos(inclo);
		cosio2 = cosio * cosio;

		/* ------------------ un-kozai the mean motion ----------------- */

		ak = pow(X_KE / no_kozai, 2.0 / 3.0);
		d1 = 0.75 * EARTH_J2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
		del = d1 / (ak * ak);
		adel = ak * (1.0 - del * del - del *
			(1.0 / 3.0 + 134.0 * del * del / 81.0));
		del = d1 / (adel * adel);
		no_unkozai = no_kozai / (1.0 + del);

		// calcualting sma for the seconds time
		double sma = pow(X_KE / (no_unkozai), 2.0 / 3.0);
		sinio = sin(inclo);
		po = sma * omeosq;
		con42 = 1.0 - 5.0 * cosio2;
		con41 = -con42 - cosio2 - cosio2;
		posq = po * po;
		rp = sma * (1.0 - ecco);
		method = 'n';
	}

	/*-----------------------------------------------------------------------------
	*
	*                             procedure sgp4init
	*
	*  this procedure initializes variables for sgp4.
	*
	*  author        : david vallado                  719-573-2600   28 jun 2005
	*
	*  inputs        :
	*    opsmode     - mode of operation afspc or improved 'a', 'i'
	*    whichconst  - which set of constants to use  72, 84
	*    satn        - satellite number
	*    bstar       - sgp4 type drag coefficient              kg/m2er
	*    ecco        - eccentricity
	*    epoch       - epoch time in days from jan 0, 1950. 0 hr
	*    argpo       - argument of perigee (output if ds)
	*    inclo       - inclination
	*    mo          - mean anomaly (output if ds)
	*    no          - mean motion
	*    nodeo       - right ascension of ascending node
	*
	*  outputs       :
	*    sat      - common values for subsequent calls
	*    return code - non-zero on error.
	*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
	*                   2 - mean motion less than 0.0
	*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
	*                   4 - semi-latus rectum < 0.0
	*                   5 - epoch elements are sub-orbital
	*                   6 - satellite has decayed
	*
	*  locals        :
	*    cnodm  , snodm  , cosim  , sinim  , cosomm , sinomm
	*    cc1sq  , cc2    , cc3
	*    coef   , coef1
	*    cosio4      -
	*    day         -
	*    dndt        -
	*    em          - eccentricity
	*    emsq        - eccentricity squared
	*    eeta        -
	*    etasq       -
	*    gam         -
	*    argpm       - argument of perigee
	*    nodem       -
	*    inclm       - inclination
	*    mm          - mean anomaly
	*    nm          - mean motion
	*    perige      - perigee
	*    pinvsq      -
	*    psisq       -
	*    qzms24      -
	*    rtemsq      -
	*    s1, s2, s3, s4, s5, s6, s7          -
	*    sfour       -
	*    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
	*    sz1, sz2, sz3
	*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
	*    tc          -
	*    temp        -
	*    temp1, temp2, temp3       -
	*    tsi         -
	*    xpidot      -
	*    xhdot1      -
	*    z1, z2, z3          -
	*    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
	*
	*  coupling      :
	*    getgravconst-
	*    initl       -
	*    dscom       -
	*    dpper       -
	*    dsinit      -
	*    sgp4        -
	*
	*  references    :
	*    hoots, roehrich, norad spacetrack report #3 1980
	*    hoots, norad spacetrack report #6 1986
	*    hoots, schumacher and glover 2004
	*    vallado, crawford, hujsak, kelso  2006
	----------------------------------------------------------------------------*/
	bool sgp4init(elsetrec& sat)
	{
		double epoch_1950 = sat.jdsatepoch + sat.jdsatepochF - 2433281.5;
		/* --------------------- local variables ------------------------ */
		double con42, cosio, sinio, cosio2, eccsq,
			omeosq, posq, rp, rteosq,
			cnodm, snodm, cosim, sinim, cosomm, sinomm, cc1sq,
			cc2, cc3, coef, coef1, cosio4, day, dndt,
			em, emsq, eeta, etasq, gam, argpm, nodem,
			inclm, mm, nm, perige, pinvsq, psisq, qzms24,
			rtemsq, s1, s2, s3, s4, s5, s6,
			s7, sfour, ss1, ss2, ss3, ss4, ss5,
			ss6, ss7, sz1, sz2, sz3, sz11, sz12,
			sz13, sz21, sz22, sz23, sz31, sz32, sz33,
			tc, temp, temp1, temp2, temp3, tsi, xpidot,
			xhdot1, z1, z2, z3, z11, z12, z13,
			z21, z22, z23, z31, z32, z33,
			qzms2t, ss, r[3], v[3],
			delmotemp, qzms2ttemp, qzms24temp;

		/* ------------------------ initialization --------------------- */
		// sgp4fix divisor for divide by zero check on inclination
		// the old check used 1.0 + cos(PI-1.0e-9), but then compared it to
		// 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
		const double temp4 = 1.5e-12;

		// single averaged mean elements
		sat.am = sat.em = sat.im = sat.Om = sat.mm = sat.nm = 0.0;

		/* ------------------------ earth constants ----------------------- */
		// sgp4fix identify constants and allow alternate values no longer needed
		// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
		ss = 78.0 / EARTH_RAD + 1.0;

		// sgp4fix use multiply for speed instead of pow
		qzms2ttemp = (120.0 - 78.0) / EARTH_RAD;
		qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;

		sat.init = 'y';
		sat.t = 0.0;

		initl
		(sat.ecco, epoch_1950, sat.inclo, sat.no_kozai,
			sat.method, sat.con41, con42, cosio, cosio2, eccsq, omeosq,
			posq, rp, rteosq, sinio, sat.no_unkozai);

		// a, alta and altp are in the unitis of earths radii
		sat.a = pow(sat.no_unkozai * TUMIN, (-2.0 / 3.0));
		sat.alta = sat.a * (1.0 + sat.ecco) - 1.0;
		sat.altp = sat.a * (1.0 - sat.ecco) - 1.0;

		std::cout << sat.a << std::endl;

		if ((omeosq >= 0.0) || (sat.no_unkozai >= 0.0))
		{
			sat.isimp = 0;
			if (rp < (220.0 / EARTH_RAD + 1.0))
				sat.isimp = 1;
			sfour = ss;
			qzms24 = qzms2t;
			perige = (rp - 1.0) * EARTH_RAD;

			/* - for perigees below 156 km, s and qoms2t are altered - */
			if (perige < 156.0)
			{
				sfour = perige - 78.0;
				if (perige < 98.0)
					sfour = 20.0;
				// sgp4fix use multiply for speed instead of pow
				qzms24temp = (120.0 - sfour) / EARTH_RAD;
				qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
				sfour = sfour / EARTH_RAD + 1.0;
			}
			pinvsq = 1.0 / posq;

			tsi = 1.0 / (sat.a - sfour);
			sat.eta = sat.a * sat.ecco * tsi;
			etasq = sat.eta * sat.eta;
			eeta = sat.ecco * sat.eta;
			psisq = fabs(1.0 - etasq);
			coef = qzms24 * pow(tsi, 4.0);
			coef1 = coef / pow(psisq, 3.5);
			cc2 = coef1 * sat.no_unkozai * (sat.a * (1.0 + 1.5 * etasq + eeta *
				(4.0 + etasq)) + 0.375 * EARTH_J2 * tsi / psisq * sat.con41 *
				(8.0 + 3.0 * etasq * (8.0 + etasq)));
			sat.cc1 = sat.bstar * cc2;
			cc3 = 0.0;
			if (sat.ecco > 1.0e-4)
				cc3 = -2.0 * coef * tsi * EARTH_J3_TO_J2 * sat.no_unkozai * sinio / sat.ecco;
			sat.x1mth2 = 1.0 - cosio2;
			sat.cc4 = 2.0 * sat.no_unkozai * coef1 * sat.a * omeosq *
				(sat.eta * (2.0 + 0.5 * etasq) + sat.ecco *
					(0.5 + 2.0 * etasq) - EARTH_J2 * tsi / (sat.a * psisq) *
					(-3.0 * sat.con41 * (1.0 - 2.0 * eeta + etasq *
						(1.5 - 0.5 * eeta)) + 0.75 * sat.x1mth2 *
						(2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * sat.argpo)));
			sat.cc5 = 2.0 * coef1 * sat.a * omeosq * (1.0 + 2.75 *
				(etasq + eeta) + eeta * etasq);
			cosio4 = cosio2 * cosio2;
			temp1 = 1.5 * EARTH_J2 * pinvsq * sat.no_unkozai;
			temp2 = 0.5 * temp1 * EARTH_J2 * pinvsq;
			temp3 = -0.46875 * EARTH_J4 * pinvsq * pinvsq * sat.no_unkozai;
			sat.mdot = sat.no_unkozai + 0.5 * temp1 * rteosq * sat.con41 + 0.0625 *
				temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
			sat.argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 *
				(7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
				temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
			xhdot1 = -temp1 * cosio;
			sat.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
				2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
			xpidot = sat.argpdot + sat.nodedot;
			sat.omgcof = sat.bstar * cc3 * cos(sat.argpo);
			sat.xmcof = 0.0;
			if (sat.ecco > 1.0e-4)
				sat.xmcof = -2.0 / 3.0 * coef * sat.bstar / eeta;
			sat.nodecf = 3.5 * omeosq * xhdot1 * sat.cc1;
			sat.t2cof = 1.5 * sat.cc1;

			// sgp4fix for divide by zero with xinco = 180 deg
			if (fabs(cosio + 1.0) > 1.5e-12)
				sat.xlcof = -0.25 * EARTH_J3_TO_J2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
			else
				sat.xlcof = -0.25 * EARTH_J3_TO_J2 * sinio * (3.0 + 5.0 * cosio) / temp4;
			sat.aycof = -0.5 * EARTH_J3_TO_J2 * sinio;

			// sgp4fix use multiply for speed instead of pow

			delmotemp = 1.0 + sat.eta * cos(sat.mo);

			sat.delmo = delmotemp * delmotemp * delmotemp;
			sat.sinmao = sin(sat.mo);
			sat.x7thm1 = 7.0 * cosio2 - 1.0;


			cc1sq = sat.cc1 * sat.cc1;
			sat.d2 = 4.0 * sat.a * tsi * cc1sq;
			temp = sat.d2 * tsi * sat.cc1 / 3.0;
			sat.d3 = (17.0 * sat.a + sfour) * temp;
			sat.d4 = 0.5 * temp * sat.a * tsi * (221.0 * sat.a + 31.0 * sfour) *
				sat.cc1;
			sat.t3cof = sat.d2 + 2.0 * cc1sq;
			sat.t4cof = 0.25 * (3.0 * sat.d3 + sat.cc1 *
				(12.0 * sat.d2 + 10.0 * cc1sq));
			sat.t5cof = 0.2 * (3.0 * sat.d4 +
				12.0 * sat.cc1 * sat.d3 +
				6.0 * sat.d2 * sat.d2 +
				15.0 * cc1sq * (2.0 * sat.d2 + cc1sq));
		}

		sgp4(sat, 0.0, r, v);

		sat.init = 'n';

		return true;
	}

	/*-----------------------------------------------------------------------------
	*
	*                             procedure sgp4
	*
	*  this procedure is the sgp4 prediction model from space command. this is an
	*    updated and combined version of sgp4 and sdp4, which were originally
	*    published separately in spacetrack report #3. this version follows the
	*    methodology from the aiaa paper (2006) describing the history and
	*    development of the code.
	*
	*  author        : david vallado                  719-573-2600   28 jun 2005
	*
	*  inputs        :
	*    sat	 - initialised structure from sgp4init() call.
	*    tsince	 - time since epoch (minutes)
	*
	*  outputs       :
	*    r           - position vector                     km
	*    v           - velocity                            km/sec
	*  return code - non-zero on error.
	*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
	*                   2 - mean motion less than 0.0
	*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
	*                   4 - semi-latus rectum < 0.0
	*                   5 - epoch elements are sub-orbital
	*                   6 - satellite has decayed
	*
	*  locals        :
	*    am          -
	*    axnl, aynl        -
	*    betal       -
	*    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
	*    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
	*    cosisq  , cossu   , sinsu   , cosu    , sinu
	*    delm        -
	*    delomg      -
	*    dndt        -
	*    eccm        -
	*    emsq        -
	*    ecose       -
	*    el2         -
	*    eo1         -
	*    eccp        -
	*    esine       -
	*    argpm       -
	*    argpp       -
	*    omgadf      -c
	*    pl          -
	*    r           -
	*    rtemsq      -
	*    rdotl       -
	*    rl          -
	*    rvdot       -
	*    rvdotl      -
	*    su          -
	*    t2  , t3   , t4    , tc
	*    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
	*    u   , ux   , uy    , uz     , vx     , vy     , vz
	*    inclm       - inclination
	*    mm          - mean anomaly
	*    nm          - mean motion
	*    nodem       - right asc of ascending node
	*    xinc        -
	*    xincp       -
	*    xl          -
	*    xlm         -
	*    mp          -
	*    xmdf        -
	*    xmx         -
	*    xmy         -
	*    nodedf      -
	*    xnode       -
	*    nodep       -
	*    np          -
	*
	*  coupling      :
	*    getgravconst- no longer used. Variables are conatined within sat
	*    dpper
	*    dpspace
	*
	*  references    :
	*    hoots, roehrich, norad spacetrack report #3 1980
	*    hoots, norad spacetrack report #6 1986
	*    hoots, schumacher and glover 2004
	*    vallado, crawford, hujsak, kelso  2006
	----------------------------------------------------------------------------*/

	bool sgp4(elsetrec& sat, double tsince, double r[3], double v[3]) {
		double am, axnl, aynl, betal, cosim, cnod,
			cos2u, coseo1, cosi, cosip, cosisq, cossu, cosu,
			delm, delomg, em, emsq, ecose, el2, eo1,
			ep, esine, argpm, argpp, argpdf, pl, mrt = 0.0,
			mvt, rdotl, rl, rvdot, rvdotl, sinim,
			sin2u, sineo1, sini, sinip, sinsu, sinu,
			snod, su, t2, t3, t4, tem5, temp,
			temp1, temp2, tempa, tempe, templ, u, ux,
			uy, uz, vx, vy, vz, inclm, mm,
			nm, nodem, xinc, xincp, xl, xlm, mp,
			xmdf, xmx, xmy, nodedf, xnode, nodep, tc, dndt, vkmpersec, delmtemp;
		int ktr;

		/* ------------------ set mathematical constants --------------- */
		// sgp4fix divisor for divide by zero check on inclination
		// the old check used 1.0 + cos(PI-1.0e-9), but then compared it to
		// 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
		const double temp4 = 1.5e-12;
		// sgp4fix identify constants and allow alternate values
		// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
		vkmpersec = EARTH_J2 * X_KE / 60.0;

		sat.t = tsince;

		/* ------- update for secular gravity and atmospheric drag ----- */
		xmdf = sat.mo + sat.mdot * sat.t;
		argpdf = sat.argpo + sat.argpdot * sat.t;
		nodedf = sat.nodeo + sat.nodedot * sat.t;
		argpm = argpdf;
		mm = xmdf;
		t2 = sat.t * sat.t;
		nodem = nodedf + sat.nodecf * t2;
		tempa = 1.0 - sat.cc1 * sat.t;
		tempe = sat.bstar * sat.cc4 * sat.t;
		templ = sat.t2cof * t2;

		if (sat.isimp != 1)
		{
			delomg = sat.omgcof * sat.t;
			// sgp4fix use mutliply for speed instead of pow
			delmtemp = 1.0 + sat.eta * cos(xmdf);
			delm = sat.xmcof *
				(delmtemp * delmtemp * delmtemp -
					sat.delmo);
			temp = delomg + delm;
			mm = xmdf + temp;
			argpm = argpdf - temp;
			t3 = t2 * sat.t;
			t4 = t3 * sat.t;
			tempa = tempa - sat.d2 * t2 - sat.d3 * t3 -
				sat.d4 * t4;
			tempe = tempe + sat.bstar * sat.cc5 * (sin(mm) -
				sat.sinmao);
			templ = templ + sat.t3cof * t3 + t4 * (sat.t4cof +
				sat.t * sat.t5cof);
		}

		nm = sat.no_unkozai;
		em = sat.ecco;
		inclm = sat.inclo;

		if (nm <= 0.0) {
			sat.error = 2;
			return false;
		}

		am = pow((X_KE / nm), 2.0 / 3.0) * tempa * tempa;
		nm = X_KE / pow(am, 1.5);
		em = em - tempe;

		// fix tolerance for error recognition
		// sgp4fix am is fixed from the previous nm check
		if ((em >= 1.0) || (em < -0.001)/* || (am < 0.95)*/) {
			sat.error = 1;
			return false;
		}

		// sgp4fix fix tolerance to avoid a divide by zero
		if (em < 1.0e-6)
			em = 1.0e-6;

		mm = mm + sat.no_unkozai * templ;
		xlm = mm + argpm + nodem;
		emsq = em * em;
		temp = 1.0 - emsq;

		nodem = fmod(nodem, TWO_PI);
		argpm = fmod(argpm, TWO_PI);
		xlm = fmod(xlm, TWO_PI);
		mm = fmod(xlm - argpm - nodem, TWO_PI);

		// sgp4fix recover singly averaged mean elements
		sat.am = am;
		sat.em = em;
		sat.im = inclm;
		sat.Om = nodem;
		sat.om = argpm;
		sat.mm = mm;
		sat.nm = nm;

		/* ----------------- compute extra mean quantities ------------- */
		sinim = sin(inclm);
		cosim = cos(inclm);

		/* -------------------- add lunar-solar periodics -------------- */
		ep = em;
		xincp = inclm;
		argpp = argpm;
		nodep = nodem;
		mp = mm;
		sinip = sinim;
		cosip = cosim;

		/* -------------------- long period periodics ------------------ */
		axnl = ep * cos(argpp);
		temp = 1.0 / (am * (1.0 - ep * ep));
		aynl = ep * sin(argpp) + temp * sat.aycof;
		xl = mp + argpp + nodep + temp * sat.xlcof * axnl;

		/* --------------------- solve kepler's equation --------------- */
		u = fmod(xl - nodep, TWO_PI);
		eo1 = u;
		tem5 = 9999.9;
		ktr = 1;
		//   sgp4fix for kepler iteration
		//   the following iteration needs better limits on corrections
		while ((fabs(tem5) >= 1.0e-12) && (ktr <= 10))
		{
			sineo1 = sin(eo1);
			coseo1 = cos(eo1);
			tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
			tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
			if (fabs(tem5) >= 0.95)
				tem5 = tem5 > 0.0 ? 0.95 : -0.95;
			eo1 = eo1 + tem5;
			ktr = ktr + 1;
		}

		/* ------------- short period preliminary quantities ----------- */
		ecose = axnl * coseo1 + aynl * sineo1;
		esine = axnl * sineo1 - aynl * coseo1;
		el2 = axnl * axnl + aynl * aynl;
		pl = am * (1.0 - el2);
		if (pl < 0.0)
		{
			sat.error = 4;
			return false;
		}

		rl = am * (1.0 - ecose);
		rdotl = sqrt(am) * esine / rl;
		rvdotl = sqrt(pl) / rl;
		betal = sqrt(1.0 - el2);
		temp = esine / (1.0 + betal);
		sinu = am / rl * (sineo1 - aynl - axnl * temp);
		cosu = am / rl * (coseo1 - axnl + aynl * temp);
		su = atan2(sinu, cosu);
		sin2u = (cosu + cosu) * sinu;
		cos2u = 1.0 - 2.0 * sinu * sinu;
		temp = 1.0 / pl;
		temp1 = 0.5 * EARTH_J2 * temp;
		temp2 = temp1 * temp;

		mrt = rl * (1.0 - 1.5 * temp2 * betal * sat.con41) +
			0.5 * temp1 * sat.x1mth2 * cos2u;
		su = su - 0.25 * temp2 * sat.x7thm1 * sin2u;
		xnode = nodep + 1.5 * temp2 * cosip * sin2u;
		xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
		mvt = rdotl - nm * temp1 * sat.x1mth2 * sin2u / X_KE;
		rvdot = rvdotl + nm * temp1 * (sat.x1mth2 * cos2u +
			1.5 * sat.con41) / X_KE;

		/* --------------------- orientation vectors ------------------- */
		sinsu = sin(su);
		cossu = cos(su);
		snod = sin(xnode);
		cnod = cos(xnode);
		sini = sin(xinc);
		cosi = cos(xinc);
		xmx = -snod * cosi;
		xmy = cnod * cosi;
		ux = xmx * sinsu + cnod * cossu;
		uy = xmy * sinsu + snod * cossu;
		uz = sini * sinsu;
		vx = xmx * cossu - cnod * sinsu;
		vy = xmy * cossu - snod * sinsu;
		vz = sini * cossu;

		/* --------- position and velocity (in km and km/sec) ---------- */
		r[0] = (mrt * ux) * EARTH_RAD;
		r[1] = (mrt * uy) * EARTH_RAD;
		r[2] = (mrt * uz) * EARTH_RAD;
		v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
		v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
		v[2] = (mvt * uz + rvdot * vz) * vkmpersec;


		// orbit decayed?
		if (mrt < 1.0) {
			sat.error = 6;
			return false;
		}

		return true;
	}

	// older sgp4io methods
	/* -----------------------------------------------------------------------------
	*
	*                           function twoline2rv
	*
	*  this function converts the two line element set character string data to
	*    variables and initializes the sgp4 variables. several intermediate varaibles
	*    and quantities are determined. note that the result is a structure so multiple
	*    satellites can be processed simaltaneously without having to reinitialize. the
	*    verification mode is an important option that permits quick checks of any
	*    changes to the underlying technical theory. this option works using a
	*    modified tle file in which the start, stop, and delta time values are
	*    included at the end of the second line of data. this only works with the
	*    verification mode. the catalog mode simply propagates from -1440 to 1440 min
	*    from epoch and is useful when performing entire catalog runs.
	*    update for alpha 5 numbering system. 4 mar 2021.
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs        :
	*    longstr1    - first line of the tle
	*    longstr2    - second line of the tle
	*    typerun     - type of run                    verification 'v', catalog 'c',
	*                                                 manual 'm'
	*    typeinput   - type of manual input           mfe 'm', epoch 'e', dayofyr 'd'
	*    opsmode     - mode of operation afspc or improved 'a', 'i'
	*    whichconst  - which set of constants to use  72, 84
	*
	*  outputs       :
	*    sat      - structure containing all the sgp4 satellite information
	*
	*  coupling      :
	*    getgravconst-
	*    days2mdhms  - conversion of days to month, day, hour, minute, second
	*    jday        - convert day month year hour minute second into julian date
	*    sgp4init    - initialize the sgp4 variables
	*
	*  references    :
	*    norad spacetrack report #3
	*    vallado, crawford, hujsak, kelso  2006
	--------------------------------------------------------------------------- */

	elsetrec twoline2rv(const sgp4::tle_set& set)
	{
		elsetrec sat;

		sat.classification = set.classification == sgp4::tle_set::classification_type::classified ? 'C' : (
			set.classification == sgp4::tle_set::classification_type::secret ? 'S' : 'U'
			);

		strcpy(sat.intldesg, (
			std::to_string(set.int_designator.year)
			+ std::to_string(set.int_designator.launch_number)
			+ std::string(set.int_designator.piece)).c_str()
		);
		sat.elnum = set.set_num;

		strcpy(sat.satnum, std::to_string(set.catalog_number).c_str());
		sat.revnum = set.rev_num;

		double julian_epoch = sgp4::time_utils::to_julian(set.epoch);
		double julian_since_midnight = julian_epoch - trunc(julian_epoch) - 0.5f;
		if (julian_since_midnight < 0.0)
			julian_since_midnight += 1.0;

		sat.jdsatepoch = julian_epoch;
		sat.jdsatepochF = julian_since_midnight;


		sat.no_kozai = set.mean_motion;
		sat.ndot = set.d_mean_motion;
		sat.nddot = set.dd_mean_motion;

		sat.inclo = set.inclination;
		sat.nodeo = set.right_ascension;
		sat.argpo = set.arg_of_perigee;
		sat.mo = set.mean_anomaly;

		sat.ecco = set.eccentricity;

		sat.bstar = set.rad_press_coef;

		// ------------ perform complete catalog evaluation, -+ 1 day -----------  why??
		double startmfe = -1440.0;
		double stopmfe = 1440.0;
		double deltamin = 10.0;

		// ---------------- initialize the orbit at sgp4epoch -------------------
		sgp4init(sat);

		return sat;
	} // twoline2rv
}