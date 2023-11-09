#include "SGP4_mod_ds.h"

#include <cmath>

static const double PI = 3.14159265358979323846;
static const double TWO_PI = 2 * PI;
static const double DEG_TO_RAD = PI / 180.0;

// gravitational parameter of the Earth [km^3/s^2]
static const double EARTH_GRV = 398600.5;

// radius of the Earth [km]
static const double EARTH_RAD = 6378.137;

// ??????
static const double X_KE = 60.0 / sqrt(EARTH_RAD * EARTH_RAD * EARTH_RAD / EARTH_GRV);

// ???? angular speed in an orbit of radius = EARTH_RAD? [rad/min]

// earth zonal harmonic model
static const double TUMIN = 1.0 / X_KE;
static const double EARTH_J2 = 0.00108262998905;
static const double EARTH_J3 = -0.00000253215306;
static const double EARTH_J4 = -0.00000161098761;
static const double EARTH_J3_TO_J2 = EARTH_J3 / EARTH_J2;

/* -----------------------------------------------------------------------------
	*
	*                           procedure dpper
	*
	*  this procedure provides deep space long period periodic contributions
	*    to the mean elements.  by design, these periodics are zero at epoch.
	*    this used to be dscom which included initialization, but it's really a
	*    recurring function.
	*
	*  author        : david vallado                  719-573-2600   28 jun 2005
	*
	*  inputs        :
	*    e3          -
	*    ee2         -
	*    peo         -
	*    pgho        -
	*    pho         -
	*    pinco       -
	*    plo         -
	*    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
	*    t           -
	*    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
	*    zmol        -
	*    zmos        -
	*    ep          - eccentricity                           0.0 - 1.0
	*    inclo       - inclination - needed for lyddane modification
	*    nodep       - right ascension of ascending node
	*    argpp       - argument of perigee
	*    mp          - mean anomaly
	*
	*  outputs       :
	*    ep          - eccentricity                           0.0 - 1.0
	*    inclp       - inclination
	*    nodep        - right ascension of ascending node
	*    argpp       - argument of perigee
	*    mp          - mean anomaly
	*
	*  locals        :
	*    alfdp       -
	*    betdp       -
	*    cosip  , sinip  , cosop  , sinop  ,
	*    dalf        -
	*    dbet        -
	*    dls         -
	*    f2, f3      -
	*    pe          -
	*    pgh         -
	*    ph          -
	*    pinc        -
	*    pl          -
	*    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
	*    sll   , sls
	*    xls         -
	*    xnoh        -
	*    zf          -
	*    zm          -
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    hoots, roehrich, norad spacetrack report #3 1980
	*    hoots, norad spacetrack report #6 1986
	*    hoots, schumacher and glover 2004
	*    vallado, crawford, hujsak, kelso  2006
	----------------------------------------------------------------------------*/

void SGP4Funcs_mod::dpper
(
	double e3, double ee2, double peo, double pgho, double pho,
	double pinco, double plo, double se2, double se3, double sgh2,
	double sgh3, double sgh4, double sh2, double sh3, double si2,
	double si3, double sl2, double sl3, double sl4, double t,
	double xgh2, double xgh3, double xgh4, double xh2, double xh3,
	double xi2, double xi3, double xl2, double xl3, double xl4,
	double zmol, double zmos, double inclo,
	char init,
	double& ep, double& inclp, double& nodep, double& argpp, double& mp
)
{
	/* --------------------- local variables ------------------------ */
	double alfdp, betdp, cosip, cosop, dalf, dbet, dls,
		f2, f3, pe, pgh, ph, pinc, pl,
		sel, ses, sghl, sghs, shll, shs, sil,
		sinip, sinop, sinzf, sis, sll, sls, xls,
		xnoh, zf, zm, zel, zes, znl, zns;

	/* ---------------------- constants ----------------------------- */
	zns = 1.19459e-5;
	zes = 0.01675;
	znl = 1.5835218e-4;
	zel = 0.05490;

	/* --------------- calculate time varying periodics ----------- */
	zm = zmos + zns * t;
	// be sure that the initial call has time set to zero
	if (init == 'y')
		zm = zmos;
	zf = zm + 2.0 * zes * sin(zm);
	sinzf = sin(zf);
	f2 = 0.5 * sinzf * sinzf - 0.25;
	f3 = -0.5 * sinzf * cos(zf);
	ses = se2 * f2 + se3 * f3;
	sis = si2 * f2 + si3 * f3;
	sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
	sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
	shs = sh2 * f2 + sh3 * f3;
	zm = zmol + znl * t;
	if (init == 'y')
		zm = zmol;
	zf = zm + 2.0 * zel * sin(zm);
	sinzf = sin(zf);
	f2 = 0.5 * sinzf * sinzf - 0.25;
	f3 = -0.5 * sinzf * cos(zf);
	sel = ee2 * f2 + e3 * f3;
	sil = xi2 * f2 + xi3 * f3;
	sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
	sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
	shll = xh2 * f2 + xh3 * f3;
	pe = ses + sel;
	pinc = sis + sil;
	pl = sls + sll;
	pgh = sghs + sghl;
	ph = shs + shll;

	if (init == 'n')
	{
		pe = pe - peo;
		pinc = pinc - pinco;
		pl = pl - plo;
		pgh = pgh - pgho;
		ph = ph - pho;
		inclp = inclp + pinc;
		ep = ep + pe;
		sinip = sin(inclp);
		cosip = cos(inclp);

		/* ----------------- apply periodics directly ------------ */
		//  sgp4fix for lyddane choice
		//  strn3 used original inclination - this is technically feasible
		//  gsfc used perturbed inclination - also technically feasible
		//  probably best to readjust the 0.2 limit value and limit discontinuity
		//  0.2 rad = 11.45916 deg
		//  use next line for original strn3 approach and original inclination
		//  if (inclo >= 0.2)
		//  use next line for gsfc version and perturbed inclination
		if (inclp >= 0.2)
		{
			ph = ph / sinip;
			pgh = pgh - cosip * ph;
			argpp = argpp + pgh;
			nodep = nodep + ph;
			mp = mp + pl;
		}
		else
		{
			/* ---- apply periodics with lyddane modification ---- */
			sinop = sin(nodep);
			cosop = cos(nodep);
			alfdp = sinip * sinop;
			betdp = sinip * cosop;
			dalf = ph * cosop + pinc * cosip * sinop;
			dbet = -ph * sinop + pinc * cosip * cosop;
			alfdp = alfdp + dalf;
			betdp = betdp + dbet;
			nodep = fmod(nodep, TWO_PI);
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			xls = mp + argpp + cosip * nodep;
			dls = pl + pgh - pinc * nodep * sinip;
			xls = xls + dls;
			xnoh = nodep;
			nodep = atan2(alfdp, betdp);
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if (fabs(xnoh - nodep) > PI)
				if (nodep < xnoh)
					nodep = nodep + TWO_PI;
				else
					nodep = nodep - TWO_PI;
			mp = mp + pl;
			argpp = xls - mp - cosip * nodep;
		}
	}   // if init == 'n'

	//#include "debug1.cpp"
}  // dpper

/*-----------------------------------------------------------------------------
*
*                           procedure dscom
*
*  this procedure provides deep space common items used by both the secular
*    and periodics subroutines.  input is provided as shown. this routine
*    used to be called dpper, but the functions inside weren't well organized.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    epoch       -
*    ep          - eccentricity
*    argpp       - argument of perigee
*    tc          -
*    inclp       - inclination
*    nodep       - right ascension of ascending node
*    np          - mean motion
*
*  outputs       :
*    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
*    day         -
*    e3          -
*    ee2         -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    gam         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    rtemsq      -
*    se2, se3         -
*    sgh2, sgh3, sgh4        -
*    sh2, sh3, si2, si3, sl2, sl3, sl4         -
*    s1, s2, s3, s4, s5, s6, s7          -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
*    nm          - mean motion
*    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*    zmol        -
*    zmos        -
*
*  locals        :
*    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
*    betasq      -
*    cc          -
*    ctem, stem        -
*    x1, x2, x3, x4, x5, x6, x7, x8          -
*    xnodce      -
*    xnoi        -
*    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
*    zcosi  , zsini  , zcosil , zsinil ,
*    zx          -
*    zy          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

 void SGP4Funcs_mod::dscom
(
	double epoch, double ep, double argpp, double tc, double inclp,
	double nodep, double np,
	double& snodm, double& cnodm, double& sinim, double& cosim, double& sinomm,
	double& cosomm, double& day, double& e3, double& ee2, double& em,
	double& emsq, double& gam, double& peo, double& pgho, double& pho,
	double& pinco, double& plo, double& rtemsq, double& se2, double& se3,
	double& sgh2, double& sgh3, double& sgh4, double& sh2, double& sh3,
	double& si2, double& si3, double& sl2, double& sl3, double& sl4,
	double& s1, double& s2, double& s3, double& s4, double& s5,
	double& s6, double& s7, double& ss1, double& ss2, double& ss3,
	double& ss4, double& ss5, double& ss6, double& ss7, double& sz1,
	double& sz2, double& sz3, double& sz11, double& sz12, double& sz13,
	double& sz21, double& sz22, double& sz23, double& sz31, double& sz32,
	double& sz33, double& xgh2, double& xgh3, double& xgh4, double& xh2,
	double& xh3, double& xi2, double& xi3, double& xl2, double& xl3,
	double& xl4, double& nm, double& z1, double& z2, double& z3,
	double& z11, double& z12, double& z13, double& z21, double& z22,
	double& z23, double& z31, double& z32, double& z33, double& zmol,
	double& zmos
)
{
	/* -------------------------- constants ------------------------- */
	const double zes = 0.01675;
	const double zel = 0.05490;
	const double c1ss = 2.9864797e-6;
	const double c1l = 4.7968065e-7;
	const double zsinis = 0.39785416;
	const double zcosis = 0.91744867;
	const double zcosgs = 0.1945905;
	const double zsings = -0.98088458;

	/* --------------------- local variables ------------------------ */
	int lsflg;
	double a1, a2, a3, a4, a5, a6, a7,
		a8, a9, a10, betasq, cc, ctem, stem,
		x1, x2, x3, x4, x5, x6, x7,
		x8, xnodce, xnoi, zcosg, zcosgl, zcosh, zcoshl,
		zcosi, zcosil, zsing, zsingl, zsinh, zsinhl, zsini,
		zsinil, zx, zy;

	nm = np;
	em = ep;
	snodm = sin(nodep);
	cnodm = cos(nodep);
	sinomm = sin(argpp);
	cosomm = cos(argpp);
	sinim = sin(inclp);
	cosim = cos(inclp);
	emsq = em * em;
	betasq = 1.0 - emsq;
	rtemsq = sqrt(betasq);

	/* ----------------- initialize lunar solar terms --------------- */
	peo = 0.0;
	pinco = 0.0;
	plo = 0.0;
	pgho = 0.0;
	pho = 0.0;
	day = epoch + 18261.5 + tc / 1440.0;
	xnodce = fmod(4.5236020 - 9.2422029e-4 * day, TWO_PI);
	stem = sin(xnodce);
	ctem = cos(xnodce);
	zcosil = 0.91375164 - 0.03568096 * ctem;
	zsinil = sqrt(1.0 - zcosil * zcosil);
	zsinhl = 0.089683511 * stem / zsinil;
	zcoshl = sqrt(1.0 - zsinhl * zsinhl);
	gam = 5.8351514 + 0.0019443680 * day;
	zx = 0.39785416 * stem / zsinil;
	zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
	zx = atan2(zx, zy);
	zx = gam + zx - xnodce;
	zcosgl = cos(zx);
	zsingl = sin(zx);

	/* ------------------------- do solar terms --------------------- */
	zcosg = zcosgs;
	zsing = zsings;
	zcosi = zcosis;
	zsini = zsinis;
	zcosh = cnodm;
	zsinh = snodm;
	cc = c1ss;
	xnoi = 1.0 / nm;

	for (lsflg = 1; lsflg <= 2; lsflg++)
	{
		a1 = zcosg * zcosh + zsing * zcosi * zsinh;
		a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
		a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
		a8 = zsing * zsini;
		a9 = zsing * zsinh + zcosg * zcosi * zcosh;
		a10 = zcosg * zsini;
		a2 = cosim * a7 + sinim * a8;
		a4 = cosim * a9 + sinim * a10;
		a5 = -sinim * a7 + cosim * a8;
		a6 = -sinim * a9 + cosim * a10;

		x1 = a1 * cosomm + a2 * sinomm;
		x2 = a3 * cosomm + a4 * sinomm;
		x3 = -a1 * sinomm + a2 * cosomm;
		x4 = -a3 * sinomm + a4 * cosomm;
		x5 = a5 * sinomm;
		x6 = a6 * sinomm;
		x7 = a5 * cosomm;
		x8 = a6 * cosomm;

		z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
		z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
		z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
		z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * emsq;
		z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * emsq;
		z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * emsq;
		z11 = -6.0 * a1 * a5 + emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
		z12 = -6.0 * (a1 * a6 + a3 * a5) + emsq *
			(-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
		z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
		z21 = 6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
		z22 = 6.0 * (a4 * a5 + a2 * a6) + emsq *
			(24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
		z23 = 6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
		z1 = z1 + z1 + betasq * z31;
		z2 = z2 + z2 + betasq * z32;
		z3 = z3 + z3 + betasq * z33;
		s3 = cc * xnoi;
		s2 = -0.5 * s3 / rtemsq;
		s4 = s3 * rtemsq;
		s1 = -15.0 * em * s4;
		s5 = x1 * x3 + x2 * x4;
		s6 = x2 * x3 + x1 * x4;
		s7 = x2 * x4 - x1 * x3;

		/* ----------------------- do lunar terms ------------------- */
		if (lsflg == 1)
		{
			ss1 = s1;
			ss2 = s2;
			ss3 = s3;
			ss4 = s4;
			ss5 = s5;
			ss6 = s6;
			ss7 = s7;
			sz1 = z1;
			sz2 = z2;
			sz3 = z3;
			sz11 = z11;
			sz12 = z12;
			sz13 = z13;
			sz21 = z21;
			sz22 = z22;
			sz23 = z23;
			sz31 = z31;
			sz32 = z32;
			sz33 = z33;
			zcosg = zcosgl;
			zsing = zsingl;
			zcosi = zcosil;
			zsini = zsinil;
			zcosh = zcoshl * cnodm + zsinhl * snodm;
			zsinh = snodm * zcoshl - cnodm * zsinhl;
			cc = c1l;
		}
	}

	zmol = fmod(4.7199672 + 0.22997150 * day - gam, TWO_PI);
	zmos = fmod(6.2565837 + 0.017201977 * day, TWO_PI);

	/* ------------------------ do solar terms ---------------------- */
	se2 = 2.0 * ss1 * ss6;
	se3 = 2.0 * ss1 * ss7;
	si2 = 2.0 * ss2 * sz12;
	si3 = 2.0 * ss2 * (sz13 - sz11);
	sl2 = -2.0 * ss3 * sz2;
	sl3 = -2.0 * ss3 * (sz3 - sz1);
	sl4 = -2.0 * ss3 * (-21.0 - 9.0 * emsq) * zes;
	sgh2 = 2.0 * ss4 * sz32;
	sgh3 = 2.0 * ss4 * (sz33 - sz31);
	sgh4 = -18.0 * ss4 * zes;
	sh2 = -2.0 * ss2 * sz22;
	sh3 = -2.0 * ss2 * (sz23 - sz21);

	/* ------------------------ do lunar terms ---------------------- */
	ee2 = 2.0 * s1 * s6;
	e3 = 2.0 * s1 * s7;
	xi2 = 2.0 * s2 * z12;
	xi3 = 2.0 * s2 * (z13 - z11);
	xl2 = -2.0 * s3 * z2;
	xl3 = -2.0 * s3 * (z3 - z1);
	xl4 = -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel;
	xgh2 = 2.0 * s4 * z32;
	xgh3 = 2.0 * s4 * (z33 - z31);
	xgh4 = -18.0 * s4 * zel;
	xh2 = -2.0 * s2 * z22;
	xh3 = -2.0 * s2 * (z23 - z21);

	//#include "debug2.cpp"
}  // dscom

/*-----------------------------------------------------------------------------
*
*                           procedure dsinit
*
*  this procedure provides deep space contributions to mean motion dot due
*    to geopotential resonance with half day and one day orbits.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    xke         - reciprocal of tumin
*    cosim, sinim-
*    emsq        - eccentricity squared
*    argpo       - argument of perigee
*    s1, s2, s3, s4, s5      -
*    ss1, ss2, ss3, ss4, ss5 -
*    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
*    t           - time
*    tc          -
*    gsto        - greenwich sidereal time                   rad
*    mo          - mean anomaly
*    mdot        - mean anomaly dot (rate)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*    nodedot     - right ascension of ascending node dot (rate)
*    xpidot      -
*    z1, z3, z11, z13, z21, z23, z31, z33 -
*    eccm        - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    xn          - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right ascension of ascending node
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    atime       -
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
*    dedt        -
*    didt        -
*    dmdt        -
*    dndt        -
*    dnodt       -
*    domdt       -
*    del1, del2, del3        -
*    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
*    theta       -
*    xfact       -
*    xlamo       -
*    xli         -
*    xni
*
*  locals        :
*    ainv2       -
*    aonv        -
*    cosisq      -
*    eoc         -
*    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
*    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
*    sini2       -
*    temp        -
*    temp1       -
*    theta       -
*    xno2        -
*
*  coupling      :
*    getgravconst- no longer used
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

void SGP4Funcs_mod::dsinit
(
	// sgp4fix just send in xke as a constant and eliminate getgravconst call
	// gravconsttype whichconst, 
	double xke,
	double cosim, double emsq, double argpo, double s1, double s2,
	double s3, double s4, double s5, double sinim, double ss1,
	double ss2, double ss3, double ss4, double ss5, double sz1,
	double sz3, double sz11, double sz13, double sz21, double sz23,
	double sz31, double sz33, double t, double tc, double gsto,
	double mo, double mdot, double no, double nodeo, double nodedot,
	double xpidot, double z1, double z3, double z11, double z13,
	double z21, double z23, double z31, double z33, double ecco,
	double eccsq, double& em, double& argpm, double& inclm, double& mm,
	double& nm, double& nodem,
	int& irez,
	double& atime, double& d2201, double& d2211, double& d3210, double& d3222,
	double& d4410, double& d4422, double& d5220, double& d5232, double& d5421,
	double& d5433, double& dedt, double& didt, double& dmdt, double& dndt,
	double& dnodt, double& domdt, double& del1, double& del2, double& del3,
	double& xfact, double& xlamo, double& xli, double& xni
)
{
	/* --------------------- local variables ------------------------ */
	double ainv2, aonv = 0.0, cosisq, eoc, f220, f221, f311,
		f321, f322, f330, f441, f442, f522, f523,
		f542, f543, g200, g201, g211, g300, g310,
		g322, g410, g422, g520, g521, g532, g533,
		ses, sgs, sghl, sghs, shs, shll, sis,
		sini2, sls, temp, temp1, theta, xno2, q22,
		q31, q33, root22, root44, root54, rptim, root32,
		root52, x2o3, znl, emo, zns, emsqo;

	q22 = 1.7891679e-6;
	q31 = 2.1460748e-6;
	q33 = 2.2123015e-7;
	root22 = 1.7891679e-6;
	root44 = 7.3636953e-9;
	root54 = 2.1765803e-9;
	rptim = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
	root32 = 3.7393792e-7;
	root52 = 1.1428639e-7;
	x2o3 = 2.0 / 3.0;
	znl = 1.5835218e-4;
	zns = 1.19459e-5;

	// sgp4fix identify constants and allow alternate values
	// just xke is used here so pass it in rather than have multiple calls
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

	/* -------------------- deep space initialization ------------ */
	irez = 0;
	if ((nm < 0.0052359877) && (nm > 0.0034906585))
		irez = 1;
	if ((nm >= 8.26e-3) && (nm <= 9.24e-3) && (em >= 0.5))
		irez = 2;

	/* ------------------------ do solar terms ------------------- */
	ses = ss1 * zns * ss5;
	sis = ss2 * zns * (sz11 + sz13);
	sls = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
	sghs = ss4 * zns * (sz31 + sz33 - 6.0);
	shs = -zns * ss2 * (sz21 + sz23);
	// sgp4fix for 180 deg incl
	if ((inclm < 5.2359877e-2) || (inclm > PI - 5.2359877e-2))
		shs = 0.0;
	if (sinim != 0.0)
		shs = shs / sinim;
	sgs = sghs - cosim * shs;

	/* ------------------------- do lunar terms ------------------ */
	dedt = ses + s1 * znl * s5;
	didt = sis + s2 * znl * (z11 + z13);
	dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
	sghl = s4 * znl * (z31 + z33 - 6.0);
	shll = -znl * s2 * (z21 + z23);
	// sgp4fix for 180 deg incl
	if ((inclm < 5.2359877e-2) || (inclm > PI - 5.2359877e-2))
		shll = 0.0;
	domdt = sgs + sghl;
	dnodt = shs;
	if (sinim != 0.0)
	{
		domdt = domdt - cosim / sinim * shll;
		dnodt = dnodt + shll / sinim;
	}

	/* ----------- calculate deep space resonance effects -------- */
	dndt = 0.0;
	theta = fmod(gsto + tc * rptim, TWO_PI);
	em = em + dedt * t;
	inclm = inclm + didt * t;
	argpm = argpm + domdt * t;
	nodem = nodem + dnodt * t;
	mm = mm + dmdt * t;
	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//if (inclm < 0.0)
	//  {
	//    inclm  = -inclm;
	//    argpm  = argpm - PI;
	//    nodem = nodem + PI;
	//  }

	/* -------------- initialize the resonance terms ------------- */
	if (irez != 0)
	{
		aonv = pow(nm / xke, x2o3);

		/* ---------- geopotential resonance for 12 hour orbits ------ */
		if (irez == 2)
		{
			cosisq = cosim * cosim;
			emo = em;
			em = ecco;
			emsqo = emsq;
			emsq = eccsq;
			eoc = em * emsq;
			g201 = -0.306 - (em - 0.64) * 0.440;

			if (em <= 0.65)
			{
				g211 = 3.616 - 13.2470 * em + 16.2900 * emsq;
				g310 = -19.302 + 117.3900 * em - 228.4190 * emsq + 156.5910 * eoc;
				g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq + 146.5816 * eoc;
				g410 = -41.122 + 242.6940 * em - 471.0940 * emsq + 313.9530 * eoc;
				g422 = -146.407 + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
				g520 = -532.114 + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
			}
			else
			{
				g211 = -72.099 + 331.819 * em - 508.738 * emsq + 266.724 * eoc;
				g310 = -346.844 + 1582.851 * em - 2415.925 * emsq + 1246.113 * eoc;
				g322 = -342.585 + 1554.908 * em - 2366.899 * emsq + 1215.972 * eoc;
				g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq + 3651.957 * eoc;
				g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
				if (em > 0.715)
					g520 = -5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
				else
					g520 = 1464.74 - 4664.75 * em + 3763.64 * emsq;
			}
			if (em < 0.7)
			{
				g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21 * eoc;
				g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
				g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4 * eoc;
			}
			else
			{
				g533 = -37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
				g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
				g532 = -40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
			}

			sini2 = sinim * sinim;
			f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq);
			f221 = 1.5 * sini2;
			f321 = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
			f322 = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
			f441 = 35.0 * sini2 * f220;
			f442 = 39.3750 * sini2 * sini2;
			f522 = 9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq) +
				0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
			f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim +
				10.0 * cosisq) + 6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
			f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim + cosisq *
				(-12.0 + 8.0 * cosim + 10.0 * cosisq));
			f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim + cosisq *
				(12.0 + 8.0 * cosim - 10.0 * cosisq));
			xno2 = nm * nm;
			ainv2 = aonv * aonv;
			temp1 = 3.0 * xno2 * ainv2;
			temp = temp1 * root22;
			d2201 = temp * f220 * g201;
			d2211 = temp * f221 * g211;
			temp1 = temp1 * aonv;
			temp = temp1 * root32;
			d3210 = temp * f321 * g310;
			d3222 = temp * f322 * g322;
			temp1 = temp1 * aonv;
			temp = 2.0 * temp1 * root44;
			d4410 = temp * f441 * g410;
			d4422 = temp * f442 * g422;
			temp1 = temp1 * aonv;
			temp = temp1 * root52;
			d5220 = temp * f522 * g520;
			d5232 = temp * f523 * g532;
			temp = 2.0 * temp1 * root54;
			d5421 = temp * f542 * g521;
			d5433 = temp * f543 * g533;
			xlamo = fmod(mo + nodeo + nodeo - theta - theta, TWO_PI);
			xfact = mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no;
			em = emo;
			emsq = emsqo;
		}

		/* ---------------- synchronous resonance terms -------------- */
		if (irez == 1)
		{
			g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
			g310 = 1.0 + 2.0 * emsq;
			g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
			f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
			f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
			f330 = 1.0 + cosim;
			f330 = 1.875 * f330 * f330 * f330;
			del1 = 3.0 * nm * nm * aonv * aonv;
			del2 = 2.0 * del1 * f220 * g200 * q22;
			del3 = 3.0 * del1 * f330 * g300 * q33 * aonv;
			del1 = del1 * f311 * g310 * q31 * aonv;
			xlamo = fmod(mo + nodeo + argpo - theta, TWO_PI);
			xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no;
		}

		/* ------------ for sgp4, initialize the integrator ---------- */
		xli = xlamo;
		xni = no;
		atime = 0.0;
		nm = no + dndt;
	}

	//#include "debug3.cpp"
}  // dsinit

/*-----------------------------------------------------------------------------
*
*                           procedure dspace
*
*  this procedure provides deep space contributions to mean elements for
*    perturbing third body.  these effects have been averaged over one
*    revolution of the sun and moon.  for earth resonance effects, the
*    effects have been averaged over no revolutions of the satellite.
*    (mean motion)
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
*    dedt        -
*    del1, del2, del3  -
*    didt        -
*    dmdt        -
*    dnodt       -
*    domdt       -
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    argpo       - argument of perigee
*    argpdot     - argument of perigee dot (rate)
*    t           - time
*    tc          -
*    gsto        - gst
*    xfact       -
*    xlamo       -
*    no          - mean motion
*    atime       -
*    em          - eccentricity
*    ft          -
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    atime       -
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         -
*    nodem       - right ascension of ascending node
*    dndt        -
*    nm          - mean motion
*
*  locals        :
*    delt        -
*    ft          -
*    theta       -
*    x2li        -
*    x2omi       -
*    xl          -
*    xldot       -
*    xnddt       -
*    xndt        -
*    xomi        -
*
*  coupling      :
*    none        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

void SGP4Funcs_mod::dspace
(
	int irez,
	double d2201, double d2211, double d3210, double d3222, double d4410,
	double d4422, double d5220, double d5232, double d5421, double d5433,
	double dedt, double del1, double del2, double del3, double didt,
	double dmdt, double dnodt, double domdt, double argpo, double argpdot,
	double t, double tc, double gsto, double xfact, double xlamo,
	double no,
	double& atime, double& em, double& argpm, double& inclm, double& xli,
	double& mm, double& xni, double& nodem, double& dndt, double& nm
)
{
	int iretn, iret;
	double delt, ft, theta, x2li, x2omi, xl, xldot, xnddt, xndt, xomi, g22, g32,
		g44, g52, g54, fasx2, fasx4, fasx6, rptim, step2, stepn, stepp;

	fasx2 = 0.13130908;
	fasx4 = 2.8843198;
	fasx6 = 0.37448087;
	g22 = 5.7686396;
	g32 = 0.95240898;
	g44 = 1.8014998;
	g52 = 1.0508330;
	g54 = 4.4108898;
	rptim = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
	stepp = 720.0;
	stepn = -720.0;
	step2 = 259200.0;

	/* ----------- calculate deep space resonance effects ----------- */
	dndt = 0.0;
	theta = fmod(gsto + tc * rptim, TWO_PI);
	em = em + dedt * t;

	inclm = inclm + didt * t;
	argpm = argpm + domdt * t;
	nodem = nodem + dnodt * t;
	mm = mm + dmdt * t;

	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//  if (inclm < 0.0)
	// {
	//    inclm = -inclm;
	//    argpm = argpm - PI;
	//    nodem = nodem + PI;
	//  }

	/* - update resonances : numerical (euler-maclaurin) integration - */
	/* ------------------------- epoch restart ----------------------  */
	//   sgp4fix for propagator problems
	//   the following integration works for negative time steps and periods
	//   the specific changes are unknown because the original code was so convoluted

	// sgp4fix take out atime = 0.0 and fix for faster operation
	ft = 0.0;
	if (irez != 0)
	{
		// sgp4fix streamline check
		if ((atime == 0.0) || (t * atime <= 0.0) || (fabs(t) < fabs(atime)))
		{
			atime = 0.0;
			xni = no;
			xli = xlamo;
		}
		// sgp4fix move check outside loop
		if (t > 0.0)
			delt = stepp;
		else
			delt = stepn;

		iretn = 381; // added for do loop
		iret = 0; // added for loop
		while (iretn == 381)
		{
			/* ------------------- dot terms calculated ------------- */
			/* ----------- near - synchronous resonance terms ------- */
			if (irez != 2)
			{
				xndt = del1 * sin(xli - fasx2) + del2 * sin(2.0 * (xli - fasx4)) +
					del3 * sin(3.0 * (xli - fasx6));
				xldot = xni + xfact;
				xnddt = del1 * cos(xli - fasx2) +
					2.0 * del2 * cos(2.0 * (xli - fasx4)) +
					3.0 * del3 * cos(3.0 * (xli - fasx6));
				xnddt = xnddt * xldot;
			}
			else
			{
				/* --------- near - half-day resonance terms -------- */
				xomi = argpo + argpdot * atime;
				x2omi = xomi + xomi;
				x2li = xli + xli;
				xndt = d2201 * sin(x2omi + xli - g22) + d2211 * sin(xli - g22) +
					d3210 * sin(xomi + xli - g32) + d3222 * sin(-xomi + xli - g32) +
					d4410 * sin(x2omi + x2li - g44) + d4422 * sin(x2li - g44) +
					d5220 * sin(xomi + xli - g52) + d5232 * sin(-xomi + xli - g52) +
					d5421 * sin(xomi + x2li - g54) + d5433 * sin(-xomi + x2li - g54);
				xldot = xni + xfact;
				xnddt = d2201 * cos(x2omi + xli - g22) + d2211 * cos(xli - g22) +
					d3210 * cos(xomi + xli - g32) + d3222 * cos(-xomi + xli - g32) +
					d5220 * cos(xomi + xli - g52) + d5232 * cos(-xomi + xli - g52) +
					2.0 * (d4410 * cos(x2omi + x2li - g44) +
						d4422 * cos(x2li - g44) + d5421 * cos(xomi + x2li - g54) +
						d5433 * cos(-xomi + x2li - g54));
				xnddt = xnddt * xldot;
			}

			/* ----------------------- integrator ------------------- */
			// sgp4fix move end checks to end of routine
			if (fabs(t - atime) >= stepp)
			{
				iret = 0;
				iretn = 381;
			}
			else // exit here
			{
				ft = t - atime;
				iretn = 0;
			}

			if (iretn == 381)
			{
				xli = xli + xldot * delt + xndt * step2;
				xni = xni + xndt * delt + xnddt * step2;
				atime = atime + delt;
			}
		}  // while iretn = 381

		nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
		xl = xli + xldot * ft + xndt * ft * ft * 0.5;
		if (irez != 1)
		{
			mm = xl - 2.0 * nodem + 2.0 * theta;
			dndt = nm - no;
		}
		else
		{
			mm = xl - nodem - argpm + theta;
			dndt = nm - no;
		}
		nm = no + dndt;
	}

	//#include "debug4.cpp"
}  // dsspace