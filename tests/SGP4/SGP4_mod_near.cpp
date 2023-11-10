#include "SGP4_mod_near.h"

#include <iostream>
#include <iomanip>
#include <utility>

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

// other SGP4 parameters

// default value of the s constant
static const double S_PARAM_DEF = 78.0 / EARTH_RAD + 1.0;
static const double Q0_PARAM_DEF = 120.0 / EARTH_RAD + 1.0;

namespace SGP4Funcs_mod_near
{
	std::pair<double, double> find_org_sma_and_mean_mot(
		double tle_mean_motion, double cos2_incl, double rt_one_min_ecc_pow2
	) {
		double a_1 = pow(X_KE / tle_mean_motion, 2.0 / 3.0);
		double delta_comm =
			0.75 * EARTH_J2 * (3.0 * cos2_incl - 1.0)
			/ (rt_one_min_ecc_pow2 * rt_one_min_ecc_pow2 * rt_one_min_ecc_pow2);

		double delta_1 = delta_comm / (a_1 * a_1);

		double a_0 = a_1 * (1.0 - delta_1 * delta_1 - delta_1 *
			(1.0 / 3.0 + 134.0 * delta_1 * delta_1 / 81.0));
		double delta_0 = delta_comm / (a_0 * a_0);
		
		double mean_motion = tle_mean_motion / (1.0 + delta_0);
		double sma = a_0 / (1.0 - delta_0);

		return std::make_pair(sma, mean_motion);
	}

	bool sgp4init(elsetrec& sat)
	{
		double epoch_1950 = sat.jdsatepoch + sat.jdsatepochF - 2433281.5;
		/* --------------------- local variables ------------------------ */
		double cosim, sinim, sinomm, cc1sq,
			dndt, temp, temp1, temp2, temp3, xpidot,
			xhdot1, r[3], v[3], delmotemp;

		double ecc_pow2 = sat.ecco * sat.ecco;
		double one_min_ecc_pow2 = 1.0 - ecc_pow2;
		double rt_one_min_ecc_pow2 = sqrt(one_min_ecc_pow2);
		double cos_incl = cos(sat.inclo);
		double sin_incl = sin(sat.inclo);
		double cos2_incl = cos_incl * cos_incl;

		double one_min_5_cos2_incl = 1.0 - 5.0 * cos2_incl;
		
		sat.con41 = -one_min_5_cos2_incl - 2 * cos2_incl;

		auto [sma, mean_motion] = find_org_sma_and_mean_mot(sat.no_kozai, cos2_incl, rt_one_min_ecc_pow2);
		sat.no_unkozai = mean_motion;
		sat.a = sma;

		// why?
		double posq = sma * one_min_ecc_pow2 * sma * one_min_ecc_pow2;

		// a, alta and altp are in the unitis of earths radii
		sat.a = pow(sat.no_unkozai * TUMIN, (-2.0 / 3.0));
		sat.alta = sat.a * (1.0 + sat.ecco) - 1.0;
		sat.altp = sat.a * (1.0 - sat.ecco) - 1.0;

		if (sat.ecco < 0.0 || sat.no_unkozai < 0.0) {
			std::cerr << "There seems to be some error" << std::endl;
			return false;
		}

		// no idea what it does
		double perigee = sat.altp * EARTH_RAD;
		
		sat.isimp = 0;
		if (perigee < 220.0)
			sat.isimp = 1;

			
		// adjust the s constant
		double s_const;
		if (perigee >= 156.0)
			s_const = S_PARAM_DEF;
		else if (perigee < 98.0)
			s_const = 20.0 / EARTH_RAD + 1.0;
		else
			s_const = sat.a * (1 - sat.ecco) - S_PARAM_DEF + 1.0;

		double q0_min_s = Q0_PARAM_DEF - s_const;
		double q0_min_s_pow4 = q0_min_s * q0_min_s * q0_min_s * q0_min_s;

		// page 11 of the report
		double tsi = 1.0 / (sat.a - s_const);
		
		sat.eta = sat.a * sat.ecco * tsi;
		double eta_pow2 = sat.eta * sat.eta;
		double ecc_eta = sat.ecco * sat.eta;

		double psisq = fabs(1.0 - eta_pow2);
		double coef = q0_min_s_pow4 * pow(tsi, 4.0);
		double coef1 = coef / pow(psisq, 3.5);
		
		double cc2 = coef1 * sat.no_unkozai * (sat.a * (1.0 + 1.5 * eta_pow2 + ecc_eta *
			(4.0 + eta_pow2)) + 0.375 * EARTH_J2 * tsi / psisq * sat.con41 *
			(8.0 + 3.0 * eta_pow2 * (8.0 + eta_pow2)));

		sat.cc1 = sat.bstar * cc2;
		double cc3 = 0.0;
		if (sat.ecco > 1.0e-4)
			cc3 = -2.0 * coef * tsi * EARTH_J3_TO_J2 * sat.no_unkozai
			* sin_incl / sat.ecco;
		
		sat.x1mth2 = 1.0 - cos2_incl;
		sat.cc4 = 2.0 * sat.no_unkozai * coef1 * sat.a * one_min_ecc_pow2 *
			(sat.eta * (2.0 + 0.5 * eta_pow2) + sat.ecco *
				(0.5 + 2.0 * eta_pow2) - EARTH_J2 * tsi / (sat.a * psisq) *
				(-3.0 * sat.con41 * (1.0 - 2.0 * ecc_eta + eta_pow2 *
					(1.5 - 0.5 * ecc_eta)) + 0.75 * sat.x1mth2 *
					(2.0 * eta_pow2 - ecc_eta * (1.0 + eta_pow2)) * cos(2.0 * sat.argpo)));
		
		sat.cc5 = 2.0 * coef1 * sat.a * one_min_ecc_pow2 * (1.0 + 2.75 *
			(eta_pow2 + ecc_eta) + ecc_eta * eta_pow2);
		
		double cos4_incl = cos2_incl * cos2_incl;
		temp1 = 1.5 * EARTH_J2 / posq * sat.no_unkozai;
		temp2 = 0.5 * temp1 * EARTH_J2 / posq;
		temp3 = -0.46875 * EARTH_J4 / posq / posq * sat.no_unkozai;
		
		sat.mdot = sat.no_unkozai + 0.5 * temp1 * rt_one_min_ecc_pow2 * sat.con41
			+ 0.0625 * temp2 * rt_one_min_ecc_pow2
			* (13.0 - 78.0 * cos2_incl + 137.0 * cos4_incl);
		
		sat.argpdot = -0.5 * temp1 * one_min_5_cos2_incl + 0.0625 * temp2 *
			(7.0 - 114.0 * cos2_incl + 395.0 * cos4_incl) +
			temp3 * (3.0 - 36.0 * cos2_incl + 49.0 * cos4_incl);
		
		xhdot1 = -temp1 * cos_incl;
		sat.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cos2_incl) +
			2.0 * temp3 * (3.0 - 7.0 * cos2_incl)) * cos_incl;

		xpidot = sat.argpdot + sat.nodedot;
		sat.omgcof = sat.bstar * cc3 * cos(sat.argpo);
		sat.xmcof = 0.0;
		
		if (sat.ecco > 1.0e-4)
			sat.xmcof = -2.0 / 3.0 * coef * sat.bstar / ecc_eta;
		sat.nodecf = 3.5 * one_min_ecc_pow2 * xhdot1 * sat.cc1;
		sat.t2cof = 1.5 * sat.cc1;

		const double num_very_close_to_0 = 1.5e-12;
		if (fabs(cos_incl + 1.0) > num_very_close_to_0)
			sat.xlcof = -0.25 * EARTH_J3_TO_J2 * sin_incl * (3.0 + 5.0 * cos_incl) / (1.0 + cos_incl);
		else
			sat.xlcof = -0.25 * EARTH_J3_TO_J2 * sin_incl * (3.0 + 5.0 * cos_incl) / num_very_close_to_0;
		
		sat.aycof = -0.5 * EARTH_J3_TO_J2 * sin_incl;

		delmotemp = 1.0 + sat.eta * cos(sat.mo);

		sat.delmo = delmotemp * delmotemp * delmotemp;
		sat.sinmao = sin(sat.mo);
		sat.x7thm1 = 7.0 * cos2_incl - 1.0;

		cc1sq = sat.cc1 * sat.cc1;
		sat.d2 = 4.0 * sat.a * tsi * cc1sq;
		
		temp = sat.d2 * tsi * sat.cc1 / 3.0;
		
		sat.d3 = (17.0 * sat.a + s_const) * temp;
		sat.d4 = 0.5 * temp * sat.a * tsi * (221.0 * sat.a + 31.0 * s_const) *
			sat.cc1;
		
		sat.t3cof = sat.d2 + 2.0 * cc1sq;
		sat.t4cof = 0.25 * (3.0 * sat.d3 + sat.cc1 *
			(12.0 * sat.d2 + 10.0 * cc1sq));

		sat.t5cof = 0.2 * (3.0 * sat.d4 +
			12.0 * sat.cc1 * sat.d3 +
			6.0 * sat.d2 * sat.d2 +
			15.0 * cc1sq * (2.0 * sat.d2 + cc1sq));

		sgp4(sat, 0.0, r, v);

		return true;
	}

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
	}
}