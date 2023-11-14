#include "SGP4_mod_near.h"

#include <iostream>
#include <iomanip>
#include <utility>
#include <tuple>
#include <stdexcept>

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
	static double get_julian_since_midnight(double julian) {
		double julian_since_midnight = julian - trunc(julian) - 0.5f;
		if (julian_since_midnight < 0.0)
			julian_since_midnight += 1.0;
		return julian_since_midnight;
	}

	static std::pair<double, double> find_org_sma_and_mean_mot(
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
	// perigee in km
	static double get_s_constant(double perigee) {
		double s_const;
		if (perigee >= 156.0)
			s_const = S_PARAM_DEF;
		else if (perigee < 98.0)
			s_const = 20.0 / EARTH_RAD + 1.0;
		else
			s_const = perigee / EARTH_RAD - S_PARAM_DEF + 2.0;

		return s_const;
	}
	
	static void calculate_c_coeffs(propagation_coeffs& sat, double s_const) {
		double q0_min_s = Q0_PARAM_DEF - s_const;
		double q0_min_s_pow4 = q0_min_s * q0_min_s * q0_min_s * q0_min_s;

		double eta_pow2 = sat.eta * sat.eta;
		double ecc_eta = sat.ecc_0 * sat.eta;

		double q0_min_s_tsi_pow4 = q0_min_s_pow4 * sat.tsi * sat.tsi * sat.tsi * sat.tsi;

		sat.c_coef_common = q0_min_s_tsi_pow4 / pow(1.0 - eta_pow2, 3.5);

		double sin_incl = sin(sat.incl_0);
		double cos_incl = cos(sat.incl_0);
		double three_cos2_incl_min_1 = 3 * cos_incl * cos_incl - 1.0;

		sat.c2 = sat.c_coef_common * sat.no_unkozai * (sat.sma_0 * (1.0 + 1.5 * eta_pow2 + ecc_eta *
			(4.0 + eta_pow2)) + 0.375 * EARTH_J2 * sat.tsi / (1.0 - eta_pow2) * three_cos2_incl_min_1 *
			(8.0 + 3.0 * eta_pow2 * (8.0 + eta_pow2)));

		sat.c1 = sat.bstar * sat.c2;

		sat.c3 = 0.0;
		if (sat.ecc_0 > 1.0e-4)
			sat.c3 = -2.0 * sat.c_coef_common * sat.tsi * EARTH_J3_TO_J2 * sat.no_unkozai
			* sin_incl / sat.ecc_0;

		double one_min_ecc_pow2 = 1.0 - sat.ecc_0 * sat.ecc_0;
		sat.c4 = 2.0 * sat.no_unkozai * sat.c_coef_common * sat.sma_0 * one_min_ecc_pow2 *
			(sat.eta * (2.0 + 0.5 * eta_pow2) + sat.ecc_0 *
				(0.5 + 2.0 * eta_pow2) - EARTH_J2 * sat.tsi / (sat.sma_0 * (1.0 - eta_pow2)) *
				(-3.0 * three_cos2_incl_min_1 * (1.0 - 2.0 * ecc_eta + eta_pow2 *
					(1.5 - 0.5 * ecc_eta)) + 0.75 * (1.0 - cos_incl * cos_incl) *
					(2.0 * eta_pow2 - ecc_eta * (1.0 + eta_pow2)) * cos(2.0 * sat.arg_p_0)));

		sat.c5 = 2.0 * sat.c_coef_common * sat.sma_0 * one_min_ecc_pow2 * (1.0 + 2.75 *
			(eta_pow2 + ecc_eta) + ecc_eta * eta_pow2);
	}

	static double get_pre_L_coeff(propagation_coeffs& sat) {
		double cos_incl = cos(sat.incl_0);
		double sin_incl = sin(sat.incl_0);

		const double num_very_close_to_0 = 1.5e-12;

		if (fabs(cos_incl + 1.0) > num_very_close_to_0)
			return -0.25 * EARTH_J3_TO_J2 * sin_incl * (3.0 + 5.0 * cos_incl) / (1.0 + cos_incl);
		
		return -0.25 * EARTH_J3_TO_J2 * sin_incl * (3.0 + 5.0 * cos_incl) / num_very_close_to_0;
	}

	static void calculate_angle_coeffs(propagation_coeffs& sat) {
		double cos_incl = cos(sat.incl_0);
		double cos2_incl = cos_incl * cos_incl;
		double cos4_incl = cos2_incl * cos2_incl;

		double one_min_ecc_pow2 = 1.0 - sat.ecc_0 * sat.ecc_0;
		double beta_0 = sqrt(one_min_ecc_pow2);

		double posq = sat.sma_0 * one_min_ecc_pow2 * sat.sma_0 * one_min_ecc_pow2;
		double temp1 = 1.5 * EARTH_J2 / posq * sat.no_unkozai;
		double temp2 = 0.5 * temp1 * EARTH_J2 / posq;
		double temp3 = -0.46875 * EARTH_J4 / posq / posq * sat.no_unkozai;

		// M_DF without (t - t_0)
		sat.mean_mot_dot = sat.no_unkozai + 0.5 * temp1 * beta_0 * (3 * cos2_incl - 1.0)
			+ 0.0625 * temp2 * beta_0
			* (13.0 - 78.0 * cos2_incl + 137.0 * cos4_incl);

		// (small) omega_DF
		sat.arg_p_dot = -0.5 * temp1 * (1.0 - 5.0 * cos2_incl) + 0.0625 * temp2 *
			(7.0 - 114.0 * cos2_incl + 395.0 * cos4_incl) +
			temp3 * (3.0 - 36.0 * cos2_incl + 49.0 * cos4_incl);


		double xhdot1 = -temp1 * cos_incl;
		// page 12
		// (capital) omega_DF
		sat.node_dot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cos2_incl) +
			2.0 * temp3 * (3.0 - 7.0 * cos2_incl)) * cos_incl;

		// delta small omega (again without t - t_0)
		sat.omgcof = sat.bstar * sat.c3 * cos(sat.arg_p_0);

		// delta M without a lot of terms
		double ecc_eta = sat.ecc_0 * sat.eta;
		sat.xmcof = 0.0;
		if (sat.ecc_0 > 1.0e-4)
			sat.xmcof = -2.0 / 3.0 * sat.c_coef_common * sat.bstar / ecc_eta;

		// this might be the capital omega (without t - t_0) and omega_DF
		sat.nodecf = 3.5 * one_min_ecc_pow2 * xhdot1 * sat.c1;
	}

	static void calculate_d_coeffs(propagation_coeffs& sat, double s_const) {
		double c1_pow2 = sat.c1 * sat.c1;
		sat.d2 = 4.0 * sat.sma_0 * sat.tsi * c1_pow2;

		double d3_d4_common = sat.d2 * sat.tsi * sat.c1 / 3.0;

		sat.d3 = (17.0 * sat.sma_0 + s_const) * d3_d4_common;
		sat.d4 = 0.5 * d3_d4_common * sat.sma_0 * sat.tsi
			* (221.0 * sat.sma_0 + 31.0 * s_const) * sat.c1;
	}

	static void calculate_t_coeffs(propagation_coeffs& sat) {
		double c1_pow2 = sat.c1 * sat.c1;
		sat.t2 = 1.5 * sat.c1;
		
		sat.t3 = sat.d2 + 2.0 * c1_pow2;

		sat.t4 = 0.25 * (3.0 * sat.d3 + sat.c1 *
			(12.0 * sat.d2 + 10.0 * c1_pow2));

		sat.t5 = 0.2 * (3.0 * sat.d4 +
			12.0 * sat.c1 * sat.d3 +
			6.0 * sat.d2 * sat.d2 +
			15.0 * c1_pow2 * (2.0 * sat.d2 + c1_pow2));
	}

	propagation_coeffs sgp4init(const sgp4::tle_set& set) {
		propagation_coeffs sat;
		
		/* Store time differently
		double julian_epoch = sgp4::time_utils::to_julian(set.epoch);
		sat.jdsatepoch = julian_epoch;
		sat.jdsatepochF = get_julian_since_midnight(julian_epoch);
		*/
		
		sat.no_kozai = set.mean_motion;
		sat.ndot = set.d_mean_motion;
		sat.nddot = set.dd_mean_motion;

		sat.incl_0 = set.inclination;
		sat.node_0 = set.right_ascension;
		sat.arg_p_0 = set.arg_of_perigee;
		sat.mean_anom_0 = set.mean_anomaly;

		sat.ecc_0 = set.eccentricity;

		sat.bstar = set.rad_press_coef;

		double ecc_pow2 = sat.ecc_0 * sat.ecc_0;

		double one_min_ecc_pow2 = 1.0 - ecc_pow2;
		double beta_0 = sqrt(one_min_ecc_pow2);
		// put these in elsetrec?
		double cos_incl = cos(sat.incl_0);
		double sin_incl = sin(sat.incl_0);
		double cos2_incl = cos_incl * cos_incl;

		auto [sma, mean_motion] = find_org_sma_and_mean_mot(sat.no_kozai, cos2_incl, beta_0);
		sat.no_unkozai = mean_motion;
		sat.sma_0 = sma;

		double perigee = (sat.sma_0 * (1.0 - sat.ecc_0) - 1.0) * EARTH_RAD;

		double s_const = get_s_constant(perigee);
		
		sat.tsi = 1.0 / (sat.sma_0 - s_const);
		sat.eta = sat.sma_0 * sat.ecc_0 * sat.tsi;

		calculate_c_coeffs(sat, s_const);

		calculate_angle_coeffs(sat);
				
		sat.pre_L_coeff = get_pre_L_coeff(sat);

		double one_p_eta_cos_mo = 1.0 + sat.eta * cos(sat.mean_anom_0);
		sat.delmo = one_p_eta_cos_mo * one_p_eta_cos_mo * one_p_eta_cos_mo;

		calculate_d_coeffs(sat, s_const);
		
		calculate_t_coeffs(sat);
		
		return sat;
	}

	static double solve_keplers_equation(double u, double a_xn, double a_yn) {
		u = fmod(u, TWO_PI);
		double eo = u;

		double delta_eo = 1.0;

		//   Vallando's comment:
		//   the following iteration needs better limits on corrections

		const int MAX_KEPLER_ITERATIONS = 10;
		for (int i = 0; i < MAX_KEPLER_ITERATIONS; i++) {
			double sin_eo = sin(eo);
			double cos_eo = cos(eo);

			delta_eo = (u - a_yn * cos_eo + a_xn * sin_eo - eo)
				/ (1.0 - cos_eo * a_xn - sin_eo * a_yn);

			if (delta_eo >= 0.95)
				delta_eo = 0.95;
			else if (delta_eo <= -0.95)
				delta_eo = -0.95;

			eo += delta_eo;

			if (fabs(delta_eo) < 1.0e-12)
				return eo;
		}
		return eo;
	}

	state_vecs sgp4(propagation_coeffs& sat, double time) {
		double temp, temp1, temp2;

		/* ------------------ set mathematical constants --------------- */
		// sgp4fix divisor for divide by zero check on inclination
		// the old check used 1.0 + cos(PI-1.0e-9), but then compared it to
		// 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
		
		/* ------- update for secular gravity and atmospheric drag ----- */
		// M_DF
		double mean_anom_df = sat.mean_anom_0 + sat.mean_mot_dot * time;
		// arg_DF
		double arg_p_df = sat.arg_p_0 + sat.arg_p_dot * time;
		// right_ascencion_DF
		double node_df = sat.node_0 + sat.node_dot * time;
		
		double t_pow2 = time * time;	
		
		// capital omega
		double node = node_df + sat.nodecf * t_pow2;
		
		// the whole if else block must be removed
		// params that depend on perigee
		// M_p, arg, ecc, sma, L (mm is both M_p and L?)
		double mean_anom, arg_p, sma, ecc, l;
		
		double perigee = (sat.sma_0 * (1.0 - sat.ecc_0) - 1.0) * EARTH_RAD;
		if (perigee > 220.0) {
			double del_arg = sat.omgcof * time;

			double temp_del_m = 1.0 + sat.eta * cos(mean_anom_df);
			double del_m = sat.xmcof *
				(temp_del_m * temp_del_m * temp_del_m - sat.delmo);

			// M
			mean_anom = mean_anom_df + del_arg + del_m;
			arg_p = arg_p_df - del_arg - del_m;
			
			double t_pow3 = t_pow2 * time;
			double t_pow4 = t_pow3 * time;

			double temp_a = 1.0 - sat.c1 * time - sat.d2 * t_pow2 - sat.d3 * t_pow3 -
				sat.d4 * t_pow4;
			sma = sat.sma_0 * temp_a * temp_a;

			double temp_e = sat.bstar * sat.c4 * time + sat.bstar * sat.c5
				* (sin(mean_anom) - sin(sat.mean_anom_0));
			ecc = sat.ecc_0 - temp_e;

			double temp_l = sat.t2 * t_pow2 + sat.t3 * t_pow3 + t_pow4 * (sat.t4 +
				time * sat.t5);
		
			l = mean_anom + arg_p + node + sat.no_unkozai * temp_l;
		} 
		else 
		{

			double temp_a = 1.0 - sat.c1 * time;
			sma = sat.sma_0 * temp_a * temp_a;
		
			ecc = sat.ecc_0 - sat.bstar * sat.c4 * time;
			
			mean_anom = mean_anom_df;
			arg_p = arg_p_df;

			l = mean_anom + arg_p + node + sat.no_unkozai * sat.t2 * t_pow2;
		}
		double nm = X_KE / pow(sma, 1.5);

		/*

		if (sat.no_unkozai <= 0.0) {
		//	sat.error = 2;
			return false;
		}

		if (ecc >= 1.0 || ecc < -0.001) {
			sat.error = 1;
			return false;
		}
		*/
		

		// sgp4fix fix tolerance to avoid a divide by zero
		if (ecc < 1.0e-6)
			ecc = 1.0e-6;

		double sin_incl = sin(sat.incl_0);
		double cos_incl = cos(sat.incl_0);

		/* -------------------- long period periodics ------------------ */
		double a_xn = ecc * cos(arg_p);

		temp = 1.0 / (sma * (1.0 - ecc * ecc));
		double a_yn = ecc * sin(arg_p) + temp * -0.5 * EARTH_J3_TO_J2 * sin_incl;

		double l_t = l + temp * sat.pre_L_coeff * a_xn;

		double ecc_anom = solve_keplers_equation(l_t - node, a_xn, a_yn);

		double sin_ecc_anom = sin(ecc_anom), cos_ecc_anom = cos(ecc_anom);


		/* ------------- short period preliminary quantities ----------- */
		double ecc_cos_e = a_xn * cos_ecc_anom + a_yn * sin_ecc_anom;
		double ecc_sin_e = a_xn * sin_ecc_anom - a_yn * cos_ecc_anom;
		double pl = sma * (1.0 - a_xn * a_xn + a_yn * a_yn);
		
		/*
		if (pl < 0.0) {
			sat.error = 4;
			return false;
		}*/

		double rl = sma * (1.0 - ecc_cos_e);
		double rdotl = sqrt(sma) * ecc_sin_e / rl;
		double rvdotl = sqrt(pl) / rl;
		double betal = sqrt(1.0 - a_xn * a_xn + a_yn * a_yn);
		temp = ecc_sin_e / (1.0 + betal);
		double sin_u = sma / rl * (sin_ecc_anom - a_yn - a_xn * temp);
		double cos_u = sma / rl * (cos_ecc_anom - a_xn + a_yn * temp);
		double sin_2u =  2 * cos_u * sin_u;
		double cos_2u = 1.0 - 2.0 * sin_u * sin_u;
		temp = 1.0 / pl;
		temp1 = 0.5 * EARTH_J2 * temp;
		temp2 = temp1 * temp;

		double cos2_incl = cos_incl * cos_incl;
		double three_cos2_incl_min_1 = 3 * cos2_incl - 1.0;
		
		double mrt = rl * (1.0 - 1.5 * temp2 * betal * three_cos2_incl_min_1) +
			0.5 * temp1 * (1.0 - cos2_incl) * cos_2u;
		double u_k = atan2(sin_u, cos_u) - 0.25 * temp2 * (7.0 * cos2_incl - 1.0) * sin_2u;
		
		double node_k = node + 1.5 * temp2 * cos_incl * sin_2u;
		
		double incl_k = sat.incl_0 + 1.5 * temp2 * cos_incl * sin_incl * cos_2u;


		double mvt = rdotl - nm * temp1 * (1.0 - cos2_incl) * sin_2u / X_KE;
		double rvdot = rvdotl + nm * temp1 * ((1.0 - cos2_incl) * cos_2u +
			1.5 * three_cos2_incl_min_1) / X_KE;

		/* --------------------- orientation vectors ------------------- */
		
		double sin_u_k = sin(u_k);
		double cos_u_k = cos(u_k);
		
		double sin_node_k = sin(node_k),  cos_node_k = cos(node_k);
		double sin_i_k = sin(incl_k), cos_incl_k = cos(incl_k);
		
		double xmx = -sin_node_k * cos_incl_k;
		double xmy = cos_node_k * cos_incl_k;

		sgp4::vec3 u = {
			xmx * sin_u_k + cos_node_k * cos_u_k,
			xmy * sin_u_k + sin_node_k * cos_u_k,
			sin_i_k * sin_u_k
		};
		sgp4::vec3 v = {
			xmx * cos_u_k - cos_node_k * sin_u_k,
			xmy * cos_u_k - sin_node_k * sin_u_k,
			sin_i_k * cos_u_k
		};
		/* --------- position and velocity (in km and km/sec) ---------- */
		/*
		// orbit decayed?
		if (mrt < 1.0) {
			sat.error = 6;
			return false;
		}*/
		const double vkmpersec = EARTH_J2 * X_KE / 60.0;
		state_vecs vecs;
		vecs.position = u * mrt * EARTH_RAD; 
		vecs.velocity = (u * mvt + v * rvdot) * vkmpersec;
		
		return vecs;
	}
}