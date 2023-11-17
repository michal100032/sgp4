#include "sgp4/propagator.hpp"

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

// KE constatn
static const double EARTH_KE = 60.0 / sqrt(EARTH_RAD * EARTH_RAD * EARTH_RAD / EARTH_GRV);

// earth zonal harmonic model
static const double EARTH_J2 = 0.00108262998905;
static const double EARTH_J3 = -0.00000253215306;
static const double EARTH_J4 = -0.00000161098761;
static const double EARTH_J3_TO_J2 = EARTH_J3 / EARTH_J2;

// other SGP4 parameters

// default value of the s constant
static const double S_PARAM_DEF = 78.0 / EARTH_RAD + 1.0;
static const double Q0_PARAM_DEF = 120.0 / EARTH_RAD + 1.0;

sgp4::propagator::propagator(std::string_view tle_string) {
	tle_set tle = sgp4::parse_tle_entry(tle_string);
	load(tle);
}

sgp4::propagator::propagator(const tle_set& set) {
	load(set);
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

void sgp4::propagator::load(const tle_set& set) {
	/* Store time differently
	double julian_epoch = sgp4::time_utils::to_julian(set.epoch);
	jdsatepoch = julian_epoch;
	jdsatepochF = get_julian_since_midnight(julian_epoch);
	*/

	no_kozai = set.mean_motion;
	ndot = set.d_mean_motion;
	nddot = set.dd_mean_motion;

	incl_0 = set.inclination;
	node_0 = set.right_ascension;
	arg_p_0 = set.arg_of_perigee;
	mean_anom_0 = set.mean_anomaly;

	ecc_0 = set.eccentricity;

	bstar = set.rad_press_coef;

	double ecc_pow2 = ecc_0 * ecc_0;

	double one_min_ecc_pow2 = 1.0 - ecc_pow2;
	double beta_0 = sqrt(one_min_ecc_pow2);
	// store these as well?
	double cos_incl = cos(incl_0);
	double sin_incl = sin(incl_0);
	double cos2_incl = cos_incl * cos_incl;

	auto [sma, mean_motion] = find_org_sma_and_mean_mot(no_kozai, cos2_incl, beta_0);
	no_unkozai = mean_motion;
	sma_0 = sma;

	double perigee = (sma_0 * (1.0 - ecc_0) - 1.0) * EARTH_RAD;

	s_const = get_s_constant(perigee);

	tsi = 1.0 / (sma_0 - s_const);
	eta = sma_0 * ecc_0 * tsi;

	calculate_c_coeffs();

	calculate_angle_coeffs();

	calculate_pre_L_coeff();

	double one_p_eta_cos_mo = 1.0 + eta * cos(mean_anom_0);
	delmo = one_p_eta_cos_mo * one_p_eta_cos_mo * one_p_eta_cos_mo;

	calculate_d_coeffs();

	calculate_t_coeffs();
}

static double get_julian_since_midnight(double julian) {
	double julian_since_midnight = julian - trunc(julian) - 0.5f;
	if (julian_since_midnight < 0.0)
		julian_since_midnight += 1.0;
	return julian_since_midnight;
}

std::pair<double, double> sgp4::propagator::find_org_sma_and_mean_mot(
	double tle_mean_motion, double cos2_incl, double rt_one_min_ecc_pow2
) {
	double a_1 = pow(EARTH_KE / tle_mean_motion, 2.0 / 3.0);
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

void sgp4::propagator::calculate_c_coeffs() {
	double q0_min_s = Q0_PARAM_DEF - s_const;
	double q0_min_s_pow4 = q0_min_s * q0_min_s * q0_min_s * q0_min_s;

	double eta_pow2 = eta * eta;
	double ecc_eta = ecc_0 * eta;

	double q0_min_s_tsi_pow4 = q0_min_s_pow4 * tsi * tsi * tsi * tsi;

	c_coef_common = q0_min_s_tsi_pow4 / pow(1.0 - eta_pow2, 3.5);

	double sin_incl = sin(incl_0), cos_incl = cos(incl_0);
	double three_cos2_incl_min_1 = 3 * cos_incl * cos_incl - 1.0;

	c2 = c_coef_common * no_unkozai * (sma_0 * (1.0 + 1.5 * eta_pow2 + ecc_eta *
		(4.0 + eta_pow2)) + 0.375 * EARTH_J2 * tsi / (1.0 - eta_pow2) * three_cos2_incl_min_1 *
		(8.0 + 3.0 * eta_pow2 * (8.0 + eta_pow2)));

	c1 = bstar * c2;

	c3 = 0.0;
	if (ecc_0 > 1.0e-4)
		c3 = -2.0 * c_coef_common * tsi * EARTH_J3_TO_J2 * no_unkozai
		* sin_incl / ecc_0;

	double one_min_ecc_pow2 = 1.0 - ecc_0 * ecc_0;
	c4 = 2.0 * no_unkozai * c_coef_common * sma_0 * one_min_ecc_pow2 *
		(eta * (2.0 + 0.5 * eta_pow2) + ecc_0 *
			(0.5 + 2.0 * eta_pow2) - EARTH_J2 * tsi / (sma_0 * (1.0 - eta_pow2)) *
			(-3.0 * three_cos2_incl_min_1 * (1.0 - 2.0 * ecc_eta + eta_pow2 *
				(1.5 - 0.5 * ecc_eta)) + 0.75 * (1.0 - cos_incl * cos_incl) *
				(2.0 * eta_pow2 - ecc_eta * (1.0 + eta_pow2)) * cos(2.0 * arg_p_0)));

	c5 = 2.0 * c_coef_common * sma_0 * one_min_ecc_pow2 * (1.0 + 2.75 *
		(eta_pow2 + ecc_eta) + ecc_eta * eta_pow2);
}

void sgp4::propagator::calculate_pre_L_coeff() {
	double cos_incl = cos(incl_0);
	double sin_incl = sin(incl_0);

	const double num_very_close_to_0 = 1.5e-12;

	if (fabs(cos_incl + 1.0) > num_very_close_to_0) {
		pre_L_coeff = -0.25 * EARTH_J3_TO_J2 * sin_incl
			* (3.0 + 5.0 * cos_incl) / (1.0 + cos_incl);
	} else {
		pre_L_coeff = -0.25 * EARTH_J3_TO_J2 * sin_incl
			* (3.0 + 5.0 * cos_incl) / num_very_close_to_0;
	}
}

void sgp4::propagator::calculate_angle_coeffs() {
	double cos_incl = cos(incl_0);
	double cos2_incl = cos_incl * cos_incl;
	double cos4_incl = cos2_incl * cos2_incl;

	double one_min_ecc_pow2 = 1.0 - ecc_0 * ecc_0;
	double beta_0 = sqrt(one_min_ecc_pow2);

	double posq = sma_0 * one_min_ecc_pow2 * sma_0 * one_min_ecc_pow2;
	double temp1 = 1.5 * EARTH_J2 / posq * no_unkozai;
	double temp2 = 0.5 * temp1 * EARTH_J2 / posq;
	double temp3 = -0.46875 * EARTH_J4 / posq / posq * no_unkozai;

	// M_DF without (t - t_0)
	mean_mot_dot = no_unkozai + 0.5 * temp1 * beta_0 * (3 * cos2_incl - 1.0)
		+ 0.0625 * temp2 * beta_0
		* (13.0 - 78.0 * cos2_incl + 137.0 * cos4_incl);

	// (small) omega_DF
	arg_p_dot = -0.5 * temp1 * (1.0 - 5.0 * cos2_incl) + 0.0625 * temp2 *
		(7.0 - 114.0 * cos2_incl + 395.0 * cos4_incl) +
		temp3 * (3.0 - 36.0 * cos2_incl + 49.0 * cos4_incl);


	double xhdot1 = -temp1 * cos_incl;
	// page 12
	// (capital) omega_DF
	node_dot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cos2_incl) +
		2.0 * temp3 * (3.0 - 7.0 * cos2_incl)) * cos_incl;

	// delta small omega (again without t - t_0)
	omgcof = bstar * c3 * cos(arg_p_0);

	// delta M without a lot of terms
	double ecc_eta = ecc_0 * eta;
	xmcof = 0.0;
	if (ecc_0 > 1.0e-4)
		xmcof = -2.0 / 3.0 * c_coef_common * bstar / ecc_eta;

	// this might be the capital omega (without t - t_0) and omega_DF
	nodecf = 3.5 * one_min_ecc_pow2 * xhdot1 * c1;
}

void sgp4::propagator::calculate_d_coeffs() {
	double c1_pow2 = c1 * c1;
	d2 = 4.0 * sma_0 * tsi * c1_pow2;

	double d3_d4_common = d2 * tsi * c1 / 3.0;

	d3 = (17.0 * sma_0 + s_const) * d3_d4_common;
	d4 = 0.5 * d3_d4_common * sma_0 * tsi
		* (221.0 * sma_0 + 31.0 * s_const) * c1;
}

void sgp4::propagator::calculate_t_coeffs() {
	double c1_pow2 = c1 * c1;
	t2 = 1.5 * c1;

	t3 = d2 + 2.0 * c1_pow2;

	t4 = 0.25 * (3.0 * d3 + c1 *
		(12.0 * d2 + 10.0 * c1_pow2));

	t5 = 0.2 * (3.0 * d4 +
		12.0 * c1 * d3 +
		6.0 * d2 * d2 +
		15.0 * c1_pow2 * (2.0 * d2 + c1_pow2));
}

static double solve_keplers_equation(double u, double a_xn, double a_yn) {
	u = fmod(u, TWO_PI);
	double eo = u;

	double delta_eo = 1.0;

	//   Vallando's comment:
	//   the following iteration needs better limits on corrections

	const int MAEARTH_KEPLER_ITERATIONS = 10;
	for (int i = 0; i < MAEARTH_KEPLER_ITERATIONS; i++) {
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

std::tuple<double, double, double, double, double> sgp4::propagator::get_perigee_dependent_terms
(double perigee, double time, double mean_anom_df, double arg_p_df, double node) {
	double mean_anom, arg_p, ecc, sma, l;
	double t_pow2 = time * time;
	if (perigee > 220.0) {
		double del_arg = omgcof * time;

		double temp_del_m = 1.0 + eta * cos(mean_anom_df);
		double del_m = xmcof *
			(temp_del_m * temp_del_m * temp_del_m - delmo);

		// M
		mean_anom = mean_anom_df + del_arg + del_m;
		arg_p = arg_p_df - del_arg - del_m;

		double t_pow3 = t_pow2 * time;
		double t_pow4 = t_pow3 * time;

		double temp_a = 1.0 - c1 * time - d2 * t_pow2 - d3 * t_pow3 -
			d4 * t_pow4;
		sma = sma_0 * temp_a * temp_a;

		double temp_e = bstar * c4 * time + bstar * c5
			* (sin(mean_anom) - sin(mean_anom_0));
		ecc = ecc_0 - temp_e;

		double temp_l = t2 * t_pow2 + t3 * t_pow3 + t_pow4 * (t4 +
			time * t5);

		l = mean_anom + arg_p + node + no_unkozai * temp_l;
	}
	else {
		double temp_a = 1.0 - c1 * time;
		sma = sma_0 * temp_a * temp_a;

		ecc = ecc_0 - bstar * c4 * time;

		mean_anom = mean_anom_df;
		arg_p = arg_p_df;

		l = mean_anom + arg_p + node + no_unkozai * t2 * t_pow2;
	}

	return { mean_anom, arg_p, ecc, sma, l };
}

// returns ecc_anom + arg_p, a_xn, a_yn
std::tuple<double, double, double> sgp4::propagator::do_long_term_periodics
(double arg_p, double ecc, double sma, double l, double node) {
	double a_xn = ecc * cos(arg_p);

	double temp = 1.0 / (sma * (1.0 - ecc * ecc));
	double a_yn = ecc * sin(arg_p) + temp * -0.5 * EARTH_J3_TO_J2 * sin(incl_0);

	double l_t = l + temp * pre_L_coeff * a_xn;

	double e = solve_keplers_equation(l_t - node, a_xn, a_yn);

	return { e, a_xn, a_yn };
}

static double get_time_from_epoch_mins(std::chrono::utc_clock::time_point time) {
	return 0.0;
}

sgp4::state_vecs sgp4::propagator::run(double time) {
	
	double mean_anom_df = mean_anom_0 + mean_mot_dot * time;
	double arg_p_df = arg_p_0 + arg_p_dot * time;
	double node_df = node_0 + node_dot * time;

	double t_pow2 = time * time;
	double perigee = (sma_0 * (1.0 - ecc_0) - 1.0) * EARTH_RAD;

	double node = node_df + nodecf * t_pow2;
	auto [mean_anom, arg_p, ecc, sma, l]
		= get_perigee_dependent_terms(perigee, time, mean_anom_df, arg_p_df, node);

	double mean_mot = EARTH_KE / pow(sma, 1.5);

	if (no_unkozai <= 0.0)
		throw std::logic_error("sgp4 error 2 (negative mean motion)");

	if (ecc >= 1.0 || ecc < -0.001)
		throw std::logic_error("sgp4 error 1 (bad eccentricity)");

	// fix tolerance to avoid a divide by zero
	if (ecc < 1.0e-6)
		ecc = 1.0e-6;

	/* -------------------- long period periodics ------------------ */
	auto [ecc_anom_arg, a_xn, a_yn] = do_long_term_periodics(arg_p, ecc, sma, l, node);
	double sin_ecc_anom = sin(ecc_anom_arg), cos_ecc_anom = cos(ecc_anom_arg);

	/* ------------- short period preliminary quantities ----------- */
	double ecc_cos_e = a_xn * cos_ecc_anom + a_yn * sin_ecc_anom;
	double ecc_sin_e = a_xn * sin_ecc_anom - a_yn * cos_ecc_anom;
	// semi latus rectum?
	double p_l = sma * (1.0 - a_xn * a_xn + a_yn * a_yn);

	if (p_l < 1.0)
		throw std::logic_error("sgp4 error 4 (bad p_l)");

	// 
	double r = sma * (1.0 - ecc_cos_e);
	double r_dot = sqrt(sma) * ecc_sin_e / r;
	double r_f_dot = sqrt(p_l) / r;

	double beta_l = sqrt(1.0 - a_xn * a_xn + a_yn * a_yn);
	;
	double sin_u = sma / r * (sin_ecc_anom - a_yn - a_xn * ecc_sin_e / (1.0 + beta_l));
	double cos_u = sma / r * (cos_ecc_anom - a_xn + a_yn * ecc_sin_e / (1.0 + beta_l));

	double sin_2u = 2 * cos_u * sin_u;
	double cos_2u = 1.0 - 2.0 * sin_u * sin_u;

	double half_j2_over_pl_pow2 = 0.5 * EARTH_J2 / p_l / p_l;

	double sin_incl = sin(incl_0), cos_incl = cos(incl_0);
	double cos2_incl = cos_incl * cos_incl;
	double three_cos2_incl_min_1 = 3 * cos2_incl - 1.0;

	double r_k = r * (1.0 - 1.5 * half_j2_over_pl_pow2 * beta_l * three_cos2_incl_min_1) +
		0.25 * EARTH_J2 / p_l * (1.0 - cos2_incl) * cos_2u;

	double u_k = atan2(sin_u, cos_u) - 0.25 * half_j2_over_pl_pow2 * (7.0 * cos2_incl - 1.0) * sin_2u;
	double node_k = node + 1.5 * half_j2_over_pl_pow2 * cos_incl * sin_2u;
	double incl_k = incl_0 + 1.5 * half_j2_over_pl_pow2 * cos_incl * sin_incl * cos_2u;


	// fuck
	double r_k_dot = r_dot - mean_mot * 0.5 * EARTH_J2 / p_l * (1.0 - cos2_incl) * sin_2u / EARTH_KE;
	double r_f_k_dot = r_f_dot + mean_mot * 0.5 * EARTH_J2 / p_l * ((1.0 - cos2_incl) * cos_2u +
		1.5 * three_cos2_incl_min_1) / EARTH_KE;

	/* --------------------- orientation vectors ------------------- */

	double sin_u_k = sin(u_k), cos_u_k = cos(u_k);
	double sin_node_k = sin(node_k), cos_node_k = cos(node_k);
	double sin_i_k = sin(incl_k), cos_incl_k = cos(incl_k);

	double m_x = -sin_node_k * cos_incl_k;
	double m_y = cos_node_k * cos_incl_k;

	sgp4::vec3 u = {
		m_x * sin_u_k + cos_node_k * cos_u_k,
		m_y * sin_u_k + sin_node_k * cos_u_k,
		sin_i_k * sin_u_k
	};
	sgp4::vec3 v = {
		m_x * cos_u_k - cos_node_k * sin_u_k,
		m_y * cos_u_k - sin_node_k * sin_u_k,
		sin_i_k * cos_u_k
	};
	/* --------- position and velocity (in km and km/sec) ---------- */
	if (r_k < 1.0)
		throw std::logic_error("sgp4 error 6 (orbit decayed)");

	// why EARTH_KE??
	const double vkmpersec = EARTH_RAD * EARTH_KE / 60.0;
	state_vecs vecs;
	vecs.position = u * r_k * EARTH_RAD;
	vecs.velocity = (u * r_k_dot + v * r_f_k_dot) * vkmpersec;

	return vecs;
}