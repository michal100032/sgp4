#pragma once

#include "tle.hpp"
#include "vec3.hpp"

#include <chrono>

namespace sgp4 {
	struct state_vecs {
		vec3 position;
		vec3 velocity;
	};

	// TODO: pimpl co to 
	class propagator {
	public:
		// loading tle data
		// refactored sgp4init
		propagator(std::string_view tle_string);
		propagator(const tle_set& set);

		void load(const tle_set& set);
	
		state_vecs run(std::chrono::utc_clock::time_point time);

	private: 
		
		// internal state
		// TODO: rename stuff
		// 	
		// coefficients that can be precalculated
		double c_coef_common = 0.0;
		double c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0, c5 = 0.0;
		double d2 = 0.0, d3 = 0.0, d4 = 0.0;
		double t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0;

		// coefficients that can be precalculated
		// rename!
		double	delmo = 0.0, eta = 0.0, tsi = 0.0, arg_p_dot = 0.0, omgcof = 0.0,
			mean_mot_dot = 0.0, node_dot = 0.0, pre_L_coeff = 0.0, xmcof = 0.0, nodecf = 0;


		// store time differently
		double jdsatepoch = 0.0, jdsatepochF = 0.0;
		std::chrono::utc_clock::time_point epoch;

		double sma_0 = 0.0, nddot = 0.0, ndot = 0.0,
			bstar = 0.0, incl_0 = 0.0, node_0 = 0.0, ecc_0 = 0.0,
			arg_p_0 = 0.0, mean_anom_0 = 0.0;

		double s_const;

		// rename!
		double no_unkozai = 0.0, no_kozai = 0.0;

		// private methods
		double get_epoch_offset_minutes(std::chrono::utc_clock::time_point time);

		void calculate_pre_L_coeff();

		void calculate_c_coeffs();
		void calculate_t_coeffs();
		void calculate_d_coeffs();

		void calculate_angle_coeffs();

		std::tuple<double, double, double> do_long_term_periodics
		(double arg_p, double ecc, double sma, double l, double node);


		std::tuple<double, double, double, double, double> get_perigee_dependent_terms
		(double perigee, double time, double mean_anom_df, double arg_p_df, double node);
	
		std::pair<double, double> find_org_sma_and_mean_mot
		(double tle_mean_motion, double cos2_incl, double rt_one_min_ecc_pow2);;
	};
}