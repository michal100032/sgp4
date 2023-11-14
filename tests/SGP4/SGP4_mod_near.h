#pragma once

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "sgp4/sgp4.hpp"

namespace SGP4Funcs_mod_near {
	// units: rad, er, minutes
	struct propagation_coeffs {
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

		double sma_0 = 0.0, nddot = 0.0, ndot = 0.0,
			bstar = 0.0, incl_0 = 0.0, node_0 = 0.0, ecc_0 = 0.0,
			arg_p_0 = 0.0, mean_anom_0 = 0.0;

		// rename!
		double no_unkozai = 0.0, no_kozai = 0.0;
	};

	propagation_coeffs sgp4init(const sgp4::tle_set& set);

	struct state_vecs {
		sgp4::vec3 position;
		sgp4::vec3 velocity;
	};

	state_vecs sgp4(propagation_coeffs& coeffs, double time);
}