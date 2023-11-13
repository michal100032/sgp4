#pragma once

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "sgp4/sgp4.hpp"

namespace SGP4Funcs_mod_near {
	struct elsetrec {
		// put elsewhere?
		int error = 0;

		// coefficients that can be precalculated
		double c_coef_common = 0.0;
		double c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0, c5 = 0.0;
		double d2 = 0.0, d3 = 0.0, d4 = 0.0;
		double t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0;
		
		// coefficients that can be precalculated
		// rename!
		double	delmo = 0.0, eta = 0.0, tsi = 0.0, argpdot = 0.0, omgcof = 0.0, sinmao = 0.0, 
			mean_mot_dot = 0.0, node_dot = 0.0, pre_L_coeff = 0.0, xmcof = 0.0, nodecf = 0;

		//
		// from TLE: 
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

		// a - semi-major axis, altp - perigee, alta - apogee
		// rename!
		double a = 0.0,
			jdsatepoch = 0.0, jdsatepochF = 0.0, nddot = 0.0, ndot = 0.0,
			bstar = 0.0, incl = 0.0, nodeo = 0.0, ecc = 0.0,
			argpo = 0.0, mo = 0.0, no_kozai = 0.0;

		// rename!
		double no_unkozai = 0.0;
	};

	elsetrec sgp4init(const sgp4::tle_set& set);

	bool sgp4(elsetrec& satrec, double tsince, double r[3], double v[3]);
}