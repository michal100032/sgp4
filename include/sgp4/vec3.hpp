#pragma once

namespace sgp4 {
	struct vec3 {
		double x, y, z;

		vec3();
		vec3(double x, double y, double z);

		double magnitude();
		vec3 normalized();
		vec3 operator+(const vec3& other);
		vec3 operator-(const vec3& other);
		vec3 operator*(float scalar);
		vec3 operator-();

		bool operator==(const vec3& other);
	};
}