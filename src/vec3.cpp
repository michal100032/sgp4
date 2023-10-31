#include "sgp4/vec3.hpp"

#include <cmath>

sgp4::vec3::vec3() :
	x(0.0), y(0.0), z(0.0) { }

sgp4::vec3::vec3(double x, double y, double z) :
	x(x), y(y), z(z) { }

double sgp4::vec3::magnitude() {
	return sqrt(x * x + y * y + z * z);
}

sgp4::vec3 sgp4::vec3::normalized() {
	double mag = magnitude();
	return { x / mag, y / mag, z / mag };
}

sgp4::vec3 sgp4::vec3::operator+(const vec3& other) {
	return { x + other.x, y + other.y, z + other.z };
}

sgp4::vec3 sgp4::vec3::operator-(const vec3& other) {
	return { x - other.x, y - other.y, z - other.z };;
}

sgp4::vec3 sgp4::vec3::operator*(float scalar) {
	return { x * scalar, y * scalar, z * scalar };;
}

sgp4::vec3 sgp4::vec3::operator-() {
	return { -x, -y, -z };
}

bool sgp4::vec3::operator==(const vec3& other) {
	return x == other.x && y == other.y && z == other.z;
}
