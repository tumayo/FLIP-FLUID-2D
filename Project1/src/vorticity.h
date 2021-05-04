#ifndef VORTICITY_H
#define VORTICITY_H

#include <iostream>
#include "vec.h"
#include "array2.h"

using namespace std;
const float PI = 3.14159265358979;
const float SIGMA = 0.0001;

// Function definitions

Array2f gamma_init(float dx, int ni, int nj) {
	Array2f gamma;
	gamma.resize(ni, nj);
	gamma.set_zero();
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			float x = i * dx;
			float y = j * dx;
			if (0.4 <= x && x <= 0.7 && 0.4 <= y && y <= 0.7) {
				gamma(i, j) = 2.0 * (y - x);
			}
			else
				gamma(i, j) = 0;
		}
	}
	return gamma;
}
float gauss_distr(Vec2f x) {
	return 1 / (2 * PI * sqr(SIGMA)) * exp(-(sqr(x[0]) + sqr(x[1])) / (2 * sqr(SIGMA)));
}

float q_function(float x) {
	return (1.0f - exp(-0.5 * sqr(x))) / (2.0f * PI);
}

Vec2f biot_savart(float dx, int ni, int nj, Vec2f pos, Array2f vor) {
	Vec2f vel = Vec2f(0, 0);
 	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			Vec2f diff = pos - Vec2f(i * dx, j * dx);
			diff[1] = -diff[1];
 			float distance = sqrt(pow(diff[0], 2) + pow(diff[1], 2));
			if (distance != 0) {
				vel -= diff * vor(i, j) * q_function(distance / SIGMA) / sqr(distance);
			}
		}
	}
	return vel;
}

Array2f vor_inversion_u(float dx, int ni, int nj, Array2f vor) {
	Array2f u;
	u.resize(ni, nj);
	u.set_zero();
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			u(i, j) = biot_savart(dx, ni, nj, Vec2f(i * dx, j * dx), vor)[0];
		}
	}
	return u;
}

Array2f vor_inversion_v(float dx, int ni, int nj, Array2f vor) {
	Array2f u;
	u.resize(ni, nj);
	u.set_zero();
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			u(i, j) = biot_savart(dx, ni, nj, Vec2f(i * dx, j * dx), vor)[1];
		}
	}
	return u;
}

Array2f curl_2D(Array2f u, Array2f v, int ni, int nj, float dx) {
	Array2f vor;
	vor.resize(ni, nj);
	vor.set_zero();
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			//curl : dv/dx - du/dy
			//dv/dx = v(x+dx) - v(x-dx) /2dx
			float dv_dx = 0;
			float du_dy = 0;
			if (i == 0) {
				dv_dx = (v(i + 1, j) - v(i, j)) / dx * 0.5;
			}
			else if (i == ni - 1) {
				dv_dx = (0 - v(i, j)) / dx * 0.5;
			}
			else {
				dv_dx = (v(i + 1, j) - v(i, j)) / dx * 0.5;
			}

			if (j == 0) {
				du_dy = (u(i, j + 1) - u(i, j)) / dx * 0.5;
			}
			else if (j == nj - 1) {
				du_dy = (0 - u(i, j)) / dx * 0.5;
			}
			else {
				du_dy = (u(i, j + 1) - u(i, j)) / dx * 0.5;
			}

			vor(i, j) = dv_dx - du_dy;
		}
	}
	return vor;
}
#endif