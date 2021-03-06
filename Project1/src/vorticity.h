#ifndef VORTICITY_H
#define VORTICITY_H

#include <iostream>
#include "vec.h"
#include "array2.h"

using namespace std;
const float PI = 3.14159265358979f;
const float SIGMA = 0.0001f;

// Function definitions

Array2f gamma_init(float dx, int ni, int nj) {
	Array2f gamma;
	gamma.resize(ni+1, nj+1);
	gamma.set_zero();
	for (int i = 0; i < ni+1; i++) {
		for (int j = 0; j < nj+1; j++) {
			float x = i * dx;
			float y = j * dx;
			if (0.2f <= x && x <= 0.6f && 0.2f <= y && y <= 0.6f) {
				gamma(i, j) = 10.0f * (y - x);
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
 	for (int i = 0; i < ni+1; i++) {
		for (int j = 0; j < nj+1; j++) {
			Vec2f diff = pos - Vec2f(i * dx, j * dx);
 			float distance = sqrt(pow(diff[0], 2) + pow(diff[1], 2));
			if (distance != 0) {
				// Particle Biot-Savart Method
				//diff[1] = -diff[1];
				//vel -= diff * vor(i, j) * q_function(distance / SIGMA) / sqr(distance);

				// MC Fluids Biot-Savart Method
				vel += Vec2f(diff[1], -diff[0]) * vor(i, j) / sqr(distance);
				
			}
		}
	}
	//return vel;
	return  vel / (ni * nj * 2.0f * PI);
}

Array2f vor_inversion_u(float dx, int ni, int nj, Array2f vor) {
	Array2f u;
	u.resize(ni, nj);
	u.set_zero();
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			u(i, j) = biot_savart(dx, ni, nj, Vec2f(i * dx, (j + 0.5f) * dx), vor)[0];
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
			u(i, j) = biot_savart(dx, ni, nj, Vec2f((i + 0.5f) * dx , j * dx), vor)[1];
		}
	}
	return u;
}

Array2f curl_2D(Array2f u, Array2f v, int ni, int nj, float dx) {
	Array2f vor;
	vor.resize(ni+1, nj+1);
	vor.set_zero();
	for (int i = 0; i < ni+1; i++) {
		for (int j = 0; j < nj+1; j++) {
			//curl : dv/dx - du/dy
			//dv/dx = v(x+0.5*dx) - v(x-0.5*dx) /dx
			float dv_dx = 0;
			float du_dy = 0;
			if (i == 0) {
				dv_dx = (v(i + 1, j) - v(i, j)) / dx;
			}
			else if (i == ni) {
				dv_dx = (v(i - 1, j) - v(i - 2, j)) / dx;
			}
			else {
				dv_dx = (v(i, j) - v(i - 1, j)) / dx;
			}

			if (j == 0) {
				du_dy = (u(i, j + 1) - u(i, j)) / dx;
			}
			else if (j == nj) {
				du_dy = (u(i, j - 1) - u(i, j - 2)) / dx;
			}
			else {
				du_dy = (u(i, j) - u(i, j - 1)) / dx;
			}

			vor(i, j) = -dv_dx + du_dy;
		}
	}
	return vor;
}
#endif
