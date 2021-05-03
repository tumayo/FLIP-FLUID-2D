#ifndef PARTICLES_H
#define PARTICLES_H

#include <iostream>
#include "array2.h"
#include "vec.h"
#include "vorticity.h"
#include "util.h"

using namespace std;

//FLIP particles
class Particles {

public:

    int num_particles; // number of particles
    Array2f x, y, u, v; // positions and velocities
	Array2f gamma; // vorticity of each particle
	Array2f sum;   // transfer stuff

	void initialize(float width, int ni, int nj) {
		num_particles = ni * nj; // each grid has one particle
		float dx = width / ni; // dx = width / nj
		gamma = gamma_init(dx, ni, nj);
		u.resize(ni, nj); u.set_zero();
		v.resize(ni, nj); v.set_zero();
		x.resize(ni, nj); x.set_zero();
		y.resize(ni, nj); y.set_zero();
		// Biot-Savart
		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nj; j++) {
				x(i, j) = i * dx;
				y(i, j) = j * dx;
				Vec2f vel = biot_savart(dx, ni, nj, Vec2f(i * dx, j * dx), gamma);
				u(i, j) = vel[0];
				v(i, j) = vel[1];
			}
		}
	}
	void transfer_to_grid(float width, int ni, int nj, Array2f& grid_u, Array2f& grid_v) {
		int i, j;
		float fx, fy;

		grid_u.resize(ni, nj);
		grid_v.resize(ni, nj);
		grid_u.set_zero();
		grid_v.set_zero();
		float dx = width / ni;

		sum.resize(ni, nj);
		sum.set_zero();
		for (int p_i = 0; p_i < sqrt(num_particles); ++p_i) {
			for (int p_j = 0; p_j < sqrt(num_particles); ++p_j) {

				get_barycentric(x(p_i, p_j) / dx, i, fx, 0, ni);
				get_barycentric(y(p_i, p_j) / dx - 0.5f, j, fy, 0, ni);
				accumulate(grid_u, u(p_i, p_j), i, j, fx, fy);
			}
		}
		for (j = 0; j < nj; ++j) {
			for (i = 0; i < ni; ++i) {
				if (sum(i, j) != 0)
					grid_u(i, j) /= sum(i, j);
			}
		}

		sum.set_zero();
		for (int p_i = 0; p_i < sqrt(num_particles); ++p_i) {
			for (int p_j = 0; p_j < sqrt(num_particles); ++p_j) {
				get_barycentric(x(p_i, p_j) / dx - 0.5f, i, fx, 0, ni);
				get_barycentric(y(p_i, p_j) / dx, j, fy, 0, ni);
				accumulate(grid_v, v(p_i, p_j), i, j, fx, fy);
			}
		}
		for (j = 0; j < nj; ++j) {
			for (i = 0; i < ni; ++i) {
				if (sum(i, j) != 0)
					grid_v(i, j) /= sum(i, j);
			}
		}
	}

private:
	void accumulate(Array2f &accum, float q, int i, int j, float fx, float fy) {
		
		float weight;
		weight = (1 - fx) * (1 - fy);
		accum(i, j) += weight * q;
		sum(i, j) += weight;

		weight = fx * (1 - fy);
		accum(i + 1, j) += weight * q;
		sum(i + 1, j) += weight;

		weight = (1 - fx) * fy;
		accum(i, j + 1) += weight * q;
		sum(i, j + 1) += weight;

		weight = fx * fy;
		accum(i + 1, j + 1) += weight * q;
		sum(i + 1, j + 1) += weight;
		
	}
};

#endif
