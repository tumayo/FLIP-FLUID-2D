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
	void move_particles_in_grid(float dt, Array2f grid_u, Array2f grid_v, float xmax, float ymax, float dx) {
		Vec2f midx, gu;
		for (int p_i = 0; p_i < sqrt(num_particles); ++p_i) {
			for (int p_j = 0; p_j < sqrt(num_particles); ++p_j) {
				// first stage of Runge-Kutta 2 (do a half Euler step)
				float x_pos = x(p_i, p_j) / dx;
				float y_pos = y(p_i, p_j) / dx;
				int v00 = floor(x_pos); int v10 = ceil(x_pos);
				int v01 = floor(y_pos); int v11 = ceil(y_pos);
				float fx = x_pos - v00;
				if (fx > v10 - x_pos)
					fx = v10 - x_pos;
				float fy = y_pos - v01;
				if (fy > v11 - y_pos)
					fy = v11 - y_pos;
				gu[0] = bilerp(grid_u(v00, v01), grid_u(v00, v11), grid_u(v10, v01), grid_u(v10, v11), fx, fy);
				gu[1] = bilerp(grid_v(v00, v01), grid_v(v00, v11), grid_v(v10, v01), grid_v(v10, v11), fx, fy);

				midx[0] = x_pos + 0.5 * dt * gu[0];
				midx[1] = y_pos + 0.5 * dt * gu[1];

				midx[0] = clamp(midx[0], 0.0f, xmax);
				midx[1] = clamp(midx[1], 0.0f, ymax);

				// second stage of Runge-Kutta 2
				v00 = floor(midx[0]); v10 = ceil(midx[0]);
				v01 = floor(midx[1]); v11 = ceil(midx[1]);

				fx = midx[0] - v00;
				if (fx > v10 - midx[0])
					fx = v10 - midx[0];
				fy = midx[1] - v01;
				if (fy > v11 - midx[1])
					fy = v11 - midx[1];

				gu[0] = bilerp(grid_u(v00, v01), grid_u(v00, v11), grid_u(v10, v01), grid_u(v10, v11), fx, fy);
				gu[1] = bilerp(grid_v(v00, v01), grid_v(v00, v11), grid_v(v10, v01), grid_v(v10, v11), fx, fy);

				x(p_i, p_j) = midx[0] * dx + dt * gu[0];
				y(p_i, p_j) = midx[1] * dx + dt * gu[1];

				x(p_i, p_j) = clamp(x(p_i, p_j), 0.0f, xmax);
				y(p_i, p_j) = clamp(y(p_i, p_j), 0.0f, ymax);

			}
		}
	}
	void update_from_grid(float dx, float ni, float nj, Array2f grid_u, Array2f grid_v, 
						  Array2f grid_du, Array2f grid_dv) {

		for (int p_i = 0; p_i < sqrt(num_particles); ++p_i) {
			for (int p_j = 0; p_j < sqrt(num_particles); ++p_j) {

				float x_pos = x(p_i, p_j) / dx;
				float y_pos = y(p_i, p_j) / dx;
				int v00 = floor(x_pos); int v10 = ceil(x_pos);
				int v01 = floor(y_pos); int v11 = ceil(y_pos);
				float fx = x_pos - v00;
				if (fx > v10 - x_pos)
					fx = v10 - x_pos;
				float fy = y_pos - v01;
				if (fy > v11 - y_pos)
					fy = v11 - y_pos;

				//FLIP
				u(p_i, p_j) += bilerp(grid_du(v00, v01), grid_du(v00, v11), grid_du(v10, v01), grid_du(v10, v11), fx, fy);
				v(p_i, p_j) += bilerp(grid_dv(v00, v01), grid_dv(v00, v11), grid_dv(v10, v01), grid_dv(v10, v11), fx, fy);

				//PIC
				//u(p_i, p_j) = bilerp(grid_u(v00, v01), grid_u(v00, v11), grid_u(v10, v01), grid_u(v10, v11), fx, fy);
				//v(p_i, p_j) = bilerp(grid_v(v00, v01), grid_v(v00, v11), grid_v(v10, v01), grid_v(v10, v11), fx, fy);
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
