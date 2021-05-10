#ifndef PARTICLES_H
#define PARTICLES_H

#include <iostream>
#include "array2.h"
#include "vec.h"
#include "vorticity.h"
#include "util.h"
#include "array2_utils.h"

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
				get_barycentric(y(p_i, p_j) / dx - 0.5f, j, fy, 0, nj);
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
				get_barycentric(y(p_i, p_j) / dx, j, fy, 0, nj);
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
	void move_particles_in_grid(float dt, Array2f grid_u, Array2f grid_v, float xmax, float ymax, float dx, int ni, int nj) {
		Vec2f midx, gu;
		for (int p_i = 0; p_i < sqrt(num_particles); ++p_i) {
			for (int p_j = 0; p_j < sqrt(num_particles); ++p_j) {
				// first stage of Runge-Kutta 2 (do a half Euler step)
				float x_pos = x(p_i, p_j) / dx;
				float y_pos = y(p_i, p_j) / dx - 0.5;
				gu[0] = interpolate_value(Vec2f(x_pos, y_pos), grid_u);// , ni, nj);

				x_pos = x(p_i, p_j) / dx - 0.5;
				y_pos = y(p_i, p_j) / dx;
				gu[1] = interpolate_value(Vec2f(x_pos, y_pos), grid_v);// , ni, nj);
				
				midx[0] = x(p_i, p_j) + 0.5 * dt * gu[0];
				midx[1] = y(p_i, p_j) + 0.5 * dt * gu[1];

				midx[0] = clamp(midx[0], 0.0f, xmax);
				midx[1] = clamp(midx[1], 0.0f, ymax);

				// second stage of Runge-Kutta 2
				x_pos = midx[0] / dx;
				y_pos = midx[1] / dx - 0.5;
				gu[0] = interpolate_value(Vec2f(x_pos, y_pos), grid_u);// , ni, nj);

				x_pos = midx[0] / dx - 0.5;
				y_pos = midx[1] / dx;
				gu[1] = interpolate_value(Vec2f(x_pos, y_pos), grid_v);// , ni, nj);

				x(p_i, p_j) += dt * gu[0];
				y(p_i, p_j) += dt * gu[1];

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
				float y_pos = y(p_i, p_j) / dx - 0.5;
				u(p_i, p_j) += interpolate_value(Vec2f(x_pos, y_pos), grid_du);// , ni, nj);

				x_pos = x(p_i, p_j) / dx - 0.5;
				y_pos = y(p_i, p_j) / dx;
				v(p_i, p_j) += interpolate_value(Vec2f(x_pos, y_pos), grid_dv);// , ni, nj);

				// What happens inside interpolate_value -> get_barycentric & bilerp
				
				/*int v00 = floor(x_pos); int v10 = ceil(x_pos);
				int v01 = floor(y_pos); int v11 = ceil(y_pos);
				float fx = x_pos - v00;
				if (fx > v10 - x_pos)
					fx = v10 - x_pos;
				float fy = y_pos - v01;
				if (fy > v11 - y_pos)
					fy = v11 - y_pos;*/

				//FLIP
				//u(p_i, p_j) += bilerp(grid_du(v00, v01), grid_du(v00, v11), grid_du(v10, v01), grid_du(v10, v11), fx, fy);
				//v(p_i, p_j) += bilerp(grid_dv(v00, v01), grid_dv(v00, v11), grid_dv(v10, v01), grid_dv(v10, v11), fx, fy);
				
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
	/*Vec2f interpolate_value(const Vec2f &point, const Array2f &grid, int ni, int nj) {
		int i, j;
		float fx, fy;

		get_barycentric(point[0], i, fx, 0, ni);
		get_barycentric(point[1], j, fy, 0, nj);

		return bilerp(
			grid(i, j), grid(i + 1, j),
			grid(i, j + 1), grid(i + 1, j + 1),
			fx, fy);
	}*/
};

#endif
