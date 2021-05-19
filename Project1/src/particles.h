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

	std::vector<Vec2f> particles_pos;
	std::vector<Vec2f> particles_vel;

	void initialize(float width, int ni, int nj) {
		num_particles = ni * nj; // each grid has one particle
		float dx = width / ni; // dx = width / nj
		gamma = gamma_init(dx, ni, nj);

		int rando = 21753;
		int part_per_cell = 6;
		for (int i = 0; i < ni; ++i) {
			for (int j = 0; j < nj; ++j) {
				for (int p = 0; p < part_per_cell; ++p) {
					float r1 = randhashf(++rando)-0.5;
					float r2 = randhashf(++rando)-0.5;
					Vec2f pos = Vec2f((i + 0.5f) * dx, (j + 0.5f) * dx) + dx*Vec2f(r1, r2);
					particles_pos.push_back(pos);
					particles_vel.push_back(biot_savart(dx, ni, nj, pos, gamma));
				}
			}
		}
		
	}


	void transfer_to_grid(float width, int ni, int nj, Array2f& grid_u, Array2f& grid_v) {
		int i, j;
		float fx, fy;
		float dx = width / (float)ni;

		sum.resize(ni+1, nj);
		sum.set_zero();

		grid_u.set_zero();
		for (int p = 0; p < particles_pos.size(); ++p) {
			get_barycentric(particles_pos[p][0] / dx, i, fx, 0, ni+1);
			get_barycentric(particles_pos[p][1] / dx - 0.5f, j, fy, 0, nj);
			accumulate(grid_u, particles_vel[p][0], i, j, fx, fy);
		}

		for (j = 0; j < nj; ++j) {
			for (i = 0; i < ni+1; ++i) {
				if (sum(i, j) != 0)
					grid_u(i, j) /= sum(i, j);
				else
					grid_u(i, j) = 0;
			}
		}

		sum.resize(ni, nj+1);
		sum.set_zero();
		grid_v.set_zero();
		for (int p = 0; p < particles_pos.size(); ++p) {
			get_barycentric(particles_pos[p][0] / dx - 0.5f, i, fx, 0, ni);
			get_barycentric(particles_pos[p][1] / dx, j, fy, 0, nj+1);
			accumulate(grid_v, particles_vel[p][1], i, j, fx, fy);
			
		}
		for (j = 0; j < nj+1; ++j) {
			for (i = 0; i < ni; ++i) {
				if (sum(i, j) != 0)
					grid_v(i, j) /= sum(i, j);
				else
					grid_v(i, j) = 0;
			}
		}
	}


	void move_particles_in_grid(float dt, Array2f grid_u, Array2f grid_v, float xmax, float ymax, float dx, int ni, int nj) {
		Vec2f midx, gu;
		for (int p = 0; p < particles_pos.size(); ++p) {
			
			// first stage of Runge-Kutta 2 (do a half Euler step)
			float x_pos = particles_pos[p][0] / dx;
			float y_pos = particles_pos[p][1] / dx - 0.5;
			gu[0] = interpolate_value(Vec2f(x_pos, y_pos), grid_u);

			x_pos = particles_pos[p][0] / dx - 0.5;
			y_pos = particles_pos[p][1] / dx;
			gu[1] = interpolate_value(Vec2f(x_pos, y_pos), grid_v);

			midx[0] = particles_pos[p][0] + 0.5 * dt * gu[0];
			midx[1] = particles_pos[p][1] + 0.5 * dt * gu[1];

			midx[0] = clamp(midx[0], 0.0f, xmax);
			midx[1] = clamp(midx[1], 0.0f, ymax);

			// second stage of Runge-Kutta 2
			x_pos = midx[0] / dx;
			y_pos = midx[1] / dx - 0.5;
			gu[0] = interpolate_value(Vec2f(x_pos, y_pos), grid_u);

			x_pos = midx[0] / dx - 0.5;
			y_pos = midx[1] / dx;
			gu[1] = interpolate_value(Vec2f(x_pos, y_pos), grid_v);

			particles_pos[p][0] += dt * gu[0];
			particles_pos[p][1] += dt * gu[1];

			particles_pos[p][0] = clamp(particles_pos[p][0], 0.0f, xmax);
			particles_pos[p][1] = clamp(particles_pos[p][1], 0.0f, ymax);

		}
	}
	

	void update_from_grid(float dx, float ni, float nj, const Array2f& grid_u, const Array2f& grid_v,
		const Array2f& grid_du, const Array2f& grid_dv) {

		for (int p = 0; p < particles_pos.size(); ++p) {
			
				float x_pos = particles_pos[p][0] / dx;
				float y_pos = particles_pos[p][1] / dx - 0.5;
				particles_vel[p][0] += interpolate_value(Vec2f(x_pos, y_pos), grid_du); // FLIP

				particles_vel[p][0] = 0.99 * particles_vel[p][0] + 0.01* interpolate_value(Vec2f(x_pos, y_pos), grid_u); //PIC/FLIP mix
				//particles_vel[p][0] = interpolate_value(Vec2f(x_pos, y_pos), grid_u); //PIC
		
				x_pos = particles_pos[p][0] / dx - 0.5;
				y_pos = particles_pos[p][1] / dx;
				particles_vel[p][1] += interpolate_value(Vec2f(x_pos, y_pos), grid_dv);

				particles_vel[p][1] = 0.99 * particles_vel[p][1] + 0.01 * interpolate_value(Vec2f(x_pos, y_pos), grid_v); //PIC/FLIP mix
				//particles_vel[p][1] = interpolate_value(Vec2f(x_pos, y_pos), grid_v); //PIC
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
