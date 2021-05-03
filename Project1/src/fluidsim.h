#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "array2.h"
#include "vec.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"

#include "particles.h"

#include <vector>

class FluidSim {

public:
   void initialize(float width, int ni_, int nj_);
   void set_boundary(float (*phi)(const Vec2f&));
   void advance(float dt);

   //Grid dimensions
   int ni,nj;
   float dx;
   
   //Fluid velocity
   Array2f u, v;
   Array2f temp_u, temp_v;

   //Fluid vorticity
   Array2f gamma;
   Array2f vor;
   
   //Tracer particles
   std::vector<Vec2f> particles;

   //FLIP particles
   Particles flip_particles;
   
   //Static geometry representation
   Array2f nodal_solid_phi;
   Array2f u_weights, v_weights;

   //Data arrays for extrapolation
   Array2c valid, old_valid;

   //Solver data
   PCGSolver<double> solver;
   SparseMatrixd matrix;
   std::vector<double> rhs;
   std::vector<double> pressure;

   Vec2f get_velocity(const Vec2f& position);
   void add_particle(const Vec2f& position);

private:

   Vec2f trace_rk2(const Vec2f& position, float dt);

   //tracer particle operations
   
   void advect_particles(float dt);

   //fluid velocity operations
   void advect(float dt);
   void add_force(float dt);

   void project(float dt);
   void compute_weights();
   void solve_pressure(float dt);
   
   void constrain_velocity();

};

#endif
