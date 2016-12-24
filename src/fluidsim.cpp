#include "fluidsim.h"

#include "array3_utils.h"
#include "levelset_util.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"

void extrapolate(Array3f& grid, Array3c& valid);

void FluidSim::initialize(float width, int ni_, int nj_, int nk_) {
   ni = ni_;
   nj = nj_;
   nk = nk_;
   dx = width / (float)ni;
   u.resize(ni+1,nj,nk); temp_u.resize(ni+1,nj,nk); u_weights.resize(ni+1,nj,nk); u_valid.resize(ni+1,nj,nk);
   v.resize(ni,nj+1,nk); temp_v.resize(ni,nj+1,nk); v_weights.resize(ni,nj+1,nk); v_valid.resize(ni,nj+1,nk);
   w.resize(ni,nj,nk+1); temp_w.resize(ni,nj,nk+1); w_weights.resize(ni,nj,nk+1); w_valid.resize(ni,nj,nk+1);

   particle_radius = (float)(dx * 1.01*sqrt(3.0)/2.0); 
   //make the particles large enough so they always appear on the grid

   u.set_zero();
   v.set_zero();
   w.set_zero();
   nodal_solid_phi.resize(ni+1,nj+1,nk+1);
   valid.resize(ni+1, nj+1, nk+1);
   old_valid.resize(ni+1, nj+1, nk+1);
   liquid_phi.resize(ni,nj,nk);

}

//Initialize the grid-based signed distance field that dictates the position of the solid boundary
void FluidSim::set_boundary(float (*phi)(const Vec3f&)) {

   for(int k = 0; k < nk+1; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni+1; ++i) {
      Vec3f pos(i*dx,j*dx,k*dx);
      nodal_solid_phi(i,j,k) = phi(pos);
   }

}

void FluidSim::set_liquid(float (*phi)(const Vec3f&)) {
   //surface.reset_phi(phi, dx, Vec3f(0.5f*dx,0.5f*dx,0.5f*dx), ni, nj, nk);
   
   //initialize particles
   int seed = 0;
   for(int k = 0; k < nk; ++k) {
      for(int j = 0; j < nj; ++j) { 
         for(int i = 0; i < ni; ++i) {

            for (int idx = 0; idx < 8; idx++) {
               Vec3f pos(i*dx,j*dx,k*dx);
               float a = randhashf(seed++); 
               float b = randhashf(seed++); 
               float c = randhashf(seed++);
               pos += dx * Vec3f(a,b,c);

               if(phi(pos) <= -particle_radius) {
                  float solid_phi = interpolate_value(pos/dx, nodal_solid_phi);
                  if(solid_phi >= 0)
                     particles.push_back(pos);
               }
            }

         }
      }
   }
}

//The main fluid simulation step
void FluidSim::advance(float dt) {
   float t = 0;

   while(t < dt) {
      float substep = cfl();   
      if(t + substep > dt)
         substep = dt - t;
      printf("Taking substep of size %f (to %0.3f%% of the frame)\n", substep, 100 * (t+substep)/dt);
      
      printf(" Surface (particle) advection\n");
      advect_particles(substep);

      printf(" Velocity advection\n");
      //Advance the velocity
      advect(substep);
      add_force(substep);

      printf(" Pressure projection\n");
      project(substep); 
       
      //Pressure projection only produces valid velocities in faces with non-zero associated face area.
      //Because the advection step may interpolate from these invalid faces, 
      //we must extrapolate velocities from the fluid domain into these invalid faces.
      printf(" Extrapolation\n");
      extrapolate(u, u_valid);
      extrapolate(v, v_valid);
      extrapolate(w, w_valid);
    
      //For extrapolated velocities, replace the normal component with
      //that of the object.
      printf(" Constrain boundary velocities\n");
      constrain_velocity();

      t+=substep;
   }
}


float FluidSim::cfl() {

   float maxvel = 0;
   for(unsigned int i = 0; i < u.a.size(); ++i)
      maxvel = fmax(maxvel, fabs(u.a[i]));
   for(unsigned int i = 0; i < v.a.size(); ++i)
      maxvel = fmax(maxvel, fabs(v.a[i]));
   for(unsigned int i = 0; i < w.a.size(); ++i)
      maxvel = fmax(maxvel, fabs(w.a[i]));
   
   return 5*dx / maxvel;
}

void FluidSim::add_particle(const Vec3f& pos) {
   particles.push_back(pos);
}

void FluidSim::add_force(float dt) {

   //gravity
   for(int k = 0;k < nk; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
      v(i,j,k) -= 9.81f * dt;
   }

}



//For extrapolated points, replace the normal component
//of velocity with the object velocity (in this case zero).
void FluidSim::constrain_velocity() {
   temp_u = u;
   temp_v = v;
   temp_w = w;

   //(At lower grid resolutions, the normal estimate from the signed
   //distance function can be poor, so it doesn't work quite as well.
   //An exact normal would do better if we had it for the geometry.)

   //constrain u
   for(int k = 0; k < u.nk;++k) for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i) {
      if(u_weights(i,j,k) == 0) {
         //apply constraint
         Vec3f pos(i*dx, (j+0.5f)*dx, (k+0.5f)*dx);
         Vec3f vel = get_velocity(pos);
         Vec3f normal(0,0,0);
         interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
         normalize(normal);
         float perp_component = dot(vel, normal);
         vel -= perp_component*normal;
         temp_u(i,j,k) = vel[0];
      }
   }

   //constrain v
   for(int k = 0; k < v.nk;++k) for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i) {
      if(v_weights(i,j,k) == 0) {
         //apply constraint
         Vec3f pos((i+0.5f)*dx, j*dx, (k+0.5f)*dx);
         Vec3f vel = get_velocity(pos);
         Vec3f normal(0,0,0);
         interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
         normalize(normal);
         float perp_component = dot(vel, normal);
         vel -= perp_component*normal;
         temp_v(i,j,k) = vel[1];
      }
   }

   //constrain w
   for(int k = 0; k < w.nk;++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i) {
      if(w_weights(i,j,k) == 0) {
         //apply constraint
         Vec3f pos((i+0.5f)*dx, (j+0.5f)*dx, k*dx);
         Vec3f vel = get_velocity(pos);
         Vec3f normal(0,0,0);
         interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
         normalize(normal);
         float perp_component = dot(vel, normal);
         vel -= perp_component*normal;
         temp_w(i,j,k) = vel[2];
      }
   }

   //update
   u = temp_u;
   v = temp_v;
   w = temp_w;

}

void FluidSim::advect_particles(float dt) { 
   for(unsigned int p = 0; p < particles.size(); ++p) {
      particles[p] = trace_rk2(particles[p], dt);
   
      //check boundaries and project exterior particles back in
      float phi_val = interpolate_value(particles[p]/dx, nodal_solid_phi); 
      if(phi_val < 0) {
         Vec3f grad;
         interpolate_gradient(grad, particles[p]/dx, nodal_solid_phi);
         if(mag(grad) > 0)
            normalize(grad);
         particles[p] -= phi_val * grad;
      }
   }
   

}

//Basic first order semi-Lagrangian advection of velocities
void FluidSim::advect(float dt) {

   temp_u.assign(0);
   temp_v.assign(0);
   temp_w.assign(0);

   //semi-Lagrangian advection on u-component of velocity
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni+1; ++i) {
      Vec3f pos(i*dx, (j+0.5f)*dx, (k+0.5f)*dx);
      pos = trace_rk2(pos, -dt);
      temp_u(i,j,k) = get_velocity(pos)[0];  
   }

   //semi-Lagrangian advection on v-component of velocity
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
      Vec3f pos((i+0.5f)*dx, j*dx, (k+0.5f)*dx);
      pos = trace_rk2(pos, -dt);
      temp_v(i,j,k) = get_velocity(pos)[1];
   }

   //semi-Lagrangian advection on w-component of velocity
   for(int k = 0; k < nk+1; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      Vec3f pos((i+0.5f)*dx, (j+0.5f)*dx, k*dx);
      pos = trace_rk2(pos, -dt);
      temp_w(i,j,k) = get_velocity(pos)[2];
   }

   //move update velocities into u/v vectors
   u = temp_u;
   v = temp_v;
   w = temp_w;
}

void FluidSim::compute_phi() {
   
   //grab from particles
   liquid_phi.assign(3*dx);
   for(unsigned int p = 0; p < particles.size(); ++p) {
      Vec3i cell_ind(particles[p] / dx);
      for(int k = max(0,cell_ind[2] - 1); k <= min(cell_ind[2]+1,nk-1); ++k) {
         for(int j = max(0,cell_ind[1] - 1); j <= min(cell_ind[1]+1,nj-1); ++j) {
            for(int i = max(0,cell_ind[0] - 1); i <= min(cell_ind[0]+1,ni-1); ++i) {
               Vec3f sample_pos((i+0.5f)*dx, (j+0.5f)*dx,(k+0.5f)*dx);
               float test_val = dist(sample_pos, particles[p]) - particle_radius;
               if(test_val < liquid_phi(i,j,k))
                  liquid_phi(i,j,k) = test_val;
            }
         }
      }
   }
   
   //extend phi slightly into solids (this is a simple, naive approach, but works reasonably well)
   Array3f phi_temp = liquid_phi;
   for(int k = 0; k < nk; ++k) {
      for(int j = 0; j < nj; ++j) {
         for(int i = 0; i < ni; ++i) {
            if(liquid_phi(i,j,k) < 0.5*dx) {
               float solid_phi_val = 0.125f*(nodal_solid_phi(i,j,k) + nodal_solid_phi(i+1,j,k) + nodal_solid_phi(i,j+1,k) + nodal_solid_phi(i+1,j+1,k)
                  + nodal_solid_phi(i,j,k+1) + nodal_solid_phi(i+1,j,k+1) + nodal_solid_phi(i,j+1,k+1) + nodal_solid_phi(i+1,j+1,k+1));
               if(solid_phi_val < 0)
                  phi_temp(i,j,k) = -0.5f*dx;
            }
         }
      }
   }
   liquid_phi = phi_temp;
   

}



void FluidSim::project(float dt) {

   //Estimate the liquid signed distance
   compute_phi();
   
   //Compute finite-volume type face area weight for each velocity sample.
   compute_weights();

   //Set up and solve the variational pressure solve.
   solve_pressure(dt);
   
}


//Apply RK2 to advect a point in the domain.
Vec3f FluidSim::trace_rk2(const Vec3f& position, float dt) {
   Vec3f input = position;
   Vec3f velocity = get_velocity(input);
   velocity = get_velocity(input + 0.5f*dt*velocity);
   input += dt*velocity;
   return input;
}

//Interpolate velocity from the MAC grid.
Vec3f FluidSim::get_velocity(const Vec3f& position) {

   //Interpolate the velocity from the u and v grids
   float u_value = interpolate_value(position / dx - Vec3f(0, 0.5f, 0.5f), u);
   float v_value = interpolate_value(position / dx - Vec3f(0.5f, 0, 0.5f), v);
   float w_value = interpolate_value(position / dx - Vec3f(0.5f, 0.5f, 0), w);

   return Vec3f(u_value, v_value, w_value);
}



//Compute finite-volume style face-weights for fluid from nodal signed distances
void FluidSim::compute_weights() {

   //Compute face area fractions (using marching squares cases).
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni+1; ++i) {
      u_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,j,  k),
                                             nodal_solid_phi(i,j+1,k),
                                             nodal_solid_phi(i,j,  k+1),
                                             nodal_solid_phi(i,j+1,k+1));
      u_weights(i,j,k) = clamp(u_weights(i,j,k),0.0f,1.0f);
   }
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
      v_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,  j,k),
                                             nodal_solid_phi(i,  j,k+1),
                                             nodal_solid_phi(i+1,j,k),
                                             nodal_solid_phi(i+1,j,k+1));
      v_weights(i,j,k) = clamp(v_weights(i,j,k),0.0f,1.0f);
   }
   for(int k = 0; k < nk+1; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      w_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,  j,  k),
                                             nodal_solid_phi(i,  j+1,k),
                                             nodal_solid_phi(i+1,j,  k),
                                             nodal_solid_phi(i+1,j+1,k));
      w_weights(i,j,k) = clamp(w_weights(i,j,k),0.0f,1.0f);
   }


}

//An implementation of the variational pressure projection solve for static geometry
void FluidSim::solve_pressure(float dt) {


   int ni = v.ni;
   int nj = u.nj;
   int nk = u.nk;

   int system_size = ni*nj*nk;
   if((int)rhs.size() != system_size) {
      rhs.resize(system_size);
      pressure.resize(system_size);
      matrix.resize(system_size);
   }
   
   matrix.zero();
   rhs.assign(rhs.size(), 0);
   pressure.assign(pressure.size(), 0);

   //Build the linear system for pressure
   for(int k = 1; k < nk-1; ++k) {
      for(int j = 1; j < nj-1; ++j) {
         for(int i = 1; i < ni-1; ++i) {
            int index = i + ni*j + ni*nj*k;

            rhs[index] = 0;
            pressure[index] = 0;
            float centre_phi = liquid_phi(i,j,k);
            if(centre_phi < 0) {

               //right neighbour
               float term = u_weights(i+1,j,k) * dt / sqr(dx);
               float right_phi = liquid_phi(i+1,j,k);
               if(right_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index + 1, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, right_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] -= u_weights(i+1,j,k)*u(i+1,j,k) / dx;

               //left neighbour
               term = u_weights(i,j,k) * dt / sqr(dx);
               float left_phi = liquid_phi(i-1,j,k);
               if(left_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index - 1, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, left_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] += u_weights(i,j,k)*u(i,j,k) / dx;

               //top neighbour
               term = v_weights(i,j+1,k) * dt / sqr(dx);
               float top_phi = liquid_phi(i,j+1,k);
               if(top_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index + ni, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, top_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] -= v_weights(i,j+1,k)*v(i,j+1,k) / dx;

               //bottom neighbour
               term = v_weights(i,j,k) * dt / sqr(dx);
               float bot_phi = liquid_phi(i,j-1,k);
               if(bot_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index - ni, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, bot_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] += v_weights(i,j,k)*v(i,j,k) / dx;


               //far neighbour
               term = w_weights(i,j,k+1) * dt / sqr(dx);
               float far_phi = liquid_phi(i,j,k+1);
               if(far_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index + ni*nj, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, far_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] -= w_weights(i,j,k+1)*w(i,j,k+1) / dx;

               //near neighbour
               term = w_weights(i,j,k) * dt / sqr(dx);
               float near_phi = liquid_phi(i,j,k-1);
               if(near_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index - ni*nj, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, near_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] += w_weights(i,j,k)*w(i,j,k) / dx;

               /*
               //far neighbour
               term = w_weights(i,j,k+1) * dt / sqr(dx);
               float far_phi = liquid_phi(i,j,k+1);
               if(far_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index + ni*nj, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, far_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] -= w_weights(i,j,k+1)*w(i,j,k+1) / dx;

               //near neighbour
               term = w_weights(i,j,k) * dt / sqr(dx);
               float near_phi = liquid_phi(i,j,k-1);
               if(near_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index - ni*nj, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, near_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] += w_weights(i,j,k)*w(i,j,k) / dx;   
               */

            }
         }
      }
   }

   //Solve the system using Robert Bridson's incomplete Cholesky PCG solver

   double tolerance;
   int iterations;
   solver.set_solver_parameters(1e-18, 1000);
   bool success = solver.solve(matrix, rhs, pressure, tolerance, iterations);
   //printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
   if(!success) {
      printf("WARNING: Pressure solve failed!************************************************\n");
   }

   //Apply the velocity update
   u_valid.assign(0);
   for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 1; i < u.ni-1; ++i) {
      int index = i + j*ni + k*ni*nj;
      if(u_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i-1,j,k) < 0)) {
         float theta = 1;
         if(liquid_phi(i,j,k) >= 0 || liquid_phi(i-1,j,k) >= 0)
            theta = fraction_inside(liquid_phi(i-1,j,k), liquid_phi(i,j,k));
         if(theta < 0.01f) theta = 0.01f;
         u(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-1]) / dx / theta; 
         u_valid(i,j,k) = 1;
      }
   }
   
   v_valid.assign(0);
   for(int k = 0; k < v.nk; ++k) for(int j = 1; j < v.nj-1; ++j) for(int i = 0; i < v.ni; ++i) {
      int index = i + j*ni + k*ni*nj;
      if(v_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i,j-1,k) < 0)) {
         float theta = 1;
         if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j-1,k) >= 0)
            theta = fraction_inside(liquid_phi(i,j-1,k), liquid_phi(i,j,k));
         if(theta < 0.01f) theta = 0.01f;
         v(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-ni]) / dx / theta; 
         v_valid(i,j,k) = 1;
      }
   }

   w_valid.assign(0);
   for(int k = 0; k < w.nk; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 1; i < w.ni-1; ++i) {
      int index = i + j*ni + k*ni*nj;
      if(w_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i,j,k-1) < 0)) {
         float theta = 1;
         if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j,k-1) >= 0)
            theta = fraction_inside(liquid_phi(i,j,k-1), liquid_phi(i,j,k));
         if(theta < 0.01f) theta = 0.01f;
         w(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-ni*nj]) / dx / theta; 
         w_valid(i,j,k) = 1;
      }
   }
 
   for(unsigned int i = 0; i < u_valid.a.size(); ++i)
      if(u_valid.a[i] == 0)
         u.a[i] = 0;
   for(unsigned int i = 0; i < v_valid.a.size(); ++i)
      if(v_valid.a[i] == 0)
         v.a[i] = 0;
   for(unsigned int i = 0; i < w_valid.a.size(); ++i)
      if(w_valid.a[i] == 0)
         w.a[i] = 0;
}


//Apply several iterations of a very simple propagation of valid velocity data in all directions
void extrapolate(Array3f& grid, Array3c& valid) {

   Array3f temp_grid = grid;
   Array3c old_valid(valid.ni,valid.nj,valid.nk);
   for(int layers = 0; layers < 10; ++layers) {
      old_valid = valid;
      for(int k = 1; k < grid.nk-1; ++k) for(int j = 1; j < grid.nj-1; ++j) for(int i = 1; i < grid.ni-1; ++i) {
         float sum = 0;
         int count = 0;

         if(!old_valid(i,j,k)) {

            if(old_valid(i+1,j,k)) {
               sum += grid(i+1,j,k);
               ++count;
            }
            if(old_valid(i-1,j,k)) {
               sum += grid(i-1,j,k);
               ++count;
            }
            if(old_valid(i,j+1,k)) {
               sum += grid(i,j+1,k);
               ++count;
            }
            if(old_valid(i,j-1,k)) {
               sum += grid(i,j-1,k);
               ++count;
            }
            if(old_valid(i,j,k+1)) {
               sum += grid(i,j,k+1);
               ++count;
            }
            if(old_valid(i,j,k-1)) {
               sum += grid(i,j,k-1);
               ++count;
            }

            //If any of neighbour cells were valid, 
            //assign the cell their average value and tag it as valid
            if(count > 0) {
               temp_grid(i,j,k) = sum /(float)count;
               valid(i,j,k) = 1;
            }

         }
      }
      grid = temp_grid;

   }

}
