#include "fluidsim.h"

#include "array3_utils.h"
#include "levelset_util.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"

void FluidSim::initialize(int i, int j, int k, float width) {
    _isize = i;
    _jsize = j;
    _ksize = k;
    _dx = width / (float)_isize;

    _u.resize(_isize + 1, _jsize, _ksize); 
    _v.resize(_isize, _jsize + 1, _ksize); 
    _w.resize(_isize, _jsize, _ksize + 1);
    _u.set_zero();
    _v.set_zero();
    _w.set_zero();

    _temp_u.resize(_isize + 1, _jsize, _ksize); 
    _temp_v.resize(_isize, _jsize + 1, _ksize); 
    _temp_w.resize(_isize, _jsize, _ksize + 1); 

    _u_weights = Array3d<float>(_isize + 1, _jsize, _ksize); 
    _v_weights = Array3d<float>(_isize, _jsize + 1, _ksize); 
    _w_weights = Array3d<float>(_isize, _jsize, _ksize + 1); 

    _u_valid = Array3d<bool>(_isize + 1, _jsize, _ksize);
    _v_valid = Array3d<bool>(_isize, _jsize + 1, _ksize); 
    _w_valid = Array3d<bool>(_isize, _jsize, _ksize + 1);

    _particle_radius = (float)(_dx * 1.01*sqrt(3.0)/2.0); 
    //make the particles large enough so they always appear on the grid

    _nodal_solid_phi = Array3d<float>(_isize + 1, _jsize + 1, _ksize + 1);
    _liquid_phi = Array3d<float>(_isize, _jsize, _ksize);
}

//Initialize the grid-based signed distance field that dictates the position of the solid boundary
void FluidSim::set_boundary(float (*phi)(vmath::vec3)) {

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                vmath::vec3 pos(i * _dx, j * _dx, k * _dx);
                _nodal_solid_phi.set(i, j, k, phi(pos));
            }
        }
    }

}

void FluidSim::set_liquid(float (*phi)(vmath::vec3)) {
    //surface.reset_phi(phi, _dx, Vec3f(0.5f*_dx,0.5f*_dx,0.5f*_dx), ni, _jsize, _ksize);
    
    //initialize particles
    int seed = 0;
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) { 
            for(int i = 0; i < _isize; i++) {

                for (int i_dx = 0; i_dx < 8; i_dx++) {
                    vmath::vec3 pos(i*_dx, j*_dx, k*_dx);
                    float a = randhashf(seed++); 
                    float b = randhashf(seed++); 
                    float c = randhashf(seed++);
                    pos += _dx * vmath::vec3(a, b, c);

                    if(phi(pos) <= -_particle_radius) {
                        float solid_phi = interpolate_value(pos/_dx, _nodal_solid_phi);
                        if(solid_phi >= 0) {
                            particles.push_back(pos);
                        }
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
        float substep = _cfl();   
        if(t + substep > dt) {
            substep = dt - t;
        }
        printf("Taking substep of size %f (to %0.3f%% of the frame)\n", substep, 100 * (t+substep)/dt);
        
        printf(" Surface (particle) advection\n");
        _advect_particles(substep);

        printf(" Velocity advection\n");
        //Advance the velocity
        _advect(substep);
        _add_force(substep);

        printf(" _pressure projection\n");
        _project(substep); 
         
        //_pressure projection only produces valid velocities in faces with non-zero associated face area.
        //Because the advection step may interpolate from these invalid faces, 
        //we must extrapolate velocities from the fluid domain into these invalid faces.
        printf(" Extrapolation\n");
        _extrapolate(_u, _u_valid);
        _extrapolate(_v, _v_valid);
        _extrapolate(_w, _w_valid);
     
        //For extrapolated velocities, replace the normal component with
        //that of the object.
        printf(" Constrain boundary velocities\n");
        _constrain_velocity();

        t += substep;
    }
}


float FluidSim::_cfl() {

    float maxvel = 0;
    for(unsigned int i = 0; i < _u.a.size(); i++) {
        maxvel = fmax(maxvel, fabs(_u.a[i]));
    }
    for(unsigned int i = 0; i < _v.a.size(); i++) {
        maxvel = fmax(maxvel, fabs(_v.a[i]));
    }
    for(unsigned int i = 0; i < _w.a.size(); i++) {
        maxvel = fmax(maxvel, fabs(_w.a[i]));
    }
    
    return 5*_dx / maxvel;
}

void FluidSim::add_particle(vmath::vec3 pos) {
    particles.push_back(pos);
}

void FluidSim::_add_force(float dt) {

    //gravity
    for(int k = 0;k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                _v(i, j, k) -= 9.81f * dt;
            }
        }
    }

}



//For extrapolated points, replace the normal component
//of velocity with the object velocity (in this case zero).
void FluidSim::_constrain_velocity() {
    _temp_u = _u;
    _temp_v = _v;
    _temp_w = _w;

    //(At lower grid resolutions, the normal estimate from the signed
    //distance function can be poor, so it doesn't work quite as well.
    //An exact normal would do better if we had it for the geometry.)

    //constrain u
    for(int k = 0; k < _u.nk; k++) {
        for(int j = 0; j < _u.nj; j++) { 
            for(int i = 0; i < _u.ni; i++) {
                if(_u_weights(i, j, k) == 0) {
                    //apply constraint
                    vmath::vec3 pos(i*_dx, (j+0.5f)*_dx, (k+0.5f)*_dx);
                    vmath::vec3 vel = _get_velocity(pos);
                    vmath::vec3 normal(0,0,0);
                    interpolate_gradient(normal, pos/_dx, _nodal_solid_phi); 
                    normal = vmath::normalize(normal);
                    float perp_component = vmath::dot(vel, normal);
                    vel -= perp_component*normal;
                    _temp_u(i, j, k) = vel[0];
                }
            }
        }
    }

    //constrain v
    for(int k = 0; k < _v.nk; k++) {
        for(int j = 0; j < _v.nj; j++) { 
            for(int i = 0; i < _v.ni; i++) {
                if(_v_weights(i, j, k) == 0) {
                    //apply constraint
                    vmath::vec3 pos((i+0.5f)*_dx, j*_dx, (k+0.5f)*_dx);
                    vmath::vec3 vel = _get_velocity(pos);
                    vmath::vec3 normal(0,0,0);
                    interpolate_gradient(normal, pos/_dx, _nodal_solid_phi); 
                    normal = vmath::normalize(normal);
                    float perp_component = vmath::dot(vel, normal);
                    vel -= perp_component*normal;
                    _temp_v(i, j, k) = vel[1];
                }
            }
        }
    }

    //constrain w
    for(int k = 0; k < _w.nk; k++) {
        for(int j = 0; j < _w.nj; j++) {
            for(int i = 0; i < _w.ni; i++) {
                if(_w_weights(i, j, k) == 0) {
                    //apply constraint
                    vmath::vec3 pos((i+0.5f)*_dx, (j+0.5f)*_dx, k*_dx);
                    vmath::vec3 vel = _get_velocity(pos);
                    vmath::vec3 normal(0,0,0);
                    interpolate_gradient(normal, pos/_dx, _nodal_solid_phi); 
                    normal = vmath::normalize(normal);
                    float perp_component = vmath::dot(vel, normal);
                    vel -= perp_component*normal;
                    _temp_w(i, j, k) = vel[2];
                }
            }
        }
    }

    //update
    _u = _temp_u;
    _v = _temp_v;
    _w = _temp_w;

}

void FluidSim::_advect_particles(float dt) { 

    for(unsigned int p = 0; p < particles.size(); p++) {
        particles[p] = _trace_rk2(particles[p], dt);
    
        //check boundaries and project exterior particles back in
        float phi_val = interpolate_value(particles[p]/_dx, _nodal_solid_phi); 
        if(phi_val < 0) {
            vmath::vec3 grad;
            interpolate_gradient(grad, particles[p]/_dx, _nodal_solid_phi);
            if(vmath::lengthsq(grad) > 0) {
                grad = vmath::normalize(grad);
            }
            particles[p] -= phi_val * grad;
        }
    }
    
}

//Basic first order semi-Lagrangian advection of velocities
void FluidSim::_advect(float dt) {

    _temp_u.assign(0);
    _temp_v.assign(0);
    _temp_w.assign(0);

    //semi-Lagrangian advection on u-component of velocity
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                vmath::vec3 pos(i*_dx, (j+0.5f)*_dx, (k+0.5f)*_dx);
                pos = _trace_rk2(pos, -dt);
                _temp_u(i, j, k) = _get_velocity(pos)[0];  
            }
        }
    }

    //semi-Lagrangian advection on v-component of velocity
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                vmath::vec3 pos((i+0.5f)*_dx, j*_dx, (k+0.5f)*_dx);
                pos = _trace_rk2(pos, -dt);
                _temp_v(i, j, k) = _get_velocity(pos)[1];
            }
        }
    }

    //semi-Lagrangian advection on w-component of velocity
    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) { 
            for(int i = 0; i < _isize; i++) {
                vmath::vec3 pos((i+0.5f)*_dx, (j+0.5f)*_dx, k*_dx);
                pos = _trace_rk2(pos, -dt);
                _temp_w(i, j, k) = _get_velocity(pos)[2];
            }
        }
    }

    //move update velocities into u/v vectors
    _u = _temp_u;
    _v = _temp_v;
    _w = _temp_w;
}

void FluidSim::_compute_phi() {
    
    //grab from particles
    _liquid_phi.fill(3*_dx);
    for(unsigned int p = 0; p < particles.size(); ++p) {

        GridIndex cell_ind = Grid3d::positionToGridIndex(particles[p], _dx);
        for(int k = max(0,cell_ind[2] - 1); k <= min(cell_ind[2]+1,_ksize-1); k++) {
            for(int j = max(0,cell_ind[1] - 1); j <= min(cell_ind[1]+1,_jsize-1); j++) {
                for(int i = max(0,cell_ind[0] - 1); i <= min(cell_ind[0]+1,_isize-1); i++) {
                    vmath::vec3 sample_pos((i+0.5f)*_dx, (j+0.5f)*_dx,(k+0.5f)*_dx);
                    float test_val = vmath::length(sample_pos - particles[p]) - _particle_radius;
                    if(test_val < _liquid_phi(i, j, k)) {
                        _liquid_phi.set(i, j, k, test_val);
                    }
                }
            }
        }

    }
    
    //extend phi slightly into solids (this is a simple, naive approach, but works reasonably well)
    Array3d<float> phi_temp = _liquid_phi;
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if(_liquid_phi(i,j,k) < 0.5*_dx) {
                    float solid_phi_val = 0.125f*(_nodal_solid_phi(i, j, k) + 
                                                  _nodal_solid_phi(i + 1, j, k) + 
                                                  _nodal_solid_phi(i, j + 1, k) + 
                                                  _nodal_solid_phi(i + 1, j + 1, k) +
                                                  _nodal_solid_phi(i, j, k + 1) + 
                                                  _nodal_solid_phi(i + 1, j, k + 1) + 
                                                  _nodal_solid_phi(i, j + 1, k + 1) + 
                                                  _nodal_solid_phi(i + 1, j + 1, k + 1));
                    if(solid_phi_val < 0) {
                        phi_temp.set(i, j, k, -0.5f * _dx);
                    }
                }
            }
        }
    }
    _liquid_phi = phi_temp;
    
}



void FluidSim::_project(float dt) {

    //Estimate the liquid signed distance
    _compute_phi();
    
    //Compute finite-volume type face area weight for each velocity sample.
    _compute_weights();

    //Set up and solve the variational _pressure solve.
    _solve_pressure(dt);
    
}


//Apply RK2 to advect a point in the domain.
vmath::vec3 FluidSim::_trace_rk2(vmath::vec3 position, float dt) {
    vmath::vec3 input = position;
    vmath::vec3 velocity = _get_velocity(input);
    velocity = _get_velocity(input + 0.5f * dt * velocity);
    input += dt * velocity;
    return input;
}

//Interpolate velocity from the MAC grid.
vmath::vec3 FluidSim::_get_velocity(vmath::vec3 position) {

    //Interpolate the velocity from the u and v grids
    float u_value = interpolate_value(position / _dx - vmath::vec3(0, 0.5f, 0.5f), _u);
    float v_value = interpolate_value(position / _dx - vmath::vec3(0.5f, 0, 0.5f), _v);
    float w_value = interpolate_value(position / _dx - vmath::vec3(0.5f, 0.5f, 0), _w);

    return vmath::vec3(u_value, v_value, w_value);
}



//Compute finite-volume style face-weights for fluid from nodal signed distances
void FluidSim::_compute_weights() {

    //Compute face area fractions (using marching squares cases).
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                float weight = 1 - fraction_inside(_nodal_solid_phi(i, j, k),
                                                   _nodal_solid_phi(i, j + 1, k),
                                                   _nodal_solid_phi(i, j, k + 1),
                                                   _nodal_solid_phi(i, j + 1, k + 1));
                _u_weights.set(i, j, k, clamp(weight, 0.0f, 1.0f));
            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                float weight = 1 - fraction_inside(_nodal_solid_phi(i, j, k),
                                                   _nodal_solid_phi(i, j, k + 1),
                                                   _nodal_solid_phi(i + 1, j, k),
                                                   _nodal_solid_phi(i + 1, j, k + 1));
                _v_weights.set(i, j, k, clamp(weight, 0.0f, 1.0f));
            }
        }
    }

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                float weight = 1 - fraction_inside(_nodal_solid_phi(i, j, k),
                                                   _nodal_solid_phi(i, j + 1, k),
                                                   _nodal_solid_phi(i + 1, j, k),
                                                   _nodal_solid_phi(i + 1, j + 1, k));
                _w_weights.set(i, j, k, clamp(weight, 0.0f, 1.0f));
            }
        }
    }

}

//An implementation of the variational _pressure projection solve for static geometry
void FluidSim::_solve_pressure(float dt) {

    int system_size = _isize * _jsize * _ksize;
    if((int)_rhs.size() != system_size) {
        _rhs.resize(system_size);
        _pressure.resize(system_size);
        _matrix.resize(system_size);
    }
    
    _matrix.zero();
    _rhs.assign(_rhs.size(), 0);
    _pressure.assign(_pressure.size(), 0);

    //Build the linear system for _pressure
    for(int k = 1; k < _ksize - 1; k++) {
        for(int j = 1; j < _jsize - 1; j++) {
            for(int i = 1; i < _isize - 1; i++) {
                int index = i + _isize*j + _isize*_jsize*k;

                _rhs[index] = 0;
                _pressure[index] = 0;
                float centre_phi = _liquid_phi(i,j,k);
                if(centre_phi >= 0) {
                    continue;
                }

                //right neighbour
                float term = _u_weights(i + 1, j, k) * dt / sqr(_dx);
                float right_phi = _liquid_phi(i + 1, j, k);
                if(right_phi < 0) {
                    _matrix.add_to_element(index, index, term);
                    _matrix.add_to_element(index, index + 1, -term);
                } else {
                    float theta = fraction_inside(centre_phi, right_phi);
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    _matrix.add_to_element(index, index, term / theta);
                }
                _rhs[index] -= _u_weights(i + 1, j, k) * _u(i + 1, j, k) / _dx;

                //left neighbour
                term = _u_weights(i, j, k) * dt / sqr(_dx);
                float left_phi = _liquid_phi(i - 1, j, k);
                if(left_phi < 0) {
                    _matrix.add_to_element(index, index, term);
                    _matrix.add_to_element(index, index - 1, -term);
                } else {
                    float theta = fraction_inside(centre_phi, left_phi);
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    _matrix.add_to_element(index, index, term / theta);
                }
                _rhs[index] += _u_weights(i, j, k) * _u(i, j, k) / _dx;

                //top neighbour
                term = _v_weights(i, j + 1, k) * dt / sqr(_dx);
                float top_phi = _liquid_phi(i, j + 1, k);
                if(top_phi < 0) {
                    _matrix.add_to_element(index, index, term);
                    _matrix.add_to_element(index, index + _isize, -term);
                } else {
                    float theta = fraction_inside(centre_phi, top_phi);
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    _matrix.add_to_element(index, index, term/theta);
                }
                _rhs[index] -= _v_weights(i, j + 1, k) * _v(i, j + 1, k) / _dx;

                //bottom neighbour
                term = _v_weights(i, j, k) * dt / sqr(_dx);
                float bot_phi = _liquid_phi(i, j - 1, k);
                if(bot_phi < 0) {
                    _matrix.add_to_element(index, index, term);
                    _matrix.add_to_element(index, index - _isize, -term);
                } else {
                    float theta = fraction_inside(centre_phi, bot_phi);
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    _matrix.add_to_element(index, index, term / theta);
                }
                _rhs[index] += _v_weights(i, j, k) * _v(i, j, k) / _dx;


                //far neighbour
                term = _w_weights(i, j, k + 1) * dt / sqr(_dx);
                float far_phi = _liquid_phi(i, j, k + 1);
                if(far_phi < 0) {
                    _matrix.add_to_element(index, index, term);
                    _matrix.add_to_element(index, index + _isize*_jsize, -term);
                } else {
                    float theta = fraction_inside(centre_phi, far_phi);
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    _matrix.add_to_element(index, index, term / theta);
                }
                _rhs[index] -= _w_weights(i, j, k + 1) * _w(i, j, k + 1) / _dx;

                //near neighbour
                term = _w_weights(i, j, k) * dt / sqr(_dx);
                float near_phi = _liquid_phi(i, j, k - 1);
                if(near_phi < 0) {
                    _matrix.add_to_element(index, index, term);
                    _matrix.add_to_element(index, index - _isize*_jsize, -term);
                } else {
                    float theta = fraction_inside(centre_phi, near_phi);
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    _matrix.add_to_element(index, index, term / theta);
                }
                _rhs[index] += _w_weights(i, j, k) * _w(i, j, k) / _dx;
            }
        }
    }

    //Solve the system using Robert Bridson's incomplete Cholesky PCG solver

    double tolerance;
    int iterations;
    _solver.set_solver_parameters(1e-18, 1000);
    bool success = _solver.solve(_matrix, _rhs, _pressure, tolerance, iterations);
    printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
    if(!success) {
        printf("WARNING: _pressure solve failed!************************************************\n");
    }

    //Apply the velocity update
    _u_valid.fill(false);
    for(int k = 0; k < _u.nk; k++) {
        for(int j = 0; j < _u.nj; j++) {
            for(int i = 1; i < _u.ni - 1; i++) {

                int index = i + j*_isize + k*_isize*_jsize;
                if(_u_weights(i, j, k) > 0 && (_liquid_phi(i, j, k) < 0 || _liquid_phi(i - 1, j, k) < 0)) {
                    float theta = 1;
                    if(_liquid_phi(i, j, k) >= 0 || _liquid_phi(i - 1, j, k) >= 0) {
                        theta = fraction_inside(_liquid_phi(i-1,j,k), _liquid_phi(i,j,k));
                    }
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    _u(i, j, k) -= dt  * (float)(_pressure[index] - _pressure[index-1]) / _dx / theta; 
                    _u_valid.set(i, j, k, true);
                }

            }
        }
    }
    
    _v_valid.fill(false);
    for(int k = 0; k < _v.nk; k++) {
        for(int j = 1; j < _v.nj - 1; j++) {
            for(int i = 0; i < _v.ni; i++) {

                int index = i + j*_isize + k*_isize*_jsize;
                if(_v_weights(i, j, k) > 0 && (_liquid_phi(i, j, k) < 0 || _liquid_phi(i, j - 1, k) < 0)) {
                    float theta = 1;
                    if(_liquid_phi(i, j, k) >= 0 || _liquid_phi(i, j - 1, k) >= 0) {
                        theta = fraction_inside(_liquid_phi(i, j - 1, k), _liquid_phi(i, j, k));
                    }
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    _v(i, j, k) -= dt  * (float)(_pressure[index] - _pressure[index-_isize]) / _dx / theta; 
                    _v_valid.set(i, j, k, true);
                }

            }
        }
    }

    _w_valid.fill(false);
    for(int k = 0; k < _w.nk; ++k) {
        for(int j = 0; j < _w.nj; ++j) {
            for(int i = 1; i < _w.ni-1; ++i) {

                int index = i + j*_isize + k*_isize*_jsize;
                if(_w_weights(i, j, k) > 0 && (_liquid_phi(i, j, k) < 0 || _liquid_phi(i, j, k - 1) < 0)) {
                    float theta = 1;
                    if(_liquid_phi(i, j, k) >= 0 || _liquid_phi(i, j, k - 1) >= 0) {
                        theta = fraction_inside(_liquid_phi(i, j, k - 1), _liquid_phi(i, j, k));
                    }
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    _w(i, j, k) -= dt  * (float)(_pressure[index] - _pressure[index-_isize*_jsize]) / _dx / theta; 
                    _w_valid.set(i, j, k, true);
                }

            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                if(!_u_valid(i, j, k)) {
                    _u(i, j, k) = 0;
                }
            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                if(!_v_valid(i, j, k)) {
                    _v(i, j, k) = 0;
                }
            }
        }
    }

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if(!_w_valid(i, j, k)) {
                    _w(i, j, k) = 0;
                }
            }
        }
    }
}


//Apply several iterations of a very simple propagation of valid velocity data in all directions
void FluidSim::_extrapolate(Array3f& grid, Array3d<bool> &valid) {

    Array3f temp_grid = grid;
    Array3d<bool> old_valid;
    for(int layers = 0; layers < 10; layers++) {

        old_valid = valid;
        for(int k = 1; k < grid.nk - 1; k++) {
            for(int j = 1; j < grid.nj - 1; j++) {
                for(int i = 1; i < grid.ni - 1; i++) {

                    if(old_valid(i,j,k)) {
                        continue;
                    }

                    float sum = 0;
                    int count = 0;
                    if(old_valid(i + 1, j, k)) {
                        sum += grid(i + 1, j, k);
                        count++;
                    }
                    if(old_valid(i - 1, j, k)) {
                        sum += grid(i - 1, j, k);
                        count++;
                    }
                    if(old_valid(i, j + 1, k)) {
                        sum += grid(i, j + 1, k);
                        count++;
                    }
                    if(old_valid(i, j - 1, k)) {
                        sum += grid(i, j - 1, k);
                        count++;
                    }
                    if(old_valid(i, j, k + 1)) {
                        sum += grid(i, j, k + 1);
                        count++;
                    }
                    if(old_valid(i, j, k - 1)) {
                        sum += grid(i, j, k - 1);
                        count++;
                    }

                    //If any of neighbour cells were valid, 
                    //assign the cell their average value and tag it as valid
                    if(count > 0) {
                        temp_grid(i, j, k) = sum /(float)count;
                        valid.set(i, j, k, true);
                    }

                }
            }
        }
        grid = temp_grid;

    }

}
