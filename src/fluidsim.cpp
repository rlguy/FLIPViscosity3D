#include "fluidsim.h"

#include "levelset_util.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"

void FluidSim::initialize(int i, int j, int k, float width) {
    _isize = i;
    _jsize = j;
    _ksize = k;
    _dx = width / (float)_isize;

    _MACVelocity = MACVelocityField(_isize, _jsize, _ksize, _dx);
    _tempMACVelocity = MACVelocityField(_isize, _jsize, _ksize, _dx);

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
    //initialize particles
    int seed = 0;
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) { 
            for(int i = 0; i < _isize; i++) {
                vmath::vec3 gpos = Grid3d::GridIndexToPosition(i, j, k, _dx);

                for (int i_dx = 0; i_dx < 8; i_dx++) {
                    float a = randhashf(seed++); 
                    float b = randhashf(seed++); 
                    float c = randhashf(seed++);
                    vmath::vec3 jitter = _dx * vmath::vec3(a, b, c);
                    vmath::vec3 pos = gpos + jitter;

                    if(phi(pos) <= -_particle_radius) {
                        float solid_phi = Interpolation::trilinearInterpolate(pos, _dx, _nodal_solid_phi);
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
        Array3d<float> *ugrid = _MACVelocity.getArray3dU();
        Array3d<float> *vgrid = _MACVelocity.getArray3dV();
        Array3d<float> *wgrid = _MACVelocity.getArray3dW();
        _extrapolate(ugrid, _u_valid);
        _extrapolate(vgrid, _v_valid);
        _extrapolate(wgrid, _w_valid);
     
        //For extrapolated velocities, replace the normal component with
        //that of the object.
        printf(" Constrain boundary velocities\n");
        _constrain_velocity();

        t += substep;
    }
}


float FluidSim::_cfl() {

    float maxvel = 0;
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                maxvel = fmax(maxvel, fabs(_MACVelocity.U(i, j, k)));
            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                maxvel = fmax(maxvel, fabs(_MACVelocity.V(i, j, k)));
            }
        }
    }

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) { 
            for(int i = 0; i < _isize; i++) {
                maxvel = fmax(maxvel, fabs(_MACVelocity.W(i, j, k)));
            }
        }
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
                double v = _MACVelocity.V(i, j, k);
                _MACVelocity.setV(i, j, k, v - 9.81f * dt);
            }
        }
    }

}



//For extrapolated points, replace the normal component
//of velocity with the object velocity (in this case zero).
void FluidSim::_constrain_velocity() {
    _tempMACVelocity.set(_MACVelocity);

    //(At lower grid resolutions, the normal estimate from the signed
    //distance function can be poor, so it doesn't work quite as well.
    //An exact normal would do better if we had it for the geometry.)

    //constrain u
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                if(_u_weights(i, j, k) == 0) {
                    //apply constraint
                    vmath::vec3 pos = Grid3d::FaceIndexToPositionU(i, j, k, _dx);
                    vmath::vec3 vel = _get_velocity(pos);
                    vmath::vec3 normal(0,0,0);
                    Interpolation::trilinearInterpolateGradient(pos, _dx, _nodal_solid_phi, &normal);
                    normal = vmath::normalize(normal);
                    float perp_component = vmath::dot(vel, normal);
                    vel -= perp_component*normal;
                    _tempMACVelocity.setU(i, j, k, vel.x);
                }
            }
        }
    }

    //constrain v
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                if(_v_weights(i, j, k) == 0) {
                    //apply constraint
                    vmath::vec3 pos = Grid3d::FaceIndexToPositionV(i, j, k, _dx);
                    vmath::vec3 vel = _get_velocity(pos);
                    vmath::vec3 normal(0,0,0);
                    Interpolation::trilinearInterpolateGradient(pos, _dx, _nodal_solid_phi, &normal);
                    normal = vmath::normalize(normal);
                    float perp_component = vmath::dot(vel, normal);
                    vel -= perp_component*normal;
                    _tempMACVelocity.setV(i, j, k, vel.y);
                }
            }
        }
    }

    //constrain w
    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) { 
            for(int i = 0; i < _isize; i++) {
                if(_w_weights(i, j, k) == 0) {
                    //apply constraint
                    vmath::vec3 pos = Grid3d::FaceIndexToPositionW(i, j, k, _dx);
                    vmath::vec3 vel = _get_velocity(pos);
                    vmath::vec3 normal(0,0,0);
                    Interpolation::trilinearInterpolateGradient(pos, _dx, _nodal_solid_phi, &normal);
                    normal = vmath::normalize(normal);
                    float perp_component = vmath::dot(vel, normal);
                    vel -= perp_component*normal;
                    _tempMACVelocity.setW(i, j, k, vel.z);
                }
            }
        }
    }

    _MACVelocity.set(_tempMACVelocity);

}

void FluidSim::_advect_particles(float dt) { 

    for(unsigned int p = 0; p < particles.size(); p++) {
        particles[p] = _trace_rk2(particles[p], dt);
    
        //check boundaries and project exterior particles back in
        float phi_val = Interpolation::trilinearInterpolate(particles[p], _dx, _nodal_solid_phi); 
        if(phi_val < 0) {
            vmath::vec3 grad;
            Interpolation::trilinearInterpolateGradient(particles[p], _dx, _nodal_solid_phi, &grad);
            if(vmath::lengthsq(grad) > 0) {
                grad = vmath::normalize(grad);
            }
            particles[p] -= phi_val * grad;
        }
    }
    
}

//Basic first order semi-Lagrangian advection of velocities
void FluidSim::_advect(float dt) {

    _tempMACVelocity.clear();

    //semi-Lagrangian advection on u-component of velocity
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                vmath::vec3 pos = Grid3d::FaceIndexToPositionU(i, j, k, _dx);
                pos = _trace_rk2(pos, -dt);
                _tempMACVelocity.setU(i, j, k, _get_velocity(pos).x);  
            }
        }
    }

    //semi-Lagrangian advection on v-component of velocity
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                vmath::vec3 pos = Grid3d::FaceIndexToPositionV(i, j, k, _dx);
                pos = _trace_rk2(pos, -dt);
                _tempMACVelocity.setV(i, j, k, _get_velocity(pos).y);
            }
        }
    }

    //semi-Lagrangian advection on w-component of velocity
    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) { 
            for(int i = 0; i < _isize; i++) {
                vmath::vec3 pos = Grid3d::FaceIndexToPositionW(i, j, k, _dx);
                pos = _trace_rk2(pos, -dt);
                _tempMACVelocity.setW(i, j, k, _get_velocity(pos).z);
            }
        }
    }

    //move update velocities into u/v vectors
    _MACVelocity.set(_tempMACVelocity);
}

void FluidSim::_compute_phi() {
    
    //grab from particles
    _liquid_phi.fill(3*_dx);
    GridIndex g, gmin, gmax;
    for(unsigned int p = 0; p < particles.size(); ++p) {

        g = Grid3d::positionToGridIndex(particles[p], _dx);
        gmin = GridIndex(max(0, g.i - 1), max(0, g.j - 1), max(0, g.k - 1));
        gmax = GridIndex(min(g.i + 1, _isize - 1), 
                         min(g.j + 1, _jsize - 1), 
                         min(g.k + 1, _ksize - 1));

        for(int k = gmin.k; k <= gmax.k; k++) {
            for(int j = gmin.j; j <= gmax.j; j++) {
                for(int i = gmin.i; i <= gmax.i; i++) {
                    vmath::vec3 sample_pos = Grid3d::GridIndexToCellCenter(i, j, k, _dx);
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
    float u_value = Interpolation::trilinearInterpolate(position - vmath::vec3(0, 0.5*_dx, 0.5*_dx), _dx, *(_MACVelocity.getArray3dU()));
    float v_value = Interpolation::trilinearInterpolate(position - vmath::vec3(0.5*_dx, 0, 0.5*_dx), _dx, *(_MACVelocity.getArray3dV()));
    float w_value = Interpolation::trilinearInterpolate(position - vmath::vec3(0.5*_dx, 0.5*_dx, 0), _dx, *(_MACVelocity.getArray3dW()));

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
                int index = Grid3d::getFlatIndex(i, j, k, _isize, _jsize);

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
                _rhs[index] -= _u_weights(i + 1, j, k) * _MACVelocity.U(i + 1, j, k) / _dx;

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
                _rhs[index] += _u_weights(i, j, k) * _MACVelocity.U(i, j, k) / _dx;

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
                _rhs[index] -= _v_weights(i, j + 1, k) * _MACVelocity.V(i, j + 1, k) / _dx;

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
                _rhs[index] += _v_weights(i, j, k) * _MACVelocity.V(i, j, k) / _dx;


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
                _rhs[index] -= _w_weights(i, j, k + 1) * _MACVelocity.W(i, j, k + 1) / _dx;

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
                _rhs[index] += _w_weights(i, j, k) * _MACVelocity.W(i, j, k) / _dx;
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
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 1; i < _isize; i++) {

                int index = Grid3d::getFlatIndex(i, j, k, _isize, _jsize);
                if(_u_weights(i, j, k) > 0 && (_liquid_phi(i, j, k) < 0 || _liquid_phi(i - 1, j, k) < 0)) {
                    float theta = 1;
                    if(_liquid_phi(i, j, k) >= 0 || _liquid_phi(i - 1, j, k) >= 0) {
                        theta = fraction_inside(_liquid_phi(i-1,j,k), _liquid_phi(i,j,k));
                    }
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    double v = _MACVelocity.U(i, j, k) - dt  * (float)(_pressure[index] - _pressure[index-1]) / _dx / theta;
                    _MACVelocity.setU(i, j, k, v);
                    _u_valid.set(i, j, k, true);
                }

            }
        }
    }
    
    _v_valid.fill(false);
    for(int k = 0; k < _ksize; k++) {
        for(int j = 1; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {

                int index = Grid3d::getFlatIndex(i, j, k, _isize, _jsize);
                if(_v_weights(i, j, k) > 0 && (_liquid_phi(i, j, k) < 0 || _liquid_phi(i, j - 1, k) < 0)) {
                    float theta = 1;
                    if(_liquid_phi(i, j, k) >= 0 || _liquid_phi(i, j - 1, k) >= 0) {
                        theta = fraction_inside(_liquid_phi(i, j - 1, k), _liquid_phi(i, j, k));
                    }
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    double v = _MACVelocity.V(i, j, k) - dt  * (float)(_pressure[index] - _pressure[index-_isize]) / _dx / theta;
                    _MACVelocity.setV(i, j, k, v);
                    _v_valid.set(i, j, k, true);
                }

            }
        }
    }

    _w_valid.fill(false);
    for(int k = 0; k < _ksize; ++k) {
        for(int j = 0; j < _jsize; ++j) {
            for(int i = 1; i < _isize; ++i) {

                int index = Grid3d::getFlatIndex(i, j, k, _isize, _jsize);
                if(_w_weights(i, j, k) > 0 && (_liquid_phi(i, j, k) < 0 || _liquid_phi(i, j, k - 1) < 0)) {
                    float theta = 1;
                    if(_liquid_phi(i, j, k) >= 0 || _liquid_phi(i, j, k - 1) >= 0) {
                        theta = fraction_inside(_liquid_phi(i, j, k - 1), _liquid_phi(i, j, k));
                    }
                    if(theta < 0.01f) {
                        theta = 0.01f;
                    }
                    double v = _MACVelocity.W(i, j, k) - dt  * (float)(_pressure[index] - _pressure[index-_isize*_jsize]) / _dx / theta;
                    _MACVelocity.setW(i, j, k, v);
                    _w_valid.set(i, j, k, true);
                }

            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                if(!_u_valid(i, j, k)) {
                    _MACVelocity.setU(i, j, k, 0.0);
                }
            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                if(!_v_valid(i, j, k)) {
                    _MACVelocity.setV(i, j, k, 0.0);
                }
            }
        }
    }

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if(!_w_valid(i, j, k)) {
                    _MACVelocity.setW(i, j, k, 0.0);
                }
            }
        }
    }
}


//Apply several iterations of a very simple propagation of valid velocity data in all directions
void FluidSim::_extrapolate(Array3d<float> *grid, Array3d<bool> &valid) {

    Array3d<float> temp_grid(grid->width, grid->height, grid->depth);
    for(int k = 0; k < grid->depth; k++) {
        for(int j = 0; j < grid->height; j++) {
            for(int i = 0; i < grid->width; i++) {
                temp_grid.set(i, j, k, grid->get(i, j, k));
            }
        }
    }

    Array3d<bool> old_valid;
    for(int layers = 0; layers < 10; layers++) {

        old_valid = valid;
        for(int k = 1; k < grid->depth - 1; k++) {
            for(int j = 1; j < grid->height - 1; j++) {
                for(int i = 1; i < grid->width - 1; i++) {

                    if(old_valid(i,j,k)) {
                        continue;
                    }

                    float sum = 0;
                    int count = 0;
                    if(old_valid(i + 1, j, k)) {
                        sum += grid->get(i + 1, j, k);
                        count++;
                    }
                    if(old_valid(i - 1, j, k)) {
                        sum += grid->get(i - 1, j, k);
                        count++;
                    }
                    if(old_valid(i, j + 1, k)) {
                        sum += grid->get(i, j + 1, k);
                        count++;
                    }
                    if(old_valid(i, j - 1, k)) {
                        sum += grid->get(i, j - 1, k);
                        count++;
                    }
                    if(old_valid(i, j, k + 1)) {
                        sum += grid->get(i, j, k + 1);
                        count++;
                    }
                    if(old_valid(i, j, k - 1)) {
                        sum += grid->get(i, j, k - 1);
                        count++;
                    }

                    //If any of neighbour cells were valid, 
                    //assign the cell their average value and tag it as valid
                    if(count > 0) {
                        temp_grid.set(i, j, k, sum /(float)count);
                        valid.set(i, j, k, true);
                    }

                }
            }
        }

        for(int k = 0; k < grid->depth; k++) {
            for(int j = 0; j < grid->height; j++) {
                for(int i = 0; i < grid->width; i++) {
                    grid->set(i, j, k, temp_grid(i, j, k));
                }
            }
        }

    }

}
