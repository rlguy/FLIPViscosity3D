#include "fluidsim.h"

#include "levelsetutils.h"

void FluidSim::initialize(int i, int j, int k, float width) {
    _isize = i;
    _jsize = j;
    _ksize = k;
    _dx = width / (float)_isize;

    _MACVelocity = MACVelocityField(_isize, _jsize, _ksize, _dx);
    _tempMACVelocity = MACVelocityField(_isize, _jsize, _ksize, _dx);

    _u_valid = Array3d<bool>(_isize + 1, _jsize, _ksize);
    _v_valid = Array3d<bool>(_isize, _jsize + 1, _ksize); 
    _w_valid = Array3d<bool>(_isize, _jsize, _ksize + 1);

    _particle_radius = (float)(_dx * 1.01*sqrt(3.0)/2.0); 
    //make the particles large enough so they always appear on the grid

    _liquidSDF = ParticleLevelSet(_isize, _jsize, _ksize, _dx);
    _weightGrid = WeightGrid(_isize, _jsize, _ksize);
}

//Initialize the grid-based signed distance field that dictates the position of the solid boundary
void FluidSim::set_boundary(MeshLevelSet &boundary) {
    int bi, bj, bk;
    boundary.getGridDimensions(&bi, &bj, &bk);
    FLUIDSIM_ASSERT(bi == _isize && bj == _jsize && bk == _ksize);
    _solidSDF = boundary;
}

void FluidSim::set_liquid(float (*phi)(vmath::vec3)) {
    //initialize particles
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) { 
            for(int i = 0; i < _isize; i++) {
                vmath::vec3 gpos = Grid3d::GridIndexToPosition(i, j, k, _dx);

                for (int i_dx = 0; i_dx < 8; i_dx++) {
                    float a = _randomDouble(0.0, _dx); 
                    float b = _randomDouble(0.0, _dx); 
                    float c = _randomDouble(0.0, _dx);
                    vmath::vec3 jitter = vmath::vec3(a, b, c);
                    vmath::vec3 pos = gpos + jitter;

                    if(phi(pos) <= -_particle_radius) {
                        float solid_phi = _solidSDF.trilinearInterpolate(pos);
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
                if(_weightGrid.U(i, j, k) == 0) {
                    //apply constraint
                    vmath::vec3 pos = Grid3d::FaceIndexToPositionU(i, j, k, _dx);
                    vmath::vec3 vel = _get_velocity(pos);
                    vmath::vec3 normal = _solidSDF.trilinearInterpolateGradient(pos);
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
                if(_weightGrid.V(i, j, k) == 0) {
                    //apply constraint
                    vmath::vec3 pos = Grid3d::FaceIndexToPositionV(i, j, k, _dx);
                    vmath::vec3 vel = _get_velocity(pos);
                    vmath::vec3 normal = _solidSDF.trilinearInterpolateGradient(pos);
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
                if(_weightGrid.W(i, j, k) == 0) {
                    //apply constraint
                    vmath::vec3 pos = Grid3d::FaceIndexToPositionW(i, j, k, _dx);
                    vmath::vec3 vel = _get_velocity(pos);
                    vmath::vec3 normal = _solidSDF.trilinearInterpolateGradient(pos);
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
        float phi_val = _solidSDF.trilinearInterpolate(particles[p]);; 
        if(phi_val < 0) {
            vmath::vec3 grad = _solidSDF.trilinearInterpolateGradient(particles[p]);
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
    _liquidSDF.calculateSignedDistanceField(particles, _particle_radius, _solidSDF);
}



void FluidSim::_project(float dt) {

    //Estimate the liquid signed distance
    _compute_phi();
    
    //Compute finite-volume type face area weight for each velocity sample.
    _compute_weights();

    //Set up and solve the variational _pressure solve.
    Array3d<float> pressureGrid = _solve_pressure(dt);
    _applyPressure(dt, pressureGrid);
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
                float weight = 1 - LevelsetUtils::fractionInside(
                                        _solidSDF(i, j, k), 
                                        _solidSDF(i, j + 1, k),
                                        _solidSDF(i, j, k + 1), 
                                        _solidSDF(i, j + 1, k + 1));
                weight = fmax(weight, 0.0);
                weight = fmin(weight, 1.0);
                _weightGrid.U.set(i, j, k, weight);
            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                float weight = 1 - LevelsetUtils::fractionInside(
                                       _solidSDF(i, j, k),
                                       _solidSDF(i, j, k + 1),
                                       _solidSDF(i + 1, j, k),
                                       _solidSDF(i + 1, j, k + 1));
                weight = fmax(weight, 0.0);
                weight = fmin(weight, 1.0);
                _weightGrid.V.set(i, j, k, weight);
            }
        }
    }

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                float weight = 1 - LevelsetUtils::fractionInside(
                                        _solidSDF(i, j, k),
                                        _solidSDF(i, j + 1, k),
                                        _solidSDF(i + 1, j, k),
                                        _solidSDF(i + 1, j + 1, k));
                weight = fmax(weight, 0.0);
                weight = fmin(weight, 1.0);
                _weightGrid.W.set(i, j, k, weight);
            }
        }
    }

}

//An implementation of the variational _pressure projection solve for static geometry
Array3d<float> FluidSim::_solve_pressure(float dt) {
    GridIndexVector pressureCells(_isize, _jsize, _ksize);
    for(int k = 1; k < _ksize - 1; k++) {
        for(int j = 1; j < _jsize - 1; j++) {
            for(int i = 1; i < _isize - 1; i++) {
                if(_liquidSDF(i, j, k) < 0) {
                    pressureCells.push_back(i, j, k);
                }
            }
        }
    }

    PressureSolverParameters params;
    params.cellwidth = _dx;
    params.density = 1.0;
    params.deltaTime = dt;
    params.pressureCells = &pressureCells;
    params.velocityField = &_MACVelocity;
    params.liquidSDF = &_liquidSDF;
    params.weightGrid = &_weightGrid;

    PressureSolver solver;
    return solver.solve(params);
}

void FluidSim::_applyPressure(float dt, Array3d<float> &pressureGrid) {
        //Apply the velocity update
    _u_valid.fill(false);
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 1; i < _isize; i++) {

                //int index = Grid3d::getFlatIndex(i, j, k, _isize, _jsize);
                if(_weightGrid.U(i, j, k) > 0 && (_liquidSDF(i, j, k) < 0 || _liquidSDF(i - 1, j, k) < 0)) {
                    float theta = 1;
                    if(_liquidSDF(i, j, k) >= 0 || _liquidSDF(i - 1, j, k) >= 0) {
                        theta = LevelsetUtils::fractionInside(_liquidSDF(i - 1, j, k), _liquidSDF(i, j, k));
                        theta = fmax(theta, _minfrac);
                    }
                    double v = _MACVelocity.U(i, j, k) - dt  * (float)(pressureGrid(i, j, k) - pressureGrid(i-1, j, k)) / _dx / theta;
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

                //int index = Grid3d::getFlatIndex(i, j, k, _isize, _jsize);
                if(_weightGrid.V(i, j, k) > 0 && (_liquidSDF(i, j, k) < 0 || _liquidSDF(i, j - 1, k) < 0)) {
                    float theta = 1.0;
                    if(_liquidSDF(i, j, k) >= 0 || _liquidSDF(i, j - 1, k) >= 0) {
                        theta = LevelsetUtils::fractionInside(_liquidSDF(i, j - 1, k), _liquidSDF(i, j, k));
                        theta = fmax(theta, _minfrac);
                    }
                    double v = _MACVelocity.V(i, j, k) - dt  * (float)(pressureGrid(i, j, k) - pressureGrid(i, j-1, k)) / _dx / theta;
                    _MACVelocity.setV(i, j, k, v);
                    _v_valid.set(i, j, k, true);
                }

            }
        }
    }

    _w_valid.fill(false);
    for(int k = 1; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {

                //int index = Grid3d::getFlatIndex(i, j, k, _isize, _jsize);
                if(_weightGrid.W(i, j, k) > 0 && (_liquidSDF(i, j, k) < 0 || _liquidSDF(i, j, k - 1) < 0)) {
                    float theta = 1.0;
                    if(_liquidSDF(i, j, k) >= 0 || _liquidSDF(i, j, k - 1) >= 0) {
                        theta = LevelsetUtils::fractionInside(_liquidSDF(i, j, k - 1), _liquidSDF(i, j, k));
                        theta = fmax(theta, _minfrac);
                    }
                    double v = _MACVelocity.W(i, j, k) - dt  * (float)(pressureGrid(i, j, k) - pressureGrid(i, j, k-1)) / _dx / theta;
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
