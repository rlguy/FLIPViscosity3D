#include "fluidsim.h"

void FluidSim::initialize(int i, int j, int k, float width) {
    _isize = i;
    _jsize = j;
    _ksize = k;
    _dx = width / (float)_isize;

    _MACVelocity = MACVelocityField(_isize, _jsize, _ksize, _dx);
    _tempMACVelocity = MACVelocityField(_isize, _jsize, _ksize, _dx);
    _validVelocities = ValidVelocityComponentGrid(_isize, _jsize, _ksize);

    //make the particles large enough so they always appear on the grid
    _particleRadius = (float)(_dx * 1.01*sqrt(3.0)/2.0); 
    _liquidSDF = ParticleLevelSet(_isize, _jsize, _ksize, _dx);
    _weightGrid = WeightGrid(_isize, _jsize, _ksize);

    _initializeBoundary();
}

void FluidSim::addBoundary(TriangleMesh &boundary, bool isInverted) {
    AABB domain(0.0, 0.0, 0.0, _isize * _dx, _jsize * _dx, _ksize *_dx);
    AABB bbox(boundary.vertices);
    FLUIDSIM_ASSERT(domain.isPointInside(bbox.getMinPoint()) &&
                    domain.isPointInside(bbox.getMaxPoint()));

    MeshLevelSet boundarySDF(_isize, _jsize, _ksize, _dx);
    boundarySDF.calculateSignedDistanceField(boundary, _meshLevelSetExactBand);
    if (isInverted) {
        boundarySDF.negate();
    }

    _solidSDF.calculateUnion(boundarySDF);
}

void FluidSim::resetBoundary() {
    _initializeBoundary();
}

void FluidSim::addLiquid(TriangleMesh &mesh) {
    AABB domain(0.0, 0.0, 0.0, _isize * _dx, _jsize * _dx, _ksize *_dx);
    AABB bbox(mesh.vertices);
    FLUIDSIM_ASSERT(domain.isPointInside(bbox.getMinPoint()) &&
                    domain.isPointInside(bbox.getMaxPoint()));

    MeshLevelSet meshSDF(_isize, _jsize, _ksize, _dx);
    meshSDF.calculateSignedDistanceField(mesh, _meshLevelSetExactBand);

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

                    if(meshSDF.trilinearInterpolate(pos) <= -_particleRadius) {
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
        _advectParticles(substep);

        printf(" Compute liquid signed distance field\n");
        _updateLiquidSDF();

        printf(" Velocity advection\n");
        //Advance the velocity
        _advectVelocityField(substep);
        _addBodyForce(substep);

        printf(" Pressure projection\n");
        _project(substep); 
         
        //_pressure projection only produces valid velocities in faces with non-zero associated face area.
        //Because the advection step may interpolate from these invalid faces, 
        //we must extrapolate velocities from the fluid domain into these invalid faces.
        printf(" Extrapolation\n");
        _extrapolateVelocityField();
     
        //For extrapolated velocities, replace the normal component with
        //that of the object.
        printf(" Constrain boundary velocities\n");
        _constrainVelocityField();

        t += substep;
    }
}

TriangleMesh FluidSim::_getTriangleMeshFromAABB(AABB bbox) {
    vmath::vec3 p = bbox.position;
    std::vector<vmath::vec3> verts{
        vmath::vec3(p.x, p.y, p.z),
        vmath::vec3(p.x + bbox.width, p.y, p.z),
        vmath::vec3(p.x + bbox.width, p.y, p.z + bbox.depth),
        vmath::vec3(p.x, p.y, p.z + bbox.depth),
        vmath::vec3(p.x, p.y + bbox.height, p.z),
        vmath::vec3(p.x + bbox.width, p.y + bbox.height, p.z),
        vmath::vec3(p.x + bbox.width, p.y + bbox.height, p.z + bbox.depth),
        vmath::vec3(p.x, p.y + bbox.height, p.z + bbox.depth)
    };

    std::vector<Triangle> tris{
        Triangle(0, 1, 2), Triangle(0, 2, 3), Triangle(4, 7, 6), Triangle(4, 6, 5),
        Triangle(0, 3, 7), Triangle(0, 7, 4), Triangle(1, 5, 6), Triangle(1, 6, 2),
        Triangle(0, 4, 5), Triangle(0, 5, 1), Triangle(3, 2, 6), Triangle(3, 6, 7)
    };


    TriangleMesh m;
    m.vertices = verts;
    m.triangles = tris;

    return m;
}

TriangleMesh FluidSim::_getBoundaryTriangleMesh() {
    double eps = 1e-6;
    AABB domainAABB(0.0, 0.0, 0.0, _isize * _dx, _jsize * _dx, _ksize * _dx);
    AABB outerAABB = domainAABB;
    outerAABB.expand(-eps);
    AABB innerAABB = domainAABB;
    innerAABB.expand(-2 * _dx - eps);

    TriangleMesh domainMesh = _getTriangleMeshFromAABB(outerAABB);
    TriangleMesh innerMesh = _getTriangleMeshFromAABB(innerAABB);
    int indexOffset = domainMesh.vertices.size();
    domainMesh.vertices.insert(domainMesh.vertices.end(), 
                               innerMesh.vertices.begin(), innerMesh.vertices.end());
    for (size_t i = 0; i < innerMesh.triangles.size(); i++) {
        Triangle t = innerMesh.triangles[i];
        t.tri[0] += indexOffset;
        t.tri[1] += indexOffset;
        t.tri[2] += indexOffset;
        domainMesh.triangles.push_back(t);
    }

    return domainMesh;
}

void FluidSim::_initializeBoundary() {
    TriangleMesh boundaryMesh = _getBoundaryTriangleMesh();
    _solidSDF = MeshLevelSet(_isize, _jsize, _ksize, _dx);
    _solidSDF.calculateSignedDistanceField(boundaryMesh, _meshLevelSetExactBand);
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

void FluidSim::_addBodyForce(float dt) {
    FluidMaterialGrid mgrid(_isize, _jsize, _ksize);
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if (_liquidSDF(i, j, k) < 0.0) {
                    mgrid.setFluid(i, j, k);
                }
            }
        }
    }

    for(int k = 0;k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                if (mgrid.isFaceBorderingFluidV(i, j, k)) {
                    _MACVelocity.addV(i, j, k, -9.81f * dt);
                }
            }
        }
    }
}


void FluidSim::_advectParticles(float dt) {

    AABB boundary(0.0, 0.0, 0.0, _isize * _dx, _jsize *_dx, _ksize * _dx);
    boundary.expand(-2 * _dx - 1e-4);

    for(unsigned int p = 0; p < particles.size(); p++) {
        particles[p] = _traceRK2(particles[p], dt);
    
        //check boundaries and project exterior particles back in
        float phi_val = _solidSDF.trilinearInterpolate(particles[p]);; 
        if(phi_val < 0) {
            vmath::vec3 grad = _solidSDF.trilinearInterpolateGradient(particles[p]);
            if(vmath::lengthsq(grad) > 0) {
                grad = vmath::normalize(grad);
            }
            particles[p] -= phi_val * grad;
        }

        if (!boundary.isPointInside(particles[p])) {
            particles[p] = boundary.getNearestPointInsideAABB(particles[p]);
        }
    }
    
}

void FluidSim::_updateLiquidSDF() {
    _liquidSDF.calculateSignedDistanceField(particles, _particleRadius, _solidSDF);
}

//Basic first order semi-Lagrangian advection of velocities
void FluidSim::_advectVelocityField(float dt) {

    FluidMaterialGrid mgrid(_isize, _jsize, _ksize);
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if (_liquidSDF(i, j, k) < 0.0) {
                    mgrid.setFluid(i, j, k);
                }
            }
        }
    }
    _tempMACVelocity.clear();

    //semi-Lagrangian advection on u-component of velocity
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                if (mgrid.isFaceBorderingFluidU(i, j, k)) {
                    vmath::vec3 pos = Grid3d::FaceIndexToPositionU(i, j, k, _dx);
                    pos = _traceRK2(pos, -dt);
                    _tempMACVelocity.setU(i, j, k, _getVelocity(pos).x); 
                } 
            }
        }
    }

    //semi-Lagrangian advection on v-component of velocity
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                if (mgrid.isFaceBorderingFluidV(i, j, k)) {
                    vmath::vec3 pos = Grid3d::FaceIndexToPositionV(i, j, k, _dx);
                    pos = _traceRK2(pos, -dt);
                    _tempMACVelocity.setV(i, j, k, _getVelocity(pos).y);
                }
            }
        }
    }

    //semi-Lagrangian advection on w-component of velocity
    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) { 
            for(int i = 0; i < _isize; i++) {
                if (mgrid.isFaceBorderingFluidW(i, j, k)) {
                    vmath::vec3 pos = Grid3d::FaceIndexToPositionW(i, j, k, _dx);
                    pos = _traceRK2(pos, -dt);
                    _tempMACVelocity.setW(i, j, k, _getVelocity(pos).z);
                }
            }
        }
    }

    //move update velocities into u/v vectors
    _MACVelocity.set(_tempMACVelocity);
}


void FluidSim::_project(float dt) {
    //Compute finite-volume type face area weight for each velocity sample.
    _computeWeights();

    //Set up and solve the variational _pressure solve.
    Array3d<float> pressureGrid = _solvePressure(dt);
    _applyPressure(dt, pressureGrid);
}


//Apply RK2 to advect a point in the domain.
vmath::vec3 FluidSim::_traceRK2(vmath::vec3 position, float dt) {
    vmath::vec3 input = position;
    vmath::vec3 velocity = _getVelocity(input);
    velocity = _getVelocity(input + 0.5f * dt * velocity);
    input += dt * velocity;
    return input;
}

//Interpolate velocity from the MAC grid.
vmath::vec3 FluidSim::_getVelocity(vmath::vec3 position) {
    float u_value = Interpolation::trilinearInterpolate(position - vmath::vec3(0, 0.5*_dx, 0.5*_dx), _dx, *(_MACVelocity.getArray3dU()));
    float v_value = Interpolation::trilinearInterpolate(position - vmath::vec3(0.5*_dx, 0, 0.5*_dx), _dx, *(_MACVelocity.getArray3dV()));
    float w_value = Interpolation::trilinearInterpolate(position - vmath::vec3(0.5*_dx, 0.5*_dx, 0), _dx, *(_MACVelocity.getArray3dW()));

    return vmath::vec3(u_value, v_value, w_value);
}



//Compute finite-volume style face-weights for fluid from nodal signed distances
void FluidSim::_computeWeights() {

    //Compute face area fractions (using marching squares cases).
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                float weight = 1.0f - _solidSDF.getFaceWeightU(i, j, k);
                weight = _clamp(weight, 0.0f, 1.0f);
                _weightGrid.U.set(i, j, k, weight);
            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                float weight = 1.0f - _solidSDF.getFaceWeightV(i, j, k);
                weight = _clamp(weight, 0.0f, 1.0f);
                _weightGrid.V.set(i, j, k, weight);
            }
        }
    }

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                float weight = 1.0f - _solidSDF.getFaceWeightW(i, j, k);
                weight = _clamp(weight, 0.0f, 1.0f);
                _weightGrid.W.set(i, j, k, weight);
            }
        }
    }

}

//An implementation of the variational _pressure projection solve for static geometry
Array3d<float> FluidSim::_solvePressure(float dt) {
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
    FluidMaterialGrid mgrid(_isize, _jsize, _ksize);
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if (_liquidSDF(i, j, k) < 0.0) {
                    mgrid.setFluid(i, j, k);
                }
            }
        }
    }

    _validVelocities.reset();
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 1; i < _isize; i++) {

                if (_weightGrid.U(i, j, k) > 0 && mgrid.isFaceBorderingFluidU(i, j, k)) {
                    float p0 = pressureGrid(i-1, j, k);
                    float p1 = pressureGrid(i, j, k);
                    float theta = fmax(_liquidSDF.getFaceWeightU(i, j, k), _minfrac);
                    _MACVelocity.addU(i, j, k, -dt * (p1 - p0) / (_dx * theta));
                    _validVelocities.validU.set(i, j, k, true);
                }

            }
        }
    }
    
    for(int k = 0; k < _ksize; k++) {
        for(int j = 1; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {

                if (_weightGrid.V(i, j, k) > 0 && mgrid.isFaceBorderingFluidV(i, j, k)) {
                    float p0 = pressureGrid(i, j - 1, k);
                    float p1 = pressureGrid(i, j, k);
                    float theta = fmax(_liquidSDF.getFaceWeightV(i, j, k), _minfrac);
                    _MACVelocity.addV(i, j, k, -dt * (p1 - p0) / (_dx * theta));
                    _validVelocities.validV.set(i, j, k, true);
                }

            }
        }
    }

    for(int k = 1; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {

                if (_weightGrid.W(i, j, k) > 0 && mgrid.isFaceBorderingFluidW(i, j, k)) {
                    float p0 = pressureGrid(i, j, k - 1);
                    float p1 = pressureGrid(i, j, k);
                    float theta = fmax(_liquidSDF.getFaceWeightW(i, j, k), _minfrac);
                    _MACVelocity.addW(i, j, k, -dt * (p1 - p0) / (_dx * theta));
                    _validVelocities.validW.set(i, j, k, true);
                }

            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                if(!_validVelocities.validU(i, j, k)) {
                    _MACVelocity.setU(i, j, k, 0.0);
                }
            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                if(!_validVelocities.validV(i, j, k)) {
                    _MACVelocity.setV(i, j, k, 0.0);
                }
            }
        }
    }

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if(!_validVelocities.validW(i, j, k)) {
                    _MACVelocity.setW(i, j, k, 0.0);
                }
            }
        }
    }
}

void FluidSim::_extrapolateVelocityField() {
    _MACVelocity.extrapolateVelocityField(_validVelocities, _numExtrapolationLayers);
}

//For extrapolated points, replace the normal component
//of velocity with the object velocity (in this case zero).
void FluidSim::_constrainVelocityField() {
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
                    vmath::vec3 vel = _getVelocity(pos);
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
                    vmath::vec3 vel = _getVelocity(pos);
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
                    vmath::vec3 vel = _getVelocity(pos);
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

