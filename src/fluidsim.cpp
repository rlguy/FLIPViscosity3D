#include "fluidsim.h"

void FluidSim::initialize(int i, int j, int k, float dx) {
    _isize = i;
    _jsize = j;
    _ksize = k;
    _dx = dx;

    _MACVelocity = MACVelocityField(_isize, _jsize, _ksize, _dx);
    _validVelocities = ValidVelocityComponentGrid(_isize, _jsize, _ksize);

    //make the particles large enough so they always appear on the grid
    _particleRadius = (float)(_dx * 1.01*sqrt(3.0)/2.0); 
    _liquidSDF = ParticleLevelSet(_isize, _jsize, _ksize, _dx);
    _weightGrid = WeightGrid(_isize, _jsize, _ksize);
    _viscosity = Array3d<float>(_isize + 1, _jsize + 1, _ksize + 1, 1.0);

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

                    if(meshSDF.trilinearInterpolate(pos) < 0.0) {
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

void FluidSim::setViscosity(float value) {
    FLUIDSIM_ASSERT(value >= 0.0);
    for (int k = 0; k < _viscosity.depth; k++) {
        for (int j = 0; j < _viscosity.height; j++) {
            for (int i = 0; i < _viscosity.width; i++) {
                _viscosity.set(i, j, k, value);
            }
        }
    }
}

void FluidSim::setViscosity(Array3d<float> &vgrid) {
    FLUIDSIM_ASSERT(vgrid.width == _isize + 1 && 
                    vgrid.height == _jsize + 1 && 
                    vgrid.depth == _ksize + 1);

    for (int k = 0; k < _viscosity.depth; k++) {
        for (int j = 0; j < _viscosity.height; j++) {
            for (int i = 0; i < _viscosity.width; i++) {
                float v = vgrid(i, j, k);
                FLUIDSIM_ASSERT(v >= 0.0);
                _viscosity.set(i, j, k, vgrid(i, j, k));
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

        std::cout << "Compute liquid signed distance field" << std::endl;
        _updateLiquidSDF();

        std::cout << "Velocity advection" << std::endl;
        _advectVelocityField(substep);

        std::cout << "Add body force" << std::endl;
        _addBodyForce(substep);

        std::cout << "Apply viscosity" << std::endl;
        _applyViscosity(substep);

        std::cout << "Pressure projection" << std::endl;
        _project(substep); 
     
        std::cout << "Constrain velocity field" << std::endl;
        _constrainVelocityField();

        std::cout << "Advect fluid particles" << std::endl;
        _advectFluidParticles(substep);

        t += substep;
    }
}

void FluidSim::_applyViscosity(float dt) {
    bool isViscosityNonZero = false;
    for (int k = 0; k < _viscosity.depth; k++) {
        for (int j = 0; j < _viscosity.height; j++) {
            for (int i = 0; i < _viscosity.width; i++) {
                if (_viscosity(i, j, k) > 0.0) {
                    isViscosityNonZero = true;
                }
            }
        }
    }

    if (!isViscosityNonZero) {
        return;
    }

    ViscositySolverParameters params;
    params.cellwidth = _dx;
    params.deltaTime = dt;
    params.velocityField = &_MACVelocity;
    params.liquidSDF = &_liquidSDF;
    params.solidSDF = &_solidSDF;
    params.viscosity = &_viscosity;

    ViscositySolver vsolver;
    vsolver.applyViscosityToVelocityField(params);
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
    innerAABB.expand(-3 * _dx - eps);

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
    _solidSDF.calculateSignedDistanceField(boundaryMesh, 50);

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize + 1; j++) { 
            for(int i = 0; i < _isize + 1; i++) {
                vmath::vec3 p = Grid3d::GridIndexToPosition(i, j, k, _dx);
            }
        }
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
    
    return (_CFLConditionNumber * _dx) / maxvel;
}

void FluidSim::_addBodyForce(float dt) {
    Array3d<bool> fgrid(_isize, _jsize, _ksize, false);
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if (_liquidSDF(i, j, k) < 0.0) {
                    fgrid.set(i, j, k, true);
                }
            }
        }
    }

    for(int k = 0;k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                if (Grid3d::isFaceBorderingValueV(i, j, k, true, fgrid)) {
                    _MACVelocity.addV(i, j, k, -9.81f * dt);
                }
            }
        }
    }
}


void FluidSim::_advectFluidParticles(float dt) {

    _updateFluidParticleVelocities();

    AABB boundary(0.0, 0.0, 0.0, _isize * _dx, _jsize *_dx, _ksize * _dx);
    boundary.expand(-2 * _dx - 1e-4);

    for(unsigned int p = 0; p < particles.size(); p++) {
        particles[p].position = _traceRK2(particles[p].position, dt);
    
        //check boundaries and project exterior particles back in
        float phi_val = _solidSDF.trilinearInterpolate(particles[p].position); 
        if(phi_val < 0) {
            vmath::vec3 grad = _solidSDF.trilinearInterpolateGradient(particles[p].position);
            if(vmath::lengthsq(grad) > 0) {
                grad = vmath::normalize(grad);
            }
            particles[p].position -= phi_val * grad;
        }

        if (!boundary.isPointInside(particles[p].position)) {
            particles[p].position = boundary.getNearestPointInsideAABB(particles[p].position);
        }
    }
}

void FluidSim::_updateFluidParticleVelocities() {

    for (size_t i = 0; i < particles.size(); i++) {
        vmath::vec3 p = particles[i].position;
        vmath::vec3 vnew = _MACVelocity.evaluateVelocityAtPositionLinear(p);
        vmath::vec3 vold = _savedVelocityField.evaluateVelocityAtPositionLinear(p);

        vmath::vec3 vPIC = vnew;
        vmath::vec3 vFLIP = particles[i].velocity + vnew - vold;
        particles[i].velocity = _ratioPICtoFLIP * vPIC + (1.0f - _ratioPICtoFLIP) * vFLIP;
    }
}

void FluidSim::_updateLiquidSDF() {
    std::vector<vmath::vec3> points;
    points.reserve(particles.size());
    for (size_t i = 0; i < particles.size(); i++) {
        points.push_back(particles[i].position);
    }

    _liquidSDF.calculateSignedDistanceField(points, _particleRadius, _solidSDF);
}

void FluidSim::_computeVelocityScalarField(Array3d<float> &field, 
                                           Array3d<bool> &isValueSet,
                                           int dir) {
    int U = 0; int V = 1; int W = 2;

    vmath::vec3 offset;
    if (dir == U) {
        offset = vmath::vec3(0.0, 0.5*_dx, 0.5*_dx);
    } else if (dir == V) {
        offset = vmath::vec3(0.5*_dx, 0.0, 0.5*_dx);
    } else if (dir == W) {
        offset = vmath::vec3(0.5*_dx, 0.5*_dx, 0.0);
    } else {
        return;
    }

    Array3d<float> weights(field.width, field.height, field.depth, 0.0);

    // coefficients for Wyvill kernel
    float r = _dx;
    float rsq = r*r;
    float coef1 = (4.0f / 9.0f) * (1.0f / (r*r*r*r*r*r));
    float coef2 = (17.0f / 9.0f) * (1.0f / (r*r*r*r));
    float coef3 = (22.0f / 9.0f) * (1.0f / (r*r));

    // transfer particle velocity component to grid
    for (size_t pidx = 0; pidx < particles.size(); pidx++) {
        vmath::vec3 p = particles[pidx].position - offset;
        double velocityComponent = particles[pidx].velocity[dir];

        GridIndex g = Grid3d::positionToGridIndex(p, _dx);
        GridIndex gmin(fmax(g.i - 1, 0), 
                       fmax(g.j - 1, 0), 
                       fmax(g.k - 1, 0));
        GridIndex gmax(fmin(g.i + 1, field.width - 1), 
                       fmin(g.j + 1, field.height - 1), 
                       fmin(g.k + 1, field.depth - 1));

        for (int k = gmin.k; k <= gmax.k; k++) {
            for (int j = gmin.j; j <= gmax.j; j++) {
                for (int i = gmin.i; i <= gmax.i; i++) {

                    vmath::vec3 gpos = Grid3d::GridIndexToPosition(i, j, k, _dx);
                    vmath::vec3 v = gpos - p;
                    double distsq = vmath::dot(v, v);
                    if (distsq < rsq) {
                        float weight = 1.0f - coef1 * distsq * distsq * distsq + 
                                              coef2 * distsq * distsq - 
                                              coef3 * distsq;
                        field.add(i, j, k, weight * velocityComponent);
                        weights.add(i, j, k, weight);
                    }
                }
            }
        }
    }

    // Divide field values by weights
    double eps = 1e-9;
    for (int k = 0; k < field.depth; k++) {
        for (int j = 0; j < field.height; j++) {
            for (int i = 0; i < field.width; i++) {
                float value = field(i, j, k);
                float weight = weights(i, j, k);

                if (weight < eps) {
                    continue;
                }
                field.set(i, j, k, value / weight);
                isValueSet.set(i, j, k, true);
            }
        }
    }
}

void FluidSim::_advectVelocityFieldU(Array3d<bool> &fluidCellGrid) {
    Array3d<float> ugrid = Array3d<float>(_isize + 1, _jsize, _ksize, 0.0f);
    Array3d<bool> isValueSet = Array3d<bool>(_isize + 1, _jsize, _ksize, false);
    _computeVelocityScalarField(ugrid, isValueSet, 0);

    _MACVelocity.clearU();
    for (int k = 0; k < ugrid.depth; k++) {
        for (int j = 0; j < ugrid.height; j++) {
            for (int i = 0; i < ugrid.width; i++) {
                if (Grid3d::isFaceBorderingValueU(i, j, k, true, fluidCellGrid)) {
                    if (isValueSet(i, j, k)) {
                        _MACVelocity.setU(i, j, k, ugrid(i, j, k));
                        _validVelocities.validU.set(i, j, k, true);
                    }
                }
            }
        }
    }
}

void FluidSim::_advectVelocityFieldV(Array3d<bool> &fluidCellGrid) {
    Array3d<float> vgrid = Array3d<float>(_isize, _jsize + 1, _ksize, 0.0f);
    Array3d<bool> isValueSet = Array3d<bool>(_isize, _jsize + 1, _ksize, false);
    _computeVelocityScalarField(vgrid, isValueSet, 1);
    
    _MACVelocity.clearV();
    for (int k = 0; k < vgrid.depth; k++) {
        for (int j = 0; j < vgrid.height; j++) {
            for (int i = 0; i < vgrid.width; i++) {
                if (Grid3d::isFaceBorderingValueV(i, j, k, true, fluidCellGrid)) {
                    if (isValueSet(i, j, k)) {
                        _MACVelocity.setV(i, j, k, vgrid(i, j, k));
                        _validVelocities.validV.set(i, j, k, true);
                    }
                }
            }
        }
    }
}

void FluidSim::_advectVelocityFieldW(Array3d<bool> &fluidCellGrid) {
    Array3d<float> wgrid = Array3d<float>(_isize, _jsize, _ksize + 1, 0.0f);
    Array3d<bool> isValueSet = Array3d<bool>(_isize, _jsize, _ksize + 1, 0.0f);
    _computeVelocityScalarField(wgrid, isValueSet, 2);
    
    _MACVelocity.clearW();
    for (int k = 0; k < wgrid.depth; k++) {
        for (int j = 0; j < wgrid.height; j++) {
            for (int i = 0; i < wgrid.width; i++) {
                if (Grid3d::isFaceBorderingValueW(i, j, k, true, fluidCellGrid)) {
                    if (isValueSet(i, j, k)) {
                        _MACVelocity.setW(i, j, k, wgrid(i, j, k));
                        _validVelocities.validW.set(i, j, k, true);
                    }
                }
            }
        }
    }
}

void FluidSim::_advectVelocityField(float dt) {
    Array3d<bool> fluidCellGrid(_isize, _jsize, _ksize, false);
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if (_liquidSDF(i, j, k) < 0.0) {
                    fluidCellGrid.set(i, j, k, true);
                }
            }
        }
    }

    _validVelocities.reset();
    _advectVelocityFieldU(fluidCellGrid);
    _advectVelocityFieldV(fluidCellGrid);
    _advectVelocityFieldW(fluidCellGrid);

    _extrapolateVelocityField(_MACVelocity, _validVelocities);
    _savedVelocityField = _MACVelocity;
}


void FluidSim::_project(float dt) {
    //Compute finite-volume type face area weight for each velocity sample.
    _computeWeights();

    //Set up and solve the variational _pressure solve.
    Array3d<float> pressureGrid = _solvePressure(dt);
    _applyPressure(dt, pressureGrid);

    _extrapolateVelocityField(_MACVelocity, _validVelocities);
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
    return _MACVelocity.evaluateVelocityAtPositionLinear(position);
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
    PressureSolverParameters params;
    params.cellwidth = _dx;
    params.density = 1.0;
    params.deltaTime = dt;
    params.velocityField = &_MACVelocity;
    params.liquidSDF = &_liquidSDF;
    params.weightGrid = &_weightGrid;

    PressureSolver solver;
    return solver.solve(params);
}

void FluidSim::_applyPressure(float dt, Array3d<float> &pressureGrid) {
    Array3d<bool> fgrid(_isize, _jsize, _ksize, false);
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if (_liquidSDF(i, j, k) < 0.0) {
                    fgrid.set(i, j, k, true);
                }
            }
        }
    }

    _validVelocities.reset();
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 1; i < _isize; i++) {

                if (_weightGrid.U(i, j, k) > 0 && Grid3d::isFaceBorderingValueU(i, j, k, true, fgrid)) {
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

                if (_weightGrid.V(i, j, k) > 0 && Grid3d::isFaceBorderingValueV(i, j, k, true, fgrid)) {
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

                if (_weightGrid.W(i, j, k) > 0 && Grid3d::isFaceBorderingValueW(i, j, k, true, fgrid)) {
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

void FluidSim::_extrapolateVelocityField(MACVelocityField &vfield, 
                                         ValidVelocityComponentGrid &valid) {
    int numLayers = ceil(_CFLConditionNumber) + 2;
    vfield.extrapolateVelocityField(valid, numLayers);
}

void FluidSim::_constrainVelocityField() {
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                if(_weightGrid.U(i, j, k) == 0) {
                    _MACVelocity.setU(i, j, k, 0.0);
                    _savedVelocityField.setU(i, j, k, 0.0);
                }
            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                if(_weightGrid.V(i, j, k) == 0) {
                    _MACVelocity.setV(i, j, k, 0.0);
                    _savedVelocityField.setV(i, j, k, 0.0);
                }
            }
        }
    }

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) { 
            for(int i = 0; i < _isize; i++) {
                if(_weightGrid.W(i, j, k) == 0) {
                    _MACVelocity.setW(i, j, k, 0.0);
                    _savedVelocityField.setW(i, j, k, 0.0);
                }
            }
        }
    }
}

