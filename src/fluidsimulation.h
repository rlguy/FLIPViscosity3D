/*
The MIT License (MIT)

Copyright (c) 2017, Ryan L. Guy

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#ifndef FLUIDSIMULATION_H
#define FLUIDSIMULATION_H

#include "vmath.h"
#include "array3d.h"
#include "grid3d.h"
#include "macvelocityfield.h"
#include "particlelevelset.h"
#include "interpolation.h"
#include "pressuresolver.h"
#include "meshlevelset.h"
#include "fluidsimassert.h"
#include "viscositysolver.h"

#include <vector>

struct FluidParticle {
    vmath::vec3 position;
    vmath::vec3 velocity;

    FluidParticle() {}
    FluidParticle(vmath::vec3 p) : position(p) {}
    FluidParticle(vmath::vec3 p, vmath::vec3 v) : 
                                  position(p),
                                  velocity(v) {}
};

class FluidSimulation {

public:
    void initialize(int i, int j, int k, float dx);
    void addBoundary(TriangleMesh &boundary, bool isInverted = false);
    void resetBoundary();
    void addLiquid(TriangleMesh &mesh);
    void setViscosity(float value);
    void setViscosity(Array3d<float> &vgrid);
    void setGravity(vmath::vec3 gravity);
    void setGravity(float gx, float gy, float gz);
    void advance(float dt);

    std::vector<FluidParticle> particles;

private:
    TriangleMesh _getTriangleMeshFromAABB(AABB bbox);
    TriangleMesh _getBoundaryTriangleMesh();
    void _initializeBoundary();
    vmath::vec3 _traceRK2(vmath::vec3 position, float dt);

    float _cfl();
    vmath::vec3 _getVelocity(vmath::vec3 position);
    void _updateLiquidSDF();

    void _advectVelocityField();
    void _advectVelocityFieldU(Array3d<bool> &fluidCellGrid);
    void _advectVelocityFieldV(Array3d<bool> &fluidCellGrid);
    void _advectVelocityFieldW(Array3d<bool> &fluidCellGrid);
    void _computeVelocityScalarField(Array3d<float> &field, 
                                     Array3d<bool> &isValueSet, 
                                     int dir);

    void _addBodyForce(float dt);
    void _project(float dt);
    void _extrapolateVelocityField(MACVelocityField &vfield, 
                                   ValidVelocityComponentGrid &valid);
    void _constrainVelocityField();

    // viscosity
    void _applyViscosity(float dt);

    //helpers for pressure projection
    void _computeWeights();
    Array3d<float> _solvePressure(float dt);
    void _applyPressure(float dt, Array3d<float> &pressureGrid);

    void _advectFluidParticles(float dt);
    void _updateFluidParticleVelocities();

    inline double _randomDouble(double min, double max) {
        return min + (double)rand() / ((double)RAND_MAX / (max - min));
    }

    template <typename T>
    T _clamp(const T& n, const T& lower, const T& upper) {
        return std::max(lower, std::min(n, upper));
    }

    //Grid dimensions
    int _isize;
    int _jsize;
    int _ksize;
    float _dx;
    
    //Fluid velocity
    MACVelocityField _MACVelocity;
    MACVelocityField _savedVelocityField;
    ValidVelocityComponentGrid _validVelocities;

    MeshLevelSet _solidSDF;
    int _meshLevelSetExactBand = 3;

    ParticleLevelSet _liquidSDF;
    float _particleRadius;

    WeightGrid _weightGrid;

    float _CFLConditionNumber = 5.0;
    float _minfrac = 0.01f;
    float _ratioPICtoFLIP = 0.05f;

    Array3d<float> _viscosity;
    vmath::vec3 _gravity;
};


#endif