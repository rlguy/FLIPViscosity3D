#ifndef FLUIDSIM_H
#define FLUIDSIM_H

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

class FluidSim {

public:
    void initialize(int i, int j, int k, float dx);
    void addBoundary(TriangleMesh &boundary, bool isInverted = false);
    void resetBoundary();
    void addLiquid(TriangleMesh &mesh);
    void setViscosity(float value);
    void setViscosity(Array3d<float> &vgrid);
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

    void _advectVelocityField(float dt);
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
    double _meshLevelSetExactBand = 3;

    ParticleLevelSet _liquidSDF;
    float _particleRadius;

    WeightGrid _weightGrid;

    double _CFLConditionNumber = 5.0;
    double _minfrac = 0.01f;
    double _ratioPICtoFLIP = 0.05f;

    Array3d<float> _viscosity;
};


#endif