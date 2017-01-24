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
#include "stopwatch.h"
#include "viscositysolver.h"

#include <vector>

class FluidSim {

public:
    void initialize(int i, int j, int k, float dx);
    void addBoundary(TriangleMesh &boundary, bool isInverted = false);
    void resetBoundary();
    void addLiquid(TriangleMesh &mesh);
    void setViscosity(float value);
    void setViscosity(Array3d<float> &vgrid);
    void advance(float dt);

    std::vector<vmath::vec3> particles;

private:
    TriangleMesh _getTriangleMeshFromAABB(AABB bbox);
    TriangleMesh _getBoundaryTriangleMesh();
    void _initializeBoundary();
    vmath::vec3 _traceRK2(vmath::vec3 position, float dt);

    float _cfl();
    vmath::vec3 _getVelocity(vmath::vec3 position);

    void _advectParticles(float dt);
    void _updateLiquidSDF();
    void _advectVelocityField(float dt);
    void _addBodyForce(float dt);
    void _project(float dt);
    void _extrapolateVelocityField();
    void _constrainVelocityField();

    // viscosity
    void _applyViscosity(float dt);

    //helpers for pressure projection
    void _computeWeights();
    Array3d<float> _solvePressure(float dt);
    void _applyPressure(float dt, Array3d<float> &pressureGrid);

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
    MACVelocityField _tempMACVelocity;
    
    ValidVelocityComponentGrid _validVelocities;
    int _numExtrapolationLayers = 10;

    MeshLevelSet _solidSDF;
    double _meshLevelSetExactBand = 3;

    ParticleLevelSet _liquidSDF;
    float _particleRadius;

    WeightGrid _weightGrid;

    double _minfrac = 0.01f;

    Array3d<float> _viscosity;
};


#endif