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

#include <vector>

class FluidSim {

public:
    void initialize(int i, int j, int k, float width);
    void addBoundary(TriangleMesh &boundary, bool isInverted = false);
    void resetBoundary();
    void set_liquid(float (*phi)(vmath::vec3));
    void add_particle(vmath::vec3 pos);

    void advance(float dt);

    std::vector<vmath::vec3> particles;

private:
    TriangleMesh _getTriangleMeshFromAABB(AABB bbox);
    TriangleMesh _getBoundaryTriangleMesh();
    void _initializeBoundary();
    vmath::vec3 _trace_rk2(vmath::vec3 position, float dt);

    float _cfl();
    vmath::vec3 _get_velocity(vmath::vec3 position);

    void _advect_particles(float dt);
    void _updateLiquidSDF();
    void _advect(float dt);
    void _add_force(float dt);
    void _project(float dt);
    void _constrain_velocity();
    void _extrapolate(Array3d<float> *grid, Array3d<bool> &valid);

    //helpers for pressure projection
    void _compute_weights();
    Array3d<float> _solve_pressure(float dt);
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
    
    Array3d<bool> _u_valid, _v_valid, _w_valid;

    MeshLevelSet _solidSDF;
    double _solidLevelSetExactBand = 3;

    ParticleLevelSet _liquidSDF;
    float _particle_radius;

    WeightGrid _weightGrid;

    double _minfrac = 0.01f;

};


#endif