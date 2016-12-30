#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "vmath.h"
#include "array3d.h"
#include "grid3d.h"
#include "macvelocityfield.h"
#include "interpolation.h"
#include "pressuresolver.h"

#include <vector>

class FluidSim {

public:
    void initialize(int i, int j, int k, float width);
    void set_boundary(float (*phi)(vmath::vec3));
    void set_liquid(float (*phi)(vmath::vec3));
    void add_particle(vmath::vec3 pos);

    void advance(float dt);

    std::vector<vmath::vec3> particles;

private:

    vmath::vec3 _trace_rk2(vmath::vec3 position, float dt);

    float _cfl();
    vmath::vec3 _get_velocity(vmath::vec3 position);

    void _advect_particles(float dt);
    void _advect(float dt);
    void _add_force(float dt);
    void _project(float dt);
    void _constrain_velocity();
    void _extrapolate(Array3d<float> *grid, Array3d<bool> &valid);

    //helpers for pressure projection
    void _compute_weights();
    void _solve_pressure(float dt);
    void _compute_phi();

    inline double _randomDouble(double min, double max) {
        return min + (double)rand() / ((double)RAND_MAX / (max - min));
    }

    //Grid dimensions
    int _isize;
    int _jsize;
    int _ksize;
    float _dx;
    
    //Fluid velocity
    MACVelocityField _MACVelocity;
    MACVelocityField _tempMACVelocity;
    
    //Static geometry representation
    Array3d<float> _nodal_solid_phi;
    Array3d<float> _u_weights, _v_weights, _w_weights;
    Array3d<bool> _u_valid, _v_valid, _w_valid;

    float _particle_radius;

    Array3d<float> _liquid_phi;

    Array3d<float> _pressureGrid;

};


#endif