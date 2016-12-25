#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "array3.h"
#include "vec.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"

#include <vector>

class FluidSim {

public:
    void initialize(int i, int j, int k, float width);
    void set_boundary(float (*phi)(const Vec3f&));
    void set_liquid(float (*phi)(const Vec3f&));
    void add_particle(const Vec3f& pos);

    void advance(float dt);

    std::vector<Vec3f> particles;

private:

    Vec3f _trace_rk2(const Vec3f& position, float dt);

    float _cfl();
    Vec3f _get_velocity(const Vec3f& position);

    void _advect_particles(float dt);
    void _advect(float dt);
    void _add_force(float dt);
    void _project(float dt);
    void _constrain_velocity();
    void _extrapolate(Array3f& grid, Array3c& valid);

    //helpers for pressure projection
    void _compute_weights();
    void _solve_pressure(float dt);
    void _compute_phi();

    //Grid dimensions
    int _isize;
    int _jsize;
    int _ksize;
    float _dx;
    
    //Fluid velocity
    Array3f _u, _v, _w;
    Array3f _temp_u, _temp_v, _temp_w;
    
    //Static geometry representation
    Array3f _nodal_solid_phi;
    Array3f _u_weights, _v_weights, _w_weights;
    Array3c _u_valid, _v_valid, _w_valid;

    float _particle_radius;

    Array3f _liquid_phi;

    //Solver data
    PCGSolver<double> _solver;
    SparseMatrixd _matrix;
    std::vector<double> _rhs;
    std::vector<double> _pressure;

};


#endif