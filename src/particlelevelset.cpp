/*
Copyright (c) 2016 Ryan L. Guy

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgement in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#include "particlelevelset.h"

ParticleLevelSet::ParticleLevelSet() {
}

ParticleLevelSet::ParticleLevelSet(int i, int j, int k, double dx) : 
                    _isize(i), _jsize(j), _ksize(k), _dx(dx) {
    _phi = Array3d<float>(i, j, k, _getMaxDistance());
}

ParticleLevelSet::~ParticleLevelSet() {
}

float ParticleLevelSet::operator()(int i, int j, int k) {
    return get(i, j, k);
}

float ParticleLevelSet::operator()(GridIndex g) {
    return get(g);
}

float ParticleLevelSet::get(int i, int j, int k) {
    FLUIDSIM_ASSERT(Grid3d::isGridIndexInRange(i, j, k, _isize, _jsize, _ksize));
    return _phi(i, j, k);
}

float ParticleLevelSet::get(GridIndex g) {
    FLUIDSIM_ASSERT(Grid3d::isGridIndexInRange(g, _isize, _jsize, _ksize));
    return _phi(g);
}

void ParticleLevelSet::calculateSignedDistanceField(std::vector<vmath::vec3> &particles, 
                                                    double radius,
                                                    MeshLevelSet &solidPhi) {
    int si, sj, sk;
    solidPhi.getGridDimensions(&si, &sj, &sk);
    FLUIDSIM_ASSERT(si == _isize && sj == _jsize && sk == _ksize);

    _computeSignedDistanceFromParticles(particles, radius);
    _extrapolateSignedDistanceIntoSolids(solidPhi);
}

float ParticleLevelSet::_getMaxDistance() {
    return 3.0 * _dx;
}

void ParticleLevelSet::_computeSignedDistanceFromParticles(std::vector<vmath::vec3> &particles, 
                                                           double radius) {
    _phi.fill(_getMaxDistance());

    GridIndex g, gmin, gmax;
    vmath::vec3 p;
    for(size_t pidx = 0; pidx < particles.size(); pidx++) {
        p = particles[pidx];
        g = Grid3d::positionToGridIndex(particles[pidx], _dx);
        gmin = GridIndex(fmax(0, g.i - 1), fmax(0, g.j - 1), fmax(0, g.k - 1));
        gmax = GridIndex(fmin(g.i + 1, _isize - 1), 
                         fmin(g.j + 1, _jsize - 1), 
                         fmin(g.k + 1, _ksize - 1));

        for(int k = gmin.k; k <= gmax.k; k++) {
            for(int j = gmin.j; j <= gmax.j; j++) {
                for(int i = gmin.i; i <= gmax.i; i++) {
                    vmath::vec3 cpos = Grid3d::GridIndexToCellCenter(i, j, k, _dx);
                    float dist = vmath::length(cpos - p) - radius;
                    if(dist < _phi(i, j, k)) {
                        _phi.set(i, j, k, dist);
                    }
                }
            }
        }

    }
}

void ParticleLevelSet::_extrapolateSignedDistanceIntoSolids(MeshLevelSet &solidPhi) {
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                if(_phi(i, j, k) < 0.5 * _dx) {
                    if(solidPhi.getDistanceAtCellCenter(i, j, k) < 0) {
                        _phi.set(i, j, k, -0.5f * _dx);
                    }
                }
            }
        }
    }
}