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
#ifndef PARTICLELEVELSET_H
#define PARTICLELEVELSET_H

#include <vector>

#include "array3d.h"
#include "grid3d.h"
#include "levelsetutils.h"
#include "meshlevelset.h"
#include "fluidsimassert.h"

class ParticleLevelSet {

public:
    ParticleLevelSet();
    ParticleLevelSet(int i, int j, int k, double dx);
    ~ParticleLevelSet();

    float operator()(int i, int j, int k);
    float operator()(GridIndex g);
    float get(int i, int j, int k);
    float get(GridIndex g);
    float getFaceWeightU(int i, int j, int k);
    float getFaceWeightU(GridIndex g);
    float getFaceWeightV(int i, int j, int k);
    float getFaceWeightV(GridIndex g);
    float getFaceWeightW(int i, int j, int k);
    float getFaceWeightW(GridIndex g);

    void calculateSignedDistanceField(std::vector<vmath::vec3> &particles, 
                                      double radius,
                                      MeshLevelSet &solidPhi);
    float trilinearInterpolate(vmath::vec3 pos);
private:

    float _getMaxDistance();
    void _computeSignedDistanceFromParticles(std::vector<vmath::vec3> &particles, 
                                             double radius);
    void _extrapolateSignedDistanceIntoSolids(MeshLevelSet &solidPhi);
    
    int _isize = 0;
    int _jsize = 0;
    int _ksize = 0;
    double _dx = 0.0;
    Array3d<float> _phi;

};

#endif
