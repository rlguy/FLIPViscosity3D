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

/*
Part of this levelset implementation was adapted from Christopher Batty's 
signed distance field generator: https://github.com/christopherbatty/SDFGen

The MIT License (MIT)

Copyright (c) 2015, Christopher Batty

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), to 
deal in the Software without restriction, including without limitation the 
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
sell copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO 
EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE 
USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#ifndef MESHLEVELSET_H
#define MESHLEVELSET_H

#include "vmath.h"
#include "grid3d.h"
#include "array3d.h"
#include "trianglemesh.h"
#include "triangle.h"
#include "interpolation.h"
#include "levelsetutils.h"

class MeshLevelSet {

public:
    MeshLevelSet();
    MeshLevelSet(int isize, int jsize, int ksize, double dx);
    ~MeshLevelSet();

    float operator()(int i, int j, int k);
    float operator()(GridIndex g);
    float get(int i, int j, int k);
    float get(GridIndex g);
    float getClosestTriangleIndex(int i, int j, int k);
    float getClosestTriangleIndex(GridIndex g);
    float getDistanceAtCellCenter(int i, int j, int k);
    float getDistanceAtCellCenter(GridIndex g);
    float trilinearInterpolate(vmath::vec3 pos);
    vmath::vec3 trilinearInterpolateGradient(vmath::vec3 pos);
    float getFaceWeightU(int i, int j, int k);
    float getFaceWeightU(GridIndex g);
    float getFaceWeightV(int i, int j, int k);
    float getFaceWeightV(GridIndex g);
    float getFaceWeightW(int i, int j, int k);
    float getFaceWeightW(GridIndex g);

    void getGridDimensions(int *i, int *j, int *k);
    TriangleMesh *getTriangleMesh();

    void calculateSignedDistanceField(TriangleMesh &m, int bandwidth = 1);
    void calculateUnion(MeshLevelSet &levelset);
    void negate();

private:

    void _computeExactBandDistanceField(int bandwidth,
                                          Array3d<int> &intersectionCounts);
    void _propagateDistanceField();
    void _computeDistanceFieldSigns(Array3d<int> &intersectionCounts);
    float _pointToTriangleDistance(vmath::vec3 x0, vmath::vec3 x1, 
                                     vmath::vec3 x2, 
                                     vmath::vec3 x3);
    bool _getBarycentricCoordinates(
              double x0, double y0, 
              double x1, double y1, double x2, double y2, double x3, double y3,
              double *a, double *b, double *c);
    float _pointToSegmentDistance(vmath::vec3 x0, vmath::vec3 x1, vmath::vec3 x2);
    int _orientation(double x1, double y1, double x2, double y2, double *twiceSignedArea);

    template <typename T>
    T _clamp(const T& n, const T& lower, const T& upper) {
      	return std::max(lower, std::min(n, upper));
    }

    int _isize = 0;
    int _jsize = 0;
    int _ksize = 0;
    double _dx = 0.0;

    TriangleMesh _mesh;
    Array3d<float> _phi;
    Array3d<int> _closestTriangles;
};

#endif
