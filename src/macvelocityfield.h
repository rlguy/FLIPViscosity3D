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
#ifndef MACVELOCITYFIELD_H
#define MACVELOCITYFIELD_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <limits>
#include <time.h>

#include "array3d.h"
#include "grid3d.h"
#include "interpolation.h"
#include "vmath.h" 
#include "fluidsimassert.h"

struct ValidVelocityComponentGrid {
    Array3d<bool> validU;
    Array3d<bool> validV;
    Array3d<bool> validW;

    ValidVelocityComponentGrid() {}
    ValidVelocityComponentGrid(int i, int j, int k) :
                               validU(i + 1, j, k, false), 
                               validV(i, j + 1, k, false), 
                               validW(i, j, k + 1, false) {}
    void reset() {
        validU.fill(false);
        validV.fill(false);
        validW.fill(false);
    }
};

class MACVelocityField
{
public:
    MACVelocityField();
    MACVelocityField(int isize, int jsize, int ksize, double dx);
    ~MACVelocityField();

    void getGridDimensions(int *i, int *j, int *k);
    double getGridCellSize();

    float U(int i, int j, int k);
    float V(int i, int j, int k);
    float W(int i, int j, int k);
    float U(GridIndex g);
    float V(GridIndex g);
    float W(GridIndex g);

    void set(MACVelocityField &vfield);
    void setU(int i, int j, int k, double val);
    void setV(int i, int j, int k, double val);
    void setW(int i, int j, int k, double val);
    void setU(GridIndex g, double val);
    void setV(GridIndex g, double val);
    void setW(GridIndex g, double val);
    void setU(Array3d<float> &ugrid);
    void setV(Array3d<float> &vgrid);
    void setW(Array3d<float> &wgrid);
    void addU(int i, int j, int k, double val);
    void addV(int i, int j, int k, double val);
    void addW(int i, int j, int k, double val);

    Array3d<float>* getArray3dU();
    Array3d<float>* getArray3dV();
    Array3d<float>* getArray3dW();

    float* getRawArrayU();
    float* getRawArrayV();
    float* getRawArrayW();

    void clear();
    void clearU();
    void clearV();
    void clearW();

    inline bool isIndexInRangeU(int i, int j, int k) {
        return Grid3d::isGridIndexInRange(i, j, k, _isize + 1, _jsize, _ksize);
    }
    inline bool isIndexInRangeV(int i, int j, int k) {
        return Grid3d::isGridIndexInRange(i, j, k, _isize, _jsize + 1, _ksize);
    }
    inline bool isIndexInRangeW(int i, int j, int k) {
        return Grid3d::isGridIndexInRange(i, j, k, _isize, _jsize, _ksize + 1);
    }
    inline bool isIndexInRangeU(GridIndex g) {
        return Grid3d::isGridIndexInRange(g, _isize + 1, _jsize, _ksize);
    }
    inline bool isIndexInRangeV(GridIndex g) {
        return Grid3d::isGridIndexInRange(g, _isize, _jsize + 1, _ksize);
    }
    inline bool isIndexInRangeW(GridIndex g) {
        return Grid3d::isGridIndexInRange(g, _isize, _jsize, _ksize + 1);
    }

    vmath::vec3 evaluateVelocityAtCellCenter(int i, int j, int k);
    float evaluateVelocityMagnitudeAtCellCenter(int i, int j, int k);
    float evaluateVelocityMagnitudeSquaredAtCellCenter(int i, int j, int k);
    float evaluateMaximumVelocityMagnitude();

    vmath::vec3 evaluateVelocityAtFaceCenterU(int i, int j, int k);
    vmath::vec3 evaluateVelocityAtFaceCenterV(int i, int j, int k);
    vmath::vec3 evaluateVelocityAtFaceCenterW(int i, int j, int k);

    vmath::vec3 evaluateVelocityAtPosition(double x, double y, double z);
    vmath::vec3 evaluateVelocityAtPosition(vmath::vec3 pos);
    vmath::vec3 evaluateVelocityAtPositionLinear(double x, double y, double z);
    vmath::vec3 evaluateVelocityAtPositionLinear(vmath::vec3 pos);

    void extrapolateVelocityField(ValidVelocityComponentGrid &validGrid, int numLayers);

private:
    void _initializeVelocityGrids();

    float _default_out_of_range_value = 0.0f;

    inline double _randomFloat(double min, double max) {
         return min + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (max - min))); 
    }

    double _interpolateU(double x, double y, double z);
    double _interpolateV(double x, double y, double z);
    double _interpolateW(double x, double y, double z);
    double _interpolateLinearU(double x, double y, double z);
    double _interpolateLinearV(double x, double y, double z);
    double _interpolateLinearW(double x, double y, double z);

    void _extrapolateGrid(Array3d<float> &grid, Array3d<bool> &valid, int numLayers);


    int _isize = 10;
    int _jsize = 10;
    int _ksize = 10;
    double _dx = 0.1;

    Array3d<float> _u;
    Array3d<float> _v;
    Array3d<float> _w;

    int _numExtrapolationLayers = 0;
};

#endif
