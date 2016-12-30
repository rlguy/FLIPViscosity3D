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
    This pressure solver is based on:

    A Fast Variational Framework for Accurate Solid-Fluid Coupling
     - Christopher Batty, Florence Bertails, Robert Bridson
    https://hal.archives-ouvertes.fr/file/index/docid/384725/filename/variationalFluids.pdf

    Source Code: https://github.com/christopherbatty/Fluid3D
*/
#ifndef PRESSURESOLVER_H
#define PRESSURESOLVER_H

#include <vector>
#include <stdio.h>
#include <iostream>
#include <limits>
#include <algorithm>

#include "macvelocityfield.h"
#include "gridindexkeymap.h"
//#include "logfile.h"
#include "grid3d.h"
#include "array3d.h"
#include "gridindexvector.h"
#include "fluidsimassert.h"

#include "levelset_util.h"

struct WeightGrid {
    Array3d<float> U;
    Array3d<float> V;
    Array3d<float> W;

    WeightGrid() {}
    WeightGrid(int i, int j, int k) :
        U(i + 1, j, k, 0.0f),
        V(i, j + 1, k, 0.0f),
        W(i, j, k + 1, 0.0f) {}
};


struct PressureSolverParameters {
    double cellwidth;
    double density;
    double deltaTime;

    GridIndexVector *pressureCells;
    MACVelocityField *velocityField;
    Array3d<float> *liquidSDF;
    WeightGrid *weightGrid;

    //LogFile *logfile;
};

/********************************************************************************
    VectorXd
********************************************************************************/

class VectorXd
{
public:
    VectorXd();
    VectorXd(int size);
    VectorXd(int size, double fill);
    VectorXd(VectorXd &vector);
    ~VectorXd();

    const double operator [](int i) const;
    double& operator[](int i);

    inline size_t size() {
        return _vector.size();
    }

    void fill(double fill);
    double dot(VectorXd &vector);
    double absMaxCoeff();

    std::vector<double> _vector;

};

/********************************************************************************
    MatrixCoefficients
********************************************************************************/

struct MatrixCell {
    float diag;
    float plusi;
    float plusj;
    float plusk;

    MatrixCell() : diag(0.0f), plusi(0.0f), plusj(0.0f), plusk(0.0f) {}
};

class MatrixCoefficients
{
public:
    MatrixCoefficients();
    MatrixCoefficients(int size);
    ~MatrixCoefficients();

    const MatrixCell operator [](int i) const;
    MatrixCell& operator [](int i);

    inline size_t size() {
        return cells.size();
    }

    std::vector<MatrixCell> cells;
};

/********************************************************************************
    PressureSolver
********************************************************************************/

class PressureSolver
{
public:
    PressureSolver();
    ~PressureSolver();

    Array3d<float> solve(PressureSolverParameters params);

private:

    inline int _GridToVectorIndex(GridIndex g) {
        return _keymap.find(g);
    }
    inline int _GridToVectorIndex(int i, int j, int k) {
        return _keymap.find(i, j, k);
    }
    inline GridIndex _VectorToGridIndex(int i) {
        return _pressureCells->at(i);
    }
    inline int _isPressureCell(GridIndex g) {
        return _keymap.find(g) != -1;
    }
    inline int _isPressureCell(int i, int j, int k) {
        return _keymap.find(i, j, k) != -1;
    }

    void _initialize(PressureSolverParameters params);
    void _initializeGridIndexKeyMap();
    void _calculateNegativeDivergenceVector(VectorXd &b);
    void _calculateMatrixCoefficients(MatrixCoefficients &A);
    void _calculatePreconditionerVector(MatrixCoefficients &A, VectorXd &precon);
    void _solvePressureSystem(MatrixCoefficients &A, 
                              VectorXd &b, 
                              VectorXd &precon,
                              VectorXd &pressure);
    void _applyPreconditioner(MatrixCoefficients &A, 
                              VectorXd &precon,
                              VectorXd &residual,
                              VectorXd &vect);
    void _applyMatrix(MatrixCoefficients &A, VectorXd &x, VectorXd &result);
    void _addScaledVector(VectorXd &v1, VectorXd &v2, double scale);
    void _addScaledVectors(VectorXd &v1, double s1, 
                           VectorXd &v2, double s2,
                           VectorXd &result);

    int _isize = 0;
    int _jsize = 0;
    int _ksize = 0;
    double _dx = 0;
    double _density = 0;
    double _deltaTime = 0;
    int _matSize = 0;

    double _pressureSolveTolerance = 1e-18;
    int _maxCGIterations = 200;
    double _minfrac = 0.01f;

    GridIndexVector *_pressureCells;
    MACVelocityField *_vField;
    Array3d<float> *_liquidSDF;
    WeightGrid *_weightGrid;

    //LogFile *_logfile;
    GridIndexKeyMap _keymap;

};

#endif