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

#ifndef VISCOSITYSOLVER_H
#define VISCOSITYSOLVER_H

#include "macvelocityfield.h"
#include "particlelevelset.h"
#include "levelsetutils.h"
#include "meshlevelset.h"
#include "pcgsolver/pcgsolver.h"

struct ViscositySolverParameters {
    float cellwidth;
    float deltaTime;

    MACVelocityField *velocityField;
    ParticleLevelSet *liquidSDF;
    MeshLevelSet *solidSDF;
    Array3d<float> *viscosity;

    //LogFile *logfile;
};

class ViscositySolver {

public:
    ViscositySolver();
    ~ViscositySolver();

    bool applyViscosityToVelocityField(ViscositySolverParameters params);

private:

    struct ViscosityVolumeGrid {
        int isize, jsize, ksize;
        Array3d<float> center;
        Array3d<float> U;
        Array3d<float> V;
        Array3d<float> W;
        Array3d<float> edgeU;
        Array3d<float> edgeV;
        Array3d<float> edgeW;
        
        ViscosityVolumeGrid() {}
        ViscosityVolumeGrid(int i, int j, int k) :
            isize(i), jsize(j), ksize(k),
            center(isize, jsize, ksize, 0.0f),
            U(isize + 1, jsize, ksize, 0.0f),
            V(isize, jsize + 1, ksize, 0.0f),
            W(isize, jsize, ksize + 1, 0.0f),
            edgeU(isize, jsize + 1, ksize + 1, 0.0f),
            edgeV(isize + 1, jsize, ksize + 1, 0.0f),
            edgeW(isize + 1, jsize + 1, ksize, 0.0f) {}

        void destroy() {
            isize = 0;
            jsize = 0;
            ksize = 0;
            center = Array3d<float>(0, 0, 0);
            U = Array3d<float>(0, 0, 0);
            V = Array3d<float>(0, 0, 0);
            W = Array3d<float>(0, 0, 0);
            edgeU = Array3d<float>(0, 0, 0);
            edgeV = Array3d<float>(0, 0, 0);
            edgeW = Array3d<float>(0, 0, 0);
        }
    };

    enum class FaceState : char { 
        air   = 0x00, 
        fluid = 0x01, 
        solid = 0x02
    };

    struct FaceStateGrid {
        int isize, jsize, ksize;
        Array3d<FaceState> U;
        Array3d<FaceState> V;
        Array3d<FaceState> W;
        
        FaceStateGrid() {}
        FaceStateGrid(int i, int j, int k) :
            isize(i), jsize(j), ksize(k),
            U(i + 1, j, k, FaceState::air),
            V(i, j + 1, k, FaceState::air),
            W(i, j, k + 1, FaceState::air) {}
    };

    struct FaceIndexer {
        int isize, jsize, ksize;

        FaceIndexer() {}
        FaceIndexer(int i, int j, int k) : isize(i), jsize(j), ksize(k) {
            _voffset = (isize + 1) * jsize * ksize;
            _woffset = _voffset + isize * (jsize + 1) * ksize;
        }

        int U(int i, int j, int k) {
            return i + (isize + 1) * (j + k * jsize);
        }

        int V(int i, int j, int k) {
            return _voffset + i + isize * (j + k * (jsize + 1));
        }

        int W(int i, int j, int k) {
            return _woffset + i + isize * (j + k * jsize);
        }

        private:

            int _voffset;
            int _woffset;
    };

    struct MatrixIndexer {
        std::vector<int> indexTable;
        FaceIndexer faceIndexer;
        int matrixSize;

        MatrixIndexer() {}
        MatrixIndexer(int i, int j, int k, std::vector<int> matrixIndexTable) :
                indexTable(matrixIndexTable), faceIndexer(i, j, k) {

            int matsize = 0;
            for (size_t idx = 0; idx < indexTable.size(); idx++) {
                if (indexTable[idx] != -1) {
                    matsize++;
                }
            }

            matrixSize = matsize;
        }

        int U(int i, int j, int k) {
            return indexTable[faceIndexer.U(i, j, k)];
        }

        int V(int i, int j, int k) {
            return indexTable[faceIndexer.V(i, j, k)];
        }

        int W(int i, int j, int k) {
            return indexTable[faceIndexer.W(i, j, k)];
        }
    };

    void _initialize(ViscositySolverParameters params);
    void _computeFaceStateGrid();
    void _computeSolidCenterPhi(Array3d<float> &solidCenterPhi);
    void _computeVolumeGrid();
    void _estimateVolumeFractions(Array3d<float> &volumes, 
                                  vmath::vec3 centerStart, 
                                  Array3d<bool> &validCells);
    void _destroyVolumeGrid();
    void _computeMatrixIndexTable();
    void _initializeLinearSystem(SparseMatrixd &matrix, std::vector<double> &rhs);
    void _initializeLinearSystemU(SparseMatrixd &matrix, std::vector<double> &rhs);
    void _initializeLinearSystemV(SparseMatrixd &matrix, std::vector<double> &rhs);
    void _initializeLinearSystemW(SparseMatrixd &matrix, std::vector<double> &rhs);
    bool _solveLinearSystem(SparseMatrixd &matrix, std::vector<double> &rhs, 
                            std::vector<double> &soln);
    void _applySolutionToVelocityField(std::vector<double> &soln);

    int _isize;
    int _jsize;
    int _ksize;
    float _dx;
    float _deltaTime;
    MACVelocityField *_velocityField;
    ParticleLevelSet *_liquidSDF;
    MeshLevelSet *_solidSDF;
    Array3d<float> *_viscosity; 

    FaceStateGrid _state;
    ViscosityVolumeGrid _volumes;
    MatrixIndexer _matrixIndex;

    double _solverTolerance = 1e-6;
    double _acceptableTolerace = 10.0;
    int _maxSolverIterations = 700;

};


#endif

