/*
    Viscosity solver adapted from Christopher Batty's viscosity3d.cpp:
        https://github.com/christopherbatty/VariationalViscosity3D/blob/master/viscosity3d.cpp

    Accurate Viscous Free Surfaces for Buckling, Coiling, and Rotating Liquids
    C. Batty and R. Bridson
    http://www.cs.ubc.ca/nest/imager/tr/2008/Batty_ViscousFluids/viscosity.pdf
*/

#include "viscositysolver.h"

ViscositySolver::ViscositySolver() {
}

ViscositySolver::~ViscositySolver() {
}

bool ViscositySolver::applyViscosityToVelocityField(ViscositySolverParameters params) {
    _initialize(params);
    _computeFaceStateGrid();
    _computeVolumeGrid();
    _computeMatrixIndexTable();

    int matsize = _matrixIndex.matrixSize;
    SparseMatrixd matrix(matsize);
    std::vector<double> rhs(matsize, 0);
    std::vector<double> soln(matsize, 0);

    _initializeLinearSystem(matrix, rhs);
    _destroyVolumeGrid();

    bool success = _solveLinearSystem(matrix, rhs, soln);
    if (!success) {
        return false;
    }

    _applySolutionToVelocityField(soln);

    return true;
}

void ViscositySolver::_initialize(ViscositySolverParameters params) {
    int isize, jsize, ksize;
    params.velocityField->getGridDimensions(&isize, &jsize, &ksize);

    _isize = isize;
    _jsize = jsize;
    _ksize = ksize;
    _dx = params.cellwidth;
    _deltaTime = params.deltaTime;
    _velocityField = params.velocityField;
    _liquidSDF = params.liquidSDF;
    _solidSDF = params.solidSDF;
    _viscosity = params.viscosity;
}

void ViscositySolver::_computeFaceStateGrid() {
    Array3d<float> solidCenterPhi(_isize, _jsize, _ksize);
    _computeSolidCenterPhi(solidCenterPhi);

    _state = FaceStateGrid(_isize, _jsize, _ksize);
    for (int k = 0; k < _state.U.depth; k++) {
        for (int j = 0; j < _state.U.height; j++) { 
            for (int i = 0; i < _state.U.width; i++) {
                bool isEdge = i == 0 || i == _state.U.width - 1;;
                if (isEdge || solidCenterPhi(i - 1, j, k) + solidCenterPhi(i, j, k) <= 0) {
                    _state.U.set(i, j, k, FaceState::solid);
                } else { 
                    _state.U.set(i, j, k, FaceState::fluid);
                }
            }
        }
    }

    for (int k = 0; k < _state.V.depth; k++) {
        for (int j = 0; j < _state.V.height; j++) {
            for (int i = 0; i < _state.V.width; i++) {
                bool isEdge = j == 0 || j == _state.V.height - 1;
                if (isEdge || solidCenterPhi(i, j - 1, k) + solidCenterPhi(i, j, k) <= 0) {
                    _state.V.set(i, j, k, FaceState::solid);
                } else { 
                    _state.V.set(i, j, k, FaceState::fluid);
                }
            }
        }
    }

    for (int k = 0; k < _state.W.depth; k++) {
        for (int j = 0; j < _state.W.height; j++) { 
            for (int i = 0; i < _state.W.width; i++) {
                bool isEdge = k == 0 || k == _state.W.depth - 1;
                if (isEdge || solidCenterPhi(i, j, k - 1) + solidCenterPhi(i, j, k) <= 0) {
                    _state.W.set(i, j, k, FaceState::solid);
                } else { 
                    _state.W.set(i, j, k, FaceState::fluid); 
                }
            }
        }
    }
}

void ViscositySolver::_computeSolidCenterPhi(Array3d<float> &solidCenterPhi) {
    for(int k = 0; k < solidCenterPhi.depth; k++) {
        for(int j = 0; j < solidCenterPhi.height; j++) { 
            for(int i = 0; i < solidCenterPhi.width; i++) {
                solidCenterPhi.set(i, j, k, _solidSDF->getDistanceAtCellCenter(i, j, k));
            }
        }
    }
}

void ViscositySolver::_computeVolumeGrid() {
    _volumes = ViscosityVolumeGrid(_isize, _jsize, _ksize);

    Array3d<bool> validCells(_isize + 1, _jsize + 1, _ksize + 1, false);
    for (int k = 0; k < _ksize; k++) {
        for (int j = 0; j < _jsize; j++) {
            for (int i = 0; i < _isize; i++) {
                if (_liquidSDF->get(i, j, k) < 0) {
                    validCells.set(i, j, k, true);
                }
            }
        }
    }

    int layers = 2;
    for (int layer = 0; layer < layers; layer++) {
        GridIndex nbs[6];
        Array3d<bool> tempValid = validCells;
        for (int k = 0; k < _ksize + 1; k++) {
            for (int j = 0; j < _jsize + 1; j++) {
                for (int i = 0; i < _isize + 1; i++) {
                    if (validCells(i, j, k)) {
                        Grid3d::getNeighbourGridIndices6(i, j, k, nbs);
                        for (int nidx = 0; nidx < 6; nidx++) {
                            if (tempValid.isIndexInRange(nbs[nidx])) {
                                tempValid.set(nbs[nidx], true);
                            }
                        }
                    }
                }
            }
        }
        validCells = tempValid;
    }

    float hdx = 0.5 * _dx;
    _estimateVolumeFractions(_volumes.center, vmath::vec3(hdx, hdx, hdx), validCells); 
    _estimateVolumeFractions(_volumes.U,      vmath::vec3(0,   hdx, hdx), validCells); 
    _estimateVolumeFractions(_volumes.V,      vmath::vec3(hdx, 0,   hdx), validCells); 
    _estimateVolumeFractions(_volumes.W,      vmath::vec3(hdx, hdx, 0  ), validCells); 
    _estimateVolumeFractions(_volumes.edgeU,  vmath::vec3(hdx, 0,   0  ), validCells); 
    _estimateVolumeFractions(_volumes.edgeV,  vmath::vec3(0,   hdx, 0  ), validCells); 
    _estimateVolumeFractions(_volumes.edgeW,  vmath::vec3(0,   0,   hdx), validCells);
}

void ViscositySolver::_estimateVolumeFractions(Array3d<float> &volumes, 
                                               vmath::vec3 centerStart, 
                                               Array3d<bool> &validCells) {

    Array3d<float> nodalPhi(volumes.width + 1, volumes.height + 1, volumes.depth + 1);
    Array3d<bool> isNodalSet(volumes.width + 1, volumes.height + 1, volumes.depth + 1, false);

    volumes.fill(0);
    float hdx = 0.5f * _dx;
    for(int k = 0; k < volumes.depth; k++) {
        for(int j = 0; j < volumes.height; j++) { 
            for(int i = 0; i < volumes.width; i++) {
                if (!validCells(i, j, k)) {
                    continue;
                }

                vmath::vec3 centre = centerStart + Grid3d::GridIndexToCellCenter(i, j, k, _dx);

                if (!isNodalSet(i, j, k)) {
                    float n = _liquidSDF->trilinearInterpolate(centre + vmath::vec3(-hdx, -hdx, -hdx));
                    nodalPhi.set(i, j, k, n);
                    isNodalSet.set(i, j, k, true);
                }
                float phi000 = nodalPhi(i, j, k);

                if (!isNodalSet(i, j, k + 1)) {
                    float n = _liquidSDF->trilinearInterpolate(centre + vmath::vec3(-hdx, -hdx, hdx));
                    nodalPhi.set(i, j, k + 1, n);
                    isNodalSet.set(i, j, k + 1, true);
                }
                float phi001 = nodalPhi(i, j, k + 1);

                if (!isNodalSet(i, j + 1, k)) {
                    float n = _liquidSDF->trilinearInterpolate(centre + vmath::vec3(-hdx, hdx, -hdx));
                    nodalPhi.set(i, j + 1, k, n);
                    isNodalSet.set(i, j + 1, k, true);
                }
                float phi010 = nodalPhi(i, j + 1, k);

                if (!isNodalSet(i, j + 1, k + 1)) {
                    float n = _liquidSDF->trilinearInterpolate(centre + vmath::vec3(-hdx, hdx, hdx));
                    nodalPhi.set(i, j + 1, k + 1, n);
                    isNodalSet.set(i, j + 1, k + 1, true);
                }
                float phi011 = nodalPhi(i, j + 1, k + 1);

                if (!isNodalSet(i + 1, j, k)) {
                    float n = _liquidSDF->trilinearInterpolate(centre + vmath::vec3(hdx, -hdx, -hdx));
                    nodalPhi.set(i + 1, j, k, n);
                    isNodalSet.set(i + 1, j, k, true);
                }
                float phi100 = nodalPhi(i + 1, j, k);

                if (!isNodalSet(i + 1, j, k + 1)) {
                    float n = _liquidSDF->trilinearInterpolate(centre + vmath::vec3(hdx, -hdx, hdx));
                    nodalPhi.set(i + 1, j, k + 1, n);
                    isNodalSet.set(i + 1, j, k + 1, true);
                }
                float phi101 = nodalPhi(i + 1, j, k + 1);

                if (!isNodalSet(i + 1, j + 1, k)) {
                    float n = _liquidSDF->trilinearInterpolate(centre + vmath::vec3(hdx, hdx, -hdx));
                    nodalPhi.set(i + 1, j + 1, k, n);
                    isNodalSet.set(i + 1, j + 1, k, true);
                }
                float phi110 = nodalPhi(i + 1, j + 1, k);

                if (!isNodalSet(i + 1, j + 1, k + 1)) {
                    float n = _liquidSDF->trilinearInterpolate(centre + vmath::vec3(hdx, hdx, hdx));
                    nodalPhi.set(i + 1, j + 1, k + 1, n);
                    isNodalSet.set(i + 1, j + 1, k + 1, true);
                }
                float phi111 = nodalPhi(i + 1, j + 1, k + 1);

                if (phi000 < 0 && phi001 < 0 && phi010 < 0 && phi011 < 0 &&
                    phi100 < 0 && phi101 < 0 && phi110 < 0 && phi111 < 0) {
                    volumes.set(i, j, k, 1.0);
                } else if (phi000 >= 0 && phi001 >= 0 && phi010 >= 0 && phi011 >= 0 &&
                    phi100 >= 0 && phi101 >= 0 && phi110 >= 0 && phi111 >= 0) {
                    volumes.set(i, j, k, 0.0);
                } else {
                    volumes.set(i, j, k, LevelsetUtils::volumeFraction(
                            phi000, phi100, phi010, phi110, phi001, phi101, phi011, phi111
                    ));
                }
            }
        }

    }

}

void ViscositySolver::_destroyVolumeGrid() {
    _volumes.destroy();
}

void ViscositySolver::_computeMatrixIndexTable() {

    int dim = (_isize + 1) * _jsize * _ksize + 
              _isize * (_jsize + 1) * _ksize + 
              _isize * _jsize * (_ksize + 1);
    FaceIndexer fidx(_isize, _jsize, _ksize);

    std::vector<bool> isIndexInMatrix(dim, false);
    for (int k = 1; k < _ksize; k++) {
        for (int j = 1; j < _jsize; j++) {
            for (int i = 1; i < _isize; i++) {
                if (_state.U(i, j, k) != FaceState::fluid) {
                    continue;
                }

                float v = _volumes.U(i, j, k);
                float vRight = _volumes.center(i, j, k);
                float vLeft = _volumes.center(i - 1, j, k);
                float vTop = _volumes.edgeW(i, j + 1, k);
                float vBottom = _volumes.edgeW(i, j, k);
                float vFront = _volumes.edgeV(i, j, k + 1);
                float vBack = _volumes.edgeV(i, j, k);

                if (v > 0.0 || vRight > 0.0 || vLeft > 0.0 || vTop > 0.0 || 
                        vBottom > 0.0 || vFront > 0.0 || vBack > 0.0) {
                    int index = fidx.U(i, j, k);
                    isIndexInMatrix[index] = true;
                }
            }
        }
    }

    for (int k = 1; k < _ksize; k++) {
        for (int j = 1; j < _jsize; j++) {
            for (int i = 1; i < _isize; i++) {
                if (_state.V(i, j, k) != FaceState::fluid) {
                    continue;
                }

                float v = _volumes.V(i, j, k);
                float vRight = _volumes.edgeW(i + 1, j, k);
                float vLeft = _volumes.edgeW(i, j, k);
                float vTop = _volumes.center(i, j, k);
                float vBottom = _volumes.center(i, j - 1, k);
                float vFront = _volumes.edgeU(i, j, k + 1);
                float vBack = _volumes.edgeU(i, j, k);

                if (v > 0.0 || vRight > 0.0 || vLeft > 0.0 || vTop > 0.0 || 
                        vBottom > 0.0 || vFront > 0.0 || vBack > 0.0) {
                    int index = fidx.V(i, j, k);
                    isIndexInMatrix[index] = true;
                }
            }
        }
    }

    for (int k = 1; k < _ksize; k++) {
        for (int j = 1; j < _jsize; j++) {
            for (int i = 1; i < _isize; i++) {
                if (_state.W(i, j, k) != FaceState::fluid) {
                    continue;
                }

                float v = _volumes.W(i, j, k);
                float vRight = _volumes.edgeV(i + 1, j, k);
                float vLeft = _volumes.edgeV(i, j, k);
                float vTop = _volumes.edgeU(i, j + 1, k);
                float vBottom = _volumes.edgeU(i, j, k);
                float vFront = _volumes.center(i, j, k);
                float vBack = _volumes.center(i, j, k - 1);

                if (v > 0.0 || vRight > 0.0 || vLeft > 0.0 || vTop > 0.0 || 
                        vBottom > 0.0 || vFront > 0.0 || vBack > 0.0) {
                    int index = fidx.W(i, j, k);
                    isIndexInMatrix[index] = true;
                }
            }
        }
    }

    std::vector<int> gridToMatrixIndex(dim, -1);
    int matrixindex = 0;
    for (size_t i = 0; i < isIndexInMatrix.size(); i++) {
        if (isIndexInMatrix[i]) {
            gridToMatrixIndex[i] = matrixindex;
            matrixindex++;
        }
    }

    _matrixIndex = MatrixIndexer(_isize, _jsize, _ksize, gridToMatrixIndex);
}

void ViscositySolver::_initializeLinearSystem(SparseMatrixd &matrix, std::vector<double> &rhs) {
    _initializeLinearSystemU(matrix, rhs);
    _initializeLinearSystemV(matrix, rhs);
    _initializeLinearSystemW(matrix, rhs);
}

void ViscositySolver::_initializeLinearSystemU(SparseMatrixd &matrix, std::vector<double> &rhs) {
    MatrixIndexer &mj = _matrixIndex;
    FaceState FLUID = FaceState::fluid;
    FaceState SOLID = FaceState::solid;

    float invdx = 1.0f / _dx;
    float factor = _deltaTime * invdx * invdx;
    for (int k = 1; k < _ksize; k++) {
        for (int j = 1; j < _jsize; j++) {
            for (int i = 1; i < _isize; i++) {
        
                if(_state.U(i, j, k) != FaceState::fluid) {
                    continue;
                }

                int row = _matrixIndex.U(i, j, k);
                if (row == -1) {
                    continue;
                }

                float viscRight = _viscosity->get(i, j, k);
                float viscLeft = _viscosity->get(i - 1, j, k);

                float viscTop    = 0.25f * (_viscosity->get(i - 1, j + 1, k) + 
                                            _viscosity->get(i - 1, j,     k) + 
                                            _viscosity->get(i,     j + 1, k) + 
                                            _viscosity->get(i,     j,     k));
                float viscBottom = 0.25f * (_viscosity->get(i - 1, j,     k) + 
                                            _viscosity->get(i - 1, j - 1, k) + 
                                            _viscosity->get(i,     j,     k) + 
                                            _viscosity->get(i,     j - 1, k));

                float viscFront = 0.25f * (_viscosity->get(i - 1, j, k + 1) + 
                                           _viscosity->get(i - 1, j, k    ) + 
                                           _viscosity->get(i,     j, k + 1) + 
                                           _viscosity->get(i,     j, k    ));
                float viscBack  = 0.25f * (_viscosity->get(i - 1, j, k    ) + 
                                           _viscosity->get(i - 1, j, k - 1) + 
                                           _viscosity->get(i,     j, k    ) + 
                                           _viscosity->get(i,     j, k - 1));

                float volRight = _volumes.center(i, j, k);
                float volLeft = _volumes.center(i-1, j, k);
                float volTop = _volumes.edgeW(i, j + 1, k);
                float volBottom = _volumes.edgeW(i, j, k);
                float volFront = _volumes.edgeV(i, j, k + 1);
                float volBack = _volumes.edgeV(i, j, k);

                float factorRight  = 2 * factor * viscRight * volRight;
                float factorLeft   = 2 * factor * viscLeft * volLeft;
                float factorTop    = factor * viscTop * volTop;
                float factorBottom = factor * viscBottom * volBottom;
                float factorFront  = factor * viscFront * volFront;
                float factorBack   = factor * viscBack * volBack;

                float diag = _volumes.U(i, j, k) + factorRight + factorLeft + factorTop + factorBottom + factorFront + factorBack;
                matrix.set(row, row, diag);
                if (_state.U(i + 1, j,     k    ) == FLUID) { matrix.add(row, mj.U(i + 1, j,     k    ), -factorRight ); }
                if (_state.U(i - 1, j,     k    ) == FLUID) { matrix.add(row, mj.U(i - 1, j,     k    ), -factorLeft  ); }
                if (_state.U(i,     j + 1, k    ) == FLUID) { matrix.add(row, mj.U(i,     j + 1, k    ), -factorTop   ); }
                if (_state.U(i,     j - 1, k    ) == FLUID) { matrix.add(row, mj.U(i,     j - 1, k    ), -factorBottom); }
                if (_state.U(i,     j,     k + 1) == FLUID) { matrix.add(row, mj.U(i,     j,     k + 1), -factorFront ); }
                if (_state.U(i,     j,     k - 1) == FLUID) { matrix.add(row, mj.U(i,     j,     k - 1), -factorBack  ); }

                if (_state.V(i,     j + 1, k    ) == FLUID) { matrix.add(row, mj.V(i,     j + 1, k    ), -factorTop   ); }
                if (_state.V(i - 1, j + 1, k    ) == FLUID) { matrix.add(row, mj.V(i - 1, j + 1, k    ),  factorTop   ); }
                if (_state.V(i,     j,     k    ) == FLUID) { matrix.add(row, mj.V(i,     j,     k    ),  factorBottom); }
                if (_state.V(i - 1, j,     k    ) == FLUID) { matrix.add(row, mj.V(i - 1, j,     k    ), -factorBottom); }
                
                if (_state.W(i,     j,     k + 1) == FLUID) { matrix.add(row, mj.W(i,     j,     k + 1), -factorFront ); }
                if (_state.W(i - 1, j,     k + 1) == FLUID) { matrix.add(row, mj.W(i - 1, j,     k + 1),  factorFront ); }
                if (_state.W(i,     j,     k    ) == FLUID) { matrix.add(row, mj.W(i,     j,     k    ),  factorBack  ); }
                if (_state.W(i - 1, j,     k    ) == FLUID) { matrix.add(row, mj.W(i - 1, j,     k    ), -factorBack  ); }

                float rval = _volumes.U(i, j, k) * _velocityField->U(i, j, k);
                if (_state.U(i + 1, j,     k)     == SOLID) { rval -= -factorRight  * _velocityField->U(i + 1, j,     k    ); }
                if (_state.U(i - 1, j,     k)     == SOLID) { rval -= -factorLeft   * _velocityField->U(i - 1, j,     k    ); }
                if (_state.U(i,     j + 1, k)     == SOLID) { rval -= -factorTop    * _velocityField->U(i,     j + 1, k    ); }
                if (_state.U(i,     j - 1, k)     == SOLID) { rval -= -factorBottom * _velocityField->U(i,     j - 1, k    ); }
                if (_state.U(i,     j,     k + 1) == SOLID) { rval -= -factorFront  * _velocityField->U(i,     j,     k + 1); }
                if (_state.U(i,     j,     k - 1) == SOLID) { rval -= -factorBack   * _velocityField->U(i,     j,     k - 1); }

                if (_state.V(i,     j + 1, k)     == SOLID) { rval -= -factorTop    * _velocityField->V(i,     j + 1, k    ); }
                if (_state.V(i - 1, j + 1, k)     == SOLID) { rval -=  factorTop    * _velocityField->V(i - 1, j + 1, k    ); }
                if (_state.V(i,     j,     k)     == SOLID) { rval -=  factorBottom * _velocityField->V(i,     j,     k    ); }
                if (_state.V(i - 1, j,     k)     == SOLID) { rval -= -factorBottom * _velocityField->V(i - 1, j,     k    ); }

                if (_state.W(i,     j,     k + 1) == SOLID) { rval -= -factorFront  * _velocityField->W(i,     j,     k + 1); } 
                if (_state.W(i - 1, j,     k + 1) == SOLID) { rval -=  factorFront  * _velocityField->W(i - 1, j,     k + 1); } 
                if (_state.W(i,     j,     k)     == SOLID) { rval -=  factorBack   * _velocityField->W(i,     j,     k    ); } 
                if (_state.W(i - 1, j,     k)     == SOLID) { rval -= -factorBack   * _velocityField->W(i - 1, j,     k    ); } 
                rhs[row] = rval;

            }
        }
    }
}

void ViscositySolver::_initializeLinearSystemV(SparseMatrixd &matrix, std::vector<double> &rhs) {
    MatrixIndexer &mj = _matrixIndex;
    FaceState FLUID = FaceState::fluid;
    FaceState SOLID = FaceState::solid;

    float invdx = 1.0f / _dx;
    float factor = _deltaTime * invdx * invdx;
    for (int k = 1; k < _ksize; k++) {
        for (int j = 1; j < _jsize; j++) {
            for (int i = 1; i < _isize; i++) {

                if (_state.V(i, j, k) != FaceState::fluid) {
                    continue;
                }

                int row = _matrixIndex.V(i, j, k);
                if (row == -1) {
                    continue;
                }   

                float viscRight = 0.25f * (_viscosity->get(i,     j - 1, k) + 
                                           _viscosity->get(i + 1, j - 1, k) + 
                                           _viscosity->get(i,     j,     k) + 
                                           _viscosity->get(i + 1, j,     k));
                float viscLeft  = 0.25f * (_viscosity->get(i,     j - 1, k) + 
                                           _viscosity->get(i - 1, j - 1, k) + 
                                           _viscosity->get(i,     j,     k) + 
                                           _viscosity->get(i - 1, j,     k));
                
                float viscTop = _viscosity->get(i, j, k);
                float viscBottom = _viscosity->get(i, j - 1, k);
                
                float viscFront = 0.25f * (_viscosity->get(i, j - 1, k    ) + 
                                           _viscosity->get(i, j - 1, k + 1) + 
                                           _viscosity->get(i, j,     k    ) + 
                                           _viscosity->get(i, j,     k + 1));
                float viscBack  = 0.25f * (_viscosity->get(i, j - 1, k    ) + 
                                           _viscosity->get(i, j - 1, k - 1) + 
                                           _viscosity->get(i, j,     k    ) + 
                                           _viscosity->get(i, j,     k - 1));

                float volRight = _volumes.edgeW(i + 1, j, k);
                float volLeft = _volumes.edgeW(i, j, k);
                float volTop = _volumes.center(i, j, k);
                float volBottom = _volumes.center(i, j - 1, k);
                float volFront = _volumes.edgeU(i, j, k + 1);
                float volBack = _volumes.edgeU(i, j, k);

                float factorRight  = factor * viscRight * volRight;
                float factorLeft   = factor * viscLeft * volLeft;
                float factorTop    = 2 * factor * viscTop * volTop;
                float factorBottom = 2 * factor * viscBottom * volBottom;
                float factorFront  = factor * viscFront * volFront;
                float factorBack   = factor * viscBack*volBack;

                float diag = _volumes.V(i, j, k) + factorRight + factorLeft + factorTop + factorBottom + factorFront + factorBack;
                matrix.set(row, row, diag);
                if (_state.V(i + 1, j,     k    ) == FLUID) { matrix.add(row, mj.V(i + 1, j,     k    ), -factorRight ); }
                if (_state.V(i - 1, j,     k    ) == FLUID) { matrix.add(row, mj.V(i - 1, j,     k    ), -factorLeft  ); }
                if (_state.V(i,     j + 1, k    ) == FLUID) { matrix.add(row, mj.V(i,     j + 1, k    ), -factorTop   ); }
                if (_state.V(i,     j - 1, k    ) == FLUID) { matrix.add(row, mj.V(i,     j - 1, k    ), -factorBottom); }
                if (_state.V(i,     j,     k + 1) == FLUID) { matrix.add(row, mj.V(i,     j,     k + 1), -factorFront ); }
                if (_state.V(i,     j,     k - 1) == FLUID) { matrix.add(row, mj.V(i,     j,     k - 1), -factorBack  ); }

                if (_state.U(i + 1, j,     k    ) == FLUID) { matrix.add(row, mj.U(i + 1, j,     k    ), -factorRight ); }
                if (_state.U(i + 1, j - 1, k    ) == FLUID) { matrix.add(row, mj.U(i + 1, j - 1, k    ),  factorRight ); }
                if (_state.U(i,     j,     k    ) == FLUID) { matrix.add(row, mj.U(i,     j,     k    ),  factorLeft  ); }
                if (_state.U(i,     j - 1, k    ) == FLUID) { matrix.add(row, mj.U(i,     j - 1, k    ), -factorLeft  ); }
            
                if (_state.W(i,     j,     k + 1) == FLUID) { matrix.add(row, mj.W(i,     j,     k + 1), -factorFront ); }
                if (_state.W(i,     j - 1, k + 1) == FLUID) { matrix.add(row, mj.W(i,     j - 1, k + 1),  factorFront ); }
                if (_state.W(i,     j,     k    ) == FLUID) { matrix.add(row, mj.W(i,     j,     k    ),  factorBack  ); }
                if (_state.W(i,     j - 1, k    ) == FLUID) { matrix.add(row, mj.W(i,     j - 1, k    ), -factorBack  ); }

                float rval = _volumes.V(i, j, k) * _velocityField->V(i, j, k);
                if (_state.V(i + 1, j,     k)     == SOLID) { rval -= -factorRight  * _velocityField->V(i + 1, j,     k    ); }
                if (_state.V(i - 1, j,     k)     == SOLID) { rval -= -factorLeft   * _velocityField->V(i - 1, j,     k    ); }
                if (_state.V(i,     j + 1, k)     == SOLID) { rval -= -factorTop    * _velocityField->V(i,     j + 1, k    ); }
                if (_state.V(i ,    j - 1, k)     == SOLID) { rval -= -factorBottom * _velocityField->V(i,     j - 1, k    ); }
                if (_state.V(i,     j,     k + 1) == SOLID) { rval -= -factorFront  * _velocityField->V(i,     j,     k + 1); }
                if (_state.V(i,     j,     k - 1) == SOLID) { rval -= -factorBack   * _velocityField->V(i,     j,     k - 1); }

                if (_state.U(i + 1, j,     k)     == SOLID) { rval -= -factorRight  * _velocityField->U(i + 1, j,     k    ); }
                if (_state.U(i + 1, j - 1, k)     == SOLID) { rval -=  factorRight  * _velocityField->U(i + 1, j - 1, k    ); }
                if (_state.U(i,     j,     k)     == SOLID) { rval -=  factorLeft   * _velocityField->U(i,     j,     k    ); }
                if (_state.U(i,     j - 1, k)     == SOLID) { rval -= -factorLeft   * _velocityField->U(i,     j - 1, k    ); }

                if (_state.W(i,     j,     k + 1) == SOLID) { rval -= -factorFront  * _velocityField->W(i,     j,     k + 1); }
                if (_state.W(i,     j - 1, k + 1) == SOLID) { rval -=  factorFront  * _velocityField->W(i,     j - 1, k + 1); }
                if (_state.W(i,     j,     k)     == SOLID) { rval -=  factorBack   * _velocityField->W(i,     j,     k    ); }
                if (_state.W(i,     j - 1, k)     == SOLID) { rval -= -factorBack   * _velocityField->W(i,     j - 1, k    ); }
                rhs[row] = rval;

            }
        }
    }
}

void ViscositySolver::_initializeLinearSystemW(SparseMatrixd &matrix, std::vector<double> &rhs) {
    MatrixIndexer &mj = _matrixIndex;
    FaceState FLUID = FaceState::fluid;
    FaceState SOLID = FaceState::solid;

    float invdx = 1.0f / _dx;
    float factor = _deltaTime * invdx * invdx;
    for (int k = 1; k < _ksize; k++) {
        for (int j = 1; j < _jsize; j++) {
            for (int i = 1; i < _isize; i++) {

                if (_state.W(i, j, k) != FaceState::fluid) {
                    continue;
                }

                int row = _matrixIndex.W(i, j, k);
                if (row == -1) {
                    continue;
                }  

                float viscRight = 0.25f * (_viscosity->get(i,     j, k    ) + 
                                           _viscosity->get(i,     j, k - 1) + 
                                           _viscosity->get(i + 1, j, k    ) + 
                                           _viscosity->get(i + 1, j, k - 1));
                float viscLeft  = 0.25f * (_viscosity->get(i,     j, k    ) + 
                                           _viscosity->get(i,     j, k - 1) + 
                                           _viscosity->get(i - 1, j, k    ) + 
                                           _viscosity->get(i - 1, j, k - 1));

                float viscTop    = 0.25f * (_viscosity->get(i, j,     k    ) + 
                                            _viscosity->get(i, j,     k - 1) + 
                                            _viscosity->get(i, j + 1, k    ) + 
                                            _viscosity->get(i, j + 1, k - 1));
                float viscBottom = 0.25f * (_viscosity->get(i, j,     k    ) + 
                                            _viscosity->get(i, j,     k - 1) + 
                                            _viscosity->get(i, j - 1, k    ) + 
                                            _viscosity->get(i, j - 1, k - 1));

                float viscFront = _viscosity->get(i, j, k);   
                float viscBack = _viscosity->get(i, j, k - 1); 

                float volRight = _volumes.edgeV(i + 1, j, k);
                float volLeft = _volumes.edgeV(i, j, k);
                float volTop = _volumes.edgeU(i, j + 1, k);
                float volBottom = _volumes.edgeU(i, j, k);
                float volFront = _volumes.center(i, j, k);
                float volBack = _volumes.center(i, j, k - 1);

                float factorRight  = factor * viscRight * volRight;
                float factorLeft   = factor * viscLeft * volLeft;
                float factorTop    = factor * viscTop * volTop;
                float factorBottom = factor * viscBottom * volBottom;
                float factorFront  = 2 * factor * viscFront * volFront;
                float factorBack   = 2 * factor * viscBack*volBack;

                float diag = _volumes.W(i, j, k) + factorRight + factorLeft + factorTop + factorBottom + factorFront + factorBack;
                matrix.set(row, row, diag);
                if (_state.W(i + 1, j,     k    ) == FLUID) { matrix.add(row, mj.W(i + 1, j,     k    ), -factorRight ); }
                if (_state.W(i - 1, j,     k    ) == FLUID) { matrix.add(row, mj.W(i - 1, j,     k    ), -factorLeft  ); }
                if (_state.W(i,     j + 1, k    ) == FLUID) { matrix.add(row, mj.W(i,     j + 1, k    ), -factorTop   ); }
                if (_state.W(i,     j - 1, k    ) == FLUID) { matrix.add(row, mj.W(i,     j - 1, k    ), -factorBottom); }
                if (_state.W(i,     j,     k + 1) == FLUID) { matrix.add(row, mj.W(i,     j,     k + 1), -factorFront ); }
                if (_state.W(i,     j,     k - 1) == FLUID) { matrix.add(row, mj.W(i,     j,     k - 1), -factorBack  ); }

                if (_state.U(i + 1, j,     k    ) == FLUID) { matrix.add(row, mj.U(i + 1, j,     k    ), -factorRight ); } 
                if (_state.U(i + 1, j,     k - 1) == FLUID) { matrix.add(row, mj.U(i + 1, j,     k - 1),  factorRight ); }
                if (_state.U(i,     j,     k    ) == FLUID) { matrix.add(row, mj.U(i,     j,     k    ),  factorLeft  ); }
                if (_state.U(i,     j,     k - 1) == FLUID) { matrix.add(row, mj.U(i,     j,     k - 1), -factorLeft  ); }
                
                if (_state.V(i,     j + 1, k    ) == FLUID) { matrix.add(row, mj.V(i,     j + 1, k    ), -factorTop   ); }
                if (_state.V(i,     j + 1, k - 1) == FLUID) { matrix.add(row, mj.V(i,     j + 1, k - 1),  factorTop   ); }
                if (_state.V(i,     j,     k    ) == FLUID) { matrix.add(row, mj.V(i,     j,     k    ),  factorBottom); }
                if (_state.V(i,     j,     k - 1) == FLUID) { matrix.add(row, mj.V(i,     j,     k - 1), -factorBottom); }

                float rval = _volumes.W(i, j, k) * _velocityField->W(i, j, k);
                if (_state.W(i + 1, j,     k)     == SOLID) { rval -= -factorRight  * _velocityField->W(i + 1, j,     k    ); }
                if (_state.W(i - 1, j,     k)     == SOLID) { rval -= -factorLeft   * _velocityField->W(i - 1, j,     k    ); }
                if (_state.W(i,     j + 1, k)     == SOLID) { rval -= -factorTop    * _velocityField->W(i,     j + 1, k    ); }
                if (_state.W(i,     j - 1, k)     == SOLID) { rval -= -factorBottom * _velocityField->W(i,     j - 1, k    ); }
                if (_state.W(i,     j,     k + 1) == SOLID) { rval -= -factorFront  * _velocityField->W(i,     j,     k + 1); }
                if (_state.W(i,     j,     k - 1) == SOLID) { rval -= -factorBack   * _velocityField->W(i,     j,     k - 1); }
                if (_state.U(i + 1, j,     k)     == SOLID) { rval -= -factorRight  * _velocityField->U(i + 1, j,     k    ); }
                if (_state.U(i + 1, j,     k - 1) == SOLID) { rval -=  factorRight  * _velocityField->U(i + 1, j,     k - 1); }
                if (_state.U(i,     j,     k)     == SOLID) { rval -=  factorLeft   * _velocityField->U(i,     j,     k    ); }
                if (_state.U(i,     j,     k - 1) == SOLID) { rval -= -factorLeft   * _velocityField->U(i,     j,     k - 1); }
                if (_state.V(i,     j + 1, k)     == SOLID) { rval -= -factorTop    * _velocityField->V(i,     j + 1, k    ); }
                if (_state.V(i,     j + 1, k - 1) == SOLID) { rval -=  factorTop    * _velocityField->V(i,     j + 1, k - 1); }
                if (_state.V(i,     j,     k)     == SOLID) { rval -=  factorBottom * _velocityField->V(i,     j,     k    ); }
                if (_state.V(i,     j,     k - 1) == SOLID) { rval -= -factorBottom * _velocityField->V(i,     j,     k - 1); }
                rhs[row] = rval;

            }
        }
    }
}

bool ViscositySolver::_solveLinearSystem(SparseMatrixd &matrix, std::vector<double> &rhs, 
                                         std::vector<double> &soln) {

    PCGSolver<double> solver;
    solver.setSolverParameters(_solverTolerance, _maxSolverIterations);

    double estimatedError;
    int numIterations;
    bool success = solver.solve(matrix, rhs, soln, estimatedError, numIterations);

    if (success) {
        std::cout << "\n\tViscosity Solver Iterations: " << numIterations <<
                     "\n\tEstimated Error: " << estimatedError << "\n\n";
        return true;
    } else if (numIterations == _maxSolverIterations && estimatedError < _acceptableTolerace) {
        std::cout << "\n\tViscosity Solver Iterations: " << numIterations <<
                     "\n\tEstimated Error: " << estimatedError << "\n\n";
        return true;
    } else {
        std::cout << "\n\t***Viscosity Solver FAILED" <<
              "\n\tViscosity Solver Iterations: " << numIterations <<
              "\n\tEstimated Error: " << estimatedError << "\n\n";
        return false;
    }
}

void ViscositySolver::_applySolutionToVelocityField(std::vector<double> &soln) {
    _velocityField->clear();
    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize + 1; i++) {
                int matidx = _matrixIndex.U(i, j, k);
                if (matidx != -1) {
                    _velocityField->setU(i, j, k, soln[matidx]);
                }
            }
        }
    }

    for(int k = 0; k < _ksize; k++) {
        for(int j = 0; j < _jsize + 1; j++) {
            for(int i = 0; i < _isize; i++) {
                int matidx = _matrixIndex.V(i, j, k);
                if (matidx != -1) {
                    _velocityField->setV(i, j, k, soln[matidx]);
                }
            }
        }
    }

    for(int k = 0; k < _ksize + 1; k++) {
        for(int j = 0; j < _jsize; j++) {
            for(int i = 0; i < _isize; i++) {
                int matidx = _matrixIndex.W(i, j, k);
                if (matidx != -1) {
                    _velocityField->setW(i, j, k, soln[matidx]);
                }
            }
        }
    }

}
