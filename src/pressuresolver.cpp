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
#include "pressuresolver.h"

/********************************************************************************
    VectorXd
********************************************************************************/

VectorXd::VectorXd() {
}

VectorXd::VectorXd(int size) : _vector(size, 0.0) {
}

VectorXd::VectorXd(int size, double fill) : _vector(size, fill) {
}

VectorXd::VectorXd(VectorXd &vector) {
    _vector.reserve(vector.size());
    for (unsigned int i = 0; i < vector.size(); i++) {
        _vector.push_back(vector[i]);
    }
}

VectorXd::~VectorXd() {
}

const double VectorXd::operator[](int i) const {
    FLUIDSIM_ASSERT(i >= 0 && i < (int)_vector.size());
    return _vector[i];
}

double& VectorXd::operator[](int i) {
    FLUIDSIM_ASSERT(i >= 0 && i < (int)_vector.size());
    return _vector[i];
}

void VectorXd::fill(double fill) {
    for (unsigned int i = 0; i < _vector.size(); i++) {
        _vector[i] = fill;
    }
}

double VectorXd::dot(VectorXd &vector) {
    FLUIDSIM_ASSERT(_vector.size() == vector._vector.size());

    double sum = 0.0;
    for (unsigned int i = 0; i < _vector.size(); i++) {
        sum += _vector[i] * vector._vector[i];
    }

    return sum;
}

double VectorXd::absMaxCoeff() {
    double max = -std::numeric_limits<double>::infinity();
    for (unsigned int i = 0; i < _vector.size(); i++) {
        if (fabs(_vector[i]) > max) {
            max = fabs(_vector[i]);
        }
    }

    return max;
}

/********************************************************************************
    MatrixCoefficients
********************************************************************************/

MatrixCoefficients::MatrixCoefficients() {
}

MatrixCoefficients::MatrixCoefficients(int size) : cells(size, MatrixCell()) {
}

MatrixCoefficients::~MatrixCoefficients() {
}

const MatrixCell MatrixCoefficients::operator[](int i) const {
    FLUIDSIM_ASSERT(i >= 0 && i < (int)cells.size());
    return cells[i];
}

MatrixCell& MatrixCoefficients::operator[](int i) {
    FLUIDSIM_ASSERT(i >= 0 && i < (int)cells.size());
    return cells[i];
}

/********************************************************************************
    PressureSolver
********************************************************************************/

PressureSolver::PressureSolver() {
}

PressureSolver::~PressureSolver() {
}

Array3d<float> PressureSolver::solve(PressureSolverParameters params) {
    _initialize(params);

    _initializeGridIndexKeyMap();

    VectorXd b(_matSize);
    _calculateNegativeDivergenceVector(b);
    if (b.absMaxCoeff() < _pressureSolveTolerance) {
        return Array3d<float>(_isize, _jsize, _ksize, 0.0);
    }

    MatrixCoefficients A(_matSize);
    _calculateMatrixCoefficients(A);

    VectorXd precon(_matSize);
    _calculatePreconditionerVector(A, precon);


    VectorXd pressure(_matSize);
    _solvePressureSystem(A, b, precon, pressure);

    Array3d<float> pressureGrid(_isize, _jsize, _ksize, 0.0);
    for (size_t i = 0; i < _pressureCells->size(); i++) {
        GridIndex g = _pressureCells->get(i);
        pressureGrid.set(g, pressure[i]);
    }

    return pressureGrid;
}

void PressureSolver::_initialize(PressureSolverParameters params) {
    params.velocityField->getGridDimensions(&_isize, &_jsize, &_ksize);
    _dx = params.cellwidth;
    _density = params.density;
    _deltaTime = params.deltaTime;

    _pressureCells = params.pressureCells;
    _vField = params.velocityField;
    _liquidSDF = params.liquidSDF;
    _weightGrid = params.weightGrid;
    //_logfile = params.logfile;
    _matSize = (int)_pressureCells->size();
}

void PressureSolver::_initializeGridIndexKeyMap() {
    _keymap = GridIndexKeyMap(_isize, _jsize, _ksize);
    for (unsigned int idx = 0; idx < _pressureCells->size(); idx++) {
        _keymap.insert(_pressureCells->at(idx), idx);
    }
}

void PressureSolver::_calculateNegativeDivergenceVector(VectorXd &b) {

    for (size_t idx = 0; idx < _pressureCells->size(); idx++) {
        int i = _pressureCells->at(idx).i;
        int j = _pressureCells->at(idx).j;
        int k = _pressureCells->at(idx).k;

        double divergence = 0.0;
        divergence -= _weightGrid->U(i + 1, j, k) * _vField->U(i + 1, j, k);
        divergence += _weightGrid->U(i, j, k)     * _vField->U(i, j, k);
        divergence -= _weightGrid->V(i, j + 1, k) * _vField->V(i, j + 1, k);
        divergence += _weightGrid->V(i, j, k)     * _vField->V(i, j, k);
        divergence -= _weightGrid->W(i, j, k + 1) * _vField->W(i, j, k + 1);
        divergence += _weightGrid->W(i, j, k)     * _vField->W(i, j, k);
        divergence /= _dx;

        b[_GridToVectorIndex(i, j, k)] = divergence;
    }
}

void PressureSolver::_calculateMatrixCoefficients(MatrixCoefficients &A) {

    double scale = _deltaTime / (_dx * _dx);
    for (size_t idx = 0; idx < _pressureCells->size(); idx++) {
        int i = _pressureCells->at(idx).i;
        int j = _pressureCells->at(idx).j;
        int k = _pressureCells->at(idx).k;
        int index = _GridToVectorIndex(i, j, k);

        //right neighbour
        float term = _weightGrid->U(i + 1, j, k) * scale;
        float phiRight = _liquidSDF->get(i + 1, j, k);
        if(phiRight < 0) {
            A[index].diag += term;
            A[index].plusi -= term;
        } else {
            float theta = fmax(_liquidSDF->getFaceWeightU(i + 1, j, k), _minfrac);
            A[index].diag += term / theta;
        }

        //left neighbour
        term = _weightGrid->U(i, j, k) * scale;
        float phiLeft = _liquidSDF->get(i - 1, j, k);
        if(phiLeft < 0) {
            A[index].diag += term;
        } else {
            float theta = fmax(_liquidSDF->getFaceWeightU(i, j, k), _minfrac);
            A[index].diag += term / theta;
        }

        //top neighbour
        term = _weightGrid->V(i, j + 1, k) * scale;
        float phiTop = _liquidSDF->get(i, j + 1, k);
        if(phiTop < 0) {
            A[index].diag += term;
            A[index].plusj -= term;
        } else {
            float theta = fmax(_liquidSDF->getFaceWeightV(i, j + 1, k), _minfrac);
            A[index].diag += term / theta;
        }

        //bottom neighbour
        term = _weightGrid->V(i, j, k) * scale;
        float phiBot = _liquidSDF->get(i, j - 1, k);
        if(phiBot < 0) {
            A[index].diag += term;
        } else {
            float theta = fmax(_liquidSDF->getFaceWeightV(i, j, k), _minfrac);
            A[index].diag += term / theta;
        }

        //far neighbour
        term = _weightGrid->W(i, j, k + 1) * scale;
        float phiFar = _liquidSDF->get(i, j, k + 1);
        if(phiFar < 0) {
            A[index].diag += term;
            A[index].plusk -= term;
        } else {
            float theta = fmax(_liquidSDF->getFaceWeightW(i, j, k + 1), _minfrac);
            A[index].diag += term / theta;
        }

        //near neighbour
        term =_weightGrid->W(i, j, k) * scale;
        float phiNear = _liquidSDF->get(i, j, k - 1);
        if(phiNear < 0) {
            A[index].diag += term;
        } else {
            float theta = fmax(_liquidSDF->getFaceWeightW(i, j, k), _minfrac);
            A[index].diag += term / theta;
        }
    }
}

void PressureSolver::_calculatePreconditionerVector(MatrixCoefficients &A, VectorXd &precon) {
    FLUIDSIM_ASSERT(A.size() == precon.size());

    double tau = 0.97;      // Tuning constant
    double sigma = 0.25;    // safety constant
    for (unsigned int idx = 0; idx < _pressureCells->size(); idx++) {
        int i = _pressureCells->at(idx).i;
        int j = _pressureCells->at(idx).j;
        int k = _pressureCells->at(idx).k;
        int vidx = _GridToVectorIndex(i, j, k);

        int vidx_im1 = _keymap.find(i - 1, j, k);
        int vidx_jm1 = _keymap.find(i, j - 1, k);
        int vidx_km1 = _keymap.find(i, j, k - 1);

        double diag = (double)A[vidx].diag;

        double plusi_im1 = vidx_im1 != -1 ? (double)A[vidx_im1].plusi : 0.0;
        double plusi_jm1 = vidx_jm1 != -1 ? (double)A[vidx_jm1].plusi : 0.0;
        double plusi_km1 = vidx_km1 != -1 ? (double)A[vidx_km1].plusi : 0.0;

        double plusj_im1 = vidx_im1 != -1 ? (double)A[vidx_im1].plusj : 0.0;
        double plusj_jm1 = vidx_jm1 != -1 ? (double)A[vidx_jm1].plusj : 0.0;
        double plusj_km1 = vidx_km1 != -1 ? (double)A[vidx_km1].plusj : 0.0;

        double plusk_im1 = vidx_im1 != -1 ? (double)A[vidx_im1].plusk : 0.0;
        double plusk_jm1 = vidx_jm1 != -1 ? (double)A[vidx_jm1].plusk : 0.0;
        double plusk_km1 = vidx_km1 != -1 ? (double)A[vidx_km1].plusk : 0.0;

        double precon_im1 = vidx_im1 != -1 ? precon[vidx_im1] : 0.0;
        double precon_jm1 = vidx_jm1 != -1 ? precon[vidx_jm1] : 0.0;
        double precon_km1 = vidx_km1 != -1 ? precon[vidx_km1] : 0.0;

        double v1 = plusi_im1 * precon_im1;
        double v2 = plusj_jm1 * precon_jm1;
        double v3 = plusk_km1 * precon_km1;
        double v4 = precon_im1 * precon_im1;
        double v5 = precon_jm1 * precon_jm1;
        double v6 = precon_km1 * precon_km1;
        
        double e = diag - v1*v1 - v2*v2 - v3*v3 - 
            tau*(plusi_im1*(plusj_im1 + plusk_im1)*v4 +
                 plusj_jm1*(plusi_jm1 + plusk_jm1)*v5 +
                 plusk_km1*(plusi_km1 + plusj_km1)*v6);

        if (e < sigma*diag) {
            e = diag;
        }

        if (fabs(e) > 10e-9) {
            precon[vidx] = 1.0 / sqrt(e);
        }
    }
}

void PressureSolver::_applyPreconditioner(MatrixCoefficients &A, 
                                          VectorXd &precon,
                                          VectorXd &residual,
                                          VectorXd &vect) {
    // Solve A*q = residual
    VectorXd q(_matSize);
    for (unsigned int idx = 0; idx < _pressureCells->size(); idx++) {
        int i = _pressureCells->at(idx).i;
        int j = _pressureCells->at(idx).j;
        int k = _pressureCells->at(idx).k;
        int vidx = _GridToVectorIndex(i, j, k);

        int vidx_im1 = _keymap.find(i - 1, j, k);
        int vidx_jm1 = _keymap.find(i, j - 1, k);
        int vidx_km1 = _keymap.find(i, j, k - 1);

        double plusi_im1 = 0.0;
        double precon_im1 = 0.0;
        double q_im1 = 0.0;
        if (vidx_im1 != -1) {
            plusi_im1  = (double)A[vidx_im1].plusi;
            precon_im1 = precon[vidx_im1];
            q_im1      = q[vidx_im1];
        }

        double plusj_jm1 = 0.0;
        double precon_jm1 = 0.0;
        double q_jm1 = 0.0;
        if (vidx_jm1 != -1) {
            plusj_jm1  = (double)A[vidx_jm1].plusj;
            precon_jm1 = precon[vidx_jm1];
            q_jm1      = q[vidx_jm1];
        }

        double plusk_km1 = 0.0;
        double precon_km1 = 0.0;
        double q_km1 = 0.0;
        if (vidx_km1 != -1) {
            plusk_km1  = (double)A[vidx_km1].plusk;
            precon_km1 = precon[vidx_km1];
            q_km1      = q[vidx_km1];
        }

        double t = residual[vidx] - plusi_im1 * precon_im1 * q_im1 -
                                    plusj_jm1 * precon_jm1 * q_jm1 -
                                    plusk_km1 * precon_km1 * q_km1;

        t = t*precon[vidx];
        q[vidx] = t;
    }

    // Solve transpose(A)*z = q
    for (int idx = (int)_pressureCells->size() - 1; idx >= 0; idx--) {
        int i = _pressureCells->at(idx).i;
        int j = _pressureCells->at(idx).j;
        int k = _pressureCells->at(idx).k;
        int vidx = _GridToVectorIndex(i, j, k);

        int vidx_ip1 = _keymap.find(i + 1, j, k);
        int vidx_jp1 = _keymap.find(i, j + 1, k);
        int vidx_kp1 = _keymap.find(i, j, k + 1);

        double vect_ip1 = vidx_ip1 != -1 ? vect[vidx_ip1] : 0.0;
        double vect_jp1 = vidx_jp1 != -1 ? vect[vidx_jp1] : 0.0;
        double vect_kp1 = vidx_kp1 != -1 ? vect[vidx_kp1] : 0.0;

        double plusi = (double)A[vidx].plusi;
        double plusj = (double)A[vidx].plusj;
        double plusk = (double)A[vidx].plusk;

        double preconval = precon[vidx];
        double t = q[vidx] - plusi * preconval * vect_ip1 -
                             plusj * preconval * vect_jp1 -
                             plusk * preconval * vect_kp1;

        t = t*preconval;
        vect[vidx] = t;
    }
}

void PressureSolver::_applyMatrix(MatrixCoefficients &A, VectorXd &x, VectorXd &result) {
    FLUIDSIM_ASSERT(A.size() == x.size() && x.size() == result.size());

    for (unsigned int idx = 0; idx < _pressureCells->size(); idx++) {
        int i = _pressureCells->at(idx).i;
        int j = _pressureCells->at(idx).j;
        int k = _pressureCells->at(idx).k;
        int ridx = _GridToVectorIndex(i, j, k);

        // val = dot product of column vector x and idxth row of matrix A
        double val = 0.0;
        int vidx = _GridToVectorIndex(i - 1, j, k);
        if (vidx != -1) { val += x._vector[vidx] * A[vidx].plusi; }

        vidx = _GridToVectorIndex(i + 1, j, k);
        if (vidx != -1) { val += x._vector[vidx] * A[ridx].plusi; }

        vidx = _GridToVectorIndex(i, j - 1, k);
        if (vidx != -1) { val += x._vector[vidx] * A[vidx].plusj; }

        vidx = _GridToVectorIndex(i, j + 1, k);
        if (vidx != -1) { val += x._vector[vidx] * A[ridx].plusj; }

        vidx = _GridToVectorIndex(i, j, k - 1);
        if (vidx != -1) { val += x._vector[vidx] * A[vidx].plusk; }

        vidx = _GridToVectorIndex(i, j, k + 1);
        if (vidx != -1) { val += x._vector[vidx] * A[ridx].plusk; }

        val += x._vector[ridx] * A.cells[ridx].diag;

        result._vector[ridx] = val;
    }
}

// v1 += v2*scale
void PressureSolver::_addScaledVector(VectorXd &v1, VectorXd &v2, double scale) {
    FLUIDSIM_ASSERT(v1.size() == v2.size());
    for (unsigned int idx = 0; idx < v1.size(); idx++) {
        v1._vector[idx] += v2._vector[idx]*scale;
    }
}

// result = v1*s1 + v2*s2
void PressureSolver::_addScaledVectors(VectorXd &v1, double s1, 
                                       VectorXd &v2, double s2,
                                       VectorXd &result) {
    FLUIDSIM_ASSERT(v1.size() == v2.size() && v2.size() == result.size());
    for (unsigned int idx = 0; idx < v1.size(); idx++) {
        result._vector[idx] = v1._vector[idx]*s1 + v2._vector[idx]*s2;
    }
}

// Solve (A*pressure = b) with Modified Incomplete Cholesky 
// Conjugate Gradient method (MICCG(0))
void PressureSolver::_solvePressureSystem(MatrixCoefficients &A, 
                                          VectorXd &b, 
                                          VectorXd &precon,
                                          VectorXd &pressure) {

    double tol = _pressureSolveTolerance;
    if (b.absMaxCoeff() < tol) {
        return;
    }

    VectorXd residual(b);
    VectorXd auxillary(_matSize);
    _applyPreconditioner(A, precon, residual, auxillary);

    VectorXd search(auxillary);

    double alpha = 0.0;
    double beta = 0.0;
    double sigma = auxillary.dot(residual);
    double sigmaNew = 0.0;
    int iterationNumber = 0;

    while (iterationNumber < _maxCGIterations) {
        _applyMatrix(A, search, auxillary);
        alpha = sigma / auxillary.dot(search);
        _addScaledVector(pressure, search, alpha);
        _addScaledVector(residual, auxillary, -alpha);

        if (residual.absMaxCoeff() < tol) {
            //_logfile->log("CG Iterations: ", iterationNumber, 1);
            std::cout << "CG Iterations: " << iterationNumber << std::endl;
            return;
        }

        _applyPreconditioner(A, precon, residual, auxillary);
        sigmaNew = auxillary.dot(residual);
        beta = sigmaNew / sigma;
        _addScaledVectors(auxillary, 1.0, search, beta, search);
        sigma = sigmaNew;

        iterationNumber++;

        if (iterationNumber % 10 == 0) {
            std::ostringstream ss;
            ss << "\tIteration #: " << iterationNumber <<
                  "\tEstimated Error: " << residual.absMaxCoeff() << std::endl;
            //_logfile->print(ss.str());
            std::cout << ss.str();
        }
    }

    //_logfile->log("Iterations limit reached.\t Estimated error : ",
    //              residual.absMaxCoeff(), 1);
}   //