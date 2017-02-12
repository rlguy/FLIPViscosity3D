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

#include "meshlevelset.h"

MeshLevelSet::MeshLevelSet() {
}

MeshLevelSet::MeshLevelSet(int isize, int jsize, int ksize, double dx) :
                _isize(isize), _jsize(jsize), _ksize(ksize), _dx(dx),
                            _closestTriangles(isize + 1, jsize + 1, ksize + 1, -1) {

    _phi = Array3d<float>(isize + 1, jsize + 1, ksize + 1, 0.0);
}

MeshLevelSet::~MeshLevelSet() {
}

float MeshLevelSet::operator()(int i, int j, int k) {
    return get(i, j, k);
}

float MeshLevelSet::operator()(GridIndex g) {
    return get(g);
}

float MeshLevelSet::get(int i, int j, int k) {
    FLUIDSIM_ASSERT(_phi.isIndexInRange(i, j, k));
    return _phi(i, j, k);
}

float MeshLevelSet::get(GridIndex g) {
    FLUIDSIM_ASSERT(_phi.isIndexInRange(g));
    return _phi(g);
}

int MeshLevelSet::getClosestTriangleIndex(int i, int j, int k) {
    FLUIDSIM_ASSERT(_closestTriangles.isIndexInRange(i, j, k));
    return _closestTriangles(i, j, k);
}

int MeshLevelSet::getClosestTriangleIndex(GridIndex g) {
    FLUIDSIM_ASSERT(_closestTriangles.isIndexInRange(g));
    return _closestTriangles(g);
}

float MeshLevelSet::getDistanceAtCellCenter(int i, int j, int k) {
    FLUIDSIM_ASSERT(Grid3d::isGridIndexInRange(i, j, k, _isize, _jsize, _ksize));
    return 0.125f * (_phi(i, j, k) + 
                     _phi(i + 1, j, k) + 
                     _phi(i, j + 1, k) + 
                     _phi(i + 1, j + 1, k) +
                     _phi(i, j, k + 1) + 
                     _phi(i + 1, j, k + 1) + 
                     _phi(i, j + 1, k + 1) + 
                     _phi(i + 1, j + 1, k + 1));
}

float MeshLevelSet::getDistanceAtCellCenter(GridIndex g) {
    return getDistanceAtCellCenter(g.i, g.j, g.k);
}

float MeshLevelSet::trilinearInterpolate(vmath::vec3 pos) {
    return (float)Interpolation::trilinearInterpolate(pos, _dx, _phi);
}

vmath::vec3 MeshLevelSet::trilinearInterpolateGradient(vmath::vec3 pos) {
    vmath::vec3 grad;
    Interpolation::trilinearInterpolateGradient(pos, _dx, _phi, &grad);
    return grad;
}

float MeshLevelSet::getFaceWeightU(int i, int j, int k) {
    FLUIDSIM_ASSERT(Grid3d::isGridIndexInRange(i, j, k, _isize + 1, _jsize, _ksize));
    return LevelsetUtils::fractionInside(_phi(i, j, k), 
                                         _phi(i, j + 1, k),
                                         _phi(i, j, k + 1), 
                                         _phi(i, j + 1, k + 1));
}

float MeshLevelSet::getFaceWeightU(GridIndex g) {
    return getFaceWeightU(g.i, g.j, g.k);
}

float MeshLevelSet::getFaceWeightV(int i, int j, int k) {
    FLUIDSIM_ASSERT(Grid3d::isGridIndexInRange(i, j, k, _isize, _jsize + 1, _ksize));
    return LevelsetUtils::fractionInside(_phi(i, j, k),
                                         _phi(i, j, k + 1),
                                         _phi(i + 1, j, k),
                                         _phi(i + 1, j, k + 1));
}

float MeshLevelSet::getFaceWeightV(GridIndex g) {
    return getFaceWeightV(g.i, g.j, g.k);
}

float MeshLevelSet::getFaceWeightW(int i, int j, int k) {
    FLUIDSIM_ASSERT(Grid3d::isGridIndexInRange(i, j, k, _isize, _jsize, _ksize + 1));
    return LevelsetUtils::fractionInside(_phi(i, j, k),
                                         _phi(i, j + 1, k),
                                         _phi(i + 1, j, k),
                                         _phi(i + 1, j + 1, k));
}

float MeshLevelSet::getFaceWeightW(GridIndex g) {
    return getFaceWeightW(g.i, g.j, g.k);
}

void MeshLevelSet::getGridDimensions(int *i, int *j, int *k) {
    *i = _isize;
    *j = _jsize;
    *k = _ksize;
}

TriangleMesh *MeshLevelSet::getTriangleMesh() {
    return &_mesh;
}

void MeshLevelSet::calculateSignedDistanceField(TriangleMesh &m, int bandwidth) {
    _mesh = m;
    Array3d<int> intersectionCounts(_phi.width, _phi.height, _phi.depth);

    // we begin by initializing distances near the mesh, and figuring out intersection counts
    _computeExactBandDistanceField(bandwidth, intersectionCounts);

    // then propagate distances outwards to the rest of the grid
    _propagateDistanceField();

    // then figure out signs (inside/outside) from intersection counts
    _computeDistanceFieldSigns(intersectionCounts);
}

void MeshLevelSet::calculateUnion(MeshLevelSet &levelset) {
    int oi, oj, ok;
    levelset.getGridDimensions(&oi, &oj, &ok);
    FLUIDSIM_ASSERT(oi == _isize && oj == _jsize && ok == _ksize);

    TriangleMesh *omesh = levelset.getTriangleMesh();
    int indexOffset = _mesh.vertices.size();
    _mesh.vertices.insert(_mesh.vertices.end(), omesh->vertices.begin(), omesh->vertices.end());

    Triangle t;
    _mesh.triangles.reserve(_mesh.triangles.size() + omesh->triangles.size());
    for (size_t i = 0; i < omesh->triangles.size(); i++) {
        t = omesh->triangles[i];
        t.tri[0] += indexOffset;
        t.tri[1] += indexOffset;
        t.tri[2] += indexOffset;
        _mesh.triangles.push_back(t);
    }

    for(int k = 0; k < _phi.depth; k++) {
        for(int j = 0; j < _phi.height; j++) {
            for(int i = 0; i < _phi.width; i++) {
                if (levelset(i, j, k) < _phi(i, j, k)) {
                    _phi.set(i, j, k, levelset(i, j, k));

                    int tidx = levelset.getClosestTriangleIndex(i, j, k) + indexOffset;
                    _closestTriangles.set(i, j, k, tidx);
                }
            }
        }
    }

}

void MeshLevelSet::negate() {
    for(int k = 0; k < _phi.depth; k++) {
        for(int j = 0; j < _phi.height; j++) {
            for(int i = 0; i < _phi.width; i++) {
                _phi.set(i, j, k, -_phi(i, j, k));
            }
        }
    }
}

void MeshLevelSet::_computeExactBandDistanceField(int bandwidth,
                                                  Array3d<int> &intersectionCounts) {
    int isize = _phi.width;
    int jsize = _phi.height;
    int ksize = _phi.depth;
    _phi.fill((float)((isize + jsize + ksize) * _dx));
    _closestTriangles.fill(-1);
    intersectionCounts.fill(0);

    Triangle t;
    double invdx = 1.0 / _dx;
    for(size_t tidx = 0; tidx < _mesh.triangles.size(); tidx++) {
        t = _mesh.triangles[tidx];
        vmath::vec3 p = _mesh.vertices[t.tri[0]];
        vmath::vec3 q = _mesh.vertices[t.tri[1]];
        vmath::vec3 r = _mesh.vertices[t.tri[2]]; 

        double fip = (double)p.x * invdx;
        double fjp = (double)p.y * invdx; 
        double fkp = (double)p.z * invdx;

        double fiq = (double)q.x * invdx;
        double fjq = (double)q.y * invdx;
        double fkq = (double)q.z * invdx;

        double fir = (double)r.x * invdx;
        double fjr = (double)r.y * invdx;
        double fkr = (double)r.z * invdx;

        int i0 = _clamp(int(fmin(fip, fmin(fiq, fir))) - bandwidth, 0, isize - 1);
        int j0 = _clamp(int(fmin(fjp, fmin(fjq, fjr))) - bandwidth, 0, jsize - 1);
        int k0 = _clamp(int(fmin(fkp, fmin(fkq, fkr))) - bandwidth, 0, ksize - 1);

        int i1 = _clamp(int(fmax(fip, fmax(fiq, fir))) + bandwidth + 1, 0, isize - 1);
        int j1 = _clamp(int(fmax(fjp, fmax(fjq, fjr))) + bandwidth + 1, 0, jsize - 1);
        int k1 = _clamp(int(fmax(fkp, fmax(fkq, fkr))) + bandwidth + 1, 0, ksize - 1);

        for(int k = k0; k <= k1; k++) {
            for(int j = j0; j<= j1; j++) { 
                for(int i = i0; i <= i1; i++){
                    vmath::vec3 gpos = Grid3d::GridIndexToPosition(i, j, k, _dx);
                    float d = _pointToTriangleDistance(gpos, p, q, r);
                    if (d < _phi(i, j, k)) {
                        _phi.set(i, j, k, d);
                        _closestTriangles.set(i, j, k, tidx);
                    }
                }
            }
        }

        // and do intersection counts
        j0 = _clamp((int)std::ceil(fmin(fjp, fmin(fjq, fjr))), 0, jsize - 1);
        k0 = _clamp((int)std::ceil(fmin(fkp, fmin(fkq, fkr))), 0, ksize - 1);

        j1 = _clamp((int)std::floor(fmax(fjp, fmax(fjq, fjr))), 0, jsize - 1);
        k1 = _clamp((int)std::floor(fmax(fkp, fmax(fkq, fkr))), 0, ksize - 1);

        for(int k = k0; k <= k1; k++) {
            for(int j = j0; j <= j1; j++){
                double a, b, c;
                if (_getBarycentricCoordinates(j, k, fjp, fkp, fjq, fkq, fjr, fkr, &a, &b, &c)) {
                    double fi = a * fip + b * fiq + c * fir;
                    int interval = int(std::ceil(fi));
                    if (interval < 0) {
                        intersectionCounts.add(0, j, k, 1);
                    } else if (interval < isize) {
                        intersectionCounts.add(interval, j, k, 1);
                    }
                }
            }
        }
    }
}

void MeshLevelSet::_propagateDistanceField() {
    int isize = _phi.width;
    int jsize = _phi.height;
    int ksize = _phi.depth;

    std::vector<GridIndex> queue;
    queue.reserve(isize * jsize * ksize);
    Array3d<bool> searchGrid(isize, jsize, ksize, false);
    for(int k = 0; k < ksize; k++) {
        for(int j = 0; j < jsize; j++) {
            for(int i = 0; i < isize; i++) {
                if (_closestTriangles(i, j, k) != -1) {
                    searchGrid.set(i, j, k, true);
                    queue.push_back(GridIndex(i, j, k));
                }
            }
        }
    }

    int unknownidx = queue.size();
    int startidx = 0;
    GridIndex g, n, nbs[6];
    while (startidx < (int)queue.size()) {
        g = queue[startidx];
        startidx++;

        Grid3d::getNeighbourGridIndices6(g, nbs);
        for (int nidx = 0; nidx < 6; nidx++) {
            n = nbs[nidx];
            if (Grid3d::isGridIndexInRange(n, isize, jsize, ksize) && !searchGrid(n)) {
                searchGrid.set(n, true);
                queue.push_back(n);
            }
        }
    }

    vmath::vec3 gpos;
    Triangle t;
    startidx = unknownidx;
    while (startidx < (int)queue.size()) {
        g = queue[startidx];
        startidx++;

        gpos = Grid3d::GridIndexToPosition(g, _dx);
        Grid3d::getNeighbourGridIndices6(g, nbs);
        for (int nidx = 0; nidx < 6; nidx++) {
            n = nbs[nidx];
            if (Grid3d::isGridIndexInRange(n, isize, jsize, ksize) && _closestTriangles(n) != -1) {
                t = _mesh.triangles[_closestTriangles(n)];
                double dist = _pointToTriangleDistance(gpos, _mesh.vertices[t.tri[0]], 
                                                             _mesh.vertices[t.tri[1]], 
                                                             _mesh.vertices[t.tri[2]]);
                if (dist < _phi(g)) {
                    _phi.set(g, (float)dist);
                    _closestTriangles.set(g, _closestTriangles(n));
                }
            }
        }
    }
}

void MeshLevelSet::_computeDistanceFieldSigns(Array3d<int> &intersectionCounts) {
    int isize = _phi.width;
    int jsize = _phi.height;
    int ksize = _phi.depth;

    for(int k = 0; k < ksize; k++) {
        for(int j = 0; j < jsize; j++){
            int tcount = 0;
            for(int i = 0; i < isize; i++){
                tcount += intersectionCounts(i, j, k);
                if(tcount % 2 == 1) {
                    _phi.set(i, j, k, -_phi(i, j, k));
                }
            }
        }
    }
}

// find distance x0 is from triangle x1-x2-x3
float MeshLevelSet::_pointToTriangleDistance(vmath::vec3 x0, vmath::vec3 x1, 
                                                             vmath::vec3 x2, 
                                                             vmath::vec3 x3) {
    // first find barycentric coordinates of closest point on infinite plane
    vmath::vec3 x13 = x1 - x3;
    vmath::vec3 x23 = x2 - x3;
    vmath::vec3 x03 = x0 - x3;

    float m13 = vmath::lengthsq(x13);
    float m23 = vmath::lengthsq(x23); 
    float d = vmath::dot(x13, x23);
    float invdet = 1.0f / fmax(m13 * m23 - d * d, 1e-30f);
    float a = vmath::dot(x13, x03);
    float b = vmath::dot(x23, x03);

    // the barycentric coordinates themselves
    float w23 = invdet * (m23 * a - d * b);
    float w31 = invdet * (m13 * b - d * a);
    float w12 = 1 - w23 - w31;
    if (w23 >= 0 && w31 >= 0 && w12 >= 0) { // if we're inside the triangle
        return vmath::length(x0 - (w23 * x1 + w31 * x2 + w12 * x3)); 
    } else { 
        // we have to clamp to one of the edges
        if (w23 > 0) { 
            // this rules out edge 2-3 for us
            float d1 = _pointToSegmentDistance(x0, x1, x2);
            float d2 = _pointToSegmentDistance(x0, x1, x3);
            return fmin(d1, d2);
        } else if(w31>0) { 
            // this rules out edge 1-3
            float d1 = _pointToSegmentDistance(x0, x1, x2);
            float d2 = _pointToSegmentDistance(x0, x2, x3);
            return fmin(d1, d2);
        } else { 
            // w12 must be >0, ruling out edge 1-2
            float d1 = _pointToSegmentDistance(x0, x1, x3);
            float d2 = _pointToSegmentDistance(x0, x2, x3);
            return fmin(d1, d2);
        }
    }
}

// robust test of (x0,y0) in the triangle (x1,y1)-(x2,y2)-(x3,y3)
// if true is returned, the barycentric coordinates are set in a,b,c.
bool MeshLevelSet::_getBarycentricCoordinates(
            double x0, double y0, 
            double x1, double y1, double x2, double y2, double x3, double y3,
            double *a, double *b, double *c) {
    x1 -= x0; 
    x2 -= x0; 
    x3 -= x0;
    y1 -= y0; 
    y2 -= y0; 
    y3 -= y0;

    double oa;
    int signa = _orientation(x2, y2, x3, y3, &oa);
    if (signa == 0) {
        return false;
    }

    double ob;
    int signb = _orientation(x3, y3, x1, y1, &ob);
    if(signb != signa) {
        return false;
    }

    double oc;
    int signc = _orientation(x1, y1, x2, y2, &oc);
    if(signc != signa) {
        return false;
    }

    double sum = oa + ob + oc;
    assert(sum != 0); // if the SOS signs match and are nonkero, there's no way all of a, b, and c are zero.
    double invsum = 1.0 / sum;

    *a = oa * invsum;
    *b = ob * invsum;
    *c = oc * invsum;

    return true;
}

// find distance x0 is from segment x1-x2
float MeshLevelSet::_pointToSegmentDistance(vmath::vec3 x0, vmath::vec3 x1, vmath::vec3 x2) {
    vmath::vec3 dx = x2 - x1;
    double m2 = vmath::lengthsq(dx);
    // find parameter value of closest point on segment
    float s12 = (float)(vmath::dot(x2 - x0, dx) / m2);
    if (s12 < 0) {
        s12 = 0;
    } else if (s12 > 1) {
        s12 = 1;
    }

    // and find the distance
    return vmath::length(x0 - (s12 * x1 + (1 - s12) * x2));
}

// calculate twice signed area of triangle (0,0)-(x1,y1)-(x2,y2)
// return an SOS-determined sign (-1, +1, or 0 only if it's a truly degenerate triangle)
int MeshLevelSet::_orientation(double x1, double y1, double x2, double y2, 
                              double *twiceSignedArea) {
    *twiceSignedArea = y1 * x2 - x1 * y2;
    if(*twiceSignedArea > 0) {
        return 1;
    } else if (*twiceSignedArea < 0) {
        return -1;
    } else if (y2 > y1) {
        return 1;
    } else if (y2 < y1) {
        return -1;
    } else if (x1 > x2) {
        return 1;
    } else if (x1 < x2) {
        return -1; 
    } else { 
        return 0; // only true when x1==x2 and y1==y2
    }
}
