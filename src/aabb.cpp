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
#include "aabb.h"

AABB::AABB() {
}

AABB::AABB(double x, double y, double z, double w, double h, double d) : 
               position((float)x, (float)y, (float)z), width(w), height(h), depth(d) {
}

AABB::AABB(vmath::vec3 p, double w, double h, double d) : 
               position(p), width(w), height(h), depth(d) {
}

AABB::AABB(vmath::vec3 p1, vmath::vec3 p2) {
    double minx = fmin(p1.x, p2.x);
    double miny = fmin(p1.y, p2.y);
    double minz = fmin(p1.z, p2.z);
    double maxx = fmax(p1.x, p2.x);
    double maxy = fmax(p1.y, p2.y);
    double maxz = fmax(p1.z, p2.z);

    position = vmath::vec3((float)minx, (float)miny, (float)minz);
    width = maxx - minx;
    height = maxy - miny;
    depth = maxz - minz;
}

AABB::AABB(std::vector<vmath::vec3> &points) {
    if (points.size() == 0) {
        return;
    }

    double minx = points[0].x;
    double miny = points[0].y;
    double minz = points[0].z;
    double maxx = points[0].x;
    double maxy = points[0].y;
    double maxz = points[0].z;

    vmath::vec3 p;
    for (unsigned int i = 0; i < points.size(); i++) {
        p = points[i];
        minx = fmin(p.x, minx);
        miny = fmin(p.y, miny);
        minz = fmin(p.z, minz);
        maxx = fmax(p.x, maxx);
        maxy = fmax(p.y, maxy);
        maxz = fmax(p.z, maxz);
    }

    double eps = 1e-9;
    position = vmath::vec3((float)minx, (float)miny, (float)minz);
    width = maxx - minx + eps;
    height = maxy - miny + eps;
    depth = maxz - minz + eps;
}

AABB::AABB(Triangle t, std::vector<vmath::vec3> &vertices) {
    vmath::vec3 points[3] = { vertices[t.tri[0]],
                            vertices[t.tri[1]],
                            vertices[t.tri[2]] };

    double minx = points[0].x;
    double miny = points[0].y;
    double minz = points[0].z;
    double maxx = points[0].x;
    double maxy = points[0].y;
    double maxz = points[0].z;

    vmath::vec3 p;
    for (int i = 0; i < 3; i++) {
        p = points[i];
        minx = fmin(p.x, minx);
        miny = fmin(p.y, miny);
        minz = fmin(p.z, minz);
        maxx = fmax(p.x, maxx);
        maxy = fmax(p.y, maxy);
        maxz = fmax(p.z, maxz);
    }

    double eps = 1e-9;
    position = vmath::vec3((float)minx, (float)miny, (float)minz);
    width = maxx - minx + eps;
    height = maxy - miny + eps;
    depth = maxz - minz + eps;
}

AABB::AABB(GridIndex g, double dx) {
    position = vmath::vec3(g.i*(float)dx, g.j*(float)dx, g.k*(float)dx);
    width = height = depth = dx;
}

AABB::~AABB() {
}

void AABB::expand(double v) {
    double h = 0.5 * v;
    position -= vmath::vec3((float)h, (float)h, (float)h);
    width += v;
    height += v;
    depth += v;
}

bool AABB::isPointInside(vmath::vec3 p) {
    return p.x >= position.x && p.y >= position.y && p.z >= position.z &&
           p.x < position.x + width && p.y < position.y + height && p.z < position.z + depth;
}

bool AABB::isLineIntersecting(vmath::vec3 p1, vmath::vec3 p2) {

    vmath::vec3 min = position;
    vmath::vec3 max = position + vmath::vec3((float)width, (float)height, (float)depth);

    vmath::vec3 d = (p2 - p1) * 0.5f;
    vmath::vec3 e = (max - min) * 0.5f;
    vmath::vec3 c = p1 + d - (min + max) * 0.5f;
    vmath::vec3 ad = vmath::vec3(fabs(d.x), fabs(d.y), fabs(d.z));

    if (fabs(c.x) > e.x + ad.x) {
        return false;
    }
    if (fabs(c.y) > e.y + ad.y) {
        return false;
    }
    if (fabs(c.z) > e.z + ad.z) {
        return false;
    }


    double eps = 10e-9;
    if (fabs(d.y * c.z - d.z * c.y) > e.y * ad.z + e.z * ad.y + eps) {
        return false;
    }
    if (fabs(d.z * c.x - d.x * c.z) > e.z * ad.x + e.x * ad.z + eps) {
        return false;
    }
    if (fabs(d.x * c.y - d.y * c.x) > e.x * ad.y + e.y * ad.x + eps) {
        return false;
    }

    return true;
}

AABB AABB::getIntersection(AABB bbox) {
    vmath::vec3 minp1 = getMinPoint();
    vmath::vec3 minp2 = bbox.getMinPoint();
    vmath::vec3 maxp1 = getMaxPoint();
    vmath::vec3 maxp2 = bbox.getMaxPoint();

    if (minp1.x > maxp2.x || minp1.y > maxp2.y || minp1.z > maxp2.z ||
        maxp1.x < minp2.x || maxp1.y < minp2.y || maxp1.z < minp2.z) {
        return AABB();
    }

    float interminx = fmax(minp1.x, minp2.x);
    float interminy = fmax(minp1.y, minp2.y);
    float interminz = fmax(minp1.z, minp2.z);
    float intermaxx = fmin(maxp1.x, maxp2.x);
    float intermaxy = fmin(maxp1.y, maxp2.y);
    float intermaxz = fmin(maxp1.z, maxp2.z);

    return AABB(vmath::vec3(interminx, interminy, interminz), 
                vmath::vec3(intermaxx, intermaxy, intermaxz));
}

AABB AABB::getUnion(AABB bbox) {
    vmath::vec3 minp1 = getMinPoint();
    vmath::vec3 minp2 = bbox.getMinPoint();
    vmath::vec3 maxp1 = getMaxPoint();
    vmath::vec3 maxp2 = bbox.getMaxPoint();

    float unionminx = fmin(minp1.x, minp2.x);
    float unionminy = fmin(minp1.y, minp2.y);
    float unionminz = fmin(minp1.z, minp2.z);
    float unionmaxx = fmax(maxp1.x, maxp2.x);
    float unionmaxy = fmax(maxp1.y, maxp2.y);
    float unionmaxz = fmax(maxp1.z, maxp2.z);

    return AABB(vmath::vec3(unionminx, unionminy, unionminz), 
                vmath::vec3(unionmaxx, unionmaxy, unionmaxz));
}

vmath::vec3 AABB::getMinPoint() {
    return position;
}

vmath::vec3 AABB::getMaxPoint() {
    return position + vmath::vec3((float)width, (float)height, (float)depth);
}

vmath::vec3 AABB::getNearestPointInsideAABB(vmath::vec3 p) {
    return getNearestPointInsideAABB(p, 1e-6);
}

vmath::vec3 AABB::getNearestPointInsideAABB(vmath::vec3 p, double eps) {
    if (isPointInside(p)) {
        return p;
    }

    vmath::vec3 min = getMinPoint();
    vmath::vec3 max = getMaxPoint();

    p.x = fmax(p.x, min.x);
    p.y = fmax(p.y, min.y);
    p.z = fmax(p.z, min.z);

    p.x = fmin(p.x, max.x - (float)eps);
    p.y = fmin(p.y, max.y - (float)eps);
    p.z = fmin(p.z, max.z - (float)eps);

    return p;
}