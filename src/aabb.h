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
#ifndef AABB_H
#define AABB_H

#include <math.h>
#include <vector>

#include "triangle.h"
#include "array3d.h"
#include "vmath.h"

class AABB
{
public:
    AABB();
    AABB(double x, double y, double z, double width, double height, double depth);
    AABB(vmath::vec3 p, double width, double height, double depth);
    AABB(vmath::vec3 p1, vmath::vec3 p2);
    AABB(std::vector<vmath::vec3> &points);
    AABB(Triangle t, std::vector<vmath::vec3> &vertices);
    AABB(GridIndex g, double dx);
    ~AABB();

    void expand(double v);
    bool isPointInside(vmath::vec3 p);
    bool isLineIntersecting(vmath::vec3 p1, vmath::vec3 p2);
    AABB getIntersection(AABB bbox);
    AABB getUnion(AABB bbox);

    vmath::vec3 getMinPoint();
    vmath::vec3 getMaxPoint();
    vmath::vec3 getNearestPointInsideAABB(vmath::vec3 p, double eps);
    vmath::vec3 getNearestPointInsideAABB(vmath::vec3 p);

    vmath::vec3 position;
    double width = 0.0;
    double height = 0.0;
    double depth = 0.0;

};

#endif
