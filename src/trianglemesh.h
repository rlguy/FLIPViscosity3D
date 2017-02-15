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
#ifndef TRIANGLEMESH_H
#define TRIANGLEMESH_H

#include <stdlib.h>
#include <queue>
#include <vector>
#include <sstream>
#include <fstream>
#include <string.h>
#include <algorithm>

#include "triangle.h"
#include "vmath.h"

class TriangleMesh
{
public:
    TriangleMesh();
    ~TriangleMesh();

    bool loadPLY(std::string PLYFilename);
    bool loadBOBJ(std::string BOBJFilename);
    bool loadOBJ(std::string filename);
    void writeMeshToPLY(std::string filename);
    void writeMeshToBOBJ(std::string filename);
    void writeMeshToOBJ(std::string filename);

    int numVertices();
    int numFaces();
    int numTriangles() { return numFaces(); }
    void translate(vmath::vec3 t);

    std::vector<vmath::vec3> vertices;
    std::vector<vmath::vec3> vertexcolors;  // r, g, b values in range [0.0, 1.0]
    std::vector<vmath::vec3> normals;
    std::vector<Triangle> triangles;

private:
    bool _getPLYHeader(std::ifstream *file, std::string *header);
    bool _getElementNumberInPlyHeader(std::string &header, 
                                      std::string &element, int *n);
    bool _getNumVerticesInPLYHeader(std::string &header, int *n);
    bool _getNumFacesInPLYHeader(std::string &header, int *n);
    bool _isVertexColorsEnabledInPLYHeader(std::string &header);
    bool _loadPLYVertexData(std::ifstream *file, std::string &header);
    bool _loadPLYTriangleData(std::ifstream *file, std::string &header);

    int _numDigitsInInteger(int num);

    template<class T>
    std::string _toString(T item) {
        std::ostringstream sstream;
        sstream << item;

        return sstream.str();
    }
};

#endif
