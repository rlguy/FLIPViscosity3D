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

enum class TriangleMeshFormat : char { 
    ply   = 0x00, 
    bobj  = 0x01
};

class TriangleMesh
{
public:
    TriangleMesh();
    ~TriangleMesh();

    bool loadPLY(std::string PLYFilename);
    bool loadBOBJ(std::string BOBJFilename);
    void writeMeshToPLY(std::string filename);
    void writeMeshToBOBJ(std::string filename);
    static std::string getFileExtension(TriangleMeshFormat fmt);

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
