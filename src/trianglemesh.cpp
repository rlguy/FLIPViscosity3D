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
#include "trianglemesh.h"

TriangleMesh::TriangleMesh() {
}

TriangleMesh::~TriangleMesh() {
}

int TriangleMesh::numVertices() {
    return (int)vertices.size();
}

int TriangleMesh::numFaces() {
    return (int)triangles.size();
}

bool TriangleMesh::loadPLY(std::string PLYFilename) {

    std::ifstream file(PLYFilename.c_str(), std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        return false;
    }

    std::string header;
    bool success = _getPLYHeader(&file, &header);
    if (!success) {
        return false;
    }

    success = _loadPLYVertexData(&file, header);
    if (!success) {
        return false;
    }

    success = _loadPLYTriangleData(&file, header);
    if (!success) {
        return false;
    }

    return true;
}

bool TriangleMesh::loadBOBJ(std::string BOBJFilename) {
    std::ifstream file(BOBJFilename.c_str(), std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        return false;
    }

    int numverts;
    file.read((char *)&numverts, sizeof(int));
    if (!file.good() || numverts < 0) {
        return false;
    }
    
    int binsize = 3 * numverts * sizeof(float);
    std::vector<vmath::vec3> vertices(numverts);
    if (numverts > 0) {
        file.read((char *)vertices.data(), binsize);
        if (!file.good()) {
            return false;
        }
    }

    int numfaces;
    file.read((char *)&numfaces, sizeof(int));
    if (!file.good() || numfaces < 0) {
        return false;
    }

    binsize = 3 * numfaces * sizeof(int);
    std::vector<Triangle> triangles(numfaces);
    if (numfaces > 0) {
        file.read((char *)triangles.data(), binsize);
        if (!file.good()) {
            return false;
        }
    }

    this->vertices = vertices;
    this->triangles = triangles;

    return true;
}

// method of loading OBJ from:
// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-7-model-loading/
bool TriangleMesh::loadOBJ(std::string filename) {
    std::vector<vmath::vec3> temp_vertices;
    std::vector<vmath::vec3> temp_normals;
    std::vector<Triangle> temp_triangles;

    FILE * file;
    file = fopen(filename.c_str(), "rb");
    if(file == NULL){
        printf("Unable to open the OBJ file!\n");
        return false;
    }

    for (;;) {
        char lineHeader[128];
        // read the first word of the line
        int res = fscanf(file, "%s", lineHeader);
        if (res == EOF) {
            break; // EOF = End Of File. Quit the loop.
        }
        
        if (strcmp( lineHeader, "v") == 0 ){
            vmath::vec3 vertex;
            fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z );
            temp_vertices.push_back(vertex);
        } else if (strcmp( lineHeader, "vn" ) == 0) {
            vmath::vec3 normal;
            fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z );
            temp_normals.push_back(normal);
        } else if (strcmp( lineHeader, "f" ) == 0) {
            long start = ftell(file);
            unsigned int vertexIndex[3];
            unsigned int uvIndex[3];
            unsigned int normalIndex[3];
            int matches = fscanf(file, "%d %d %d\n", &vertexIndex[0], &vertexIndex[1], &vertexIndex[2]);

            if (matches != 3){
                long diff = ftell(file) - start;
                fseek (file, -diff , SEEK_CUR);
                start = ftell(file);
                matches = fscanf(file, "%d//%d %d//%d %d//%d\n", 
                                 &vertexIndex[0], &normalIndex[0], 
                                 &vertexIndex[1], &normalIndex[1],
                                 &vertexIndex[2], &normalIndex[2]);

                if (matches != 6) {
                    long diff = ftell(file) - start;
                    fseek (file, -diff , SEEK_CUR);
                    matches = fscanf(file, "%d/%d/%d %d/%d/%d %d/%d/%d\n", 
                                     &vertexIndex[0], &normalIndex[0], &uvIndex[0],
                                     &vertexIndex[1], &normalIndex[1], &uvIndex[1],
                                     &vertexIndex[2], &normalIndex[2], &uvIndex[2]);

                    if (matches != 9) {
                        printf("File can't be read by our simple parser : ( Try exporting with other options\n");
                        return false;
                    }
                }
            }

            Triangle t = Triangle(vertexIndex[0] - 1,
                                  vertexIndex[1] - 1,
                                  vertexIndex[2] - 1);
            temp_triangles.push_back(t);
        }
    }

    fclose(file);

    vertices.clear();
    triangles.clear();
    vertices.insert(vertices.end(), temp_vertices.begin(), temp_vertices.end());
    triangles.insert(triangles.end(), temp_triangles.begin(), temp_triangles.end());

    normals.clear();
    if (normals.size() == vertices.size()) {
        normals.insert(normals.end(), normals.begin(), normals.end());
    }

    return true;
}

void TriangleMesh::writeMeshToPLY(std::string filename) {
    // Header format:
    /*
        ply
        format binary_little_endian 1.0
        element vertex FILL_IN_NUMBER_OF_VERTICES
        property float x
        property float y
        property float z
        element face FILL_IN_NUMBER_OF_FACES
        property list uchar int vertex_index
        end_header
    */
    
    char header1[51] = {'p', 'l', 'y', '\n', 
                        'f', 'o', 'r', 'm', 'a', 't', ' ', 'b', 'i', 'n', 'a', 'r', 'y', '_', 'l', 
                        'i', 't', 't', 'l', 'e', '_', 'e', 'n', 'd', 'i', 'a', 'n', ' ', '1', '.', '0', '\n',
                        'e', 'l', 'e', 'm', 'e', 'n', 't', ' ', 'v', 'e', 'r', 't', 'e', 'x', ' '};
                    
    char header2[65] = {'\n', 'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'x', '\n',
                              'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'y', '\n',
                              'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'z', '\n',
                              'e', 'l', 'e', 'm', 'e', 'n', 't', ' ', 'f', 'a', 'c', 'e', ' '};

    char header2color[125] = {'\n', 'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'x', '\n',
                                    'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'y', '\n',
                                    'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'z', '\n',
                                    'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'u', 'c', 'h', 'a', 'r', ' ', 'r', 'e', 'd', '\n',
                                    'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'u', 'c', 'h', 'a', 'r', ' ', 'g', 'r', 'e', 'e', 'n', '\n',
                                    'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'u', 'c', 'h', 'a', 'r', ' ', 'b', 'l', 'u', 'e', '\n',
                                    'e', 'l', 'e', 'm', 'e', 'n', 't', ' ', 'f', 'a', 'c', 'e', ' '};
                          
    char header3[49] = {'\n', 'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'l', 'i', 's', 't', ' ', 
                              'u', 'c', 'h', 'a', 'r', ' ', 'i', 'n', 't', ' ', 
                              'v', 'e', 'r', 't', 'e', 'x', '_', 'i', 'n', 'd', 'e', 'x', '\n',
                              'e', 'n', 'd', '_', 'h', 'e', 'a', 'd', 'e', 'r', '\n'};

    bool isColorEnabled = vertices.size() == vertexcolors.size();

    std::string vertstring = _toString(vertices.size());
    std::string facestring = _toString(triangles.size());
    int vertdigits = (int)vertstring.length();
    int facedigits = (int)facestring.length();

    int offset = 0;
    int headersize;
    if (isColorEnabled) {
        headersize = 51 + vertdigits + 125 + facedigits + 49;
    } else {
        headersize = 51 + vertdigits + 65 + facedigits + 49;
    }

    int binsize;
    if (isColorEnabled) {
        binsize = headersize + 3*(sizeof(float)*(int)vertices.size() + sizeof(unsigned char)*(int)vertices.size())
                             + (sizeof(unsigned char) + 3*sizeof(int))*(int)triangles.size();
    } else {
        binsize = headersize + 3*sizeof(float)*(int)vertices.size()
                             + (sizeof(unsigned char) + 3*sizeof(int))*(int)triangles.size();
    }
    char *bin = new char[binsize];

    memcpy(bin + offset, header1, 51);
    offset += 51;
    memcpy(bin + offset, vertstring.c_str(), vertdigits*sizeof(char));
    offset += vertdigits*sizeof(char);

    if (isColorEnabled) { 
        memcpy(bin + offset, header2color, 125);
        offset += 125;
    } else {
        memcpy(bin + offset, header2, 65);
        offset += 65;
    }

    memcpy(bin + offset, facestring.c_str(), facedigits*sizeof(char));
    offset += facedigits*sizeof(char);
    memcpy(bin + offset, header3, 49);
    offset += 49;

    if (isColorEnabled) {
        float *vertdata = new float[3*vertices.size()];
        vmath::vec3 v;
        for (unsigned int i = 0; i < vertices.size(); i++) {
            v = vertices[i];
            vertdata[3*i] = v.x;
            vertdata[3*i + 1] = v.y;
            vertdata[3*i + 2] = v.z;
        }

        unsigned char *colordata = new unsigned char[3*vertexcolors.size()];
        vmath::vec3 c;
        for (unsigned int i = 0; i < vertexcolors.size(); i++) {
            c = vertexcolors[i];
            colordata[3*i] = (unsigned char)((c.x/1.0)*255.0);
            colordata[3*i + 1] = (unsigned char)((c.y/1.0)*255.0);
            colordata[3*i + 2] = (unsigned char)((c.z/1.0)*255.0);
        }

        int vertoffset = 0;
        int coloroffset = 0;
        int vertsize = 3*sizeof(float);
        int colorsize = 3*sizeof(unsigned char);
        for (unsigned int i = 0; i < vertices.size(); i++) {
            memcpy(bin + offset, vertdata + vertoffset, vertsize);
            offset += vertsize;
            vertoffset += 3;

            memcpy(bin + offset, colordata + coloroffset, colorsize);
            offset += colorsize;
            coloroffset += 3;
        }

        delete[] colordata;
        delete[] vertdata;
    } else {
        float *vertdata = new float[3*vertices.size()];
        vmath::vec3 v;
        for (unsigned int i = 0; i < vertices.size(); i++) {
            v = vertices[i];
            vertdata[3*i] = v.x;
            vertdata[3*i + 1] = v.y;
            vertdata[3*i + 2] = v.z;
        }
        memcpy(bin + offset, vertdata, 3*sizeof(float)*vertices.size());
        offset += 3*sizeof(float)*(int)vertices.size();
        delete[] vertdata;
    }

    Triangle t;
    int verts[3];
    for (unsigned int i = 0; i < triangles.size(); i++) {
        t = triangles[i];
        verts[0] = t.tri[0];
        verts[1] = t.tri[1];
        verts[2] = t.tri[2];

        bin[offset] = 0x03;
        offset += sizeof(unsigned char);

        memcpy(bin + offset, verts, 3*sizeof(int));
        offset += 3*sizeof(int);
    }

    std::ofstream erasefile;
    erasefile.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
    erasefile.close();

    std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);
    file.write(bin, binsize);
    file.close();

    delete[] bin;
}

int TriangleMesh::_numDigitsInInteger(int num) {
    if (num == 0) {
        return 1;
    }

    int count = 0;
    while(num != 0) {
        num /= 10;
        count++;
    }

    return count;
}

void TriangleMesh::writeMeshToBOBJ(std::string filename) {
    std::ofstream erasefile;
    erasefile.open(filename, std::ofstream::out | std::ofstream::trunc);
    erasefile.close();

    std::ofstream bobj(filename.c_str(), std::ios::out | std::ios::binary);

    int numVertices = (int)vertices.size();
    bobj.write((char *)&numVertices, sizeof(int));

    int binsize = 3 * numVertices * sizeof(float);
    bobj.write((char *)vertices.data(), binsize);

    int numTriangles = (int)triangles.size();
    bobj.write((char *)&numTriangles, sizeof(int));

    binsize = 3 * numTriangles * sizeof(int);
    bobj.write((char *)triangles.data(), binsize);

    bobj.close();
}

void TriangleMesh::writeMeshToOBJ(std::string filename) {
    std::ostringstream str;

    str << "# OBJ file format with ext .obj" << std::endl;
    str << "# vertex count = " << vertices.size() << std::endl;
    str << "# face count = " << triangles.size() << std::endl;

    vmath::vec3 p;
    for (unsigned int i = 0; i < vertices.size(); i++) {
        p = vertices[i];
        str << "v " << p.x << " " << p.y << " " << p.z << std::endl;
    }

    if (normals.size() == vertices.size()) {
        vmath::vec3 n;
        for (unsigned int i = 0; i < normals.size(); i++) {
            n = normals[i];
            str << "vn " << n.x << " " << n.y << " " << n.z << std::endl;
        }
    }

    Triangle t;
    int v1, v2, v3;
    for (unsigned int i = 0; i < triangles.size(); i++) {
        t = triangles[i];
        v1 = t.tri[0] + 1;
        v2 = t.tri[1] + 1;
        v3 = t.tri[2] + 1;

        str << "f " << v1 << "//" << v1 << " " <<
            v2 << "//" << v2 << " " <<
            v3 << "//" << v3 << std::endl;
    }

    std::ofstream out(filename);
    out << str.str();
    out.close();
}

void TriangleMesh::translate(vmath::vec3 t) {
    for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i] += t;
    }
}

bool TriangleMesh::_getPLYHeader(std::ifstream *file, std::string *header) {
    file->seekg(0, std::ios_base::beg);

    int maxHeaderSize = 2048;
    char headerBufferChars[2048];
    file->read(headerBufferChars, maxHeaderSize);
    std::string headerBufferString(headerBufferChars, 2048);

    std::string endHeaderString("end_header\n");

    std::size_t match = headerBufferString.find(endHeaderString);
    if (match == std::string::npos) {
        return false;
    }

    *header = headerBufferString.substr(0, match + endHeaderString.size());

    return true;
}

bool TriangleMesh::_getElementNumberInPlyHeader(std::string &header, 
                                                std::string &element, int *n) {
    std::size_t match = header.find(element);
    if (match == std::string::npos) {
        return false;
    }

    int startidx = (int)match + (int)element.size();
    int endidx = 0;
    bool numberFound = false;

    for (unsigned int i = startidx; i < header.size(); i++) {
        if (header[i] == '\n') {
            endidx = i - 1;
            numberFound = true;
            break;
        }
    }

    if (!numberFound) {
        return false;
    }

    std::string numberString = header.substr(startidx, endidx - startidx + 1);
    std::istringstream ss(numberString);
    ss >> *n;

    if (ss.fail()) {
        return false;
    }

    return true;
}

bool TriangleMesh::_getNumVerticesInPLYHeader(std::string &header, int *n) {
    std::string vertexString("element vertex ");
    bool success = _getElementNumberInPlyHeader(header, vertexString, n);

    return success;
}

bool TriangleMesh::_getNumFacesInPLYHeader(std::string &header, int *n) {
    std::string faceString("element face ");
    bool success = _getElementNumberInPlyHeader(header, faceString, n);

    return success;
}

bool TriangleMesh::_isVertexColorsEnabledInPLYHeader(std::string &header) {
    std::string colorString("property uchar red\nproperty uchar green\nproperty uchar blue\n");
    std::size_t match = header.find(colorString);
    return match != std::string::npos;
}

bool TriangleMesh::_loadPLYVertexData(std::ifstream *file, std::string &header) {
    int numVertices;
    bool success = _getNumVerticesInPLYHeader(header, &numVertices);
    if (!success) {
        return false;
    }

    if (numVertices == 0) {
        return true;
    }

    bool isColorEnabled = _isVertexColorsEnabledInPLYHeader(header);

    int vertexSize = 3*sizeof(float);
    if (isColorEnabled) {
        vertexSize = 3*sizeof(float) + 3*sizeof(char);
    }

    int vertexDataSize = numVertices*vertexSize;
    int vertexDataOffset = (int)header.size();

    file->seekg(vertexDataOffset, std::ios_base::beg);
    char *vertexData = new char[vertexDataSize];
    if (!file->read(vertexData, vertexDataSize)) {
        return false;
    }

    vertices.reserve(numVertices);
    if (isColorEnabled) {
        vertexcolors.reserve(numVertices);

        vmath::vec3 p;
        int offset = 0;
        for (int i = 0; i < numVertices; i++) {
            memcpy(&p, vertexData + offset, 3*sizeof(float));
            offset += 3*sizeof(float);

            unsigned char r = vertexData[offset + 0];
            unsigned char g = vertexData[offset + 1];
            unsigned char b = vertexData[offset + 2];
            offset += 3*sizeof(char);

            vertices.push_back(p);
            vertexcolors.push_back(vmath::vec3(r / 255.0f, g / 255.0f, b / 255.0f));
        }
    } else {
        vertices.assign((vmath::vec3*)vertexData, (vmath::vec3*)vertexData + numVertices);
    }
    delete[] vertexData;

    return true;
}

bool TriangleMesh::_loadPLYTriangleData(std::ifstream *file, std::string &header) {
    int numVertices;
    bool success = _getNumVerticesInPLYHeader(header, &numVertices);
    if (!success) {
        return false;
    }

    bool isColorEnabled = _isVertexColorsEnabledInPLYHeader(header);

    int vertexSize = 3*sizeof(float);
    if (isColorEnabled) {
        vertexSize = 3*sizeof(float) + 3*sizeof(char);
    }

    int vertexDataSize = numVertices*vertexSize;
    int vertexDataOffset = (int)header.size();

    int numFaces;
    success = _getNumFacesInPLYHeader(header, &numFaces);
    if (!success) {
        return false;
    }

    if (numFaces == 0) {
        return true;
    }

    int faceSize = sizeof(char) + 3*sizeof(int);
    int faceDataSize = numFaces*faceSize;
    int faceDataOffset = vertexDataOffset + vertexDataSize;

    file->seekg(faceDataOffset, std::ios_base::beg);
    char *faceData = new char[faceDataSize];
    if (!file->read(faceData, faceDataSize)) {
        return false;
    }

    int offset = 0;
    Triangle t;
    triangles.reserve(numFaces);
    for (int i = 0; i < numFaces; i++) {
        unsigned int faceverts = faceData[offset];
        offset += sizeof(char);

        if (faceverts != 0x03) {
            return false;
        }

        memcpy(&(t.tri), faceData + offset, 3*sizeof(int));
        offset += 3*sizeof(int);

        if (t.tri[0] < 0 || t.tri[0] >= numVertices || 
            t.tri[1] < 0 || t.tri[1] >= numVertices || 
            t.tri[2] < 0 || t.tri[2] >= numVertices) {
            return false;
        }
        triangles.push_back(t);
    }

    delete[] faceData;

    return true;
}
