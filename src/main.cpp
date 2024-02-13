#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "fluidsimulation.h"
#include "trianglemesh.h"

bool EXPORT_OBJ = true;
bool EXPORT_PLY = false;

void export_particles(int frame, std::vector<FluidParticle> &particles) {
    TriangleMesh mesh;
    float scale = 1.0f;
    for (unsigned int i = 0; i < particles.size(); i++) {
        mesh.vertices.push_back(scale * particles[i].position);
    }

    int zerofill = 4;
    std::ostringstream ss;
    ss << frame;
    std::string framestr = ss.str();
    framestr.insert(framestr.begin(), zerofill - framestr.size(), '0');

    // Write mesh as OBJ file (https://en.wikipedia.org/wiki/Wavefront_.obj_file)
    if (EXPORT_OBJ) {
        std::string filenameOBJ = framestr + ".obj";
        mesh.writeMeshToOBJ(filenameOBJ);
        std::cout << "Exporting Particles to file: " << filenameOBJ << std::endl;
    }

    // Write mesh as PLY file (http://paulbourke.net/dataformats/ply/)
    if (EXPORT_PLY) {
        std::string filenamePLY = framestr + ".ply";
        mesh.writeMeshToPLY(filenamePLY);
        std::cout << "Exporting Particles to file: " << filenamePLY << std::endl;
    }
}

int main() {

    // This simulation example will drop a mass of fluid in the shape of the
    // Stanford Bunny inside of a spherical container

    FluidSimulation fluidsim;

    std::cout << "Initializing Data" << std::endl;
    int isize = 64;
    int jsize = 64;
    int ksize = 64;
    float dx = 1.0f / std::max(std::max(isize, jsize), ksize);
    fluidsim.initialize(isize, jsize, ksize, dx);
    
    std::cout << "Initializing Boundary" << std::endl;
    TriangleMesh boundaryMesh;
    std::string boundaryMeshPath = "sample_meshes/sphere_large.ply";
    bool invertMesh = true;
    if (!boundaryMesh.loadPLY(boundaryMeshPath)) {
        std::cout << "Error loading boundary mesh: " << boundaryMeshPath << std::endl;
        return 0;
    }
    fluidsim.addBoundary(boundaryMesh, invertMesh);
    
    std::cout << "Initializing Liquid" << std::endl;
    TriangleMesh liquidMesh;
    std::string liquidMeshPath = "sample_meshes/stanford_bunny.ply";
    if (!liquidMesh.loadPLY(liquidMeshPath)) {
        std::cout << "Error loading liquid mesh: " << liquidMeshPath << std::endl;
        return 0;
    }
    fluidsim.addLiquid(liquidMesh);

    std::cout << "Initializing Settings" << std::endl;
    fluidsim.setViscosity(5.0f);
    fluidsim.setGravity(0.0f, -9.81f, 0.0f);

    int numFrames = 300;
    float timestep = 0.01f;
    for(int frameno = 0; frameno < numFrames; frameno++) {
        export_particles(frameno, fluidsim.particles);

        std::cout << "Simulating Liquid" << std::endl;
        fluidsim.advance(timestep);
        std::cout << "-------------------------------------------" << std::endl;
    }

    return 0;
}
