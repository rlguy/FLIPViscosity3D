#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "fluidsimulation.h"
#include "trianglemesh.h"

void export_particles(int frame, std::vector<FluidParticle> &particles) {
    TriangleMesh mesh;
    float scale = 10.0f;
    for(unsigned int i = 0; i < particles.size(); i++) {
        mesh.vertices.push_back(scale * particles[i].position);
    }

    std::ostringstream ss;
    ss << frame;
    std::string currentFrame = ss.str();
    currentFrame.insert(currentFrame.begin(), 6 - currentFrame.size(), '0');
    std::string filename = "C:/Users/Ryan/Downloads/variational_test_cache_flip_fluid2/bakefiles/" + currentFrame + ".bobj";
    mesh.writeMeshToBOBJ(filename);

    std::cout << "Exporting Particles: " << filename << std::endl;
}

int main() {

    int isize = 64;
    int jsize = 64;
    int ksize = 64;
    float dx = 1.0f / fmax(fmax(isize, jsize), ksize);
    FluidSimulation fluidsim;

    std::cout << "Initializing Data" << std::endl;
    fluidsim.initialize(isize, jsize, ksize, dx);
    
    std::cout << "Initializing Boundary" << std::endl;
    TriangleMesh boundaryMesh;
    bool success = boundaryMesh.loadPLY("sample_meshes/sphere_large.ply");
    FLUIDSIM_ASSERT(success);
    bool invertMesh = true;
    fluidsim.addBoundary(boundaryMesh, invertMesh);
    
    std::cout << "Initializing Liquid" << std::endl;
    TriangleMesh liquidMesh;
    success = liquidMesh.loadPLY("sample_meshes/stanford_bunny.ply");
    FLUIDSIM_ASSERT(success);
    fluidsim.addLiquid(liquidMesh);

    std::cout << "Initializing Settings" << std::endl;
    fluidsim.setViscosity(5.0f);
    fluidsim.setGravity(0.0f, -9.81f, 0.0f);

    int num_frames = 300;
    float timestep = 0.01f;
    for(int frame = 0; frame < num_frames; frame++) {
        export_particles(frame, fluidsim.particles);

        std::cout << "Simulation Liquid" << std::endl;
        fluidsim.advance(timestep);
    }

    return 0;
}
