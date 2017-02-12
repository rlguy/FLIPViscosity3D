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
    float scale = 10.0;
    for(unsigned int i = 0; i < particles.size(); i++) {
        mesh.vertices.push_back(scale * particles[i].position);
    }

    std::ostringstream ss;
    ss << frame;
    std::string currentFrame = ss.str();
    currentFrame.insert(currentFrame.begin(), 6 - currentFrame.size(), '0');
    std::string filename = "C:/Users/Ryan/Downloads/variational_test_cache_flip_fluid2/bakefiles/" + currentFrame + ".bobj";
    mesh.writeMeshToBOBJ(filename);

    std::cout << "Exporting: " << filename << std::endl;
}

int main() {
    int isize = 32;
    int jsize = 32;
    int ksize = 32;
    float dx = 1.0f / (float)isize;
    FluidSimulation fluidsim;

    printf("Initializing data\n");
    fluidsim.initialize(isize, jsize, ksize, dx);
    
    printf("Initializing boundary\n");
    TriangleMesh boundaryMesh;
    bool success = boundaryMesh.loadBOBJ("bunny.bobj");
    boundaryMesh.writeMeshToPLY("bunny.ply");
    FLUIDSIM_ASSERT(success);
    fluidsim.addBoundary(boundaryMesh, false);
    
    printf("Initializing liquid\n");
    TriangleMesh liquidMesh;
    success = liquidMesh.loadBOBJ("sphere.bobj");
    FLUIDSIM_ASSERT(success);
    fluidsim.addLiquid(liquidMesh);

    fluidsim.setViscosity(0.05f);
    fluidsim.setGravity(0.0f, -9.81f, 0.0f);

    int num_frames = 300;
    float timestep = 0.01f;
    for(int frame = 0; frame < num_frames; frame++) {
        export_particles(frame, fluidsim.particles);

        printf("Simulating liquid\n");
        fluidsim.advance(timestep);
    }

    return 0;
}
