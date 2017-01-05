#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "fluidsim.h"
#include "trianglemesh.h"
#include "meshlevelset.h"
#include "fluidmaterialgrid.h"

using namespace std;

void export_particles(int frame, std::vector<vmath::vec3> &particles) {
    TriangleMesh mesh;
    double scale = 10.0;
    for(unsigned int p = 0; p < particles.size(); ++p) {
        vmath::vec3 v(particles[p][0] * scale, particles[p][1] * scale, particles[p][2] * scale);
        mesh.vertices.push_back(v);
    }

    std::ostringstream ss;
    ss << frame;
    std::string currentFrame = ss.str();
    currentFrame.insert(currentFrame.begin(), 6 - currentFrame.size(), '0');
    std::string filename = "C:/Users/Ryan/Downloads/variational_test_cache_flip_fluid2/bakefiles/" + currentFrame + ".bobj";
    mesh.writeMeshToBOBJ(filename);

    std::cout << "Exporting: " << filename << std::endl;
}

//Main testing code
//-------------
int main(int argc, char **argv)
{  

    int grid_resolution = 32;
    float timestep = 0.01f;
    float grid_width = 1;
    FluidSim sim;

    printf("Initializing data\n");
    sim.initialize(grid_resolution, grid_resolution, grid_resolution, grid_width);
    
    printf("Initializing boundary\n");
    TriangleMesh boundaryMesh;
    bool success = boundaryMesh.loadBOBJ("bunny.bobj");
    FLUIDSIM_ASSERT(success);
    sim.addBoundary(boundaryMesh);
    
    printf("Initializing liquid\n");
    TriangleMesh liquidMesh;
    success = liquidMesh.loadBOBJ("sphere.bobj");
    FLUIDSIM_ASSERT(success);
    sim.addLiquid(liquidMesh);

    int num_frames = 300;
    for(int frame = 0; frame < num_frames; frame++) {
        export_particles(frame, sim.particles);

        //Simulate
        printf("Simulating liquid\n");
        sim.advance(timestep);
    }

    return 0;
}


