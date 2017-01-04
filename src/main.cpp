#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "fluidsim.h"
#include "trianglemesh.h"
#include "meshlevelset.h"

using namespace std;

float sphere_phi(vmath::vec3 position, vmath::vec3 centre, float radius) {
    return vmath::length(position - centre) - radius;
}

float liquid_phi(vmath::vec3 position) {
    return sphere_phi(position, vmath::vec3(0.55f, 0.55f, 0.4f), 0.23f);
}

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

TriangleMesh getMeshFromAABB(AABB bbox) {
    vmath::vec3 p = bbox.position;
    std::vector<vmath::vec3> verts{
        vmath::vec3(p.x, p.y, p.z),
        vmath::vec3(p.x, p.y, p.z),
        vmath::vec3(p.x + bbox.width, p.y, p.z),
        vmath::vec3(p.x + bbox.width, p.y, p.z + bbox.depth),
        vmath::vec3(p.x, p.y, p.z + bbox.depth),
        vmath::vec3(p.x, p.y + bbox.height, p.z),
        vmath::vec3(p.x + bbox.width, p.y + bbox.height, p.z),
        vmath::vec3(p.x + bbox.width, p.y + bbox.height, p.z + bbox.depth),
        vmath::vec3(p.x, p.y + bbox.height, p.z + bbox.depth),
    };

    std::vector<Triangle> tris{
        Triangle(0, 1, 2), Triangle(0, 2, 3), Triangle(4, 7, 6), Triangle(4, 6, 5),
        Triangle(0, 3, 7), Triangle(0, 7, 4), Triangle(1, 5, 6), Triangle(1, 6, 2),
        Triangle(0, 4, 5), Triangle(0, 5, 1), Triangle(3, 2, 6), Triangle(3, 6, 7)
    };


    TriangleMesh m;
    m.vertices = verts;
    m.triangles = tris;

    return m;
}

//Main testing code
//-------------
int main(int argc, char **argv)
{  

    int grid_resolution = 32;
    float timestep = 0.01f;
    float grid_width = 1;
    float dx = grid_width / (double)grid_resolution;
    FluidSim sim;

    TriangleMesh boundaryMesh;
    bool success = boundaryMesh.loadBOBJ("sphere.bobj");
    FLUIDSIM_ASSERT(success);
    MeshLevelSet boundary(grid_resolution, grid_resolution, grid_resolution, dx);
    boundary.calculateSignedDistanceField(boundaryMesh, 5);
    boundary.negate();
    
    /*
    AABB domain(0.0, 0.0, 0.0, grid_resolution * dx, grid_resolution * dx, grid_resolution * dx);
    AABB outer = domain;
    outer.expand(-1e-6);
    AABB inner = domain;
    inner.expand(-6*dx - 1e-6);

    TriangleMesh dmesh = getMeshFromAABB(outer);
    TriangleMesh imesh = getMeshFromAABB(inner);
    int indexOffset = dmesh.vertices.size();
    dmesh.vertices.insert(dmesh.vertices.end(), imesh.vertices.begin(), imesh.vertices.end());
    for (size_t i = 0; i < imesh.triangles.size(); i++) {
        Triangle t = imesh.triangles[i];
        t.tri[0] += indexOffset;
        t.tri[1] += indexOffset;
        t.tri[2] += indexOffset;
        dmesh.triangles.push_back(t);
    }
    MeshLevelSet boundary(grid_resolution, grid_resolution, grid_resolution, dx);
    boundary.calculateSignedDistanceField(dmesh, 5);
    */

    printf("Initializing data\n");
    sim.initialize(grid_resolution, grid_resolution, grid_resolution, grid_width);
    
    printf("Initializing boundary\n");
    sim.set_boundary(boundary);
    
    printf("Initializing liquid\n");
    sim.set_liquid(liquid_phi);

    int num_frames = 300;
    for(int frame = 0; frame < num_frames; frame++) {
        export_particles(frame, sim.particles);

        //Simulate
        printf("Simulating liquid\n");
        sim.advance(timestep);
        
    }

    return 0;
}


