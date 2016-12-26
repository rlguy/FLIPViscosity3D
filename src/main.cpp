#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "fluidsim.h"
#include "trianglemesh.h"

using namespace std;

float sphere_phi(vmath::vec3 position, vmath::vec3 centre, float radius) {
   return vmath::length(position - centre) - radius;
}

vmath::vec3 c0(0.5f,0.5f,0.5f);
float rad0 = 0.35f;

float boundary_phi(vmath::vec3 position) {
   return -sphere_phi(position, c0, rad0);
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
   std::string filename = "C:/Users/Ryan/Downloads/variational_test_cache_flip_fluid/bakefiles/" + currentFrame + ".bobj";
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
   sim.set_boundary(boundary_phi);
   
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


